"""Submit every group edge explicitly after the range generator lost an edge."""
from pathlib import Path
from types import FunctionType
import argparse
import json
import signal
import time
import numpy as np
import def_photon_line_resolved as previous

native=previous.native;ex=previous.ex;h=previous.h
OUT=previous.OUT.parent/'def-photon-explicit-grid'


def prepare():
    assert not OUT.exists();OUT.mkdir();original=json.loads((previous.OUT/'requests.json').read_text());rows=[]
    for r in original:
        row=dict(r);f=dict(row['fields']);edges=np.geomspace(float(f['egplow']),float(f['egphigh']),1000)
        name=r['name'].replace('phLine','phExact');f.update(mixname=name,egrid='specific',energies=' '.join(f'{x:.17g}' for x in edges))
        row.update(name=name,fields=f);rows.append(row)
    ex.write(OUT/'requests.json',rows);(OUT/'lanl-tops-form.html').write_bytes((native.OUT/'lanl-tops-form.html').read_bytes())
    ex.write(previous.OUT/'reader-failure.json',dict(classification='Counterexample candidate',passed=False,
        reason='Each request specified 1000 range-generated edges, but each returned Photon grid has 998 lower edges and 998 group rows, not 999. The last printed lower edge matches requested index 997. The requested egphigh is only echoed, not proof of the actual last bin upper edge. Do not silently stretch a bin or invent the missing interval.',
        original_data_preserved=True,scientific_gate_evaluated=False))
    plan=json.loads((previous.OUT/'plan.json').read_text());paths=[Path(__file__),Path(previous.__file__),previous.OUT/'plan.json',previous.OUT/'reader-failure.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    plan['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths})
    plan['reader_repair']='Submit the same 1000 boundaries explicitly instead of invoking the provider logarithmic range generator. The native group count and every printed lower edge must match. Preserve the failed range records. No opacity, physical state, source gate or target interval changes.'
    plan['budget']=dict(queries=16,pilot_queries=1,remaining_queries=15,forecast_network_seconds=36,hard_network_seconds=60,per_query_hard_seconds=20,CPU_workers=1,new_native_EOS_calls=0,new_stellar_steps=0,automatic_retry=False)
    ex.write(OUT/'plan.json',plan)


def fetch(pilot=False):
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    if not pilot:assert json.loads((OUT/'pilot.json').read_text())['passed']
    requests=json.loads((OUT/'requests.json').read_text());rows=requests[:1] if pilot else requests[1:]
    env=dict(vars(native.retrieval),OUT=OUT);env['save']=lambda name,v:ex.write(OUT/name,v)
    function=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__)
    start=time.monotonic();results=[]
    for row in rows:
        assert not list(OUT.glob(row['name']+'-*'));began=time.monotonic();signal.alarm(20)
        try:result=function(row)
        finally:signal.alarm(0)
        result['seconds']=time.monotonic()-began
        lines=[s.strip() for s in (OUT/(row['name']+'-table.txt')).read_text().splitlines() if s.strip()]
        count=int(next(s for s in lines if s.startswith('Photon grid')).split()[-2]);assert count==999
        j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
        tokens=[s.split() for s in lines[j+1:j+1000]];assert np.array(tokens).shape==(999,3)
        submitted=np.array(row['fields']['energies'].split(),float)
        assert max(native.audit.score(native.audit.F(float(x)),t[0]) for x,t in zip(submitted[:-1],tokens))<=1
        result.update(passed=not result['warnings_present'],returned_group_count=count);results.append(result)
        ex.write(OUT/'progress.json',results);print('QUERY',row['name'],result['seconds'],flush=True)
        assert time.monotonic()-start<60
    ex.write(OUT/('pilot.json' if pilot else 'retrieval.json'),dict(passed=all(r['passed'] for r in results),seconds=time.monotonic()-start,queries=len(results)))


def analyze():
    env=dict(vars(previous),OUT=OUT)
    FunctionType(previous.analyze.__code__,env,argdefs=previous.analyze.__defaults__)()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','fetch','analyze']);a=p.parse_args().action
    fetch(True) if a=='pilot' else globals()[a]()
