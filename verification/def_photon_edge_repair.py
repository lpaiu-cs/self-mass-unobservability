"""Audit the final generated group and retrieve only its missing neighbor."""
from pathlib import Path
from types import FunctionType
import argparse
import json
import signal
import time
import numpy as np
import def_photon_line_resolved as line

native=line.native;ex=line.ex;h=line.h
OUT=line.OUT.parent/'def-photon-edge-repair'


def prepare():
    assert not OUT.exists();OUT.mkdir();rows=[]
    for row in json.loads((line.OUT/'requests.json').read_text()):
        fields=dict(row['fields']);edge=np.geomspace(float(fields['egplow']),float(fields['egphigh']),1000)
        name=row['name'].replace('phLine','phEdge');fields.update(mixname=name,egrid='specific',energies=' '.join(f'{v:.17g}' for v in edge[-3:]))
        rows.append(dict(row,name=name,fields=fields,original_name=row['name']))
    ex.write(OUT/'requests.json',rows);(OUT/'lanl-tops-form.html').write_bytes((native.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(line.__file__),line.OUT/'plan.json',line.OUT/'reader-failure.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Determine the actual last range-generated bin by an overlapping two-bin explicit query and supply its missing last neighbor; reuse all other source data.',
        comparison='First explicit interval [e997,e998] must match the last returned range row in exact printed intervals for both opacity means. Then [e998,e999] is appended as a separately retrieved bin. A mismatch rejects this repair.',
        queries='Three explicitly supplied edges per request, same material and temperature/density. Full 1000-edge explicit pilot returned HTTP 403; no retry of that rejected payload.',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(pilot_queries=1,remaining_queries=15,hard_network_seconds=60,per_query_hard_seconds=20,automatic_retry=False)))


def parsed(name):
    path=OUT/(name+'-table.txt');lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
    count=int(next(s for s in lines if s.startswith('Photon grid')).split()[-2]);assert count==2
    j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
    tokens=[s.split() for s in lines[j+1:j+3]];assert np.array(tokens).shape==(2,3)
    return tokens


def fetch(pilot=False):
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    if not pilot:assert json.loads((OUT/'pilot.json').read_text())['passed']
    rows=json.loads((OUT/'requests.json').read_text());rows=rows[:1] if pilot else rows[1:]
    env=dict(vars(native.retrieval),OUT=OUT);env['save']=lambda name,v:ex.write(OUT/name,v)
    function=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__)
    start=time.monotonic();results=[];header=native.reader.namespace()['header']
    for row in rows:
        assert not list(OUT.glob(row['name']+'-*'));signal.alarm(20)
        try:result=function(row)
        except Exception as exc:
            ex.write(OUT/f"{row['name']}-failure.json",dict(error=repr(exc)));raise
        finally:signal.alarm(0)
        data=header(OUT/(row['name']+'-table.txt'),dict(row,fields=dict(row['fields'],datype='gray')));assert data['input_passed']
        tokens=parsed(row['name']);edge=np.array(row['fields']['energies'].split(),float)
        assert max(native.audit.score(native.audit.F(float(a)),b[0]) for a,b in zip(edge[:-1],tokens))<=1
        oldlast=(line.OUT/(row['original_name']+'-table.txt')).read_text().splitlines()[-1].split()
        scores=[]
        for a,b in zip(oldlast[1:],tokens[0][1:]):
            la,ua=native.audit.interval(a);lb,ub=native.audit.interval(b);scores.append(bool(max(la,lb)<=min(ua,ub)))
        result.update(passed=all(scores),overlap_matches=scores,old_last=oldlast,explicit_first=tokens[0]);results.append(result)
        ex.write(OUT/'progress.json',results);print('EDGE',result,flush=True)
        assert result['passed'],'Last-bin interpretation rejected';assert time.monotonic()-start<60
    ex.write(OUT/('pilot.json' if pilot else 'retrieval.json'),dict(passed=all(r['passed'] for r in results),seconds=time.monotonic()-start,queries=len(results)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','fetch']);a=p.parse_args().action
    fetch(True) if a=='pilot' else globals()[a]()
