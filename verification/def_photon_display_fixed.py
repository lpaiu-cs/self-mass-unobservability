"""Recover all computed groups by widening only the results display window."""
from pathlib import Path
from types import FunctionType
import argparse
import inspect
import json
import signal
import time
import numpy as np
import def_photon_line_resolved as previous
import def_photon_display_control as control

native=previous.native;ex=previous.ex;h=previous.h
OUT=previous.OUT.parent/'def-photon-display-fixed'


def prepare():
    assert not OUT.exists();OUT.mkdir();check=json.loads((control.OUT/'result.json').read_text())['checks']
    assert [r['groups'] for r in check]==[1,2] and check[0]['rows'][0]==check[1]['rows'][0]
    for name in ['requests.json','lanl-tops-form.html']:(OUT/name).write_bytes((previous.OUT/name).read_bytes())
    source=inspect.getsource(previous.analyze);before="'-results-request.json'";assert source.count(before)==1
    source=source.replace(before,"'-generated-results-request.json'");(OUT/'corrected-analyze.py').write_text(source)
    plan=json.loads((previous.OUT/'plan.json').read_text())
    paths=[Path(__file__),Path(previous.__file__),control.OUT/'plan.json',control.OUT/'result.json',OUT/'requests.json',OUT/'lanl-tops-form.html',OUT/'corrected-analyze.py']
    plan['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths})
    plan['display_repair']='A same-submit/same-cookie intervention returned one group with the original numerical-output energy limits and two with only those output limits widened. The common row was identical. The range-generation-loss hypothesis is withdrawn. Resubmit the unchanged 16 physical requests and widen only their output window; require all 999 groups and every old 998-row prefix to match bitwise as parsed tokens. The generated unmodified results form remains the requested endpoint echo.'
    plan['budget']=dict(queries=16,forecast_network_seconds=40,hard_network_seconds=60,per_query_hard_seconds=20,CPU_workers=1,new_native_EOS_calls=0,new_stellar_steps=0,automatic_retry=False)
    ex.write(OUT/'plan.json',plan)


def fetch():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();results=[]
    for row in json.loads((OUT/'requests.json').read_text()):
        assert not list(OUT.glob(row['name']+'-*'))
        def display(form):
            original=native.retrieval.fields(form);ex.write(OUT/(row['name']+'-generated-results-request.json'),original)
            return dict(original,egplow=f"{float(original['egplow'])*(1-1e-4):.17g}",egphigh=f"{float(original['egphigh'])*(1+1e-4):.17g}")
        env=dict(vars(native.retrieval),OUT=OUT,fields=display);env['save']=lambda name,v:ex.write(OUT/name,v)
        function=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__)
        began=time.monotonic();signal.alarm(20)
        try:result=function(row)
        finally:signal.alarm(0)
        def group_tokens(path):
            lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
            n=int(next(s for s in lines if s.startswith('Photon grid')).split()[-2])
            j=next(j for j,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
            return [s.split() for s in lines[j+1:j+1+n]]
        before=group_tokens(previous.OUT/(row['name']+'-table.txt'));after=group_tokens(OUT/(row['name']+'-table.txt'))
        assert len(before)==998 and len(after)==999 and before==after[:998]
        result.update(seconds=time.monotonic()-began,all_previous_groups_unchanged=True,returned_groups=len(after))
        results.append(result);ex.write(OUT/'progress.json',results);print('QUERY',row['name'],result['seconds'],flush=True)
        assert time.monotonic()-start<60
    ex.write(OUT/'retrieval.json',dict(passed=all(not r['warnings_present'] for r in results),seconds=time.monotonic()-start,queries=len(results)))


def analyze():
    env=dict(vars(previous),OUT=OUT)
    exec(compile((OUT/'corrected-analyze.py').read_text(),str(OUT/'corrected-analyze.py'),'exec'),env)
    env['analyze']()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','fetch','analyze']);globals()[p.parse_args().action]()
