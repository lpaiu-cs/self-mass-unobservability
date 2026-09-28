"""Native isentropic continuation of the saved finite-pressure material surface.

This constructs mechanical Cauchy data. It does not impose radiative/thermal
stationarity or infer an atmosphere from a Rosseland mean alone.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import def_hydrostatic_background as h

OUT=h.molecular.g.OUT/'def-native-atmosphere'
BACKGROUND=h.OUT/'absolute-shoot'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),BACKGROUND/'background-0.001.npz',BACKGROUND/'lapse.npz',
           BACKGROUND/'native-audit.json',h.OLD/'molecular-adiabats-17.npz',
           h.ROOT/'verification/gr_radiative_boundary.py']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='fd402042',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        claim='Continue the actual molecular EOS at the outer shell entropy/composition beyond the nonzero-pressure optical boundary; solve material/scalar/metric atmosphere and identify a controlled free-surface limit.',
        material='Original interior inventories are unchanged. Any added atmosphere inventory is recorded separately and must be retained in future paired backgrounds.',
        thermal='Mechanical Cauchy data only. No radiative equilibrium, physical opacity calibration or orbital stationarity is assumed.',
        pilot_log_pressure_drops=[0,1,2,4,8,12],
        budget=dict(pilot_seconds=60,table_seconds=120,atmosphere_seconds=60,maximum_new_entropy_roots=128,automatic_expansion=False),
        gates=dict(native_root='Same strict energy-unit entropy inverse as Phase47',junction=1e-8,
                   declared_atmosphere_baryon_fraction=1e-12,unresolved_tail_baryon_fraction=1e-14),
        stop='Stop on EOS/root/domain failure. Do not replace a failed molecular adiabat by a fitted polytrope or claim the mass budget bounds all dynamic boundary errors.'))


def pilot():
    assert not (OUT/'pilot.json').exists()
    plan=json.loads((OUT/'plan.json').read_text());data,tab=h.inputs()
    eos=h.molecular.model.EOS();stats=dict(calls=0,evaluations=0,maximum_score=0.,label='atmosphere-pilot')
    inverse=h.molecular.inverse(stats,OUT/'pilot-roots')
    entropy=tab['reference'][0,3];composition=data['X'][0]
    ps=float(data['boundary_logP']);previous=ps;lt=float(tab['lnT'][0]);nabla=.3;rows=[];start=time.monotonic()
    for drop in plan['pilot_log_pressure_drops']:
        lp=ps-drop;guess=max(4.,lt+nabla*(lp-previous))
        a,lt,_=inverse(eos,lp,entropy,composition,guess)
        cpT=a[10]-a[1]/a[0]*a[8];nabla=-a[1]/a[0]*a[8]/cpT
        dlnT_dlnrho=(a[1]/a[0]-a[9])/a[10]
        gamma=a[5]+a[6]*dlnT_dlnrho
        rows.append(dict(drop=drop,logP=lp,lnT=float(lt),raw=a.tolist(),gamma1=float(gamma),nabla=float(nabla)))
        h.write(OUT/'pilot-progress.json',dict(classification='Counterexample candidate',rows=rows,statistics=stats,seconds=time.monotonic()-start))
        previous=lp
    h.write(OUT/'pilot.json',dict(classification='Counterexample candidate',passed=True,rows=rows,statistics=stats,seconds=time.monotonic()-start))
    print(json.dumps(dict(seconds=time.monotonic()-start,rows=[{k:v for k,v in r.items() if k!='raw'} for r in rows],statistics=stats)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot'])
    globals()[p.parse_args().action]()
