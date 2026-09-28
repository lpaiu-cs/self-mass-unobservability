"""Resolve the composition transition using already-computed local states.

Reuse the permutation-invariant positive-measure evaluator. Promoted failed
holdouts are explicitly anchors; two fresh transition states remain withheld.
"""
from pathlib import Path
from contextlib import contextmanager
import argparse
import json
import numpy as np
import def_conduction_state_bank as state
import def_conduction_mixture_bank as mixture

model=state.model;ex=model.ex;h=model.h
OUT=state.OUT.parent/'def-conduction-resolved-bank'
PROMOTED=[2575,2588,2600,2625]
ANCHORS=np.sort(np.r_[state.ANCHORS,PROMOTED])
WITHHELD=np.array([2612,2637])
DEVELOPMENT=np.array([1786,2249,2900,3062,3312,4189,5491,1980,3031,4210,4050])
PATHS={int(i):state.previous.source(i) for i in np.r_[state.ANCHORS,state.previous.WITHHELD]}
PATHS.update({int(i):state.OUT/'direct'/f'cell-{i}.npz' for i in state.WITHHELD})
PATHS.update({int(i):mixture.OUT/'direct'/f'cell-{i}.npz' for i in mixture.WITHHELD})


@contextmanager
def configuration():
    before=(state.OUT,state.ANCHORS,state.WITHHELD,state.previous.WITHHELD,state.previous.source)
    state.OUT=OUT;state.ANCHORS=ANCHORS;state.WITHHELD=WITHHELD
    state.previous.WITHHELD=DEVELOPMENT;state.previous.source=lambda i:PATHS[int(i)]
    try:yield
    finally:state.OUT,state.ANCHORS,state.WITHHELD,state.previous.WITHHELD,state.previous.source=before


def bank(points):
    with configuration():return state.bank(points)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not set(WITHHELD)&set(PATHS)
    paths=[Path(__file__),Path(state.__file__),state.OUT/'result.json',mixture.OUT/'result.json']+list(PATHS.values())
    forecast=json.loads((mixture.OUT/'direct/result.json').read_text())['seconds']/3*2+5
    plan=dict(classification='Counterexample candidate',checkpoint='98f78085',
        claim='Resolve the narrow changing-composition layer using four existing direct states, with independent midpoint checks and no new anchor calculations.',
        reassessment='A one-coordinate wide segment and then a two-coordinate straight segment each failed the original tau gate in the H/He transition. Four direct states already exist there. Promote only those to anchors instead of fitting a new coordinate metric or calculating the whole radial mesh. The original failed verdicts remain unchanged.',
        anchors=ANCHORS.tolist(),promoted_prior_holdouts=PROMOTED,withheld=WITHHELD.tolist(),development=DEVELOPMENT.tolist(),
        state_paths={str(i):p.relative_to(h.ROOT).as_posix() for i,p in PATHS.items()},
        method='Unchanged eta-coordinate PCHIP DC trend and convex mixing of whole normalized positive measures from def_conduction_state_bank. Local transition anchors now carry the actual mixture dependence. No composition-averaged collision operator and no metric fitting.',
        gates=dict(K_relative=.01,tau_relative=.01,normalized_complex_response_absolute=.01,Drude_relative=1e-8,anchor_reproduction=1e-10),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(fresh_states=2,new_anchors_computed=0,events_per_state=4*2**17,forecast_seconds=forecast,hard_seconds=35,CPU_workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='Only two new independent transition checks; older nonanchor checks are development data. No uniform or physical error certificate. Unknown outer faces remain explicit and the whole-star thermal problem is unclosed.')
    ex.write(OUT/'plan.json',plan);ex.write(OUT/'symbolic.json',state.symbolic())
    direct=OUT/'direct';direct.mkdir();parent=json.loads((mixture.OUT/'direct/plan.json').read_text())
    parent['cells']=WITHHELD.tolist();parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'plan.json']})
    ex.write(direct/'plan.json',parent);ex.write(direct/'pilot.json',dict(forecast_seconds=forecast))
    print('FORECAST',forecast,flush=True)


def run():
    with configuration():state.run()
    assert json.loads((OUT/'result.json').read_text())['seconds']<35


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
