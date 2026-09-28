"""Include composition in whole-measure mixing on the fixed background.

The full mixture still enters each direct collision solve. Z_eff below is
only an interpolation coordinate, not an averaged-ion collision model.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import def_conduction_state_bank as state

model=state.model;ex=model.ex;h=model.h
OUT=state.OUT.parent/'def-conduction-mixture-bank'
ANCHORS=state.ANCHORS;WITHHELD=np.array([2588,2625,4050])


def bank(points):
    weights,rates,A=state.bank(points);K=weights.sum(1)
    rows=[dict(np.load(state.previous.source(i))) for i in ANCHORS]
    eta=np.array([s['eta'] for s in rows]);_,data=model.base.inputs()
    active=np.arange(ANCHORS[0],ANCHORS[-1]+1)
    eta_all=ex.model.distribution(data['T'][active],data['ne'][active],256)[0]
    z=np.interp(points,active,eta_all)
    Z=model.base.thermal.g.c.Z;effective=(data['ion']@(Z*Z))/(data['ion']@Z)
    logZ=np.log(effective);localZ=np.interp(points,np.arange(len(logZ)),logZ)
    lo=np.clip(np.searchsorted(eta,z)-1,0,len(eta)-2);hi=lo+1
    de=eta[hi]-eta[lo];dz=logZ[ANCHORS[hi]]-logZ[ANCHORS[lo]]
    # Nearest point on a segment in (eta, log Z_eff), not extrapolation.
    # ponytail: fixed-background two-coordinate surrogate; no general EOS table.
    f=np.clip(((z-eta[lo])*de+(localZ-logZ[ANCHORS[lo]])*dz)/(de*de+dz*dz),0,1)
    p=np.array([s['fine_weights']/s['fine_poles'] for s in rows]);p/=p.sum(1)[:,None]
    weights=K[:,None]*np.c_[(1-f)[:,None]*p[lo],f[:,None]*p[hi]]
    return weights,rates,A


def prepare():
    assert not OUT.exists();OUT.mkdir()
    development=np.r_[state.previous.WITHHELD,state.WITHHELD]
    assert not set(WITHHELD)&set(np.r_[ANCHORS,development])
    paths=[Path(__file__),Path(state.__file__),state.OUT/'plan.json',state.OUT/'result.json']
    paths += [state.previous.source(i) for i in ANCHORS]
    forecast=json.loads((state.OUT/'direct/result.json').read_text())['seconds']/4*3+5
    plan=dict(classification='Counterexample candidate',checkpoint='98f78085',
        claim='Include actual mixture variation in the positive response measure without adding anchors or altering the full-mixture microscopic collision equations.',
        reassessment='The eta-only bank failed tau at cell 2600 by 1.1804 percent. Between 2500 and 2650 Z_eff changes 1.5782 to 2.0444, and 2600 has Z_eff=1.8731. A single eta coordinate misses independent composition variation. This repair changes the interpolation rule, not the collision integral or acceptance thresholds.',
        method='Keep the eta PCHIP DC ratio. For whole-measure convex mixing, project the actual (eta, log Z_eff) point onto the adjacent eta-anchor segment in the Euclidean metric, clipping the mixing fraction to [0,1]. Z_eff=sum(nZ^2)/sum(nZ) is used only as a fractional elastic-strength coordinate; no fitted metric weights. The actual 28 anchor collision operators retain every charge group.',
        anchors=ANCHORS.tolist(),withheld=WITHHELD.tolist(),development=development.tolist(),
        gates=dict(K_relative=.01,tau_relative=.01,normalized_complex_response_absolute=.01,Drude_relative=1e-8,anchor_reproduction=1e-10),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(fresh_states=3,events_per_state=4*2**17,forecast_seconds=forecast,hard_seconds=40,CPU_workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='Finite fixed-background surrogate tests only; no uniform bound, no general composition equation of state, no physical outer current or new insulating boundary.')
    ex.write(OUT/'plan.json',plan);ex.write(OUT/'symbolic.json',state.symbolic())
    direct=OUT/'direct';direct.mkdir();parent=json.loads((state.OUT/'direct/plan.json').read_text())
    parent['cells']=WITHHELD.tolist();parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'plan.json']})
    ex.write(direct/'plan.json',parent);ex.write(direct/'pilot.json',dict(forecast_seconds=forecast))
    print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();oldout,oldcells,oldevents=model.OUT,model.CELLS,ex.events
    model.OUT=OUT/'direct';model.CELLS=WITHHELD.tolist();ex.events=state.old.warm.events
    try:model.run()
    finally:model.OUT,model.CELLS,ex.events=oldout,oldcells,oldevents
    checks=[]
    cohorts=[('independent',WITHHELD),('development',np.array(plan['development']))]
    for cohort,indices in cohorts:
        w,r,_=bank(indices)
        for j,i in enumerate(indices):
            if cohort=='independent':path=OUT/'direct'/f'cell-{i}.npz'
            elif i in state.WITHHELD:path=state.OUT/'direct'/f'cell-{i}.npz'
            else:path=state.previous.source(i)
            s=np.load(path);K=w[j].sum();tau=np.sum(w[j]/r[j])/K;omega=np.array([.01,1,100])/s['fine_tau']
            response=np.sum(w[j,None,:]/(1+1j*omega[:,None]/r[j]),axis=1)/s['fine_K']
            checks.append(dict(index=int(i),cohort=cohort,K_relative=float(abs(K/s['fine_K']-1)),
                tau_relative=float(abs(tau/s['fine_tau']-1)),normalized_complex_response_absolute=float(np.max(abs(response-s['fine_frequency_response'])))))
    w,r,_=bank(ANCHORS);anchor_error=max(max(abs(w[j].sum()/np.load(state.previous.source(i))['fine_K']-1),
        abs(np.sum(w[j]/r[j])/w[j].sum()/np.load(state.previous.source(i))['fine_tau']-1)) for j,i in enumerate(ANCHORS))
    faces=np.arange(ANCHORS[0]+1,ANCHORS[-1]+1);w,r,A=bank(faces-.5)
    drude_error=float(np.max(abs(np.sum(w*r,axis=1)/A-1)))
    passed=json.loads((OUT/'direct/result.json').read_text())['passed'] and anchor_error<1e-10 and drude_error<1e-8 and all(
        all(row[k]<.01 for k in ['K_relative','tau_relative','normalized_complex_response_absolute']) for row in checks)
    np.savez_compressed(OUT/'partial-bank.npz',known_faces=faces,unknown_faces=np.arange(ANCHORS[0]+1),
        mode_K_SI=w,poles_proper_s=r,Drude_weight=A,passed=passed)
    result=dict(classification='Counterexample candidate',passed=passed,checks=checks,
        anchor_reproduction_relative=float(anchor_error),Drude_relative=drude_error,
        covered_internal_faces=len(faces),seconds=time.monotonic()-start,
        uniform_interpolation_certified=False,physical_outer_flux_specified=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);assert result['seconds']<40


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
