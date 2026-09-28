"""Positive spectral interpolation with a frozen, bounded radial sampling plan.

Old failed holdouts are development data; fresh holdouts never enter the bank.
The finite-state comparisons are not a uniform or physical error certificate.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import def_conduction_envelope_bank as old

model=old.model;ex=model.ex;h=model.h
OUT=old.OUT.parent/'def-conduction-shape-bank'
REUSED=np.sort(np.r_[old.ANCHORS,old.WITHHELD])
EXTRA=np.concatenate([np.linspace(a,b,int(np.ceil((b-a)/200))+1,dtype=int)[1:-1]
    for a,b in zip(REUSED[:-1],REUSED[1:])])
ANCHORS=np.sort(np.r_[REUSED,EXTRA])
WITHHELD=np.array([1786,2249,2575,2900,3062,3312,4189,5491])


def source(i):
    if i in old.ANCHORS:
        canonical=json.loads((old.SOURCE/'canonical-result.json').read_text())
        return h.ROOT/next(r['state'] for r in canonical['checks'] if r['index']==i)
    return (old.OUT/'withheld' if i in old.WITHHELD else OUT/'direct')/f'cell-{i}.npz'


def bank(points):
    points=np.asarray(points);assert np.all((points>=ANCHORS[0])&(points<=ANCHORS[-1]))
    states=[dict(np.load(source(i))) for i in ANCHORS]
    K=np.array([s['fine_K'] for s in states]);native=np.array([s['native_K'] for s in states])
    rates=np.array([s['fine_poles'] for s in states])
    p=np.array([s['fine_weights']/s['fine_poles'] for s in states]);p/=p.sum(1)[:,None]
    A=old.drude(ANCHORS);r=rates*K[:,None]/A[:,None]
    assert np.max(abs(np.sum(p*r,axis=1)-1))<1e-8
    smooth=lambda v:np.exp(PchipInterpolator(ANCHORS,np.log(v),axis=0,extrapolate=False)(points))
    data,_=model.base.inputs();native_all=data['K'][:,1]*1e-5
    Kp=smooth(K/native)*np.exp(np.interp(points,np.arange(len(native_all)),np.log(native_all)))
    pp=smooth(p);pp/=pp.sum(1)[:,None]
    rr=smooth(r);rr/=np.sum(pp*rr,axis=1)[:,None]
    active=np.arange(ANCHORS[0],ANCHORS[-1]+1)
    Ap=np.exp(np.interp(points,active,np.log(old.drude(active))))
    mode=Kp[:,None]*pp;pole=rr*Ap[:,None]/Kp[:,None]
    assert np.min(mode)>0 and np.min(pole)>0
    return mode,pole,Ap


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not set(ANCHORS)&set(WITHHELD)
    paths=[Path(__file__),Path(old.__file__),Path(model.__file__),Path(old.warm.__file__),
        old.OUT/'plan.json',old.OUT/'result.json',old.OUT/'withheld/result.json',
        old.SOURCE/'canonical-result.json',model.base.thermal.OUT/'coefficients.npz']
    paths += [source(i) for i in REUSED]
    forecast=22*json.loads((old.OUT/'withheld/result.json').read_text())['seconds']/4+10
    plan=dict(classification='Counterexample candidate',checkpoint='98f78085',
        claim='Replace coarse index-linear pole mixing by positive shape-preserving spectral interpolation, then export a usable partial-domain bank after independent state and frequency checks.',
        reassessment='Ten anchors miss radial curvature. Saved four failed holdouts are development data, not fresh validation. Ordinary cubic or matrix interpolation on the same ten anchors did not remove the discrepancy. A fixed maximum 200-cell anchor gap resolves the known sparse sampling; no adaptive refinement or gate relaxation is allowed in this run.',
        method='PCHIP of log(K/native_K), log positive DC probabilities and log dimensionless ordered poles r=lambda*K/A. Normalize probabilities and sum(p*r)=1, restoring exact collision-independent A and interpolated K. Native K supplies a coordinate trend only. No collision-rate calibration.',
        reused=REUSED.tolist(),new_anchors=EXTRA.tolist(),anchors=ANCHORS.tolist(),withheld=WITHHELD.tolist(),
        fresh_validation='Eight fixed states across classical, transition and degenerate layers. Neither these states nor their responses enter the bank.',
        gates=dict(K_relative=.01,tau_relative=.01,normalized_complex_response_absolute=.01,Drude_relative=1e-8,anchor_reproduction=1e-10),
        frequencies_times_reference_tau=[.01,1.,100.],
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(fresh_states=22,events_per_state=4*2**17,forecast_seconds=forecast,hard_seconds=180,CPU_workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='No uniform interpolation bound, no continuum collision certificate, no physical outer flux or new insulating edge. Known internal faces only; full coupled evolution remains unclosed.')
    assert forecast<150;ex.write(OUT/'plan.json',plan)
    direct=OUT/'direct';direct.mkdir();parent=json.loads((old.OUT/'withheld/plan.json').read_text())
    parent['cells']=np.r_[EXTRA,WITHHELD].tolist()
    parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'plan.json']})
    ex.write(direct/'plan.json',parent);ex.write(direct/'pilot.json',dict(forecast_seconds=forecast))
    print('PLAN',len(EXTRA),'new anchors',len(WITHHELD),'independent states; forecast',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();oldout,oldcells,oldevents=model.OUT,model.CELLS,ex.events
    model.OUT=OUT/'direct';model.CELLS=np.r_[EXTRA,WITHHELD].tolist();ex.events=old.warm.events
    try:model.run()
    finally:model.OUT,model.CELLS,ex.events=oldout,oldcells,oldevents
    weights,rates,A=bank(WITHHELD);rows=[]
    for j,i in enumerate(WITHHELD):
        direct=np.load(source(i));K=weights[j].sum();tau=np.sum(weights[j]/rates[j])/K
        omega=np.array(plan['frequencies_times_reference_tau'])/direct['fine_tau']
        response=np.sum(weights[j,None,:]/(1+1j*omega[:,None]/rates[j]),axis=1)/direct['fine_K']
        rows.append(dict(index=int(i),K_relative=float(abs(K/direct['fine_K']-1)),
            tau_relative=float(abs(tau/direct['fine_tau']-1)),
            normalized_complex_response_absolute=float(np.max(abs(response-direct['fine_frequency_response'])))))
    aw,ar,_=bank(ANCHORS);anchor_error=max(max(abs(aw[j].sum()/np.load(source(i))['fine_K']-1),
        abs(np.sum(aw[j]/ar[j])/aw[j].sum()/np.load(source(i))['fine_tau']-1)) for j,i in enumerate(ANCHORS))
    faces=np.arange(ANCHORS[0]+1,ANCHORS[-1]+1);weights,rates,A=bank(faces-.5)
    drude_error=float(np.max(abs(np.sum(weights*rates,axis=1)/A-1)))
    direct=json.loads((OUT/'direct/result.json').read_text())
    passed=direct['passed'] and anchor_error<1e-10 and drude_error<1e-8 and all(
        all(row[k]<.01 for k in ['K_relative','tau_relative','normalized_complex_response_absolute']) for row in rows)
    np.savez_compressed(OUT/'partial-bank.npz',known_faces=faces,unknown_faces=np.arange(ANCHORS[0]+1),
        mode_K_SI=weights,poles_proper_s=rates,Drude_weight=A,passed=passed)
    result=dict(classification='Counterexample candidate',passed=passed,withheld=rows,
        anchor_reproduction_relative=float(anchor_error),Drude_relative=drude_error,
        covered_internal_faces=len(faces),new_direct_states=22,seconds=time.monotonic()-start,
        uniform_interpolation_certified=False,physical_outer_flux_specified=False,
        whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)
    assert result['seconds']<180,'bounded run budget'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
