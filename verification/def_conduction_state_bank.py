"""Background-specific interpolation in degeneracy, mixing whole measures.

No eigenmode identity is assumed between independently sampled operators.
Eta is an interpolation coordinate on this fixed background, not a claim
that general mixture transport is a function of eta alone.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import def_conduction_shape_bank as previous

model=previous.model;ex=model.ex;h=model.h;old=previous.old
OUT=previous.OUT.parent/'def-conduction-state-bank'
ANCHORS=previous.ANCHORS
WITHHELD=np.array([1980,2600,3031,4210])


def bank(points):
    points=np.asarray(points);assert np.all((points>=ANCHORS[0])&(points<=ANCHORS[-1]))
    states=[dict(np.load(previous.source(i))) for i in ANCHORS]
    coordinate=np.array([s['eta'] for s in states]);assert np.all(np.diff(coordinate)>0)
    data,physical=model.base.inputs();active=np.arange(ANCHORS[0],ANCHORS[-1]+1)
    eta=ex.model.distribution(physical['T'][active],physical['ne'][active],256)[0]
    z=np.interp(points,active,eta)
    assert np.all((z>=coordinate[0])&(z<=coordinate[-1])), 'No eta extrapolation'
    lo=np.clip(np.searchsorted(coordinate,z)-1,0,len(coordinate)-2);hi=lo+1
    f=(z-coordinate[lo])/(coordinate[hi]-coordinate[lo]);assert np.all((f>=0)&(f<=1))
    K=np.array([s['fine_K'] for s in states]);native=np.array([s['native_K'] for s in states])
    native_all=data['K'][:,1]*1e-5
    Kp=np.exp(PchipInterpolator(coordinate,np.log(K/native),extrapolate=False)(z)
        +np.interp(points,np.arange(len(native_all)),np.log(native_all)))
    p=np.array([s['fine_weights']/s['fine_poles'] for s in states]);p/=p.sum(1)[:,None]
    A=old.drude(ANCHORS);r=np.array([s['fine_poles'] for s in states])*K[:,None]/A[:,None]
    assert np.max(abs(np.sum(p*r,axis=1)-1))<1e-8
    # Whole-measure convex mixing is invariant under any per-state pole permutation.
    pp=np.c_[(1-f)[:,None]*p[lo],f[:,None]*p[hi]];rr=np.c_[r[lo],r[hi]]
    Ap=np.exp(np.interp(points,active,np.log(old.drude(active))))
    weights=Kp[:,None]*pp;rates=rr*Ap[:,None]/Kp[:,None]
    assert np.min(weights)>=0 and np.min(rates)>0
    return weights,rates,Ap


def symbolic():
    import sympy as sp
    f,K,A,p,q,r,s=sp.symbols('f K A p q r s',positive=True)
    # One arbitrary component from each endpoint; sum endpoint p*r=1.
    assert sp.simplify(K*(1-f)*p*(r*A/K)+K*f*q*(s*A/K)-A*((1-f)*p*r+f*q*s))==0
    return dict(classification='Proven',passed=True,
        DC='sum p=1 gives sum K*p=K.',
        Drude='If both endpoint measures obey sum p*r=1, their convex mixture also does; lambda=r*A/K gives sum K*p*lambda=A.',
        passivity='For K_j>=0 and lambda_j>0, Re(sum K_j/(1+s/lambda_j))>=0 for Re(s)>=0. Finite first moment is sum K_j/lambda_j / sum K_j.',
        permutation='Sums over whole measures are unchanged by independently permuting endpoint eigenmodes. No mode tracking is required.',
        limitation='These identities do not bound interpolation error or identify the physical continuous collision operator.')


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not set(WITHHELD)&set(np.r_[ANCHORS,previous.WITHHELD])
    paths=[Path(__file__),Path(previous.__file__),previous.OUT/'plan.json',previous.OUT/'result.json']
    paths += [previous.source(i) for i in ANCHORS]
    forecast=json.loads((previous.OUT/'direct/result.json').read_text())['seconds']/22*4+5
    plan=dict(classification='Counterexample candidate',checkpoint='98f78085',
        claim='Repair radial coordinate and eigenmode-matching artifacts without adding anchors, events or basis functions. Validate four new states not used to develop the method.',
        reassessment='The failed 28-anchor bank missed near-constant eta across the temperature turnover and mixed unequal spectral labels near the core. Eta-coordinate whole-measure interpolation passed all eight already inspected development comparisons; those are not independent validation.',
        anchors=ANCHORS.tolist(),withheld=WITHHELD.tolist(),development=previous.WITHHELD.tolist(),
        selection='1980 samples a small eta reversal; 2600 the thermodynamic turnover; 3031 the degeneracy transition; 4210 the changed eigenmode weights.',
        method='PCHIP log computed/native K ratio versus actual eta. Convex mix entire positive DC measures of r=lambda*K/A between neighboring eta anchors, then restore local exact susceptibility A. No pointwise pole matching, no clamping or extrapolation.',
        gates=dict(K_relative=.01,tau_relative=.01,normalized_complex_response_absolute=.01,Drude_relative=1e-8,anchor_reproduction=1e-10),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(fresh_states=4,events_per_state=4*2**17,forecast_seconds=forecast,hard_seconds=45,CPU_workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False),
        limits='Eta is not globally monotone on native cells and is not a sufficient physical state variable. This is a fixed-background surrogate with finite withheld tests; no uniform or physical certificate. Unknown outer faces remain explicit.')
    ex.write(OUT/'plan.json',plan);ex.write(OUT/'symbolic.json',symbolic())
    direct=OUT/'direct';direct.mkdir();parent=json.loads((previous.OUT/'direct/plan.json').read_text())
    parent['cells']=WITHHELD.tolist();parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'plan.json']})
    ex.write(direct/'plan.json',parent);ex.write(direct/'pilot.json',dict(forecast_seconds=forecast))
    print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();oldout,oldcells,oldevents=model.OUT,model.CELLS,ex.events
    model.OUT=OUT/'direct';model.CELLS=WITHHELD.tolist();ex.events=old.warm.events
    try:model.run()
    finally:model.OUT,model.CELLS,ex.events=oldout,oldcells,oldevents
    rows=[]
    for cohort,indices in [('independent',WITHHELD),('development',previous.WITHHELD)]:
        w,r,_=bank(indices)
        for j,i in enumerate(indices):
            path=OUT/'direct'/f'cell-{i}.npz' if cohort=='independent' else previous.source(i)
            s=np.load(path);K=w[j].sum();tau=np.sum(w[j]/r[j])/K;omega=np.array([.01,1,100])/s['fine_tau']
            response=np.sum(w[j,None,:]/(1+1j*omega[:,None]/r[j]),axis=1)/s['fine_K']
            rows.append(dict(index=int(i),cohort=cohort,K_relative=float(abs(K/s['fine_K']-1)),
                tau_relative=float(abs(tau/s['fine_tau']-1)),normalized_complex_response_absolute=float(np.max(abs(response-s['fine_frequency_response'])))))
    w,r,_=bank(ANCHORS);anchor_error=max(max(abs(w[j].sum()/np.load(previous.source(i))['fine_K']-1),
        abs(np.sum(w[j]/r[j])/w[j].sum()/np.load(previous.source(i))['fine_tau']-1)) for j,i in enumerate(ANCHORS))
    faces=np.arange(ANCHORS[0]+1,ANCHORS[-1]+1);w,r,A=bank(faces-.5)
    drude_error=float(np.max(abs(np.sum(w*r,axis=1)/A-1)))
    passed=json.loads((OUT/'direct/result.json').read_text())['passed'] and anchor_error<1e-10 and drude_error<1e-8 and all(
        all(row[k]<.01 for k in ['K_relative','tau_relative','normalized_complex_response_absolute']) for row in rows)
    np.savez_compressed(OUT/'partial-bank.npz',known_faces=faces,unknown_faces=np.arange(ANCHORS[0]+1),
        mode_K_SI=w,poles_proper_s=r,Drude_weight=A,passed=passed)
    result=dict(classification='Counterexample candidate',passed=passed,checks=rows,
        anchor_reproduction_relative=float(anchor_error),Drude_relative=drude_error,
        covered_internal_faces=len(faces),seconds=time.monotonic()-start,
        uniform_interpolation_certified=False,physical_outer_flux_specified=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);assert result['seconds']<45


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
