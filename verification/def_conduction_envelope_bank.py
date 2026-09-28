"""Positive pole interpolation preserving DC and the thermal Drude weight.

Unknown outer faces remain explicit; this bank cannot authorize a new
insulated cohort boundary in the stellar heat equation.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import def_ionic_structure_transport as model
import def_electron_warm_quadrature as warm

ex=model.ex
h=model.h
SOURCE=model.OUT
OUT=SOURCE.parent/'def-conduction-envelope-bank'
ANCHORS=np.array(model.CELLS)
WITHHELD=[1850,2650,3250,4525]


def drude(indices):
    _,data=model.base.inputs();T=data['T'][indices];ne=data['ne'][indices]
    eta,x,p,fp,dx,_,_=ex.model.distribution(T,ne,256)
    theta=model.base.k*T/(model.base.m_e*model.base.c**2);E=np.sqrt(1+p*p)
    measure=theta[:,None]*p**3/E*fp*dx/(3*np.pi**2)
    z=x-eta[:,None];mean=np.sum(measure*z,axis=1)/measure.sum(1)
    return ex.KUNIT*ex.RATE*np.sum(measure*(z-mean[:,None])**2,axis=1)


def rows():
    canonical=json.loads((SOURCE/'canonical-result.json').read_text())
    assert canonical['passed']
    return [dict(np.load(h.ROOT/r['state'])) for r in canonical['checks']]


def bank(points):
    d,_=model.base.inputs();native=d['K'][:,1]*1e-5;states=rows()
    a=np.clip(np.searchsorted(ANCHORS,points)-1,0,len(ANCHORS)-2);b=a+1
    f=(points-ANCHORS[a])/(ANCHORS[b]-ANCHORS[a]);assert np.all((f>=0)&(f<=1))
    K0=np.array([q['fine_K'] for q in states]);native0=np.array([q['native_K'] for q in states])
    K=np.exp(np.interp(points,np.arange(len(native)),np.log(native))+(1-f)*np.log(K0[a]/native0[a])+f*np.log(K0[b]/native0[b]))
    rates=np.array([q['fine_poles'] for q in states]);p=np.array([q['fine_weights']/q['fine_poles'] for q in states]);p/=p.sum(1)[:,None]
    A=drude(ANCHORS);exactA=K0*np.sum(p*rates,axis=1)
    assert np.max(abs(A/exactA-1))<1e-8,'continuum versus finite-basis susceptibility'
    dimensionless=rates*K0[:,None]/exactA[:,None]
    probability=np.c_[(1-f)[:,None]*p[a],f[:,None]*p[b]]
    dimensionless=np.c_[dimensionless[a],dimensionless[b]]
    # Evaluate the collision-independent susceptibility on the original cells
    # once; its logarithm supplies the face midpoint value.
    active=np.arange(ANCHORS[0],ANCHORS[-1]+1)
    Ap=np.exp(np.interp(points,active,np.log(drude(active))))
    poles=dimensionless*Ap[:,None]/K[:,None];modeK=K[:,None]*probability
    assert np.min(modeK)>=0 and np.min(poles)>0
    return modeK,poles,Ap,float(np.max(abs(A/exactA-1)))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(model.__file__),Path(warm.__file__),SOURCE/'manifest.json',SOURCE/'canonical-result.json',model.base.thermal.OUT/'coefficients.npz']
    plan=dict(classification='Counterexample candidate',checkpoint='98f78085',
        claim='Interpolate the accepted cross-degeneracy collision response onto the actual radial cells while preserving the positive DC measure and its collision-independent high-frequency thermal susceptibility. Test new states withheld from interpolation.',
        interpolation='Use the native conductivity only as a coordinate trend for the computed/native DC ratio, with no fitted normalization. Normalize every pole by Drude_weight/K; mix the two adjacent positive DC probabilities in original cell coordinate, then restore the exact susceptibility and interpolated K at each point.',
        boundary='Cells 0..1721 and the external thermal surface remain unspecified. Export unknown faces explicitly. No zero-flux replacement, no new stellar path.',
        anchors=ANCHORS.tolist(),withheld=WITHHELD,gates=dict(K_relative=.01,tau_relative=.01,Drude_relative=1e-8),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        budget=dict(fresh_states=4,events_per_state=4*2**17,hard_seconds=90,CPU_workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False))
    ex.write(OUT/'plan.json',plan)
    pilot=json.loads((SOURCE/'warm-quadrature/pilot.json').read_text())['seconds']
    hnc=json.loads((SOURCE/'pilot.json').read_text())
    forecast=4*(pilot+hnc['HNC_coarse_seconds']+hnc['HNC_fine_seconds'])+10
    ex.write(OUT/'pilot.json',dict(forecast_seconds=forecast,basis='Four measured warm-state collision solves plus two measured HNC grids each, and 10 seconds for whole-domain susceptibility evaluation and export; input/output overhead excluded.'))
    direct=OUT/'withheld';direct.mkdir()
    parent=json.loads((SOURCE/'plan.json').read_text());parent['cells']=WITHHELD;parent['box_controls']['cells']=[]
    parent['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'plan.json']})
    parent['proposal_override']='Use Phase57 exact defensive energy reweighting on every withheld state; the four unchanged seeds and event/basis gates are retained.'
    ex.write(direct/'plan.json',parent);ex.write(direct/'pilot.json',dict(forecast_seconds=forecast))
    print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    assert json.loads((OUT/'pilot.json').read_text())['forecast_seconds']<80
    start=time.monotonic();oldout,oldcells,oldevents=model.OUT,model.CELLS,ex.events
    model.OUT=OUT/'withheld';model.CELLS=WITHHELD;ex.events=warm.events
    try:model.run()
    finally:model.OUT,model.CELLS,ex.events=oldout,oldcells,oldevents
    w,r,A,finite=bank(np.array(WITHHELD));comparisons=[]
    for j,i in enumerate(WHELD:=WITHHELD):
        direct=np.load(OUT/'withheld'/f'cell-{i}.npz')
        comparisons.append(dict(index=i,K_relative=float(abs(w[j].sum()/direct['fine_K']-1)),
            tau_relative=float(abs(np.sum(w[j]/r[j])/w[j].sum()/direct['fine_tau']-1))))
    faces=np.arange(ANCHORS[0]+1,ANCHORS[-1]+1);w,r,A,finite=bank(faces-.5)
    np.savez_compressed(OUT/'partial-bank.npz',known_faces=faces,unknown_faces=np.arange(0,ANCHORS[0]+1),
        mode_K_SI=w,poles_proper_s=r,Drude_weight=A)
    direct=json.loads((OUT/'withheld/result.json').read_text());passed=direct['passed'] and all(x['K_relative']<.01 and x['tau_relative']<.01 for x in comparisons)
    result=dict(classification='Counterexample candidate',passed=passed,withheld=comparisons,
        susceptibility_continuum_finite_difference=finite,interpolated_Drude_relative=float(np.max(abs(np.sum(w*r,axis=1)/A-1))),
        covered_internal_faces=len(faces),seconds=time.monotonic()-start,
        physical_outer_flux_specified=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
