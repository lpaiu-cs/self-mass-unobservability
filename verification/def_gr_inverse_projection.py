"""Shift-invert Galerkin propagation of the unchanged coupled GR equation.

The mass-whitened equation w_tt+A*w=b is equivalently
R*w_tt+sigma*(I-R)*w=R*b, R=sigma*(A+sigma*I)^-1.
Project both R and the complete transformed source, never clip modes.
"""
from pathlib import Path
import argparse
import json
import time
import signal
import resource
import numpy as np
from scipy.linalg import eigh
import def_gr_hierarchical as space

task=space.task;modes=task.modes;OUT=task.OUT.parent/'def-gr-inverse-projection'
write=task.write;FIELDS=task.FIELDS


def projected(P,white):
    n=white.shape[1];sigma=P.sigma;T=P.C@white;S=np.empty_like(T)
    for lo in range(0,n,16):
        sl=slice(lo,min(lo+16,n));rhs=T[:,sl];x=P.lu.solve(rhs/P.scale[:,None])/P.scale[:,None]
        longrhs=rhs.astype(np.longdouble)
        for _ in range(2):
            defect=longrhs-P.extended@x.astype(np.longdouble)
            x+=P.lu.solve(np.asarray(defect/P.scale[:,None],float))/P.scale[:,None]
        defect=longrhs-P.extended@x.astype(np.longdouble)
        error=float(np.max(abs(defect)/(abs(P.extended)@abs(x)+abs(rhs)+1e-100)))
        P.error=max(P.error,error);assert P.error<1e-9;S[:,sl]=x
    R=sigma*(T.T@S);skew=float(np.max(abs(R-R.T))/max(np.max(abs(R)),1e-100));assert skew<1e-12
    # Both the left operator and right source receive the same inverse map.
    load=sigma*(S.T@P.model.load)
    return (R+R.T)/2,load,dict(resolvent_skew=skew,linear_residual=P.error)


def response(model,R,load,physical):
    theta,V=eigh(R);assert theta.min()>0
    k=model_sigma*(1-theta)/theta
    assert k.min()>-1.,float(k.min())
    heat=model.heat;tc=model.original.radiation.geometry.tc;lam=heat.rates*tc
    faces=(V.T@load)/theta[:,None];coefficients=[];bound=0.
    for i,ki in enumerate(k):
        amp=faces[i,heat.face_ids,None]*heat.amplitude*tc;denom=lam*lam+ki;assert np.all(denom>0)
        coefficients.append([np.sum(amp),np.sum(amp/lam),np.sum(amp/denom),np.sum(amp/(lam*denom))])
        bound=max(bound,float(np.sum(abs(amp)/denom)))
    modal=np.zeros((65,len(k)));velocity=modal.copy()
    for j,t in enumerate(modes.TIMES[1:],1):
        for i,(ki,(a,b,c,d)) in enumerate(zip(k,coefficients)):
            S,C,F,cosine=modes.functions(ki,t)
            modal[j,i]=a*F-b*C+c*S-d*cosine
            velocity[j,i]=a*C-b*S+c*cosine+d*ki*S
    assert lam.min()*modes.TIMES[1]>700
    Q=physical@V
    return Q@modal.T,Q@velocity.T,dict(minimum_stiffness=float(k.min()),maximum_stiffness=float(k.max()),
        resolvent_eigenvalue_minimum=float(theta.min()),omitted_exponential_modal_velocity_bound=bound*np.exp(-700))


model_sigma=128.**2


def series(model,R,load,Q,label):
    began=time.monotonic();u,v,meta=response(model,R,load,Q);history=[];vs=[];fs=[]
    r=model.original.native;N,a=model.original.radiation.geometry.metric(r)
    speed=a/N*r*task.fem.base.task.h.gr.C;w=model.original.weights
    for j,t in enumerate(modes.TIMES):
        flux,energy=model.heat.faces(t);q=u[:,j]-model.H@energy;qt=v[:,j]-model.H@(flux*model.original.radiation.geometry.tc)
        cv=speed*(model.nativeV[0]@qt);cf=model.nativeV[1]@q
        history.append(dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(w@(cv*cv))),scalar_mass_RMS=float(np.sqrt(w@(cf*cf))),
            **{key:float(np.sqrt(w[mask]@cv[mask]**2/w[mask].sum())) for key,mask in zip(FIELDS[2:],model.original.masks)}))
        vs.append(cv);fs.append(cf)
    np.savez_compressed(OUT/(label+'.npz'),grid=model.grid,cells=model.cells,indices=model.indices,q=q,qt=qt,
        native_radius=r,native_velocity=vs,native_scalar=fs,weights=w,masks=model.original.masks,heat_energy=energy,heat_flux=flux)
    result=dict(classification='Counterexample candidate',degree=model.degree,projection=Q.shape[1],history=history,seconds=time.monotonic()-began,
        heat_telescoping=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100)),**meta)
    write(OUT/(label+'.json'),result);print('INVERSE',label,result['seconds'],history[-1],flush=True);return result


def comparison(rows):
    result={}
    for key in FIELDS:
        a,b,c=[np.array([r[key] for r in row['history']]) for row in rows]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        result[key]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
    return result


def control():
    import sympy as s
    a,b,sigma=s.symbols('a b sigma',positive=True);A=s.Matrix([[a,b],[b,a+2*b]])
    R=sigma*(A+sigma*s.eye(2)).inv()
    assert s.simplify(R*A-sigma*(s.eye(2)-R))==s.zeros(2)
    assert s.simplify(R.inv()*R)==s.eye(2)
    # Full-space transfer must agree independently, including weak coupling.
    K=np.array([[9.,1e-11],[1e-11,25.]]);force=np.array([1.,2e-12]);R=128.**2*np.linalg.inv(K+128.**2*np.eye(2))
    theta,V=eigh(R);lam=128.**2*(1-theta)/theta;rhs=V.T@(R@force)/theta
    errors=[]
    for z in [2+3j,6+100j]:
        direct=np.linalg.solve(K+z*z*np.eye(2),force);project=V@(rhs/(lam+z*z))
        errors.append(float(np.max(abs(project-direct)/np.maximum(abs(direct),1e-100))))
    assert max(errors)<1e-9,errors
    return dict(classification='Proven',equivalence_checked=True,identity='R*A=sigma*(I-R). Transform both forcing and inertia before projection; no change to full equations.',
        numerical_classification='Counterexample candidate',two_field_transfer_errors=errors,
        scope='Full-space identity and weakly coupled control; projected GR histories still require all original gates.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Test whether bounded shift-invert projection removes scalar loss while preserving already resolved rapid velocity in the same degree4 coupled GR space.',
        method='Same mass-orthogonal rational basis and128/256/512 dimensions, but project R inertia and sigma*(I-R) stiffness together with the complete R-transformed heat source. Same exact pole propagation, horizon, spatial basis and readouts.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Pilot64, then degree4 first. Stop on any four-field propagation failure. If accepted, fixed degree2 and1 paths; only after spatial acceptance run coefficient/outer/quadrature contrasts. No larger basis, degree, time interval or input.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=120,total_cap_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(space.__file__),Path(task.__file__),Path(modes.__file__),task.OUT/'result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');P=task.Projection(model)
    setup=time.monotonic()-start;start=time.monotonic();white,Q,orth=P.basis(64);basis=time.monotonic()-start
    start=time.monotonic();R,F,meta=projected(P,white);projection=time.monotonic()-start;case=series(model,R,F,Q,'pilot-64')
    # Prior same512-dimensional basis measured13.96s; inverse projection is new.
    forecast=1.4*(setup+13.96+64*projection+8*case['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,basis64_seconds=basis,
        projection64_seconds=projection,response64_seconds=case['seconds'],first_case_forecast_seconds=forecast,
        orthogonality=orth,**meta,assumption='Reuse measured13.96s degree4/512 basis cost; project conservatively quadratically and response linearly from64 to512,40 percent margin. Subsequent cases unmeasured and remain conditional.'))
    signal.alarm(0);print('INVERSE FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));began=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');P=task.Projection(model);white,Q,orth=P.basis(512)
    R,F,meta=projected(P,white);rows=[series(model,R[:n,:n],F[:n],Q[:,:n],f'p4-{n}') for n in [128,256,512]]
    cmp=comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,orthogonality=orth,seconds=time.monotonic()-began,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('INVERSE RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
