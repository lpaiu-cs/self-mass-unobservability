"""Keep field subspaces exact, then balance reduced coordinates by input norm.

No physical source rescaling or one-way approximation: both cross blocks stay
in the reduced matrix. A similarity transform balances arithmetic before a
general eigensolve; original physical coordinates are restored for readout.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import time
import json
import signal
import resource
import numpy as np
from scipy.linalg import eig,solve
import def_gr_inverse_projection as common

task=common.task;space=common.space;modes=common.modes
OUT=common.OUT/'fields';write=common.write


def basis(model,size):
    fluid=model.indices[:model.surface_index+1,0];scalar=model.indices[:-1,1]
    assert len(fluid)+len(scalar)==model.size and not np.intersect1d(fluid,scalar).size
    Q=np.zeros((model.size,size));scales=np.empty(size);checks=[]
    for j,indices in enumerate([fluid,scalar]):
        part=SimpleNamespace(K=model.K[indices][:,indices],M=model.M[indices][:,indices],size=len(indices),
            heat=model.heat,load=model.load[indices],original=model.original)
        P=task.Projection(part);_,local,orth=P.basis(size//2)
        Q[np.ix_(indices,np.arange(j,size,2))]=local
        # Norm of the complete physical force, mass-whitened in this field.
        from scipy.sparse.linalg import spsolve_triangular
        scale=np.linalg.norm(spsolve_triangular(P.C,P.forcing,lower=True));assert scale>0
        scales[j::2]=scale;checks.append(dict(field=j,force_norm=float(scale),orthogonality=orth,linear_residual=P.error))
    K,skew=model.energy_matrix(Q);load=Q.T@model.load
    return K,load,Q,scales,dict(fields=checks,stiffness_skew=skew)


def diagonalize(K,load,Q,scales):
    balanced=K*scales[None,:]/scales[:,None];forcing=load/scales[:,None]
    k,V=eig(balanced);assert np.all(k.imag==0),float(np.max(abs(k.imag)))
    k=k.real;assert k.min()>-1.
    residual=float(np.max(abs(balanced@V-V*k))/max(np.max(abs(balanced)@abs(V)),1e-100));assert residual<1e-10
    faces=solve(V,forcing);physical=Q@(scales[:,None]*V)
    assert np.max(abs(faces.imag))==0 and np.max(abs(physical.imag))==0
    return k,faces.real,physical.real,dict(balanced_eigen_residual=residual,field_coordinate_scale_ratio=float(scales[::2][0]/scales[1::2][0]))


def response(model,K,load,packed):
    Q,scales=packed;k,faces,physical,meta=diagonalize(K,load,Q,scales)
    heat=model.heat;tc=model.original.radiation.geometry.tc;lam=heat.rates*tc;coefficients=[];bound=0.
    for i,ki in enumerate(k):
        amp=faces[i,heat.face_ids,None]*heat.amplitude*tc;denom=lam*lam+ki;assert np.all(denom>0)
        coefficients.append([np.sum(amp),np.sum(amp/lam),np.sum(amp/denom),np.sum(amp/(lam*denom))])
        bound=max(bound,float(np.sum(abs(amp)/denom)))
    modal=np.zeros((65,len(k)));velocity=modal.copy()
    for j,t in enumerate(modes.TIMES[1:],1):
        for i,(ki,(a,b,c,d)) in enumerate(zip(k,coefficients)):
            S,C,F,cosine=modes.functions(ki,t)
            modal[j,i]=a*F-b*C+c*S-d*cosine;velocity[j,i]=a*C-b*S+c*cosine+d*ki*S
    assert lam.min()*modes.TIMES[1]>700
    return physical@modal.T,physical@velocity.T,dict(minimum_stiffness=float(k.min()),maximum_stiffness=float(k.max()),
        omitted_exponential_modal_velocity_bound=bound*np.exp(-700),**meta)


def series(model,K,load,Q,scales,label):
    # Reuse the same physical native readout and history writer. Packed Q has a
    # shape solely to preserve the existing dimension metadata at that boundary.
    class Packed(tuple):
        @property
        def shape(self):return self[0].shape
    common.OUT=OUT;common.response=response
    return common.series(model,K,load,Packed((Q,scales)),label)


def control():
    import sympy as s
    x,y,a,b,c,d=s.symbols('x y a b c d',nonzero=True);D=s.diag(x,y);K=s.Matrix([[a,b],[b,c]]);v=s.Matrix([c,d])
    assert s.simplify(D*(D.inv()*K*D)*D.inv()-K)==s.zeros(2)
    assert s.simplify(D*(D.inv()*v)-v)==s.zeros(2,1)
    matrix=np.array([[9.,1e-11],[1e-11,25.]]);force=np.array([[1.],[2e-12]])
    k,f,Q,meta=diagonalize(matrix,force,np.eye(2),np.array([1.,2e-12]));errors=[]
    for z in [2+3j,6+100j]:
        actual=Q@(f[:,0]/(k+z*z));expected=np.linalg.solve(matrix+z*z*np.eye(2),force[:,0])
        errors.append(float(np.max(abs(actual-expected)/abs(expected))))
    assert max(errors)<1e-12
    return dict(classification='Proven',similarity_and_source_identity_checked=True,
        numerical_classification='Counterexample candidate',weak_component_transfer_errors=errors,
        scope='Exact full similarity and a weak-field control; both physical coupling blocks remain. Actual GR convergence is separate.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Preserve the weak scalar component with exactly field-separated mass-orthogonal trial bases and input-norm-balanced arithmetic, then evolve the unchanged coupled system.',
        method='Total128/256/512 columns, equally split fluid/scalar; independent shifted Krylov bases. Retain full cross-field weak stiffness and source. Similarity scales from original mass-whitened force norms, no fit to response. General eigensolve and analytic original heat poles.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Same original degrees4,2,1 and dimensions; degree4 first. Stop on any component failure. Only after propagation and spatial gates run coefficient, outer and quadrature contrasts. No larger space or fitted source.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=120,total_cap_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(common.__file__),Path(space.__file__),Path(task.__file__),common.OUT/'stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    start=time.monotonic();K,F,Q,S,meta=basis(model,64);build=time.monotonic()-start
    row=series(model,K,F,Q,S,'pilot-64')
    forecast=1.4*(setup+64*build+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,build64_seconds=build,response64_seconds=row['seconds'],
        first_case_forecast_seconds=forecast,assumption='Quadratic build and linear response from64 to512 with40 percent margin; first-case cap120s, remaining cases conditional.',**meta))
    signal.alarm(0);print('FIELD FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');K,F,Q,S,meta=basis(model,512)
    rows=[series(model,K[:n,:n],F[:n],Q[:,:n],S[:n],f'p4-{n}') for n in [128,256,512]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('FIELD RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
