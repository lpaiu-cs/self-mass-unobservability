"""Locate coupled GR error at full, projected and diagonalized transfers."""
from pathlib import Path
import argparse
import json
import signal
import resource
import time
import numpy as np
from scipy.linalg import eig,solve
from scipy.sparse.linalg import splu
import def_gr_allocated_projection as previous

task=previous.task;space=previous.space;gram=previous.gram;common=previous.common
OUT=task.OUT.parent/'def-gr-transfer-repair';write=task.write;sigma=128.**2


def readout(model,u,z,E):
    q=u-model.H@E;r=model.original.native;N,a=model.original.radiation.geometry.metric(r)
    speed=a/N*r*task.fem.base.task.h.gr.C
    return np.array([speed*(model.nativeV[0]@(z*q)),model.nativeV[1]@q])


def norms(model,x):
    w=model.original.weights;v,f=x
    return np.array([np.sqrt(w@abs(v)**2),np.sqrt(w@abs(f)**2),
        *[np.sqrt(w[m]@abs(v[m])**2/w[m].sum()) for m in model.original.masks]])


def source(model,z):
    heat=model.heat;tc=model.original.radiation.geometry.tc;lam=heat.rates*tc
    E=np.zeros(len(heat.edges),complex)
    E[heat.face_ids]=np.sum(heat.amplitude*tc*lam/(z*z*(z+lam)),axis=1)
    return E


def full(K,M,rhs,z):
    matrix=(K+z*z*M).tocsc();scale=np.asarray(abs(matrix).sum(1)).ravel()
    lu=splu(matrix.multiply((1/scale)[:,None]).tocsc());x=lu.solve(rhs/scale).astype(np.clongdouble)
    extended=matrix.astype(np.clongdouble);b=rhs.astype(np.clongdouble)
    for _ in range(3):x+=lu.solve(np.asarray((b-extended@x)/scale,complex))
    residual=float(np.max(abs(b-extended@x)/(abs(extended)@abs(x)+abs(b)+1e-100)))
    assert residual<1e-9
    return x,residual


def run():
    assert not OUT.exists();OUT.mkdir();signal.alarm(120)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Locate transfer error in the same failed fixed512 coupled projection before another trajectory: full sparse solve versus projected dense solve versus current eigendecomposition.',
        method='Reuse degree4 same source/bank and phase78 fixed7:1 basis. Frequencies4+i*[0,4,16,64,256,1024] are declared before inspection. Symmetric full assembly matches the Gram method. Preserve all four native readouts and both coupling blocks.',
        budget=dict(hard_seconds=120,CPU_threads=1,memory_GB=4,new_time_paths=0,new_EOS_calls=0),
        decision='Only a diagnosis. Preserve failed convergence. Use the component-wise transfer errors to choose a concrete projection or diagonalization repair; do not increase dimensions or run trajectories from a small residual alone.',
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(previous.__file__),Path(gram.__file__),Path(space.__file__),Path(task.__file__),previous.OUT/'stage-result.json',task.BANK/'fine-bank.npz']}))
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');R,F,Q,S,tail,meta=previous.basis(model,512)
    K=(model.K+model.K.T)/2;M=(model.M+model.M.T)/2;load=model.load+(K-model.K)@model.H
    cases=[];reduced={};mass=Q.T@(M@Q)
    for n in [128,256,512]:
        A=R[:n,:n]*S[None,:n]/S[:n,None];theta,V=eig(A);assert np.all(theta.imag==0)
        reduced[n]=(A,theta.real,V.real)
    for frequency in [0.,4.,16.,64.,256.,1024.]:
        z=4+1j*frequency;E=source(model,z);rhs=load@E
        direct,residual=full(K,M,rhs,z);reference=readout(model,direct,z,E);norm=norms(model,reference)
        source_projected=sigma*(tail[0].T@rhs)
        for n in [128,256,512]:
            A,theta,V=reduced[n];f=source_projected[:n]/S[:n]
            operator=z*z*A+sigma*(np.eye(n)-A)
            y=solve(operator,f);ue=Q[:,:n]@(S[:n]*(V@(solve(V,f)/(z*z*theta+sigma*(1-theta)))))
            u=Q[:,:n]@(S[:n]*y)
            # The original transformed Galerkin form has Q.T M Q, not an
            # assumed identity. Measure whether its finite arithmetic matters.
            G=mass[:n,:n]*S[None,:n]/S[:n,None]
            ym=solve(z*z*A+sigma*(G-A),f);um=Q[:,:n]@(S[:n]*ym)
            projected=readout(model,u,z,E);diagonalized=readout(model,ue,z,E);mass_corrected=readout(model,um,z,E)
            cases.append(dict(frequency=frequency,dimension=n,full_residual=residual,
                projection_relative=dict(zip(task.FIELDS,(norms(model,projected-reference)/norm).astype(float))),
                diagonalization_relative=dict(zip(task.FIELDS,(norms(model,diagonalized-projected)/norm).astype(float))),
                actual_mass_relative=dict(zip(task.FIELDS,(norms(model,mass_corrected-projected)/norm).astype(float))),
                reduced_eigenvalue_min=float(theta.min()),reduced_eigenvalue_max=float(theta.max())))
    np.savez_compressed(OUT/'reduced.npz',R=R,F=F,scales=S,mass=mass)
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,cases=cases,**meta,
        original_failure_resolved=False,propagation_passed=False,full_dynamic_charge_solved=False)
    write(OUT/'transfer.json',result);signal.alarm(0)
    for n in [128,256,512]:
        print('TRANSFER',n,json.dumps({method:{key:max(r[method][key] for r in cases if r['dimension']==n) for key in task.FIELDS}
            for method in ['projection_relative','diagonalization_relative','actual_mass_relative']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run']);globals()[p.parse_args().action]()
