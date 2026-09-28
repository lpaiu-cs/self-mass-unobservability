"""Full coupled response snapshots replace positive-real-only trial spaces."""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import signal
import resource
import time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu,spsolve_triangular
import def_gr_transfer_repair as transfer

previous=transfer.previous;gram=transfer.gram;task=transfer.task;space=transfer.space;common=transfer.common
OUT=task.OUT.parent/'def-gr-response-basis';write=task.write


def frequency(index):
    # Nested, deterministic van der Corput points on a cosine frequency grid.
    x=0.;weight=.5
    while index:x+=(index%2)*weight;index//=2;weight/=2
    return 1024*np.sin(np.pi*x/2)**2


def full(model,z):
    K=(model.K+model.K.T)/2;M=(model.M+model.M.T)/2
    load=model.load+(K-model.K)@model.H;E=transfer.source(model,z);rhs=load@E
    matrix=(K+z*z*M).tocsc();scale=np.sqrt(abs(matrix.diagonal()));D=diags(1/scale)
    lu=splu((D@matrix@D).tocsc());u=(lu.solve(rhs/scale)/scale).astype(np.clongdouble)
    extended=matrix.astype(np.clongdouble);b=rhs.astype(np.clongdouble)
    for _ in range(3):u+=(lu.solve(np.asarray((b-extended@u)/scale,complex))/scale).astype(np.clongdouble)
    residual=float(np.max(abs(b-extended@u)/(abs(extended)@abs(u)+abs(b)+1e-100)));assert residual<1e-9
    return np.asarray(u,complex),residual


def field_basis(model,size):
    assert size%4==0
    groups=[model.indices[:model.surface_index+1,0],model.indices[:-1,1]];operators=[];norms=[]
    for ids in groups:
        part=SimpleNamespace(K=model.K[ids][:,ids],M=model.M[ids][:,ids],size=len(ids),heat=model.heat,load=model.load[ids],original=model.original)
        P=task.Projection(part);operators.append(P)
        norms.append(float(np.linalg.norm(spsolve_triangular(P.C,P.forcing,lower=True))))
    whites=[np.empty((len(ids),size//2)) for ids in groups];smallest=1.;samples=[]
    for j in range(size//4):
        omega=frequency(j+1);u,error=full(model,4+1j*omega)
        for field,(ids,P,W) in enumerate(zip(groups,operators,whites)):
            v=P.C.T@u[ids]
            for phase,raw in enumerate([v.real,v.imag]):
                column=2*j+phase;norm=np.linalg.norm(raw);q=raw/norm
                for _ in range(2):q-=W[:,:column]@(W[:,:column].T@q)
                residual=np.linalg.norm(q);smallest=min(smallest,float(residual))
                assert residual>1e-13,(j,field,phase,float(residual))
                W[:,column]=q/residual
        samples.append(dict(frequency=omega,linear_residual=error))
    Q=np.zeros((model.size,size));scales=np.empty(size);checks=[]
    for field,(ids,P,W) in enumerate(zip(groups,operators,whites)):
        physical=spsolve_triangular(P.C.T.tocsr(),W,lower=False)
        columns=np.array([4*j+2*phase+field for j in range(size//4) for phase in range(2)])
        Q[np.ix_(ids,columns)]=physical;scales[columns]=norms[field]
        orth=float(np.max(abs(W.T@W-np.eye(size//2))));assert orth<1e-10
        checks.append(dict(field=field,orthogonality=orth,columns=size//2,force_norm=norms[field]))
    return None,None,Q,scales,dict(fields=checks,response_frequencies=samples,smallest_new_direction=smallest)


def basis(model,n):
    gram.multi.basis=field_basis
    return gram.basis(model,n)


def series(*args):
    gram.OUT=OUT
    return gram.series(*args)


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Repair the measured complex-frequency projection error using actual full coupled response snapshots, without growing the total subspace.',
        evidence='Full versus projected transfer dominates direct-versus-eigen and actual-mass effects on the declared6 frequencies. Positive-real trial recursions miss rapid coupled wave responses.',
        method='Nested32/64/128 full solves at z=4+i*1024*sin(pi*vdc(j)/2)^2, j starts1. Split real/imaginary responses into fluid/scalar mass-orthonormal columns, exactly total128/256/512. Each full snapshot retains both coupling directions. Same coupled Cholesky Gram and analytic heat poles. Do not fit physical source or output amplitudes.',
        holdout_frequencies=[0.,7.,23.,73.,233.,733.,1001.],
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13,pole_tail_relative=1e-12),
        decision='Pilot64. Fixed degree4 triple only within120s forecast. Reject on any original propagation gate; only after acceptance proceed to fixed degree2/1 and conditional contrasts under300s total. No alternative frequency grid,dimension or allocation after failure in this branch.',
        budget=dict(pilot_seconds=90,first_case_seconds=120,total_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(transfer.__file__),Path(gram.__file__),Path(space.__file__),Path(task.__file__),transfer.OUT/'scaled/transfer.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',gram.control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    amp=model.heat.amplitude;assert np.all((amp>=0).all(1)|(amp<=0).all(1))
    start=time.monotonic();u,error=full(model,4+512j);one=time.monotonic()-start
    start=time.monotonic();R,F,Q,S,tail,meta=basis(model,64);build=time.monotonic()-start;row=series(model,R,F,Q,S,tail,'pilot-64')
    # Prior full512 Gram timing includes fixed setup and quadratic trial work;
    # add all128 measured direct snapshot costs and pilot response allowance.
    prior=json.loads((previous.OUT/'stage-result.json').read_text())['seconds']
    forecast=1.4*(setup+2*prior+128*one+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,single_snapshot_seconds=one,build64_seconds=build,response64_seconds=row['seconds'],
        first_case_forecast_seconds=forecast,assumption='Twice measured prior512 Gram case plus128 full snapshot solves at the measured central frequency,setup and8 pilot responses,40 percent margin. Other frequencies and the larger QR cost remain assumptions.',**meta))
    signal.alarm(0);print('RESPONSE BASIS FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');amp=model.heat.amplitude;assert np.all((amp>=0).all(1)|(amp<=0).all(1))
    R,F,Q,S,tail,meta=basis(model,512)
    # Unseen frequencies check transfer coverage independently of the sampled
    # response columns. These are diagnostics, not substitutes for time gates.
    from scipy.linalg import solve
    holdouts=[];K=(model.K+model.K.T)/2;load=model.load+(K-model.K)@model.H
    A=R*S[None,:]/S[:,None]
    for omega in plan['holdout_frequencies']:
        z=4+1j*omega;E=transfer.source(model,z);u,residual=full(model,z)
        f=128.**2*(tail[0].T@(load@E))/S
        y=solve(z*z*A+128.**2*(np.eye(512)-A),f);projected=Q@(S*y)
        ref=transfer.readout(model,u,z,E);approx=transfer.readout(model,projected,z,E)
        errors=transfer.norms(model,approx-ref)/transfer.norms(model,ref)
        holdouts.append(dict(frequency=omega,relative=dict(zip(task.FIELDS,errors.astype(float)))))
    rows=[series(model,R[:n,:n],F[:n],Q[:,:n],S[:n],(tail[0][:,:n],tail[1],tail[2]),f'p4-{n}') for n in [128,256,512]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,comparisons=cmp,holdouts=holdouts,
        seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('RESPONSE BASIS RESULT',json.dumps({k:v for k,v in result.items() if k!='response_frequencies'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
