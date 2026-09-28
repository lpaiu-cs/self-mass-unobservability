"""Fixed broad-band rational bases for the unchanged two-field GR operator."""
from pathlib import Path
from types import SimpleNamespace
import argparse
import time
import signal
import json
import resource
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu,spsolve_triangular
import def_gr_field_projection as fields

task=fields.task;space=fields.space;common=fields.common;OUT=common.OUT/'multishift';write=task.write
FREQUENCIES=(4.,16.,64.,256.,1024.)


def local_basis(part,size):
    operators=[]
    for frequency in FREQUENCIES:
        P=task.Projection(part);P.sigma=frequency**2;P.shift=part.K+P.sigma*part.M
        P.scale=np.sqrt(P.shift.diagonal());D=diags(1/P.scale)
        P.lu=splu((D@P.shift@D).tocsc());P.extended=P.shift.astype(np.longdouble)
        operators.append(P)
    Q=np.empty((part.size,size),order='F');Q[:,0]=operators[0].initial
    for j in range(1,size):
        P=operators[(j-1)%len(operators)];q=P.sigma*(P.C.T@P.solve(P.C@Q[:,j-1]))
        for _ in range(2):q-=Q[:,:j]@(Q[:,:j].T@q)
        norm=np.linalg.norm(q);assert norm>1e-14,(j,norm);Q[:,j]=q/norm
    orth=float(np.max(abs(Q.T@Q-np.eye(size))));assert orth<1e-10
    P=operators[0];physical=spsolve_triangular(P.C.T.tocsr(),Q,lower=False)
    scale=np.linalg.norm(spsolve_triangular(P.C,P.forcing,lower=True))
    return physical,scale,dict(orthogonality=orth,linear_residual=max(P.error for P in operators))


def basis(model,size):
    groups=[model.indices[:model.surface_index+1,0],model.indices[:-1,1]]
    Q=np.zeros((model.size,size));scales=np.empty(size);checks=[]
    for j,indices in enumerate(groups):
        part=SimpleNamespace(K=model.K[indices][:,indices],M=model.M[indices][:,indices],size=len(indices),
            heat=model.heat,load=model.load[indices],original=model.original)
        local,scale,check=local_basis(part,size//2)
        Q[np.ix_(indices,np.arange(j,size,2))]=local;scales[j::2]=scale
        checks.append(dict(field=j,force_norm=float(scale),**check))
    K,skew=model.energy_matrix(Q)
    return K,Q.T@model.load,Q,scales,dict(fields=checks,stiffness_skew=skew)


def series(*args):
    fields.OUT=OUT
    return fields.series(*args)


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Use the same total basis budget across a declared broad frequency range to retain both rapid fluid response and weak scalar response in one coupled solve.',
        method='Same separated fluid/scalar bases and force-norm similarity. Cycle fixed positive frequencies4,16,64,256,1024 instead of repeating128. Same total128/256/512 columns, physical stiffness, source, polynomial cells and exact heat-pole propagation. Frequencies are not fit to output amplitudes.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Degree4 first, stop on any propagation failure. Only on simultaneous four-field acceptance run already fixed degree2/1, then conditional coefficient/outer/quadrature contrasts. No greater degree, dimensions or horizon.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=120,total_cap_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(fields.__file__),Path(common.__file__),Path(space.__file__),Path(task.__file__),common.OUT/'gauss/stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',fields.control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    start=time.monotonic();K,F,Q,S,meta=basis(model,64);build=time.monotonic()-start;row=series(model,K,F,Q,S,'pilot-64')
    # Direct scaling of all pilot setup charges is conservative and recorded.
    forecast=1.4*(setup+64*build+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,build64_seconds=build,
        response64_seconds=row['seconds'],first_case_forecast_seconds=forecast,
        assumption='Quadratic64-to512 scaling includes repeated shifted factorizations;40 percent allowance, first-case120s and total300s unchanged.',**meta))
    signal.alarm(0);print('MULTISHIFT FORECAST',forecast,flush=True)


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
    write(OUT/'stage-result.json',result);signal.alarm(0);print('MULTISHIFT RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
