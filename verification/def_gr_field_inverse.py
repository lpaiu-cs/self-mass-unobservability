"""Combine exact field coordinates with the full coupled inverse projection."""
from pathlib import Path
import argparse
import time
import json
import signal
import resource
import numpy as np
import def_gr_multishift as multi

fields=multi.fields;common=multi.common;task=multi.task;space=multi.space
OUT=common.OUT/'field-inverse';write=task.write
original_diagonalize=fields.diagonalize


def diagonalize(R,load,Q,scales):
    theta,forcing,physical,meta=original_diagonalize(R,load,Q,scales)
    assert theta.min()>0
    k=128.**2*(1-theta)/theta;assert k.min()>-1.
    return k,forcing/theta[:,None],physical,dict(resolvent_eigenvalue_minimum=float(theta.min()),**meta)


def basis(model,size):
    _,_,Q,scales,checks=multi.basis(model,size)
    P=task.Projection(model);white=P.C.T@Q
    orth=float(np.max(abs(white.T@white-np.eye(size))));assert orth<1e-10
    R,load,meta=common.projected(P,white)
    return R,load,Q,scales,dict(full_mass_orthogonality=orth,**checks,**meta)


def series(*args):
    fields.diagonalize=diagonalize;fields.OUT=OUT
    return fields.series(*args)


def control():
    inverse=common.control();similarity=fields.control()
    matrix=np.array([[9.,1e-11],[1e-11,25.]]);force=np.array([[1.],[2e-12]])
    R=128.**2*np.linalg.inv(matrix+128.**2*np.eye(2))
    k,f,Q,meta=diagonalize(R,R@force,np.eye(2),np.array([1.,2e-12]));errors=[]
    for z in [2+3j,6+100j]:
        actual=Q@(f[:,0]/(k+z*z));expected=np.linalg.solve(matrix+z*z*np.eye(2),force[:,0])
        errors.append(float(np.max(abs(actual-expected)/abs(expected))))
    assert max(errors)<1e-9
    return dict(classification='Proven',inverse_identity=inverse['equivalence_checked'],similarity_identity=similarity['similarity_and_source_identity_checked'],
        numerical_classification='Counterexample candidate',combined_transfer_errors=errors,
        scope='Full-space identities and two-field transfer only. All actual component, space and contrast gates remain mandatory.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Combine field-separated, input-balanced arithmetic with a full coupled inverse projection so the same transformed operator and source govern both components.',
        method='Reuse fixed4/16/64/256/1024 shifted trial directions and total128/256/512 columns. Project the full coupled R, R-transformed source and inertia, then the same input-norm similarity. No one-way approximation, eigenvalue clipping or component splicing.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Degree4 first; stop on any component failure. Only if all pass run fixed degree2/1 and conditional coefficient/outer/quadrature cases. No larger dimensions, cells, degree or horizon.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=120,total_cap_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(multi.__file__),Path(fields.__file__),Path(common.__file__),Path(space.__file__),Path(task.__file__),common.OUT/'multishift/stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    start=time.monotonic();R,F,Q,S,meta=basis(model,64);build=time.monotonic()-start;row=series(model,R,F,Q,S,'pilot-64')
    # Use actual512-dimensional multishift build; only the additional inverse
    # projection scales from its separate measured64-dimensional timing.
    prior=json.loads((common.OUT.parent/'multishift/stage-result.json').read_text()) if common.OUT.name=='field-inverse' else json.loads((common.OUT/'multishift/stage-result.json').read_text())
    inverse=json.loads((OUT.parent/'pilot-budget.json').read_text())
    forecast=1.4*(setup+prior['seconds']+64*inverse['projection64_seconds']+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,build64_seconds=build,
        response64_seconds=row['seconds'],first_case_forecast_seconds=forecast,
        assumption='Charge the entire measured21.67s multishift path plus new setup, quadratically scaled inverse projection and full response with40 percent allowance; no unmeasured later cases launched.',**meta))
    signal.alarm(0);print('FIELD INVERSE FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');R,F,Q,S,meta=basis(model,512)
    rows=[series(model,R[:n,:n],F[:n],Q[:,:n],S[:n],f'p4-{n}') for n in [128,256,512]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('FIELD INVERSE RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
