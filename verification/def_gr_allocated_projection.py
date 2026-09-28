"""Reallocate fixed trial dimension toward the unresolved fluid field."""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import signal
import time
import resource
import numpy as np
import def_gr_inverse_gram as gram

multi=gram.multi;task=gram.task;space=gram.space;common=gram.common
OUT=task.OUT.parent/'def-gr-allocated-projection';write=task.write


def field_basis(model,size):
    assert size%8==0
    groups=[model.indices[:model.surface_index+1,0],model.indices[:-1,1]]
    columns=[np.flatnonzero(np.arange(size)%8!=7),np.arange(7,size,8)]
    Q=np.zeros((model.size,size));scales=np.empty(size);checks=[]
    for j,(indices,cols) in enumerate(zip(groups,columns)):
        part=SimpleNamespace(K=model.K[indices][:,indices],M=model.M[indices][:,indices],size=len(indices),
            heat=model.heat,load=model.load[indices],original=model.original)
        local,scale,check=multi.local_basis(part,len(cols))
        Q[np.ix_(indices,cols)]=local;scales[cols]=scale
        checks.append(dict(field=j,columns=len(cols),force_norm=float(scale),**check))
    # Gram projection consumes Q and its scales, not a separately projected K.
    # Both full coupling blocks enter the original shifted Cholesky afterward.
    return None,None,Q,scales,dict(fields=checks,allocation='7 fluid : 1 scalar')


def basis(model,n):
    multi.basis=field_basis
    return gram.basis(model,n)


def series(*args):
    gram.OUT=OUT
    return gram.series(*args)


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    previous=json.loads((common.OUT/'inverse-gram/stage-result.json').read_text())
    scalar=previous['comparisons']['scalar_mass_RMS']
    forecast_scalar=scalar['last']*4**scalar['order']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='caf5c5b',
        claim='Spend the same total trial dimension on unresolved fluid dynamics while retaining the weak scalar and full coupling, instead of increasing time steps or dimensions.',
        evidence='SDIRK and refined Gauss distinguish unresolved global/old velocity from weak-readout arithmetic. Equal-allocation Gram already has scalar order2.64 and last0.0001704; reducing final scalar dimension256 to64 forecasts about0.00662 if that empirical order persists. This forecast is not acceptance.',
        predicted_scalar_difference=forecast_scalar,
        method='Nested total128/256/512 dimensions now have fluid112/224/448 and scalar16/32/64. Same4,16,64,256,1024 shifts, full coupled Cholesky Gram, source-norm coordinate balance and bounded analytic heat response. One complete solution per level, no field splicing.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13,pole_tail_relative=1e-12),
        decision='Pilot64, then degree4 first only within120s forecast. Stop any failed component; no alternative ratios or shifts in this branch. Conditional fixed degrees2/1 and original three contrasts under300s total. No larger basis,grid,time period or source.',
        budget=dict(pilot_seconds=90,first_case_seconds=120,total_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(gram.__file__),Path(multi.__file__),Path(space.__file__),Path(task.__file__),common.OUT/'inverse-gram/stage-result.json',task.OUT.parent/'def-gr-gauss-refined/arithmetic-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',gram.control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    amp=model.heat.amplitude;assert np.all((amp>=0).all(1)|(amp<=0).all(1))
    start=time.monotonic();R,F,Q,S,tail,meta=basis(model,64);build=time.monotonic()-start
    row=series(model,R,F,Q,S,tail,'pilot-64')
    # A fixed total has at most twice the old equal-field quadratic work; the
    # previous actual512 run measured all setup, responses and factorization.
    forecast=1.4*(setup+2*previous['seconds']+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,build64_seconds=build,response64_seconds=row['seconds'],
        first_case_forecast_seconds=forecast,assumption='Same512 total. At most2 times equal-split quadratic orthogonalization work; use2 times the prior measured24.06s entire Gram case plus setup,pilot response scaling and40 percent allowance. New allocation timing is unmeasured at512.',**meta))
    signal.alarm(0);print('ALLOCATION FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');amp=model.heat.amplitude;assert np.all((amp>=0).all(1)|(amp<=0).all(1))
    R,F,Q,S,tail,meta=basis(model,512)
    rows=[series(model,R[:n,:n],F[:n],Q[:,:n],S[:n],(tail[0][:,:n],tail[1],tail[2]),f'p4-{n}') for n in [128,256,512]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,comparisons=cmp,
        seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('ALLOCATION RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
