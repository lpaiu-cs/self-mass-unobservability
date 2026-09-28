"""Importance sampling the same finite-T collision integral at fixed cost.

The old logistic proposal undersamples high polynomial moments in warm
states. A defensive gamma mixture changes the proposal, not the physics.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import numpy as np
from scipy.special import expit,gammainc,gammaincinv,gammaln
import def_ionic_structure_transport as parent

ex=parent.ex
h=parent.h
OUT=parent.OUT/'warm-quadrature'
CELLS=[1722,2000,2500,2800,3000]
original_energy=ex.energy
original_events=ex.events


def proposal(u,state):
    eta=state['eta'];top=eta+40;normal=expit(eta)-expit(-40.)
    branch=u<.5;v=np.where(branch,2*u,2*u-1)
    x=original_energy(v,state)[0]
    # ponytail: one fixed defensive mixture, no learned proposal or tuning.
    # shape 8 / scale 1.5 covers the degree-six collision Gram tails.
    cutoff=gammainc(8.,top/1.5)
    x[~branch]=1.5*gammaincinv(8.,v[~branch]*cutoff)
    f=expit(eta-x);old=f*(1-f)/normal
    tail=np.exp(7*np.log(x)-x/1.5-gammaln(8.)-8*np.log(1.5))/cutoff
    density=.5*(old+tail);ratio=old/density
    assert np.all((x>=0)&(x<=top)) and np.all((ratio>0)&(ratio<=2))
    return x,f,normal,ratio


def events(sample,state,fixed_p=None):
    assert fixed_p is None,'This change is only the two-incoming-energy bracket.'
    ratios=[]
    def sampled(u,s):
        x,f,n,r=proposal(u,s);ratios.append(r);return x,f,n
    ex.energy=sampled
    try:w,delta,checks=original_events(sample,state)
    finally:ex.energy=original_energy
    assert len(ratios)==2
    return w*ratios[0]*ratios[1],delta,checks


def controls():
    from scipy.integrate import quad
    rows=[]
    for index in [1722,3000,5734]:
        s=dict(np.load(parent.OUT/f'cell-{index}.npz'));u=np.linspace(1e-8,1-1e-8,1024)
        x,f,n,r=proposal(u,s);branch=u<.5;v=np.where(branch,2*u,2*u-1)
        recovered=np.where(branch,(expit(s['eta'])-f)/n,gammainc(8.,x/1.5)/gammainc(8.,(s['eta']+40)/1.5))
        error=float(np.max(abs(v-recovered)));assert error<1e-12
        integral=quad(lambda x:expit(s['eta']-x)*expit(x-s['eta'])/n,0,s['eta']+40,epsabs=1e-12)[0]
        assert abs(integral-1)<1e-12
        rows.append(dict(index=index,inverse_CDF_absolute=error,old_density_normalization=integral,max_weight=float(r.max())))
    return dict(classification='Proven',passed=True,checks=rows,
        identity='q_mix=(q_FD+q_gamma)/2. Multiplying the old collision weight by q_FD(x1)q_FD(x2)/(q_mix(x1)q_mix(x2)) leaves the integral unchanged. Each ratio lies in (0,2], the event measure remains positive, and event kinematics/conservation are unchanged.',
        limitation='Floating controls check the implementation; quadrature convergence and finite energy/basis truncation still require their independent gates.')


def cell(index):
    start=time.monotonic();saved=dict(np.load(parent.OUT/f'cell-{index}.npz'));s=saved
    ex.amplitude=parent.previous.microscopic.leading.amplitude(s,s['phase'],s['polarization'])
    ex.events=events;fine=[];coarse=[];checks=[]
    try:
        for seed in ex.SEEDS:
            a,b,c=ex.bracket(s,17,seed,15);fine.append(a);coarse.append(b);checks.append(c)
    finally:ex.events=original_events
    EE=np.mean(fine,axis=0);smallEE=np.mean(coarse,axis=0)
    answer=ex.transfer(s,EE,7);small=ex.transfer(s,smallEE,7);low=ex.transfer(s,EE,5)
    errors={k:dict(quadrature=float(abs(small[k]/answer[k]-1)),basis=float(abs(low[k]/answer[k]-1))) for k in ['K','tau']}
    changed={k:float(answer[k]/saved['fine_'+k]-1) for k in ['K','tau']}
    saved.update(EE=EE,coarse_EE=smallEE,replicate_EE=fine,**{'fine_'+k:v for k,v in answer.items()})
    np.savez_compressed(OUT/f'cell-{index}.npz',**saved)
    row=dict(index=index,passed=all(e['quadrature']<.03 and e['basis']<.02 for e in errors.values()),
        errors=errors,K=answer['K'],tau=answer['tau'],native_ratio=answer['K']/s['native_K'],
        change_from_failed=changed,event_checks=checks,seconds=time.monotonic()-start)
    ex.write(OUT/f'cell-{index}.json',row);print('CELL',index,row,flush=True);return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(parent.__file__),parent.OUT/'plan.json',parent.OUT/'result.json']+[parent.OUT/f'cell-{i}.npz' for i in CELLS]
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Resolve the warm collision-moment quadrature failure by an exactly reweighted incoming-energy proposal. No physical kernel, EI matrices, state, event count, basis order or acceptance threshold changes.',
        proposal='Equal mixture of the original truncated FD-logistic density and Gamma(shape=8,scale=1.5) truncated to the same 0..eta+40 domain. Fixed before the new results; bounded density ratio <=2 for each incoming particle.',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},cells=CELLS,seeds=ex.SEEDS,powers=[15,17],basis=[5,7],
        gates=dict(quadrature=.03,basis=.02),
        budget=dict(pilot_cells=1,production_cells=4,hard_seconds=90,workers=1,BLAS_threads=1,new_native_calls=0,new_stellar_steps=0,automatic_expansion=False)))
    ex.write(OUT/'controls.json',controls());row=cell(CELLS[0])
    ex.write(OUT/'pilot.json',dict(seconds=row['seconds'],remaining_forecast_seconds=4*row['seconds'],
        stop='If the pilot fails, do not automatically repeat the remaining states or increase samples. Reassess the unresolved mathematical or physical obstruction.'))


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    first=json.loads((OUT/f'cell-{CELLS[0]}.json').read_text());assert first['passed'],'pilot gate'
    assert json.loads((OUT/'pilot.json').read_text())['remaining_forecast_seconds']<80
    start=time.monotonic();rows=[first]
    for index in CELLS[1:]:
        rows.append(cell(index));assert time.monotonic()-start<85
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),checks=rows,
        seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        old_failed_result_preserved=True,physics_and_thresholds_unchanged=True,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result['passed'],'SECONDS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
