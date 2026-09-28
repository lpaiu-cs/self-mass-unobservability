"""Reproduce the finite-jump coupled response and its fixed-budget comparisons."""
from pathlib import Path
import argparse
import gc
import json
import resource
import signal
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec
import def_photon_finite_jump as m

OUT=m.OUT


def prepare():
    assert not (OUT/'coupling-plan.json').exists()
    prior=json.loads((OUT/'plan.json').read_text())
    for p,sha in prior['bindings'].items():
        source=OUT/'pilot-source.py' if p==Path(m.__file__).relative_to(m.old.h.ROOT).as_posix() else m.old.h.ROOT/p
        assert m.old.h.digest(source)==sha,p
    pilot=json.loads((OUT/'pilot.json').read_text())
    assert pilot['energy_residual']<1e-9 and pilot['solver_residual']<1e-11
    paths=[Path(__file__),Path(m.__file__),OUT/'plan.json',OUT/'pilot.json',OUT/'pilot-source.py',OUT/'matvec-pilot.json']
    m.old.ex.write(OUT/'coupling-plan.json',dict(classification='Counterexample candidate',
        reason='The four-mode pilot measured 12.65 s build + 10.28 s evolution. The initial 240 s production estimate omitted the cost of repeated nonlocal multiplies. Reassess BEFORE production: four kernel builds, all predeclared comparisons, roughly 300-650 s, hard cap 720 s. No accuracy gate or physical duration changed.',
        forecast_seconds=[300,650],hard_seconds=720,memory_cap_GB=4,CPU_workers=1,
        new_EOS_calls=0,new_physical_queries=0,new_stellar_steps=0,automatic_expansion=False,
        source_change='Split real/imaginary sparse matvec avoids full matrix dtype promotion; measured vector difference zero. The continuum normalization control uses all terms of the cited Eq.16, not only its linear photon-energy truncation.',
        bindings={p.relative_to(m.old.h.ROOT).as_posix():m.old.h.digest(p) for p in paths}))
    print('PRODUCTION PLAN',json.loads((OUT/'coupling-plan.json').read_text()),flush=True)


def narrow_tail(theta):
    x,w=leggauss(48);t=(x+1)/2;mu=1-2*t*t;aw=2*w*t;rows=[]
    for u in [.1,1,5,15,50]:
        def f(z):
            v=u*np.exp(2*np.sqrt(theta)*z)
            return float(np.sum(m.kernel(u,v,mu,theta)*aw))*2*np.sqrt(theta)*v
        tail=sum(quad_vec(f,a,b,epsabs=1e-16,epsrel=1e-8)[0] for a,b in [(-9,-6),(6,9)])
        rows.append(dict(u=u,probability_outside_pair_window=float(tail)))
    return rows


def jump_rhs(op,T,E):
    z=E.reshape(-1,op.order);r=-op.diagonal*z-op.off(E).reshape(z.shape)
    r[:,0]-=op.jump_q*T
    return -op.jump_aa*T-op.jump_q@z[:,0],r.ravel()


def algebra(op,b):
    root=np.sqrt(op.Ci);e=np.zeros(op.size);e[::op.order]=root
    q=np.zeros(op.size);q[::op.order]=root/op.u
    vE=np.r_[np.sqrt(op.Cm),e];vN=np.r_[0,q];values=[]
    for v in [vE,vN]:
        a,z=jump_rhs(op,v[0],v[1:]);values.append(float(np.linalg.norm(np.r_[a,z])/(float(b['rate_e'])*np.linalg.norm(v))))
    heating=op.jump_aa*op.Cm
    continuum=4*float(b['arad'])*float(b['T'])**3*float(b['rate_C'])
    return dict(energy_number_null_relative=values,initial_heating=heating,
        continuum_heating=continuum,heating_relative=float(abs(heating/continuum-1)))


def diagnostics(op,r):
    return dict(balance=r['balance'],energy_residual=r['energy_equation_residual'],
        entropy_growth=r['entropy_growth'],solver_residual=op.max_solver_error,iterations=op.max_iterations)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    plan=json.loads((OUT/'coupling-plan.json').read_text());original=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert m.old.h.digest(m.old.h.ROOT/p)==sha,p
    signal.alarm(plan['hard_seconds']);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    b=dict(np.load(m.old.OUT/'bank.npz'));duration=float(b['H']/b['c']);rows=[]
    symbolic=m.old.symbolic();symbolic['finite_jump_extension']='For each frequency pair, angular moments obey |alpha_l|<=alpha_0 by positive quadrature and |P_l|<=1. Each l>0 two-bin block is PSD. The l=0 recoil bracket is the same positive energy/number-orthogonal outer product for ANY frequency jump. Same-bin angular loss is nonnegative. The resulting whole collision operator is PSD and streaming is skew-Hermitian.'
    m.old.ex.write(OUT/'symbolic.json',symbolic)
    moment=m.moments(float(b['theta']));tails=narrow_tail(float(b['theta']))
    m.old.ex.write(OUT/'continuum-controls.json',dict(classification='Counterexample candidate',moments=moment,tails=tails))
    op=m.Operator(b,8,1,24,1);runs=[op.evolve(duration,n) for n in [16,32,64]]
    error=[m.old.compare(a,z,op) for a,z in zip(runs[:-1],runs[1:])];order=float(np.log2(error[0]/error[1]))
    reference=runs[-1];base_algebra=algebra(op,b)
    rows.append(dict(name='base-8-24-1',coefficients=op.info,time_order=order,time_differences_initial=error,
        algebra=base_algebra,**diagnostics(op,reference)))
    np.savez_compressed(OUT/'base.npz',T=reference['T'],E=reference['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=op.Cm,duration=duration)
    m.old.ex.write(OUT/'base-result.json',rows[-1]);print('BASE',rows[-1],flush=True)
    del op;gc.collect()
    comparisons={}
    for name,L,nangle,nfreq in [('angular',12,24,1),('scattering-quadrature',8,48,1),('frequency-quadrature',8,24,2)]:
        op=m.Operator(b,L,1,nangle,nfreq);r=op.evolve(duration,64)
        E=r['E'].reshape(-1,L);R=reference['E'].reshape(-1,8)
        delta=float(np.sqrt(abs(r['T']-reference['T'])**2+np.linalg.norm(E[:,:8]-R)**2+np.linalg.norm(E[:,8:])**2)/np.sqrt(op.Cm))
        comparisons[name]=delta;alg=algebra(op,b)
        row=dict(name=name,coefficients=op.info,difference_initial=delta,algebra=alg,**diagnostics(op,r));rows.append(row)
        m.old.ex.write(OUT/(name+'.json'),row);print('COMPARISON',row,flush=True)
        if nfreq==2:
            np.savez_compressed(OUT/'frequency-reference.npz',T=r['T'],E=E,u=op.u,Ci=op.Ci,Cm=op.Cm,duration=duration)
        del op;gc.collect()
    # Same free-electron population for deciding diffusion vs finite jumps.
    controls=[]
    for name,bank,options in [
        ('free-electron-Kompaneets',dict(b,rate_s=np.full_like(b['rate_s'],float(b['rate_e']))),{}),
        ('free-electron-elastic',dict(b,rate_s=np.full_like(b['rate_s'],float(b['rate_e']))),dict(compton=False)),
        ('native-elastic',b,dict(compton=False))]:
        op=m.old.Operator(bank,8,1,**options);r=op.evolve(duration,64)
        delta=m.old.compare(r,reference,op)
        row=dict(name=name,difference_from_finite_jump_initial=delta,
            temperature_K=float((r['T']/np.sqrt(op.Cm)).real),
            temperature_difference_K=float(((r['T']-reference['T'])/np.sqrt(op.Cm)).real))
        controls.append(row)
        np.savez_compressed(OUT/(name+'.npz'),T=r['T'],E=r['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=op.Cm)
        del op
    finite_temperature=float((reference['T']/np.sqrt(float(b['Cm']))).real)
    original_path=np.load(m.old.OUT/'mode-1.npz')
    original_difference=float(np.sqrt(abs(original_path['T']-reference['T'])**2
        +np.linalg.norm(original_path['E'][:,:8]-reference['E'].reshape(-1,8))**2
        +np.linalg.norm(original_path['E'][:,8:])**2)/np.sqrt(float(b['Cm'])))
    gates=original['gates']
    passed=bool(order>=gates['time_order'] and error[-1]<.001 and max(comparisons.values())<.001
        and all(r['balance']<1e-9 and r['energy_residual']<1e-9 and r['entropy_growth']<1e-10
            and r['solver_residual']<1e-11 and max(r['algebra']['energy_number_null_relative'])<1e-10
            and r['algebra']['heating_relative']<.005 for r in rows)
        and all(r['normalization_error']<1e-7 and max(r['drift_relative'],r['diffusion_relative'])<.005 for r in moment))
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,
        comparisons=comparisons,controls=controls,final_temperature_K=finite_temperature,
        original_phase61_difference_initial=original_difference,
        maximum_moment_normalization_error=max(r['normalization_error'] for r in moment),
        maximum_cutoff_probability=max(r['probability_outside_pair_window'] for r in tails),
        seconds=time.monotonic()-start,peak_memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        free_electron_finite_jump_coupled=True,native_medium_scattering_certified=False,
        actual_temperature_opacity_certified=False,physical_atmosphere_closed=False,
        full_GR_photon_feedback_evolved=False,full_dynamic_charge_solved=False)
    m.old.ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
