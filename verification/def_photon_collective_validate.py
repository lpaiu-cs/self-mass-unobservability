"""Fixed-budget EOS-consistent collective photon/material paths."""
from pathlib import Path
import argparse
import gc
import json
import resource
import signal
import time
import numpy as np
import def_photon_collective as m
import def_photon_finite_jump_validate as checks

OUT=m.OUT


def prepare():
    assert not (OUT/'coupling-plan.json').exists()
    paths=[Path(__file__),Path(m.__file__),Path(m.previous.__file__),Path(checks.__file__),
        OUT/'inventory.npz',OUT/'inventory.json',OUT/'pilot.json',OUT/'pilot-source.py',
        OUT/'controls.json',OUT/'symbolic.json',m.old.OUT/'bank.npz',m.previous.OUT/'manifest.json']
    assert json.loads((OUT/'controls.json').read_text())['maximum_relative']<1e-4
    m.old.ex.write(OUT/'coupling-plan.json',dict(classification='Counterexample candidate',checkpoint='9b1fce346',
        claim='Evolve the same material/absorption/Fourier mode using its actual EOS electron and multi-ion collective density spectrum. Integrate resonances over cells, preserve number/energy, and compare the collective correction against the SAME independent-electron discretization.',
        primary=dict(velocity_quadrature=8,q_grid=129,angles=24,modes=8,steps=[16,32,64],kH=1),
        comparisons=[dict(name='velocity',velocity_quadrature=16,q_grid=129,angles=24,modes=8),
            dict(name='q-grid',velocity_quadrature=16,q_grid=257,angles=24,modes=8),
            dict(name='angular',velocity_quadrature=16,q_grid=257,angles=48,modes=12)],
        control='Analytic independent Maxwell Gaussian cell convolution at 32/64 steps; compare with Phase62 free-electron finite jumps and pair RPA-minus-free differences at both time steps.',
        quadrature='Six fixed positive-velocity intervals [0,.025,.1,.5,2,5,9], 8/16 Gauss-Legendre nodes EACH, 48/96 nodes per species. Ten actual ionic mass groups plus electrons. SVD of Gibbs-whitened Vlasov matrix. Nested q grids reuse common spectral states without recomputation. No rate/opacity normalization.',
        pilot='Old 4-mode, 16-step Hermite-velocity pilot: 16.005 s total, 10.996 s build. Composite velocity quadrature resolves the low phase-velocity dielectric control; independent controls took 3.813 s.',
        forecast_seconds=[300,650],hard_seconds=900,memory_cap_GB=4,CPU_workers=1,GPU=False,
        new_EOS_calls=0,new_opacity_queries=0,new_stellar_steps=0,automatic_expansion=False,
        gates=dict(time_order=1.8,time_difference_initial=.001,paired_time_difference_initial=1e-6,
            velocity_difference_initial=1e-6,q_grid_difference_initial=1e-6,angular_difference_initial=1e-6,
            free_model_difference_initial=1e-5,energy_balance=1e-9,energy_residual=1e-9,
            null_modes=1e-10,solver_residual=1e-11,entropy_growth=1e-10),
        limitations='Collisionless weak-coupling Maxwell/RPA, FDT and vacuum-photon kinematics. Native scalar opacity is NOT assumed a total differential cross section. Bound-electron dynamics, nonideal/collisional corrections, sub-plasma-frequency propagation, actual-temperature absorption and full GR remain separate requirements.',
        bindings={p.relative_to(m.old.h.ROOT).as_posix():m.old.h.digest(p) for p in paths}))
    print('PREPARED PRODUCTION',OUT,flush=True)


def distance(T,E,U,W,Cm):return float(np.sqrt(abs(T-U)**2+np.linalg.norm(E-W)**2)/np.sqrt(Cm))


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();plan=json.loads((OUT/'coupling-plan.json').read_text())
    for p,sha in plan['bindings'].items():assert m.old.h.digest(m.old.h.ROOT/p)==sha,p
    signal.alarm(plan['hard_seconds']);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    b=dict(np.load(m.old.OUT/'bank.npz'));duration=float(b['H']/b['c']);Cm=float(b['Cm']);rows=[]
    base=dict(b,velocity_order=8,q_points=129);op=m.Operator(base,8,1,24)
    runs=[op.evolve(duration,n) for n in [16,32,64]]
    errors=[m.old.compare(a,z,op) for a,z in zip(runs[:-1],runs[1:])];time_order=float(np.log2(errors[0]/errors[1]))
    reference=runs[-1];middle=runs[-2]
    row=dict(name='base',coefficients=op.info,time_order=time_order,time_differences_initial=errors,
        algebra=checks.algebra(op,b),**checks.diagnostics(op,reference));rows.append(row)
    np.savez_compressed(OUT/'base.npz',T=reference['T'],E=reference['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    m.old.ex.write(OUT/'base-result.json',row);print('BASE',row,flush=True);del op;gc.collect()
    previous=reference;comparisons={}
    for spec in plan['comparisons']:
        name=spec['name'];L=spec['modes'];bank=dict(b,velocity_order=spec['velocity_quadrature'],q_points=spec['q_grid'])
        op=m.Operator(bank,L,1,spec['angles']);r=op.evolve(duration,64)
        z=r['E'].reshape(-1,L);p=previous['E'].reshape(-1,8)
        delta=distance(r['T'],z[:,:8],previous['T'],p,Cm)
        if L>8:delta=float(np.sqrt(delta*delta+np.linalg.norm(z[:,8:])**2/Cm))
        comparisons[name]=delta
        row=dict(name=name,difference_initial=delta,coefficients=op.info,algebra=checks.algebra(op,b),**checks.diagnostics(op,r));rows.append(row)
        m.old.ex.write(OUT/(name+'.json'),row);print('COMPARISON',row,flush=True)
        np.savez_compressed(OUT/(name+'.npz'),T=r['T'],E=z,u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
        if name=='q-grid':best8=r
        previous=r;del op;gc.collect()
    free=dict(b,collective_free=True);op=m.Operator(free,8,1,24)
    f32=op.evolve(duration,32);f64=op.evolve(duration,64)
    paired=distance(reference['T']-f64['T'],reference['E']-f64['E'],middle['T']-f32['T'],middle['E']-f32['E'],Cm)
    collective_difference=m.old.compare(best8,f64,op)
    old_free=np.load(m.previous.OUT/'base.npz')
    free_difference=distance(f64['T'],f64['E'].reshape(-1,8),old_free['T'],old_free['E'],Cm)
    np.savez_compressed(OUT/'independent-electron.npz',T=f64['T'],E=f64['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    row=dict(name='independent-electron',algebra=checks.algebra(op,b),**checks.diagnostics(op,f64));rows.append(row)
    # Native scalar-rate transport-vs-total comparison is diagnostic, not a fit.
    s=np.load(OUT/'inventory.npz');mu,w=m.leggauss(128);Z=float(s['ion_strength'].sum());u=op.u
    x2=2*(1-mu)*(u[:,None]*float(s['wave_lambda']))**2;S=(x2+Z)/(x2+1+Z)
    transport=(S*(3/8*(1+mu**2)*(1-mu)*w)).sum(1)
    totals=(S*(3/8*(1+mu**2)*w)).sum(1)
    np.savez_compressed(OUT/'transport-comparison.npz',u=u,total=totals,transport=transport,native=b['rate_s'][:len(u)]/b['rate_e'])
    photon_capacity_below_plasma=float(op.Ci[u<float(s['plasma_u'])].sum()/op.Ci.sum())
    passed=bool(time_order>=1.8 and errors[-1]<.001 and max(comparisons.values())<1e-6 and paired<1e-6 and free_difference<1e-5
        and all(r['balance']<1e-9 and r['energy_residual']<1e-9 and r['solver_residual']<1e-11
            and r['entropy_growth']<1e-10 and max(r['algebra']['energy_number_null_relative'])<1e-10 for r in rows))
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,comparisons=comparisons,
        paired_time_difference_initial=paired,free_reference_difference_initial=free_difference,
        collective_effect_initial=collective_difference,collective_temperature_effect_K=float(((best8['T']-f64['T'])/np.sqrt(Cm)).real),
        below_plasma_LTE_capacity_fraction=photon_capacity_below_plasma,
        seconds=time.monotonic()-start,peak_memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        EOS_population_matched_collective_response_evolved=True,bound_electron_kernel_complete=False,
        nonideal_collisional_quantum_error_certified=False,native_scalar_kernel_identified=False,
        actual_temperature_absorption_certified=False,physical_atmosphere_closed=False,
        full_GR_feedback_evolved=False,full_dynamic_charge_solved=False)
    m.old.ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
