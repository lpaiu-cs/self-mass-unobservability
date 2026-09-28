"""Reassessed fixed-budget recovery of the rejected discrete-velocity response."""
from pathlib import Path
import argparse
import gc
import json
import resource
import signal
import time
import numpy as np
import def_photon_collective_continuum as c
import def_photon_collective_validate as previous
from def_photon_finite_jump_validate import algebra,diagnostics

m=c.m;OUT=c.OUT


def prepare():
    assert not (OUT/'coupling-plan.json').exists()
    assert json.loads((OUT/'uniform-log-pilot.json').read_text())['passed']
    paths=[Path(__file__),Path(c.__file__),Path(m.__file__),Path(previous.__file__),
        OUT/'uniform-log-pilot.json',m.OUT/'result.json',m.OUT/'independent-electron.npz',
        m.OUT/'inventory.npz',m.old.OUT/'bank.npz']
    m.old.ex.write(OUT/'coupling-plan.json',dict(classification='Counterexample candidate',
        checkpoint='9b1fce346',reason='Discrete velocity atoms cause response aliasing. Fix the representation, not the rejected thresholds. An arctan-to-linear interpolation pilot also failed; log-distance wings fix that separate interpolation defect. Preserve all failures.',
        claim='Resolve the EOS-matched collective correction in the SAME 2.316801 ms matter/absorption/photon/spatial problem with unchanged response and conservation gates.',
        reuse='Native EOS and opacity, photon cells, SDIRK2 solver, prior free-electron 64-step path; no repeated long GR history. Source-cell geometry and atomic opacity remain conditional.',
        primary=dict(mesh=128,q_grid=129,angles=24,modes=8,steps=[16,32,64]),
        comparisons=[dict(name='density-mesh',mesh=256,q_grid=129,angles=24,modes=8),
            dict(name='q-grid',mesh=256,q_grid=257,angles=24,modes=8),
            dict(name='angular',mesh=256,q_grid=257,angles=48,modes=12)],
        controls='Free Gaussian 32-step path only; reuse saved 64-step result. Pair collective-minus-free at 32/64. Preserve previous Gaussian-vs-Sazonov control.',
        gates=dict(time_order=1.8,time_difference_initial=.001,response_comparison_initial=1e-6,
            paired_time_difference_initial=1e-6,energy_balance=1e-9,energy_residual=1e-9,
            null_modes=1e-10,solver_residual=1e-11,entropy_growth=1e-10,spectral_moments=2e-4),
        measured='Discrete original: 604.928 s total, 2.279 GB; cached kernel 17.782 s at 48 angles. Continuous ten-x 32/64/128 control took 0.122 s; all-q 64/128 tables completed in the same sub-4-second command including imports. Extra arithmetic in continuous CDF queries is unmeasured at full size.',
        forecast_seconds=[400,850],hard_seconds=1200,memory_cap_GB=4,CPU_workers=1,GPU=False,
        new_native_EOS_calls=0,new_opacity_queries=0,new_stellar_steps=0,automatic_expansion=False,
        stop='Stop at bound or failed positivity/solver/spectral algebra; record response failures without changing gates. Do not expand duration or matrices automatically.',
        limitations='Classical collisionless Maxwell/RPA with specified FDT, pole residue approximation, neutral atomic mass convention and vacuum photon kinematics. Nonideal/quantum/collisional, bound-electron, actual-temperature absorption, atmosphere and GR are not certified.',
        bindings={p.relative_to(m.old.h.ROOT).as_posix():m.old.h.digest(p) for p in paths}))
    print('PREPARED RECOVERY',flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();plan=json.loads((OUT/'coupling-plan.json').read_text())
    for p,sha in plan['bindings'].items():assert m.old.h.digest(m.old.h.ROOT/p)==sha,p
    signal.alarm(plan['hard_seconds']);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    b=dict(np.load(m.old.OUT/'bank.npz'));duration=float(b['H']/b['c']);Cm=float(b['Cm']);rows=[]
    op=c.Operator(dict(b,velocity_order=128,q_points=129),8,1,24)
    runs=[op.evolve(duration,n) for n in [16,32,64]]
    errors=[m.old.compare(a,z,op) for a,z in zip(runs[:-1],runs[1:])]
    time_order=float(np.log2(errors[0]/errors[1]));middle=runs[1];reference=runs[2]
    row=dict(name='base',coefficients=op.info,time_order=time_order,time_differences_initial=errors,algebra=algebra(op,b),**diagnostics(op,reference))
    rows.append(row);m.old.ex.write(OUT/'base.json',row);print('BASE',row,flush=True)
    np.savez_compressed(OUT/'base.npz',T=reference['T'],E=reference['E'].reshape(-1,8),T32=middle['T'],E32=middle['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    del op;gc.collect();prior=reference;comparisons={}
    for spec in plan['comparisons']:
        name=spec['name'];L=spec['modes']
        op=c.Operator(dict(b,velocity_order=spec['mesh'],q_points=spec['q_grid']),L,1,spec['angles']);r=op.evolve(duration,64)
        E=r['E'].reshape(-1,L);P=prior['E'].reshape(-1,8)
        delta=previous.distance(r['T'],E[:,:8],prior['T'],P,Cm)
        if L>8:delta=float(np.sqrt(delta*delta+np.linalg.norm(E[:,8:])**2/Cm))
        comparisons[name]=delta
        row=dict(name=name,difference_initial=delta,coefficients=op.info,algebra=algebra(op,b),**diagnostics(op,r));rows.append(row)
        m.old.ex.write(OUT/(name+'.json'),row);print('COMPARISON',row,flush=True)
        np.savez_compressed(OUT/(name+'.npz'),T=r['T'],E=E,u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
        if name=='q-grid':best=r
        prior=r;del op;gc.collect()
    op=m.Operator(dict(b,collective_free=True),8,1,24);f32=op.evolve(duration,32)
    f64=np.load(m.OUT/'independent-electron.npz');F=f64['E'].ravel();U=f64['T']
    np.savez_compressed(OUT/'independent-electron-32.npz',T=f32['T'],E=f32['E'].reshape(-1,8),u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
    paired=previous.distance(reference['T']-U,reference['E']-F,middle['T']-f32['T'],middle['E']-f32['E'],Cm)
    effect=previous.distance(best['T'],best['E'],U,F,Cm)
    row=dict(name='independent-electron-32',algebra=algebra(op,b),**diagnostics(op,f32));rows.append(row)
    old=np.load(m.OUT/'q-grid.npz');replaced=previous.distance(best['T'],best['E'].reshape(-1,8),old['T'],old['E'],Cm)
    passed=bool(time_order>=1.8 and errors[-1]<.001 and max(comparisons.values())<1e-6 and paired<1e-6
        and all(r['balance']<1e-9 and r['energy_residual']<1e-9 and r['solver_residual']<1e-11
            and r['entropy_growth']<1e-10 and max(r['algebra']['energy_number_null_relative'])<1e-10 for r in rows))
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,comparisons=comparisons,
        paired_time_difference_initial=paired,collective_effect_initial=effect,
        collective_temperature_effect_K=float(((best['T']-U)/np.sqrt(Cm)).real),
        final_temperature_perturbation_K=float((best['T']/np.sqrt(Cm)).real),
        rejected_discrete_difference_initial=replaced,
        seconds=time.monotonic()-start,peak_memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_discrete_result_preserved=True,original_discrete_passed=False,
        EOS_population_matched_collective_response_evolved=True,bound_electron_kernel_complete=False,
        nonideal_collisional_quantum_error_certified=False,native_scalar_kernel_identified=False,
        actual_temperature_absorption_certified=False,physical_atmosphere_closed=False,
        full_GR_feedback_evolved=False,full_dynamic_charge_solved=False)
    m.old.ex.write(OUT/'result.json',result);print('RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
