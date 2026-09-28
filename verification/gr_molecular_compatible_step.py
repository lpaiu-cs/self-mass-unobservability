"""One budgeted native molecular-EOS step through the validated GR operator.

Counterexample candidate. This transfers the new initial state to the corrected
spatial solver; it does not add scalar forcing or certify long-time convergence.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import time

import numpy as np
import gr_molecular_shell_initial as initial
import gr_cached_material_evolution as cached
import gr_iteration_budget_recovery as repaired
from gr_implicit_coupled_evolution import CachedOnly

m, e, ld = cached.original, initial.e, initial.ld
OUT = e.g.OUT/'gr-molecular-compatible-step'


def save(name, value):
    e.write(OUT/name, value)


def initialize(pool=None):
    star = initial.initialize(CachedOnly() if pool is None else cached.CachedMaterialPool(pool))
    star.__class__ = m.CompatibleStar
    return m.attach_equilibrium(star)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel, value in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == value, rel
    initial.bindings()
    return plan


def prepare():
    assert not OUT.exists()
    assert json.loads((initial.OUT/'restriction.json').read_text())['passed']
    initial.original.worker_init()
    star = initialize()
    zero = np.zeros_like(star.base); z = star.evaluate(zero)
    assert np.all(z['dU'] == 0) and np.all(z['dBX'] == 0)
    flux = np.zeros(star.n+1,dtype=ld)
    template = json.loads((m.BASE/'production-resume271-consoleless/plan.json').read_text())
    budget = m.budget(star,z,flux,flux)[0]; cone = m.parent.cones(z)
    assert m.budget_passed(budget,template) and cone['sampled_cone_inside_light_cone']
    assert m.prior.symbolic()['passed']
    OUT.mkdir()
    files = [Path(__file__),Path(initial.__file__),Path(m.__file__),Path(cached.__file__),
        Path(repaired.__file__),Path(m.method.__file__),Path(m.prior.__file__),
        initial.OUT/'initial-manifest.json',initial.OUT/'initial.npz',
        m.BASE/'production-resume271-consoleless/execution-manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files},
        duration_seconds=str(e.TAU/128),steps=1,cells=star.n,
        finite_conservation_gates=template['finite_conservation_gates'],
        nonlinear_absolute_tolerances=m.ATOL.astype(float).tolist(),maximum_stage_iterations=24,
        method='One initial backward-Euler step with the corrected compatible pressure/gravity operator, full native molecular EOS/opacity and 26 advected species. Reuse the exact material cache and the within-24-iterations donor repair. No old accepted time state or old-model EOS auxiliary is imported.',
        decision='Can the new conservative molecular initial state actually enter the validated spherical GR evolution while satisfying the original native residual, conservation and cone gates? A passing single step closes the provider/initialization handoff, not a time-convergence test.',
        initial_state=dict(budget=budget,cone=cone,zero_conservative_increments=True),
        budget=dict(workers=8,blas_threads=1,gpu=False,hard_timeout_seconds=600,
            estimated_wall_seconds=[120,480],maximum_runs=1,
            estimate='Seven native restriction states measured 0.58 CPU seconds. The stage adds full-grid composition probes and nonlinear updates; their counts are not yet measured. Stop on failure, 24 records or 600 seconds, with no automatic step-size change or longer run.'),
        full_duration_evolution=False,scalar_forcing_included=False,
        physical_EOS_certified=False,time_convergence_verified=False,observational_closure=False))
    print('PREPARED native molecular compatible GR step',cone['maximum_local_rest_characteristic_speed_over_c'],flush=True)


def run():
    plan = bindings(); assert not (OUT/'step-0000.npz').exists()
    began = time.monotonic(); h = ld(plan['duration_seconds']); records=[]
    with ProcessPoolExecutor(max_workers=plan['budget']['workers'],initializer=initial.original.worker_init) as pool:
        star = initialize(pool); delta = np.zeros_like(star.base); z0 = star.evaluate(delta)
        previous = (delta,z0); coefficients = m.weights(h,None)
        zeros = np.zeros(star.n+1,dtype=ld)
        b0,ed0,bd0 = m.budget(star,z0,zeros,zeros)
        m.save_state(OUT/'step-0000.npz',delta,z0,ld(0),zeros,zeros,ed0,bd0)
        def log(row):
            records.append(row)
            with (OUT/'iterations.jsonl').open('a') as stream:
                stream.write(json.dumps(row)+'\n')
            print('MOLECULAR NATIVE GR STAGE',row['iteration'],row['residual_norm'],flush=True)
        try:
            delta,z = repaired.stage(star,previous,previous,h,coefficients,log)
            area=4*np.pi*star.rf**2
            energy,baryons=h*e.C*area*z['fluxes'][1],h*e.C*area*z['fluxes'][0]
            budget,ed,bd=m.budget(star,z,energy,baryons); cone=m.parent.cones(z)
            assert m.budget_passed(budget,plan) and cone['sampled_cone_inside_light_cone'],(budget,cone)
            assert len(records)<=24 and records[-1]['residual_norm']<=1
            m.save_state(OUT/'step-0001.npz',delta,z,h,energy,baryons,ed,bd)
            result=dict(classification='Counterexample candidate',completed=True,passed=True,
                duration_seconds=float(h),cells=star.n,iterations=len(records),
                native_residual_norm=records[-1]['residual_norm'],budget=budget,cone=cone,
                maximum_primitive_changes=np.max(abs(delta[:,:5]),axis=0).astype(float).tolist(),
                native_material_calls=star.pool.evaluations,native_cache_hits=star.pool.hits,
                seconds=time.monotonic()-began,full_duration_evolution=False,
                scalar_forcing_included=False,time_convergence_verified=False,
                physical_EOS_certified=False,observational_closure=False)
            save('result.json',result)
        except Exception as error:
            if hasattr(star,'last_delta'):
                np.savez_compressed(OUT/'failed-iterate.npz',delta=star.last_delta,aux=star.last_state['aux'])
            save('failure.json',dict(classification='Counterexample candidate',completed=False,error=repr(error)))
            raise
    verify()
    print('PASS native molecular compatible GR step',json.dumps(result),flush=True)


def verify():
    plan=bindings(); result=json.loads((OUT/'result.json').read_text())
    assert result['passed'] and not (OUT/'failure.json').exists()
    star=initialize(); old=np.load(OUT/'step-0000.npz'); new=np.load(OUT/'step-0001.npz')
    states=[]
    for saved in [old,new]:
        delta=saved['delta']; y=star.base+delta
        star.material_cache={e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),saved['aux'])}
        state=star.evaluate(delta)
        for key in ['m','mf','a','N','Q','aux','dU']:
            assert np.array_equal(state[key],saved[key]),key
        states.append((delta,state))
    value,_=m.residual(star,states[1][0],states[0],states[0],ld(plan['duration_seconds']),m.weights(ld(plan['duration_seconds']),None))
    norm=float(np.max(abs(value)/m.ATOL)); assert norm<=1
    b,ed,bd=m.budget(star,states[1][1],new['integrated_energy_flux'],new['integrated_baryon_flux'])
    assert m.budget_passed(b,plan)
    assert np.array_equal(ed,new['normalized_energy_defect']) and np.array_equal(bd,new['normalized_baryon_defect'])
    save('replay.json',dict(classification='Counterexample candidate',passed=True,native_residual_norm=norm,exact_saved_state_arrays=True))
    save('manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p)
        for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run','verify'])
    globals()[parser.parse_args().action]()
