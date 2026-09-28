"""Counterexample candidate: evolve the second already declared constant heat time.

Start from the same conservative initial data, at the same coordinate time nodes.
Retain the full entropy transport term and reject an acausal converged stage.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path
import argparse
import json
import subprocess
import time

import numpy as np

import gr_implicit_coupled_anderson as solver
import gr_implicit_coupled_analysis as analysis

old, e, ld = solver.old, solver.e, solver.ld
OUT = e.g.OUT/'gr-implicit-second-heat-time'
TIMES = e.g.OUT/'gr-heat-entropy-closure/constant-times.json'
FIRST_TAU = ld('0.0004113000088929766')


def second_tau():
    record = json.loads(TIMES.read_text())
    values = list(map(Fraction, record['exact_proper_seconds']))
    assert values[1]/values[0] == Fraction(3, 2)
    assert record['physical_calibration'] is False
    return ld(values[1].numerator)/ld(values[1].denominator)


def time_nodes(plan, refinement):
    edges = np.array(plan['coordinate_edges_seconds'], dtype=ld)
    return np.r_[ld(0), np.concatenate([np.linspace(a, b, refinement+1, dtype=ld)[1:]
                                      for a, b in zip(edges[:-1], edges[1:])])]


def check():
    solver.check()
    previous = json.loads((old.OUT/'plan.json').read_text())
    edges = np.array(previous['time_edges_tau'], dtype=ld)*FIRST_TAU
    plan = dict(coordinate_edges_seconds=list(map(str, edges)))
    for refinement in [1, 2, 4]:
        nodes = time_nodes(plan, refinement)
        assert np.array_equal(nodes[::refinement], edges)
        assert float(nodes[-1]) == previous['duration_seconds']
    e.TAU = second_tau()
    assert e.TAU > FIRST_TAU
    star = old.initialize(None)
    rates, state = star.rhs(np.zeros_like(star.base))
    assert np.isfinite(rates).all()
    cone = analysis.cones(state)
    assert cone['sampled_cone_inside_light_cone']
    print('PASS second constant time, identical coordinate nodes and initial cone', cone, flush=True)


def prepare():
    assert not OUT.exists(), 'Preserve every prior candidate and verdict.'
    check()
    plan = json.loads((solver.OUT/'plan.json').read_text())
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    plan.update(checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        candidate='Second predeclared spatially and temporally constant proper heat time; restart at t=0.',
        proper_tau_seconds=str(second_tau()),
        coordinate_edges_seconds=list(map(str, np.array(plan['time_edges_tau'], dtype=ld)*FIRST_TAU)),
        duration_tau=None,
        parent_plan_sha256=e.digest(solver.OUT/'plan.json'),
        cone_gate='Every converged SDIRK stage and every stored endpoint must have real local rest-frame characteristics strictly inside the light cone. Save an inadmissible stage and stop; no coefficient clipping or silent time adaptation.',
        physical_heat_time_calibrated=False)
    # The old dimensionless edges would refer to a different physical time.
    del plan['time_edges_tau']
    for path in [Path(__file__), Path(analysis.__file__), TIMES, solver.OUT/'plan.json']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    print('PREPARED SECOND HEAT TIME', plan['proper_tau_seconds'], flush=True)


def save_state(path, delta, z, t, exchange, mass_flux, defect):
    np.savez_compressed(path, delta=delta, **{k:z[k] for k in ['m', 'mf', 'a', 'N', 'Q', 'aux']},
        time_seconds=t, boundary_exchange=exchange, integrated_face_mass_flux=mass_flux,
        normalized_cell_energy_defect=defect)


def run(refinement, workers):
    plan = json.loads((OUT/'plan.json').read_text())
    assert refinement in plan['refinements'] and workers > 0
    for rel, digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    for rel, digest in plan['runtime_sha256'].items():
        assert e.digest(Path(rel)) == digest, rel
    e.TAU = second_tau()
    assert str(e.TAU) == plan['proper_tau_seconds']
    times = time_nodes(plan, refinement)
    folder = OUT/f'path-{refinement}'
    assert not folder.exists()
    folder.mkdir()
    began = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, initializer=e.worker_init) as pool:
        star = old.initialize(pool)
        delta = np.zeros_like(star.base)
        exchange, mass_flux = np.zeros(2, dtype=ld), np.zeros(star.n+1, dtype=ld)
        for step, t in enumerate(times):
            _, z = star.rhs(delta)
            defect, budget = old.energy_budget(star, delta, z, mass_flux)
            cone = analysis.cones(z)
            save_state(folder/f'step-{step:04d}.npz', delta, z, t, exchange, mass_flux, defect)
            assert cone['sampled_cone_inside_light_cone'], cone
            row = dict(classification='Counterexample candidate', step=step, time_seconds=float(t),
                maximum_changes=np.max(abs(delta[:, :4]), axis=0).astype(float).tolist(),
                elapsed_seconds=time.monotonic()-began, **budget, **cone)
            e.write(folder/'progress.json', row)
            print('SECOND HEAT TIME STEP', refinement, json.dumps(row), flush=True)
            if step == len(times)-1:
                break
            dt = times[step+1]-t

            def log(record):
                with (folder/'iterations.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(step=step+1, stage=stage_number, **record))+'\n')
                print('SECOND TIME NATIVE STAGE', refinement, step+1, stage_number,
                      record['iteration'], record['residual_norm'], flush=True)

            def stage(base, seed):
                value, rates, state = solver.stage(star, base, seed, old.GAMMA*dt, log)
                cone = analysis.cones(state)
                with (folder/'stage-cones.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(step=step+1, stage=stage_number, **cone))+'\n')
                if not cone['sampled_cone_inside_light_cone']:
                    np.savez_compressed(folder/'inadmissible-stage.npz', delta=value, aux=state['aux'],
                        time_seconds=t+(old.GAMMA if stage_number == 1 else 1)*dt)
                    raise RuntimeError(('Converged stage outside declared causal domain', cone))
                return value, rates, state

            try:
                stage_number = 1
                first, f1, z1 = stage(delta, delta)
                stage_number = 2
                delta, _, z2 = stage(delta+(1-old.GAMMA)*dt*f1, first)
            except Exception as error:
                e.write(folder/'failure.json', dict(classification='Counterexample candidate', completed=False,
                    step=step+1, stage=stage_number, reason=repr(error)))
                raise
            exchange += dt*((1-old.GAMMA)*z1['boundary_rates']+old.GAMMA*z2['boundary_rates'])
            flux = ((1-old.GAMMA)*star.faces(z1['N']*z1['S']/z1['a'], odd=True)
                    +old.GAMMA*star.faces(z2['N']*z2['S']/z2['a'], odd=True))
            mass_flux += dt*e.GRAV*e.C*4*np.pi*star.rf**2*flux
        e.write(folder/'result.json', dict(classification='Counterexample candidate', completed=True,
            steps=len(times)-1, duration_seconds=float(times[-1]), actual_native_nonlinear_path=True,
            physical_EOS_certified=False, physical_exterior_match=False, observational_closure=False))
    e.write(folder/'manifest.json', dict(plan_sha256=e.digest(OUT/'plan.json'),
        sha256={p.name:e.digest(p) for p in folder.iterdir() if p.is_file()}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check', 'prepare', 'run'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement, args.workers)
    else:
        globals()[args.command]()
