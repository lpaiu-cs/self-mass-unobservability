"""Conservative restriction of the existing reference, followed by actual evolution.

Counterexample candidate. The restriction preserves each stored high-order
baryon/energy moment. It does not make the inherited spatial operator high order.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import subprocess
import time

import numpy as np
import gr_coupled_evolution as evolution

ld = evolution.ld
ROOT = evolution.ROOT
OUT = evolution.g.OUT/'gr-conservative-evolution'
write, digest = evolution.write, evolution.digest


def restrict_cell(row):
    index, lr, lt, target_u, x = row
    for iteration in range(16):
        aux = evolution.material((lr, lt, x))
        residual = (ld(aux[2])-target_u)/ld(aux[10])
        assert np.isfinite(residual) and aux[10] > 0, (index, iteration)
        if abs(residual) < 1e-8:
            return lt, aux, float(residual), iteration+1
        lt -= np.clip(residual, ld('-.25'), ld('.25'))
    raise RuntimeError(('Unresolved conservative restriction', index, float(residual)))


def prepare(workers):
    assert not (OUT/'plan.json').exists(), 'Keep every prior attempted plan immutable.'
    OUT.mkdir(exist_ok=True)
    source = evolution.g.OUT/'gr-metric-coupled-subcell/plan.json'
    previous = json.loads(source.read_text())
    state_path = evolution.g.OUT/'initial-state-17-4.npz'
    original = dict(np.load(state_path))
    count = len(original['dm'])
    baryons, energy = np.zeros((2, count), dtype=ld)
    seen = set()
    inputs = [Path(__file__), Path(evolution.__file__), source, state_path,
              evolution.g.OUT/'gr-microphysics/auxiliaries.npz',
              evolution.g.OUT/'gr-transport/diagnostics.npz']
    for rel in sorted(set(previous['cell_sources'].values())):
        path = ROOT/rel
        path = path.with_name(path.stem+'-nodes-16.npz')
        assert digest(path) == previous['bindings'][path.relative_to(ROOT).as_posix()]
        inputs.append(path)
        data = np.load(path)
        cells = data['cells']
        assert not seen.intersection(cells)
        seen.update(cells)
        rho, u = data['eos'][:, :, 0].astype(ld), data['eos'][:, :, 2].astype(ld)
        weights = data['coordinate_weights_cm3'].astype(ld)
        baryons[cells] = np.sum(rho*data['metric_a']*weights, axis=1)
        energy[cells] = np.sum(rho*(data['C_X'].astype(ld)[:, None]*evolution.C**2+u)*weights, axis=1)
    assert seen == set(range(count))
    # Freeze the moments actually present in the reference. Never substitute the
    # requested total stellar mass or rescale quadrature weights to fit it.
    baryons, energy = baryons[::-1], energy[::-1]
    inherited = json.loads((evolution.g.OUT/'gr-coupled-evolution/radial-5735-4/plan.json').read_text())
    for rel, h in inherited['bindings'].items():
        assert digest(ROOT/rel) == h, rel
    for rel, h in previous['runtime_sha256'].items():
        assert digest(Path(rel)) == h, rel
    plan = dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        bindings=dict(inherited['bindings'], **{p.relative_to(ROOT).as_posix(): digest(p) for p in inputs}),
        runtime_sha256=previous['runtime_sha256'],
        reference_nodes=16, cells=count, steps=[8, 16], duration_seconds=float(evolution.TAU),
        restriction='At v=0, set E_i=reference_coordinate_energy_i/V_i. Compute common mass faces and the inherited constant-energy nodal metric a_i from these energies. Set rho_i=B_i/(a_i V_i), then solve native u(rho_i,T_i,X_i)=E_i/rho_i-C_X. Thus the two moments determine the primitive state; no requested global mass enters the inverse.',
        gates=dict(maximum_cell_baryon_relative_defect=1e-14,
                   maximum_cell_energy_relative_defect=1e-13,
                   maximum_native_thermal_residual=1e-8,
                   maximum_original_baryon_relative_defect=1e-11,
                   maximum_original_mass_relative_defect=2e-12),
        evolution='Unmodified existing Star.rhs and SSPRK2 on primitive increments for one proper-tau coordinate interval. No nuclear reactions; advected composition. Retain the original physical Q at t=0. Reconstruct mass/lapse at each stage. Also accumulate every common-face energy flux to measure the discrete Einstein time-constraint defect.',
        boundary='Same regular-centre and computational outer-face extension as gr_coupled_evolution. No atmosphere/exterior match.',
        limitations='Conservative restriction is a change of spatial representation, not high-order spatial evolution, an isentropic remap, a native-EOS certificate or a rigorous time/space error bound.',
        symbolic=evolution.symbolic())
    write(OUT/'plan.json', plan)
    with ProcessPoolExecutor(max_workers=workers, initializer=evolution.worker_init) as pool:
        star = evolution.Star(count, pool)
        mean_energy = energy/star.volume
        faces = np.r_[ld(0), np.cumsum(evolution.GRAV*energy)]
        m = faces[:-1]+evolution.GRAV*mean_energy*4*np.pi/3*(star.r**3-star.rf[:-1]**3)
        a = 1/np.sqrt(1-2*m/star.r)
        rho = baryons/(a*star.volume)
        lr = np.log(rho)
        # Use the identical represented density and composition rest energy as
        # Star.state, retaining long-double subtraction before the native solve.
        rest = (star.base[:, 4:]/evolution.g.c.A)@evolution.g.c.W*evolution.C**2
        target_u = mean_energy/np.exp(lr)-rest
        rows = list(pool.map(restrict_cell, zip(star.indices, lr, star.base[:, 1],
                                                target_u, star.base[:, 4:]), chunksize=4))
        base = star.base.copy()
        base[:, 0], base[:, 1] = lr, [row[0] for row in rows]
        aux = np.asarray([row[1] for row in rows], dtype=ld)
        star.base = base
        star.material_cache = {evolution.material_key(row): value for row, value in
                               zip(zip(base[:, 0], base[:, 1], base[:, 4:]), aux)}
        z = star.state(base)
        actual_b = z['a']*z['D']*star.volume
        actual_e = z['E']*star.volume
        report = dict(classification='Counterexample candidate', cells=count,
            maximum_cell_baryon_relative_defect=float(np.max(abs(actual_b/baryons-1))),
            maximum_cell_energy_relative_defect=float(np.max(abs(actual_e/energy-1))),
            maximum_native_thermal_residual=max(abs(row[2]) for row in rows),
            maximum_original_baryon_relative_defect=float(abs(actual_b.sum()/np.sum(original['dm'].astype(ld))-1)),
            maximum_original_mass_relative_defect=float(abs(evolution.GRAV*actual_e.sum()/star.original_mass-1)),
            maximum_initial_primitive_changes=np.max(abs(base[:, :2]-np.column_stack([
                original['lnd'][star.indices], original['lnT'][star.indices]])), axis=0).astype(float).tolist(),
            maximum_native_calls_per_cell=max(row[3] for row in rows))
        report['passed'] = all(report[k] < limit for k, limit in plan['gates'].items())
        write(OUT/'restriction.json', report)
        assert report['passed'], report
        np.savez_compressed(OUT/'initial.npz', base=base, aux=aux, qscale=star.qscale,
                            radius=star.r, volume=star.volume, indices=star.indices,
                            baryons=baryons, energy=energy, restriction_residual=np.array([r[2] for r in rows]))
    write(OUT/'initial-manifest.json', dict(sha256={p.name: digest(p) for p in
          [OUT/'initial.npz', OUT/'restriction.json', OUT/'plan.json']}))
    print('CONSERVATIVE INITIAL STATE', json.dumps(report), flush=True)


def verify_initial():
    for rel, h in json.loads((OUT/'initial-manifest.json').read_text())['sha256'].items():
        assert digest(OUT/rel) == h, rel
    plan = json.loads((OUT/'plan.json').read_text())
    for rel, h in plan['bindings'].items():
        assert digest(ROOT/rel) == h, rel
    for rel, h in plan['runtime_sha256'].items():
        assert digest(Path(rel)) == h, rel
    assert json.loads((OUT/'restriction.json').read_text())['passed']
    assert evolution.symbolic()['passed']
    return plan


def run(steps, workers):
    plan = verify_initial()
    assert steps in plan['steps']
    folder = OUT/f'path-{steps}'
    assert not folder.exists()
    folder.mkdir()
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, initializer=evolution.worker_init) as pool:
        star = evolution.Star(plan['cells'], pool)
        initial = np.load(OUT/'initial.npz')
        assert np.array_equal(star.r, initial['radius'])
        assert np.array_equal(star.volume, initial['volume'])
        star.base = initial['base'].copy()
        star.qscale = initial['qscale'].copy()
        star.material_cache = {evolution.material_key(row): value for row, value in zip(
            zip(star.base[:, 0], star.base[:, 1], star.base[:, 4:]), initial['aux'])}
        delta = np.zeros_like(star.base)
        exchange = np.zeros(2, dtype=ld)
        mass_flux_integral = np.zeros(star.n+1, dtype=ld)
        dt = evolution.TAU/steps
        history = []
        for step in range(steps+1):
            k1, z = star.rhs(delta)
            totals = np.array([np.sum(z['a']*z['D']*star.volume), np.sum(z['E']*star.volume)])
            if step == 0:
                baseline = totals.copy()
                mf0 = z['mf'].copy()
                thermal = np.cumsum(z['rho']*z['aux'][:, 10]*star.volume)
            mass_defect = z['mf']-mf0+mass_flux_integral
            cell_defect = np.diff(mass_defect)/(evolution.GRAV*np.diff(np.r_[ld(0), thermal]))
            row = dict(step=step, time_seconds=float(step*dt),
                maximum_changes=np.max(abs(delta[:, :4]), axis=0).astype(float).tolist(),
                relative_budget_defect=np.asarray((totals-baseline-exchange)/baseline, float).tolist(),
                maximum_face_mass_defect_over_total=float(np.max(abs(mass_defect))/mf0[-1]),
                maximum_cell_energy_defect_over_initial_heat_capacity=float(np.max(abs(cell_defect))),
                maximum_time_matrix_residual=z['time_matrix_residual'], elapsed_seconds=time.monotonic()-started)
            history.append(row)
            np.savez_compressed(folder/f'step-{step:04d}.npz', delta=delta, m=z['m'], mf=z['mf'],
                a=z['a'], N=z['N'], Q=z['Q'], totals=totals, boundary_exchange=exchange,
                integrated_face_mass_flux=mass_flux_integral, mass_time_constraint_defect=mass_defect)
            write(folder/'progress.json', dict(classification='Counterexample candidate', rows=history))
            print('CONSERVATIVE COUPLED STEP', steps, json.dumps(row), flush=True)
            if step == steps:
                break
            k2, z2 = star.rhs(delta+dt*k1)
            delta += dt*(k1+k2)/2
            exchange += dt*(z['boundary_rates']+z2['boundary_rates'])/2
            flux = sum(star.faces(q['N']*q['S']/q['a'], odd=True) for q in [z, z2])/2
            mass_flux_integral += dt*evolution.GRAV*evolution.C*4*np.pi*star.rf**2*flux
        assert np.all(np.max(abs(delta[:, :4]), axis=0) > 0)
        write(folder/'result.json', dict(classification='Counterexample candidate', completed=True,
            actual_nonlinear_finite_time_path=True, cells=star.n, steps=steps, last=history[-1],
            conservative_initial_restriction=True, high_order_spatial_evolution=False,
            exact_discrete_energy_conservation=False, physical_exterior_match=False,
            physical_EOS_certified=False, observational_closure=False))
    write(folder/'manifest.json', dict(initial_manifest_sha256=digest(OUT/'initial-manifest.json'),
          sha256={p.name: digest(p) for p in folder.iterdir() if p.is_file()}))
    print('CONSERVATIVE COUPLED PATH COMPLETE', steps, flush=True)


def compare():
    plan = verify_initial()
    assert not (OUT/'comparison.json').exists()
    arrays, records = [], []
    for steps in plan['steps']:
        folder = OUT/f'path-{steps}'
        manifest = json.loads((folder/'manifest.json').read_text())
        assert manifest['initial_manifest_sha256'] == digest(OUT/'initial-manifest.json')
        for rel, h in manifest['sha256'].items():
            assert digest(folder/rel) == h, rel
        result = json.loads((folder/'result.json').read_text())
        assert result['completed'] and result['cells'] == plan['cells']
        arrays.append(np.load(folder/f'step-{steps:04d}.npz')['delta'])
        records.append(result)
    for key in ['m', 'mf', 'a', 'N', 'Q', 'totals']:
        assert np.array_equal(*(np.load(OUT/f'path-{n}/step-0000.npz')[key] for n in plan['steps']))
    error = np.max(abs(arrays[1][:, :4]-arrays[0][:, :4]), axis=0)
    write(OUT/'comparison.json', dict(classification='Counterexample candidate', completed=True,
        same_initial_state=True, endpoint_maximum_8_16_differences=error.astype(float).tolist(),
        runs=records, restriction=json.loads((OUT/'restriction.json').read_text()),
        observed_order=None, rigorous_time_error_bound=False,
        note='Two resolutions measure sensitivity to dt; they do not determine an observed convergence order.'))
    write(OUT/'comparison-manifest.json', dict(sha256={str(p.relative_to(ROOT)): digest(p) for p in
        [OUT/'initial-manifest.json', OUT/'comparison.json']+[OUT/f'path-{n}/manifest.json' for n in plan['steps']]}))
    print('CONSERVATIVE COUPLED COMPARISON', error, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run', 'verify_initial', 'compare'])
    parser.add_argument('--steps', type=int, default=16)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.workers)
    elif args.command == 'run': run(args.steps, args.workers)
    else: globals()[args.command]()
