"""Counterexample candidate: reuse accepted GR states for a scalar readout.

This is a frozen-snapshot elliptic readout, not scalar evolution or an orbital
transfer function. No EOS call or accepted hydrodynamic state is recomputed.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import sympy as sp

import common_eos as eos
import gr_resume271_consoleless as accepted
from gr_scalar_reference import riccati

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/direct-eos-gr33/gr-driven-readout'


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def symbolic():
    # The trace is frame invariant even with radial velocity and heat flux.
    v, eps, p, q = sp.symbols('v eps p q')
    E = (eps + p*v*v + 2*q*v)/(1-v*v)
    R = (eps*v*v + p + 2*q*v)/(1-v*v)
    assert sp.cancel(-E + R + 2*p - (-eps + 3*p)) == 0
    # DEF parity: matter forcing is quadratic at zero scalar background.
    b, phi0, dphi, a, s, L, h, f, chi = sp.symbols('b phi0 dphi a s L h f chi')
    force = b*(phi0 + a*dphi)**2/2
    assert sp.diff(force, a).subs(a, 0) == b*phi0*dphi
    assert sp.diff(force, a).subs({a: 0, phi0: 0}) == 0
    du = phi0*f*dphi/(s-L)
    dq = chi*dphi + phi0*h*du
    assert sp.cancel(dq/dphi - chi - phi0**2*h*f/(s-L)) == 0
    # Closed outer-face energy flux gives constant mass in variable-step BDF2.
    r, M, c1, c2 = sp.symbols('r M c1 c2')
    c0 = (1+2*r)/(1+r)
    assert sp.cancel(c0 - (1+r) + r*r/(1+r)) == 0
    return dict(classification='Proven', passed=True,
        trace='-E+R+2P=-epsilon+3P, including radial heat flux.',
        parity='For A(phi)=exp(beta*phi^2/2), the matter source beta*phi*T*d_mu(phi) is quadratic about phi=0. Linear scalar perturbations do not excite the fluid/thermal block on the unscalarized branch.',
        mixed_response='For a smooth nonzero-background expansion with bounded stationary resolvent: deltaU=phi0*(sI-L)^(-1)*f*delta_phi+..., deltaq=chi*delta_phi+phi0*h*deltaU+..., so the fluid-mediated term in deltaq/delta_phi starts at phi0^2. This is power counting, not a magnitude bound near a pole.',
        mass='A shared-face conservative energy update with zero central/outer energy flux and constant starting mass preserves total mass by induction, up to the accepted residual.',
        limitations='No claim of stability, small resolvent norm, physical EOS, stationary evolving background, nonzero scalar backreaction, or observational closure.')


def prepare():
    assert not OUT.exists(), 'Preserve prior outcomes.'
    plan = accepted.bindings()
    verdict = json.loads((accepted.OUT / 'time-refinement.json').read_text())
    assert verdict['passed']
    sources = [Path(__file__), Path(eos.__file__), Path(riccati.__code__.co_filename),
               accepted.OUT/'plan.json', accepted.OUT/'time-refinement.json',
               accepted.OUT/'time-refinement-manifest.json']
    for refinement, steps in [(1, 70), (2, 140), (4, 280)]:
        folder = accepted.OUT / f'path-{refinement}'
        manifest = json.loads((folder/'manifest.json').read_text())
        for step in [0, steps]:
            path = folder / f'step-{step:04}.npz'
            assert accepted.e.digest(path) == manifest['sha256'][path.name]
            sources.append(path)
    OUT.mkdir()
    save('plan.json', dict(classification='Counterexample candidate',
        checkpoint='a9b8878c', bindings={p.relative_to(ROOT).as_posix():accepted.e.digest(p) for p in sources},
        beta=-4, tolerances=[1e-10, 2e-12], quadrature_nodes=[6, 12],
        absolute_alpha_readout_agreement=1e-9,
        change_over_discretization_difference_gate=10,
        readout='alpha/phi_infinity=-chi/M, exact Schwarzschild logarithmic scalar tail on each frozen snapshot. The isotropic comoving trace is rho*(rest+u)-3P; do not use Eulerian E-3P.',
        drive='Subsequent candidate: delta_phi_ext(t) on a specified small nonzero DEF background phi0, with consistent matter source beta*phi*T*d_mu(phi). This run applies no drive; it measures whether the already computed GR trajectory changes the static scalar readout.',
        comparator='For this run: constant initial susceptibility. For a driven experiment: co-fit the registered instantaneous response and justified derivative/nuisance directions; a transient difference alone is not a nonadiabatic observable.',
        decision='Use existing endpoints to assess a resolvable charge-readout change. Closed-GR mass drift is a negative control, never the candidate signal. Do not launch molecular full-duration evolution on this evidence alone.',
        resource_budget=dict(cpu_workers=1, gpu=False, maximum_wall_seconds=300,
            maximum_model_integrations=8, maximum_Riccati_integrations=6, native_EOS_calls=0,
            memory_estimate='Below 1 GiB for saved arrays and two scalar interpolants; measured peak recorded by /usr/bin/time.',
            wall_estimate='Unmeasured scalar-readout runtime: bounded pilot, not a long-run ETA. Hard external timeout 300 seconds.',
            reuse='Saved 70/140/280 endpoints, initial state and native auxiliaries.',
            stop='Stop on input mismatch, failed numerical control, or wall-time cap; preserve partial result; no automatic tolerance/refinement enlargement.'),
        limitations=['Frozen nonstationary snapshots do not define a causal time response.',
            'Finite-pressure wall is not a physical scalar-tensor exterior/atmosphere match.',
            'Two ODE tolerances and quadratures are finite comparisons, not rigorous error bounds.',
            'The new molecular EOS, finite scalar backreaction, orbital transfer and data inference are absent.']))
    save('symbolic.json', symbolic())
    print('PREPARED bounded saved-state scalar readout', flush=True)


def snapshot(star, refinement, step):
    path = accepted.OUT / f'path-{refinement}' / f'step-{step:04}.npz'
    z = dict(np.load(path))
    y = star.base + z['delta']
    data = dict(radius_faces_m=star.rf[::-1]/100,
        mass_faces_geom=z['mf'][::-1]/100, r_mid_m=star.r[::-1]/100,
        m_mid_geom=z['m'][::-1]/100, lnd=y[::-1, 0], X=y[::-1, 5:],
        logP=np.log(z['aux'][::-1, 1]), u_W=z['aux'][::-1, 2], nu=np.log(z['N'][::-1]))
    return data, z


def difference(base, end, count):
    assert np.array_equal(base['x'], end['x'])
    # Same PCHIP mesh: subtract coefficients before evaluating to reduce loss.
    x = base['x']; width = np.diff(x)
    nodes, weights = np.polynomial.legendre.leggauss(count)
    offset = width[:, None]*(nodes+1)/2
    xx = x[:-1, None] + offset
    def coefficient_difference(key):
        coefficients = end[key].c - base[key].c
        return sum(coefficients[k, :, None]*offset**(3-k) for k in range(4))
    f0, g0 = base['field'](xx.ravel())
    f1, g1 = end['field'](xx.ravel())
    values = f0*f1*coefficient_difference('v').ravel() + g0*g1*coefficient_difference('p').ravel()
    integral = np.sum(width/2 * ((xx*xx*values.reshape(xx.shape)) @ weights))
    t = (nodes+1)/2
    exterior = np.dot(weights, t/((1-2*base['mu']*t)*(1-2*end['mu']*t)))/2
    da = -integral + 2*(end['mu']-base['mu'])*base['a']*end['a']*exterior
    normalization = ((end['R']-base['R'])-base['R']/base['M']*(end['M']-base['M']))/end['M']
    return float(-(end['R']/end['M']*da + base['a']*normalization))


def run():
    began = time.monotonic()
    plan = json.loads((OUT/'plan.json').read_text())
    assert not (OUT/'result.json').exists()
    for rel, digest in plan['bindings'].items():
        assert accepted.e.digest(ROOT/rel) == digest, rel
    assert symbolic()['passed']
    star = accepted.context['initialize'](None)  # cached-only native provider
    base_data, base_z = snapshot(star, 1, 0)
    rows = []
    for tolerance in plan['tolerances']:
        base = eos.scalar_model('accepted-initial', base_data, tolerance)
        assert difference(base, base, 12) == 0
        for refinement, step in [(1, 70), (2, 140), (4, 280)]:
            assert time.monotonic()-began < plan['resource_budget']['maximum_wall_seconds']
            data, z = snapshot(star, refinement, step)
            end = eos.scalar_model('accepted-endpoint', data, tolerance)
            control = riccati(end, tolerance)
            w = control.y[0, -1]
            independent = end['R']/end['M']*w/(1-w*np.log1p(-2*end['mu'])/(2*end['mu']))
            alpha0, alpha1 = -base['chi']/base['M'], -end['chi']/end['M']
            changes = [difference(base, end, n) for n in plan['quadrature_nodes']]
            row = dict(refinement=refinement, steps=step, tolerance=tolerance,
                time_seconds=float(z['time_seconds']), initial_alpha_over_phi=float(alpha0),
                end_alpha_over_phi=float(alpha1), delta_alpha_direct=float(alpha1-alpha0),
                delta_alpha_wronskian=changes[-1], quadrature_difference=abs(changes[1]-changes[0]),
                direct_identity_difference=abs(float(alpha1-alpha0)-changes[-1]),
                independent_equation_difference=abs(float(alpha1)-float(independent)),
                conditional_weak_companion_force_coefficient=float(-4*changes[-1]),
                mass_relative_change=float((z['mf'][-1]-base_z['mf'][-1])/base_z['mf'][-1]),
                outer_integrated_energy_flux=float(z['integrated_energy_flux'][-1]),
                elapsed_seconds=time.monotonic()-began)
            rows.append(row)
            save('partial.json', dict(classification='Counterexample candidate', rows=rows))
            print('READOUT', json.dumps(row), flush=True)
    fine = [row for row in rows if row['tolerance'] == plan['tolerances'][-1]]
    delta = [row['delta_alpha_wronskian'] for row in fine]
    numerical = max(max(row[k] for k in ['quadrature_difference', 'direct_identity_difference',
        'independent_equation_difference']) for row in rows)
    tolerance_difference = max(abs(rows[k]['delta_alpha_wronskian']-rows[k+3]['delta_alpha_wronskian']) for k in range(3))
    resolution = max(numerical, tolerance_difference, abs(delta[2]-delta[1]), np.finfo(float).eps)
    controls = numerical <= plan['absolute_alpha_readout_agreement'] and tolerance_difference <= plan['absolute_alpha_readout_agreement']
    result = dict(classification='Counterexample candidate', completed=True, rows=rows,
        maximum_numerical_comparison_difference=numerical,
        maximum_tolerance_difference=tolerance_difference,
        readout_refinement_differences=[abs(delta[1]-delta[0]), abs(delta[2]-delta[1])],
        diagnostic_change_over_difference=abs(delta[-1])/resolution,
        numerical_controls_passed=controls,
        change_resolved_against_declared_finite_comparisons=controls and abs(delta[-1]) > plan['change_over_discretization_difference_gate']*resolution,
        maximum_mass_relative_change=max(abs(row['mass_relative_change']) for row in rows),
        outer_energy_flux_exactly_zero=all(row['outer_integrated_energy_flux'] == 0 for row in rows),
        elapsed_seconds=time.monotonic()-began,
        physical_EOS_certified=False, causal_scalar_evolution=False,
        driven_thermal_response=False, rigorous_error_bound=False, observational_closure=False)
    save('result.json', result)
    save('manifest.json', dict(sha256={p.relative_to(ROOT).as_posix():accepted.e.digest(p)
        for p in OUT.iterdir() if p.is_file() and p.name != 'manifest.json'}))
    assert controls, result
    print('COMPLETE saved-state readout', json.dumps({k:v for k,v in result.items() if k != 'rows'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run', 'symbolic'])
    args = parser.parse_args()
    if args.command == 'symbolic':
        print(json.dumps(symbolic(), indent=2))
    else:
        globals()[args.command]()
