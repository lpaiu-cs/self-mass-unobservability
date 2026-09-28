"""Counterexample candidate: compatible pressure work and hydrostatic balance.

This is a new spatial discretization. Historical GR sources and verdicts are
immutable. Only the equilibrium momentum quadrature is balanced; physical
heat transport, energy/baryon fluxes and all native 31 equations remain active.
"""
import os
import gr_conservative_composition_tangent as method
from gr_conservative_composition_tangent import *

BASE = e.g.OUT/'gr-compatible-equilibrium-v2'
SOURCE = e.g.OUT/'gr-conservative-composition-recovery-20260914'
OUT = BASE/'pilot'
REFINEMENT, PREFIX = None, 0


def pressure_faces(star, value):
    """Dual of the original velocity interpolation, with even pressure parity."""
    t = (star.rf[1:-1]-star.r[:-1])/np.diff(star.r)
    return np.r_[value[0], value[:-1]+(1-t)*(value[1:]-value[:-1]), value[-1]]


def acoustic_stress(star, z):
    beta = (z['P']/z['rho']-z['aux'][:,9])/z['aux'][:,10]
    sound2 = z['P']/z['w']*(z['aux'][:,5]+z['aux'][:,6]*beta)
    assert np.all((sound2 > 0)&(sound2 < 1)), 'Invalid native acoustic impedance'
    impedance = z['N']*z['w']*np.sqrt(sound2)
    # The fixed one-half acoustic Riemann impedance damps velocity jumps.
    # It vanishes at hydrostatic equilibrium and under smooth mesh refinement.
    return np.r_[ld(0), -.5*np.maximum(impedance[:-1],impedance[1:])*np.diff(z['v']), ld(0)]


def fluxes(star, z, donor=None):
    velocity = star.faces(z['N']*z['v']/z['a'], odd=True)
    direction = velocity if donor is None else donor
    fB = velocity*prior.upwind(z['B'], direction)
    fE = velocity*prior.upwind(z['E'], direction)+star.faces(
        z['N']/z['a']*(z['P']*z['v']+z['Q']), odd=True)
    fS = velocity*prior.upwind(z['AS'], direction)+pressure_faces(star, z['N']*z['P'])
    fS += acoustic_stress(star, z)
    fS += star.faces(z['N']*z['Q']*z['v'])
    fS[-1] = z['N'][-1]*z['R'][-1]
    return fB, fE, fS


def finish_fluxes(star, z, donor=None):
    z['fluxes'] = fluxes(star, z, donor)
    mass_face_rate = -e.GRAV*e.C*4*np.pi*star.rf**2*z['fluxes'][1]
    mdot = mass_face_rate[:-1]+star.fraction*np.diff(mass_face_rate)
    z['at'] = z['a']**3*mdot/star.r


class CompatibleStar(prior.ConservativeStar):
    def finish_moments(self, delta, z, dE):
        prior.ConservativeStar.finish_moments(self, delta, z, dE)
        finish_fluxes(self, z)


def attach_equilibrium(star):
    # The reference hydrostatic projection has v=Q=0. Do not subtract the
    # actual thermal drive or recompute this reference at evolving states.
    ref = star.reference
    star.equilibrium_pressure_flux = pressure_faces(star, ref['N']*ref['P'])
    star.equilibrium_gravity = ref['N']*ref['nur']*ref['E']
    star.equilibrium_pressure = ref['N']*ref['P']
    return star


def initialize(pool):
    star = method.initialize(pool)
    star.__class__ = CompatibleStar
    return attach_equilibrium(star)


def residual(star, delta, previous, older, h, coefficients):
    value, z = prior.residual(star, delta, previous, older, h, coefficients)
    c0, c1, c2 = coefficients
    momentum_increment = c0*z['dU'][:, 2]+c1*previous[1]['dU'][:, 2]+c2*older[1]['dU'][:, 2]
    force_difference = -z['at']*z['S']-e.C*(z['N']*z['nur']*z['E']-star.equilibrium_gravity)
    force_difference += e.C*(z['N']*z['P']-star.equilibrium_pressure)*star.area_difference_over_volume
    value[:, 2] = (momentum_increment+h*e.C*star.divergence(
        z['fluxes'][2]-star.equilibrium_pressure_flux)-h*force_difference)/star.momentum_scale
    return value, z


class CompatibleTangent(method.CompositionTangent):
    def finish_moments(self, delta, z, dE):
        prior.ConservativeStar.finish_moments(self, delta, z, dE)
        anchor = self.linearization
        donor = self.faces(anchor['N']*anchor['v']/anchor['a'], odd=True)
        finish_fluxes(self, z, donor)


def tangent(star, delta, z):
    model = method.tangent(star, delta, z)
    model.__class__ = CompatibleTangent
    return model


# Reuse the actual native solver and all conservation gates. Both its physical
# residual and its iteration matrix route through the corrected spatial code.
jacobian = FunctionType(prior.jacobian.__code__, globals())
stage = FunctionType(method.stage.__code__, globals())
_bindings = FunctionType(prior.bindings.__code__, globals())
_run = FunctionType(prior.run.__code__, globals())
completed = FunctionType(prior.completed.__code__, globals())
compare = FunctionType(prior.compare.__code__, globals())


def bindings():
    plan = _bindings()
    assert plan['imported_prefix_steps'] == 0
    assert plan['nonlinear_absolute_tolerances'] == ATOL.astype(float).tolist()
    assert plan['reference_balance_rows'] == [2]
    return plan


def prepare(kind):
    global OUT
    assert kind in ['pilot', 'production']
    OUT = BASE/kind
    assert not OUT.exists(), 'Preserve every attempted new path.'
    design = json.loads((BASE/'design.json').read_text())
    check = json.loads((BASE/'preflight.json').read_text())
    for relative, digest in json.loads((BASE/'preflight-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/relative) == digest, relative
    assert check['passed'] and check['implementation_sha256'] == e.digest(Path(__file__))
    if kind == 'production':
        for relative, digest in json.loads((BASE/'pilot/time-refinement-manifest.json').read_text())['sha256'].items():
            assert e.digest(e.ROOT/relative) == digest, relative
        pilot = json.loads((BASE/'pilot/time-refinement.json').read_text())
        assert pilot['passed'], 'Full-grid native pilot must pass before production.'
    template = json.loads((SOURCE/'plan.json').read_text())
    multiplier = check['selected_time_multiplier']
    edges = parent.wall.prior.time_nodes(template, multiplier)
    if kind == 'pilot':
        edges = np.linspace(ld(0), 4*e.TAU, 9, dtype=ld)
    files = [Path(__file__), BASE/'design.json', BASE/'preflight.json', BASE/'preflight-manifest.json',
             SOURCE/'plan.json', SOURCE/'time-refinement.json']
    plan = dict(template, checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        classification='Counterexample candidate', spatial_discretization=design['spatial_discretization'],
        reference_balance_rows=[2], equilibrium_scope=design['equilibrium_scope'],
        imported_prefix_steps=0, imported_prefix='No old accepted time states. Start each path from the unchanged conserved initial projection.',
        coordinate_edges_seconds=[str(t) for t in edges], duration_seconds=float(edges[-1]),
        refinements=[1, 2, 4], base_time_multiplier=multiplier,
        phase=kind, predecessor_failed_verdict_sha256=e.digest(SOURCE/'time-refinement.json'),
        continuum_equations_intended_unchanged=True, spatial_equations_modified=True,
        reference_equilibrium_error_certified=False, physical_EOS_certified=False)
    for key in ['symbolic','implementation_check','operator_check','inherited_numerical_limits']:
        plan.pop(key, None)  # Old tests remain in the bound predecessor plan.
    plan.update(candidate='Compatible pressure interpolation, fixed acoustic impedance stress, and fixed hydrostatic momentum balance with the original native EOS and 26-species evolution.',
        method='Variable-step BDF2 with initial backward Euler on the predeclared new time grid. Shared original upwind baryon, energy and isotope fluxes; dual pressure faces plus acoustic velocity-jump stress. Same finite-volume metric-rate reconstruction and full two-carrier heat equations.',
        limitations='Conditional LTE two-carrier physics and reflecting wall. First-order acoustic viscosity and material advection. No full nonlinear entropy theorem, physical EOS, atmosphere, continuum-error bound, nuclear reactions or observational closure.',
        nonlinear_operator=template['nonlinear_operator'].replace('physical model, time grid and scientific gates.', 'native physical laws and acceptance tolerances. The new momentum discretization is used in both the native residual and iteration matrix.'),
        cone_gate='Every converged BDF step and stored endpoint must satisfy the original local rest-frame light-cone gate.',
        previous_paths='The bound original 39/78/156-step experiment failed time convergence; both stopped v1 acoustic preflights are preserved separately.',
        parent_plan_sha256=e.digest(SOURCE/'plan.json'))
    plan['bindings'] = dict(template['bindings'], **{p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files})
    if kind == 'production':
        p = BASE/'pilot/time-refinement-manifest.json'
        plan['bindings'][p.relative_to(e.ROOT).as_posix()] = e.digest(p)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    bindings()


def run(refinement, workers):
    global REFINEMENT, PREFIX
    REFINEMENT, PREFIX = refinement, 0
    return _run(refinement, workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run', 'completed', 'compare', 'chain'])
    parser.add_argument('--phase', choices=['pilot', 'production'], default='pilot')
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert set(os.sched_getaffinity(0)) <= set(range(16)) and 1 <= args.workers <= 16
    OUT = BASE/args.phase
    if args.command == 'prepare':
        prepare(args.phase)
    elif args.command == 'chain':
        for r in [1, 2, 4]:
            run(r, args.workers)
        compare()
    elif args.command == 'run':
        run(args.refinement, args.workers)
    elif args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1]))
    else:
        compare()
