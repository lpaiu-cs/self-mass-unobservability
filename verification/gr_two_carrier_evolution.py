"""Counterexample candidate: separate radiative and conductive heat relaxation.

One-temperature, diagonal quadratic entropy closure. Radiation uses the grey
transport collision time; conduction retains the old uncalibrated proper time.
This is not non-LTE radiation hydrodynamics or a physical stellar exterior.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import subprocess

import numpy as np
import sympy as sp
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu

import gr_implicit_reflecting_boundary as wall
import gr_radiative_boundary as radiative

old, e, ld = wall.old, wall.e, wall.e.ld
OUT = e.g.OUT/'gr-two-carrier-evolution'
SCALE = np.array([1e-6, .1, 1e-8, 1e-10, 1e-10], dtype=ld)
ATOL = np.array([1e-13, 1e-9, 1e-14, 1e-17, 1e-17]+[1e-13]*26, dtype=ld)


def opacity_parts(model, row):
    lr, lt, x = row
    nuclei = x/e.g.c.A
    hydrogen, helium = x[e.g.c.Z == 1].sum(), x[e.g.c.Z == 2].sum()
    p = np.asarray([nuclei@e.g.c.Z/nuclei.sum(), hydrogen,
                    np.clip(1-hydrogen-helium, 0, .1), lr/np.log(10.), lt/np.log(10.)], float)
    rad = radiative.radiative(model, p)
    zbar, _, _, r, t = p
    grid = model.data['conduction-logzs']
    z = np.clip(np.log10(zbar), grid[0], grid[-1])
    j = np.clip(np.searchsorted(grid, z, side='right')-1, 0, len(grid)-2)
    u = (z-grid[j])/(grid[j+1]-grid[j])
    cr, ct = model.data['conduction-logrhos'], model.data['conduction-logts']
    coeff = model.data['conduction-f_ary']
    spline = radiative.tables.spline
    cond = (1-u)*spline(cr, ct, coeff[:, :, :, j], r, t)+u*spline(cr, ct, coeff[:, :, :, j+1], r, t)
    rc, tc = np.clip(r, cr[0], cr[-1]), np.clip(t, ct[0], ct[-1])
    cond = np.array([3*tc-rc-cond[0]-3.51937938116756,
                     (-1-cond[1])*float(rc == r), (3-cond[2])*float(tc == t)])
    # Direct conduction coefficients avoid subtracting almost equal K values.
    return np.r_[10**rad[0], rad[1:], 10**cond[0], cond[1:]]


def material(row):
    return np.r_[e.material(row), opacity_parts(e.OPACITY, row)]


def combined(y):
    return np.column_stack([y[:, :4], y[:, 5:]])


native_state = FunctionType(e.Star.state.__code__, dict(vars(e), material=material))


class TwoCarrierStar(wall.ReflectingStar):
    def state(self, y):
        self.primitives = combined(y)
        z = native_state(self, self.primitives)
        aux, rho, T = z['aux'], z['rho'], z['T']
        factor = 16*ld('5.670400e-5')*T**3/(3*rho)
        z['Krad'], z['Kcond'] = factor/aux[:, 24], factor/aux[:, 27]
        z['tau_rad'] = 1/(e.C*rho*aux[:, 24])
        z['Qrad'] = y[:, 4]*self.qscale
        z['Qcond'] = z['Q']-z['Qrad']
        assert np.max(abs((z['Krad']+z['Kcond'])/z['K']-1)) < 2e-14
        self.current = z
        return z

    def gradient(self, value, odd=False):
        z = getattr(self, 'current', None)
        if z is not None and value is z['Qrad']:
            a, N = z['a'], z['N']
            ar = a*a*(4*np.pi*e.GRAV*self.r*z['E']-z['m']/self.r**2)
            return a/N*self.divergence(self.faces(N*value/a, odd=True))-value*(z['nur']-ar+2/self.r)
        return super().gradient(value, odd=odd)

    def rhs(self, delta):
        y = self.base+delta
        z = self.state(y)
        rho,T,P,w,v,Q,W,D,E,S,R,a,N,nur,at = [z[k] for k in
            ['rho','T','P','w','v','Q','W','D','E','S','R','a','N','nur','at']]
        aux = z['aux']
        speed = e.C*N*v/a
        fB, fE, fS = self.faces(N*D*v, odd=True), self.faces(N*S/a, odd=True), self.faces(N*R)
        Bdot = -e.C*self.divergence(fB)-at*D
        Sdot = -e.C*self.divergence(fS)-at*S-N*nur*e.C*E+e.C*N*P*self.area_difference_over_volume
        matrix = np.zeros((self.n, 5, 5), dtype=ld)
        right = np.zeros((self.n, 5), dtype=ld)
        matrix[:, 0, 0], matrix[:, 0, 2] = a*D, a*rho*W**3*v
        right[:, 0] = Bdot+a*speed*self.gradient(D)
        matrix[:, 1, :4] = np.column_stack([rho*aux[:, 9]-P, rho*aux[:, 10], 2*Q*W*W, v*self.qscale])
        right[:, 1] = -e.C*N/(a*W*W)*self.gradient(Q, odd=True)-2*Q*(e.C*N*nur/a+v*at/a+e.C*N/(a*self.r))
        er, et = rho*(z['rest']+z['u']+aux[:, 9]), rho*aux[:, 10]
        pr, pt = P*aux[:, 5], P*aux[:, 6]
        matrix[:, 2, :4] = a[:, None]*np.column_stack([W*W*v*(er+pr), W*W*v*(et+pt),
                               W**4*((1+v*v)*w+4*Q*v), (1+v*v)*W*W*self.qscale])
        right[:, 2] = Sdot+a*speed*self.gradient(S, odd=True)
        force = e.C*self.gradient(y[:, 1])/(a*W)+e.C*W*nur/a+W*v*at/(N*a)
        for row, q, K, tau, kr, kt in [
            (3, z['Qcond'], z['Kcond'], e.TAU, -aux[:, 28], 5-aux[:, 29]),
            (4, z['Qrad'], z['Krad'], z['tau_rad'], 1, 5)]:
            h = K*T/e.C**2
            matrix[:, row, :3] = (W/N)[:, None]*np.column_stack(
                [-tau*q*kr/2, h*v-tau*q*kt/2, h*W*W])
            matrix[:, row, row] = W/N*tau*self.qscale
            if row == 3:
                matrix[:, row, 4] = -W/N*tau*self.qscale
            right[:, row] = -q-h*force
        scales = np.max(abs(matrix), axis=2)
        A, b = np.asarray(matrix/scales[:, :, None], float), np.asarray(right/scales, float)
        rates = np.linalg.solve(A, b[:, :, None])[:, :, 0].astype(ld)
        residual = float(np.max(abs(np.einsum('nij,nj->ni', A, rates)-b)/(1+abs(b))))
        assert residual < 1e-10, residual
        out = np.zeros_like(y)
        for j in range(5):
            spatial = (self.gradient(z['Q'] if j == 3 else z['Qrad'], odd=True)/self.qscale
                       if j >= 3 else self.gradient(y[:, j], odd=(j == 2)))
            out[:, j] = rates[:, j]-speed*spatial
        left = np.vstack([y[0, 5:]*0, np.diff(y[:, 5:], axis=0)/np.diff(self.r)[:, None]])
        right_x = np.vstack([np.diff(y[:, 5:], axis=0)/np.diff(self.r)[:, None], y[-1, 5:]*0])
        out[:, 5:] = -speed[:, None]*np.where(speed[:, None] >= 0, left, right_x)
        z['time_matrix_residual'] = residual
        z['boundary_rates'] = np.array([-e.C*4*np.pi*self.rf[-1]**2*fB[-1], -e.C*4*np.pi*self.rf[-1]**2*fE[-1]])
        return out, z


def initialize(pool):
    star = old.initialize(pool)
    model = radiative.tables.Opacity()
    rows = zip(star.base[:, 0], star.base[:, 1], star.base[:, 4:])
    parts = np.array([opacity_parts(model, row) for row in rows], dtype=ld)
    aux = np.column_stack([star.initial_aux, parts])
    rad_fraction = aux[:, 21]/aux[:, 24]
    star.base = np.column_stack([star.base[:, :4], star.base[:, 3]*rad_fraction, star.base[:, 4:]])
    star.__class__ = TwoCarrierStar
    star.initial_aux = aux
    star.material_cache = {e.material_key(row): value for row, value in zip(
        zip(star.base[:, 0], star.base[:, 1], star.base[:, 5:]), aux)}
    return star


class LocalTangent(TwoCarrierStar):
    def state(self, y):
        d, aux = y-self.anchor, self.reference['aux'].copy()
        aux[:, 1] *= np.exp(aux[:, 5]*d[:, 0]+aux[:, 6]*d[:, 1])
        aux[:, 2] += aux[:, 9]*d[:, 0]+aux[:, 10]*d[:, 1]
        for j in [21, 24, 27]:
            aux[:, j] *= np.exp(aux[:, j+1]*d[:, 0]+aux[:, j+2]*d[:, 1])
        # Enforce the exact combiner in the approximate tangent too.
        aux[:, 21] = aux[:, 24]*aux[:, 27]/(aux[:, 24]+aux[:, 27])
        self.material_cache = {e.material_key(row): value for row, value in zip(
            zip(y[:, 0], y[:, 1], y[:, 5:]), aux)}
        z = super().state(y)
        for key in ['mf', 'm', 'a', 'N', 'nur']:
            z[key] = self.reference[key]
        z['at'] = -4*np.pi*e.GRAV*e.C*self.r*z['N']*z['a']**2*z['S']
        return z


def tangent(star, delta, z):
    model = LocalTangent.__new__(LocalTangent)
    model.__dict__ = dict(star.__dict__)
    model.pool = old.CachedOnly()
    model.anchor, model.reference = star.base+delta, z
    return model


def jacobian(model, delta):
    rows, columns, values = [], [], []
    step = ld('.001')
    for field in range(5):
        for color in range(5):
            selected = np.arange(color, model.n, 5)
            change = np.zeros_like(delta)
            change[selected, field] = step*SCALE[field]
            response = (model.rhs(delta+change)[0][:, :5]-model.rhs(delta-change)[0][:, :5])/(2*step*SCALE)
            for offset in range(-2, 3):
                target = selected+offset
                valid = (target >= 0) & (target < model.n)
                target, source = target[valid], selected[valid]
                rows.extend((5*target[:, None]+np.arange(5)).ravel())
                columns.extend(np.repeat(5*source+field, 5))
                values.extend(np.asarray(response[target], float).ravel())
    return coo_matrix((values, (rows, columns)), shape=(5*model.n, 5*model.n)).tocsc()


def species(star, stage_base, h, y):
    view = SimpleNamespace(**dict(star.__dict__, base=combined(star.base)))
    return old.species(view, combined(stage_base), h, combined(y))


def energy_budget(star, delta, z, mass_flux):
    view = SimpleNamespace(**dict(star.__dict__, base=combined(star.base)))
    return old.energy_budget(view, combined(delta), z, mass_flux)


def cones(z):
    w, aux = z['w'], z['aux']
    j = z['Q']/w
    matrix, spatial = np.zeros((len(w), 5, 5)), np.zeros((len(w), 5, 5))
    matrix[:, 0, 0], spatial[:, 0, 2] = 1, 1
    matrix[:, 1, 1], matrix[:, 1, 2] = z['rho']*aux[:, 10]/w, 2*j
    spatial[:, 1, 2], spatial[:, 1, 3] = (z['P']-z['rho']*aux[:, 9])/w, 1
    matrix[:, 2, 2], matrix[:, 2, 3] = 1, 1
    spatial[:, 2, 0], spatial[:, 2, 1], spatial[:, 2, 2] = z['P']*aux[:, 5]/w, z['P']*aux[:, 6]/w, 2*j
    for row, q, K, tau, kr, kt in [
        (3, z['Qcond'], z['Kcond'], e.TAU, -aux[:, 28], 5-aux[:, 29]),
        (4, z['Qrad'], z['Krad'], z['tau_rad'], 1, 5)]:
        h = K*z['T']/(e.C**2*w*tau)
        matrix[:, row, 0], matrix[:, row, 1], matrix[:, row, 2] = -q/w*kr/2, -q/w*kt/2, h
        matrix[:, row, row] = 1
        if row == 3:
            matrix[:, row, 4] = -1
        spatial[:, row, 1] = h
    roots = np.linalg.eigvals(np.linalg.solve(matrix, spatial))
    speed, imaginary = float(abs(roots.real).max()), float(abs(roots.imag).max())
    return dict(maximum_local_rest_characteristic_speed_over_c=speed,
        maximum_characteristic_imaginary_part=imaginary,
        sampled_cone_inside_light_cone=bool(speed < 1 and imaginary < 1e-10))


def symbolic():
    rho,T,k,c,a = sp.symbols('rho T k c a', positive=True)
    tau, K = 1/(rho*k*c), 4*a*c*T**3/(3*rho*k)
    beta = sp.cancel(tau/(K*T**2))
    assert sp.simplify(beta-3/(4*a*c**2*T**5)) == 0
    assert sp.simplify(rho*sp.diff(sp.log(beta), rho)) == 0
    assert sp.simplify(T*sp.diff(sp.log(beta), T)+5) == 0
    # theta=-D ln rho; variable tau must be differentiated before using it.
    kr, kt, tr, tt, nr, nt = sp.symbols('kr kt tr tt nr nt')
    divbeta = -nr+(tr-kr)*nr+(tt-kt-2)*nt
    assert sp.expand(divbeta+(1+kr-tr)*nr+(2+kt-tt)*nt) == 0
    q, dq, force, B, conductivity, temperature, relaxation = sp.symbols('q dq force B K T tau')
    entropy = -q*force/temperature**2-relaxation*q*dq/(conductivity*temperature**2)-q*q*B/2
    law = -(relaxation*dq+q)/conductivity-temperature**2*q*B/2
    assert sp.simplify(entropy.subs(force, law)-q*q/(conductivity*temperature**2)) == 0
    return dict(classification='Proven', passed=True,
        beta='tau_rad/(K_rad*T^2)=3/(4*a_rad*c^2*T^5)',
        material_entropy_derivative='theta+D ln beta=-D ln rho_B-5 D ln T',
        scope='Algebra within a diagonal sum of two heat entropy currents and a grey transport-time constitutive choice. No LTE, opacity, electron-time or continuum certificate.')


def check():
    assert symbolic()['passed']
    solver_check()
    star = initialize(None)
    delta = np.zeros_like(star.base)
    rates, z = star.rhs(delta)
    original = old.initialize(None)
    before = original.state(original.base)
    for key in ['m', 'mf', 'a', 'N', 'Q', 'P', 'u']:
        assert np.array_equal(z[key], before[key]), key
    assert np.isfinite(rates).all() and np.array_equal(z['boundary_rates'], np.zeros(2))
    assert np.all(energy_budget(star, delta, z, np.zeros(star.n+1, dtype=ld))[0] == 0)
    initial_cone = cones(z)
    assert initial_cone['sampled_cone_inside_light_cone'], initial_cone
    opacity = radiative.tables.Opacity()
    for i in [0, 3000, 5500, 5734]:
        row = (star.base[i, 0], star.base[i, 1], star.base[i, 5:])
        values = opacity_parts(opacity, row)
        for coordinate in [0, 1]:
            pair = []
            for sign in [-1, 1]:
                shifted = list(row)
                shifted[coordinate] += sign*ld('0.000005')
                pair.append(opacity_parts(opacity, shifted)[[0, 3]])
            estimated = np.log(pair[1]/pair[0])/1e-5
            expected = values[[1+coordinate, 4+coordinate]]
            assert np.max(abs(estimated-expected)/np.maximum(1, abs(expected))) < 1e-3
    model = tangent(star, delta, z)
    matrix = jacobian(model, delta)
    rng = np.random.default_rng(33)
    direction = rng.uniform(-1, 1, size=(star.n, 5)).astype(ld)
    shift = np.zeros_like(delta)
    step = ld('.0001')
    shift[:, :5] = step*direction*SCALE
    finite = (model.rhs(delta+shift)[0][:, :5]-model.rhs(delta-shift)[0][:, :5])/(2*step*SCALE)
    predicted = matrix.dot(np.asarray(direction, float).ravel()).reshape(star.n, 5)
    # Measure each field; a very stiff radiation row must not hide fluid errors.
    errors = np.max(abs(finite-predicted), axis=0)/np.maximum(np.max(abs(finite), axis=0), 1e-30)
    assert errors.max() < 1e-4, errors
    result = dict(classification='Counterexample candidate', passed=True, symbolic=symbolic(),
        initial_geometry_total_heat_and_EOS_bitwise_unchanged=True, initial_cone=initial_cone,
        tangent_direction_relative_errors=errors.astype(float).tolist(),
        radiation_collision_time_range_seconds=[float(z['tau_rad'].min()), float(z['tau_rad'].max())],
        full_GR_evolution=False, physical_exterior_match=False)
    print('TWO CARRIER CHECK', json.dumps(result), flush=True)
    return result


def solver_check():
    # Exercise the actual fifth-field solve and backtracking, not just its Jacobian.
    from scipy.sparse import diags
    class Toy:
        n = 1
        base = np.zeros((1, 31), dtype=ld)
        def rhs(self, value):
            rates = np.zeros_like(value)
            rates[:, :5] = -20*value[:, :5]-1e6*value[:, :5]**3
            self.current = {}
            return rates, self.current
    surrogate = SimpleNamespace(ATOL=np.full(31, 1e-14, dtype=ld), SCALE=np.ones(5, dtype=ld),
        tangent=lambda *args: None, jacobian=lambda *args: diags(np.full(5, -.1)),
        species=lambda *args: np.zeros((1, 26), dtype=ld))
    tested = FunctionType(stage.__code__, dict(globals(), engine=surrogate))
    base = np.zeros((1, 31), dtype=ld)
    base[0, :5] = [.2, .1, .001, .0001, .15]
    logs = []
    value, rates, _ = tested(Toy(), base, np.zeros_like(base), ld(1), logs.append)
    assert np.max(abs(value-base-rates)/(surrogate.ATOL+ld('1e-7')*np.maximum(abs(value), abs(base)))) <= 1
    assert value[0, 4] != 0
    assert any(t.get('fraction', 1) < 1 for row in logs for t in row['preceding_trials'])


def prepare():
    assert not OUT.exists(), 'Preserve the historical single-heat trajectories.'
    checked = check()
    plan = json.loads((wall.OUT/'plan.json').read_text())
    for rel, value in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == value, rel
    plan.update(checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        candidate='Total heat plus independently relaxing radiative heat; restart from the same initial star.',
        variables='lnrho_B,lnT,v/c,Qtotal/w0,Qrad/w0,X26; Qcond=Qtotal-Qrad',
        radiation_time='1/(c*rho_B*kappa_rad), with the full material derivative of tau/(K*T^2).',
        conduction_time=str(e.TAU), radiation_time_physically_calibrated=False,
        conduction_time_physically_calibrated=False,
        initial_split='Qrad=Qtotal*Krad/(Krad+Kcond); preserve Qtotal and every initial fluid/metric quantity.',
        nonlinear_absolute_tolerances=ATOL.astype(float).tolist(), implementation_check=checked,
        imported_sources=[dict(classification='Imported from prior work', url='https://arxiv.org/abs/2412.00275',
                              use='Optically thick grey radiative conductivity and collision-time relation only.'),
                          dict(classification='Imported from prior work', url='https://arxiv.org/abs/astro-ph/9609119',
                              use='State-dependent full quadratic heat entropy transport, equation 2.22.')],
        limitations='One common temperature, diagonal heat entropy with no cross coefficients or anisotropic radiation stress, original prescribed electron time, original table mass/composition conventions, and a reflecting wall. Rosseland time is a grey constitutive approximation, not a measured absorption time. Does not solve the known optically thin outer layer, dense-plasma photon EOS, exterior, reactions or observations.')
    for path in [Path(__file__), Path(radiative.__file__), wall.OUT/'plan.json']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    print('PREPARED TWO CARRIER NATIVE EVOLUTION', flush=True)


# The proven current-inverse secant helper is unchanged.
direction = wall.prior.solver.direction


def stage(star, stage_base, seed, h, log):
    delta, history, factor = seed.copy(), [], None
    rates, z = star.rhs(delta)
    evaluations, trials = 1, []
    for iteration in range(24):
        residual = delta-stage_base-h*rates
        tolerance = engine.ATOL+ld('1e-7')*np.maximum(abs(delta), abs(stage_base))
        norm = float(np.max(abs(residual)/tolerance))
        log(dict(iteration=iteration, residual_norm=norm, native_evaluations=evaluations,
            maximum_absolute_residual=np.max(abs(residual), axis=0).astype(float).tolist(),
            preceding_trials=trials))
        if norm <= 1:
            return delta, rates, z
        if iteration == 23:
            break
        if factor is None or iteration % 4 == 0:
            J = engine.jacobian(engine.tangent(star, delta, z), delta)
            factor = splu(eye(5*star.n, format='csc')-float(h)*J)
        x = (delta[:, :5]/engine.SCALE).ravel().copy()
        f = (residual[:, :5]/engine.SCALE).ravel().copy()
        mixed, raw = direction(x, f, factor, history)
        trials, accepted = [], None
        proposals = [mixed, raw] if history else [raw]
        for proposal in proposals:
            correction = proposal.reshape(star.n, 5).astype(ld)*engine.SCALE
            fraction = min(1., .1/max(float(np.max(abs(correction[:, :2]))), 1e-300),
                           .01/max(float(np.max(abs(correction[:, 2]))), 1e-300))
            for backtrack in range(8):
                candidate = delta.copy()
                candidate[:, :5] += fraction*correction
                star.current = z
                candidate[:, 5:] = engine.species(star, stage_base, h, star.base+candidate)
                evaluations += 1
                try:
                    candidate_rates, candidate_z = star.rhs(candidate)
                    candidate_residual = candidate-stage_base-h*candidate_rates
                    candidate_tolerance = engine.ATOL+ld('1e-7')*np.maximum(abs(candidate), abs(stage_base))
                    score = float(np.max(abs(candidate_residual)/candidate_tolerance))
                    trials.append(dict(fraction=fraction, norm=score))
                    if score <= 1 or score < norm*(1-1e-4*fraction):
                        accepted = candidate, candidate_rates, candidate_z
                        break
                except (AssertionError, np.linalg.LinAlgError, FloatingPointError) as error:
                    trials.append(dict(fraction=fraction, invalid_trial=repr(error)))
                fraction /= 2
            if accepted is not None:
                break
        if accepted is None:
            raise RuntimeError(('Native residual line search failed', norm, trials))
        history = (history+[(x, f)])[-6:]
        delta, rates, z = accepted
    raise RuntimeError(('Native accelerated implicit stage did not converge', norm, trials))



engine = SimpleNamespace(**dict(vars(old), initialize=initialize, tangent=tangent,
    jacobian=jacobian, species=species, energy_budget=energy_budget, SCALE=SCALE, ATOL=ATOL))
solver = SimpleNamespace(**dict(vars(wall.prior.solver), stage=stage))
namespace = dict(vars(wall.prior), OUT=OUT, old=engine, solver=solver,
                 analysis=SimpleNamespace(cones=cones), second_tau=lambda: wall.prior.FIRST_TAU)
run = FunctionType(wall.prior.run.__code__, namespace)


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
