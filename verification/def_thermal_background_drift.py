"""Evolve the same affine heat/GR system to test its frozen-background premise.

Counterexample candidate: zero initial currents, frozen coefficients, internal
transport only. This is a drift test, not an orbit-long nonlinear star or a
physical radiative atmosphere. Native reaction directions are reported, not
silently replaced by external heat sources.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu
import def_orbital_radiative_transport as old

OUT = old.OUT.parent/'def-thermal-background-drift'
write = old.write


def symbolic():
    import sympy as s
    h, lam, k, G, H, M, K, F = s.symbols('h lam k G H M K F', nonzero=True)
    q, w, E, f, qr, wr, Er, fr = s.symbols('q w E f qr wr Er fr')
    f1 = (fr+h*lam*k*G)/(1+h*lam)
    gain = h*h*lam*k/(1+h*lam)
    assert s.simplify(Er+h*f1-(Er+h*fr/(1+h*lam)+gain*G)) == 0
    w1 = (q-qr+H*(E-Er))/h
    gr = s.expand(M*(w1-wr)/h+K*q-F*E)
    eliminated = (K+M/h**2)*q-(F-M*H/h**2)*E-M*(qr+h*wr+H*Er)/h**2
    assert s.simplify(gr-eliminated) == 0
    return dict(classification='Proven', passed=True,
        identity='SDIRK stage elimination retains qdot=w-H*Edot, M*wdot=-K*q+F*E, Edot=sum f and fdot=lambda*(k*theta-f), including the affine background theta.',
        scope='Algebraic identity of the declared frozen-coefficient model; no nonlinear/atmosphere or actual-background error certificate.')


class Problem(old.Problem):
    def __init__(self, degree=4):
        super().__init__(degree)
        m = self.model
        self.kappa = self.conductance*m.heat.geometry.tc
        self.theta0 = self.background_difference.astype(np.longdouble)
        self.H = self.H.real.astype(np.longdouble)
        self.Gq = self.Gq.real.astype(np.longdouble)
        self.GEraw = self.GEraw.real.astype(np.longdouble)
        self.Kr = self.Kx.real.astype(np.longdouble)
        self.Mr = self.Mx.real.astype(np.longdouble)
        self.F = self.original_force.real.astype(np.longdouble)
        self.max_error = 0.
        self.max_heat_error = 0.

    def stage(self, h):
        h = np.longdouble(h)
        m = self.model
        G = self.Kr+self.Mr/(h*h)
        B = self.F-(self.Mr@self.H)/(h*h)
        gain = h*h*np.sum(self.lam*self.kappa/(1+h*self.lam), axis=1)
        C = diags(1/gain)-self.GEraw
        block = bmat([[G, -B@diags(self.energy_scale)],
                      [-self.Gq, C@diags(self.energy_scale)]], format='csc')
        perm = self.permutation
        A = block[perm, :][:, perm].astype(float).tocsc()
        row = np.asarray(abs(A).max(axis=1).toarray()).ravel()
        scaled = diags(1/row)@A
        col = np.asarray(abs(scaled).max(axis=0).toarray()).ravel()
        factor = splu((scaled@diags(1/col)).tocsc(), permc_spec='NATURAL')
        scale = np.sqrt(abs(G.diagonal()))
        D = diags(1/scale)
        grfactor = splu((D@G@D).astype(float).tocsc(), permc_spec='NATURAL')

        def invert(rhs):
            y = factor.solve(np.asarray(rhs[perm]/row, float))/col
            answer = np.empty_like(y)
            answer[perm] = y
            return answer.astype(np.longdouble)

        def step(state):
            qr, wr, Er, fr = state
            rhs = self.Mr@(qr/(h*h)+wr/h+self.H@Er/(h*h))
            Ebase = Er+h*np.sum(fr/(1+h*self.lam), axis=1)
            bvec = np.r_[rhs, Ebase/gain+self.theta0]
            answer = invert(bvec)
            for _ in range(3):
                answer += invert(bvec-block@answer)
            E = answer[m.size:]*self.energy_scale
            # Preserve tiny scalar components separately from fluid/heat units.
            target = rhs+B@E
            q = (grfactor.solve(np.asarray(target/scale, float))/scale).astype(np.longdouble)
            for _ in range(3):
                defect = target-self.Kr@q-self.Mr@q/(h*h)
                q += (grfactor.solve(np.asarray(defect/scale, float))/scale).astype(np.longdouble)
            answer[:m.size] = q
            defect = bvec-block@answer
            error = float(np.max(abs(defect)/(abs(bvec)+abs(block)@abs(answer)+1e-100)))
            self.max_error = max(self.max_error, error)
            assert error < 1e-9, error
            theta = self.theta0+self.Gq@q+self.GEraw@E
            f = (fr+h*self.lam*self.kappa*theta[:, None])/(1+h*self.lam)
            heat_error = float(max(abs(E-Er-h*f.sum(1)))/max(max(abs(E)), max(abs(Er)), 1e-100))
            self.max_heat_error = max(self.max_heat_error, heat_error)
            assert heat_error < 1e-9, heat_error
            w = (q-qr+self.H@(E-Er))/h
            return q, w, E, f
        return step

    def read(self, state, seconds):
        q, w, E, f = state
        m = self.model
        full = np.zeros(len(m.heat.edges), np.longdouble)
        full[self.active] = E
        temp = (self.Tq@q+self.TE@full).real
        velocity = self.speed*(m.nativeV[0]@(w-self.H@f.sum(1)))
        scalar = m.nativeV[1]@q
        index = int(np.argmax(abs(temp)))
        native_index = len(temp)-1-index
        weights = m.original.weights
        balance = float(abs(np.sum(-np.diff(full), dtype=np.longdouble))/max(max(abs(E)), 1e-100))
        return dict(seconds=float(seconds), maximum_delta_lnT=float(max(abs(temp))),
            hottest_change_native_cell=native_index, delta_lnT_at_extremum=float(temp[index]),
            outer_delta_lnT=float(temp[-1]), cell132_delta_lnT=float(temp[len(temp)-1-132]),
            temperature_mass_RMS=float(np.sqrt(weights@temp**2)),
            velocity_mass_RMS_m_s=float(np.sqrt(weights@velocity**2)),
            scalar_mass_RMS=float(np.sqrt(weights@scalar**2)),
            maximum_face_energy_erg=float(max(abs(E))), internal_energy_balance=balance), temp, velocity

    def evolve(self, seconds, steps, label):
        start = time.monotonic()
        m = self.model
        dt = np.longdouble(seconds/m.heat.geometry.tc/steps)
        gamma = 1-1/np.sqrt(np.longdouble(2))
        step = self.stage(gamma*dt)
        state = (np.zeros(m.size, np.longdouble), np.zeros(m.size, np.longdouble),
                 np.zeros(len(self.active), np.longdouble), np.zeros_like(self.lam))
        history = [self.read(state, 0)[0]]
        temperatures = [np.zeros(len(self.theta))]
        for i in range(steps):
            first = step(state)
            base = tuple(a+(1-gamma)/gamma*(b-a) for a, b in zip(state, first))
            state = step(base)
            row, temp, velocity = self.read(state, seconds*(i+1)/steps)
            history.append(row)
            temperatures.append(temp)
            assert row['maximum_delta_lnT'] < .05, 'Declared tangent window exceeded five percent'
        np.savez_compressed(OUT/(label+'.npz'), q=state[0], w=state[1], E=state[2],
                            currents=state[3], temperature=temperatures,
                            native_radius=m.original.native, active_faces=self.active,
                            native_velocity=velocity, grid=m.grid, indices=m.indices)
        result = dict(classification='Counterexample candidate', history=history,
            steps=steps, degree=m.degree, seconds=time.monotonic()-start,
            maximum_linear_residual=self.max_error, maximum_heat_residual=self.max_heat_error,
            maximum_energy_balance=max(x['internal_energy_balance'] for x in history),
            full_nonlinear_background=False, physical_atmosphere=False)
        write(OUT/(label+'.json'), result)
        print('THERMAL DRIFT', label, result['seconds'], history[-1], flush=True)
        return result


def prepare():
    assert not OUT.exists()
    OUT.mkdir()
    write(OUT/'plan.json', dict(classification='Counterexample candidate', checkpoint='5c754a94c',
        claim='Actually evolve the same full-radius affine radiative/conductive GR equations from zero heat currents and the stored temperature gradient, to test whether their frozen-background premise can persist toward an orbit.',
        decision='Use a short time window selected from the stored constitutive startup, not an orbit integration. Record the first one-percent temperature-drift crossing. A crossing invalidates a one-percent frozen-state premise for this declared transport-only closed-native-boundary candidate, not every physical atmosphere. Do not claim a physical nonlinear drift bound from this linear test.',
        horizon_rule='Find the first time at which the no-motion constitutive startup gives maximum reconstructed |delta lnT|=0.02, capped at the orbital period. Use the same horizon for every path. Stop if actual tangent drift exceeds0.05. Do not extend after a failed gate.',
        equations='qdot=w-H*Edot; M*wdot=-K*q+F_source*E; Edot=sum f_l; fdot_l=lambda_l*(tc*conductance_l*(theta0_difference+Gq*q+GE*E)-f_l). Same native thermal capacity, mass constraint and heat momentum; no scalar orbital forcing during this secular-background test.',
        scope='Frozen grey single-temperature transport plus saved microscopic conduction. Zero omitted outer flux is the declared test boundary, not a certified isolated-star atmosphere. Frozen nuclear composition and neutrino source directions are reported separately, not inserted as nonconservative heating. No scalar charge inference from sub-light-crossing boundary samples.',
        paths=['p4-16', 'p4-32', 'p4-64', 'p2-64'],
        gates=dict(temperature_time=.01, temperature_order=1.5, temperature_space=.01,
                   GR_linear=1e-9, heat_linear=1e-9, balance=2e-13),
        diagnostic='Velocity and scalar time/space contrasts are retained separately; only passing temperature convergence supports the bounded frozen-temperature conclusion, never full four-component convergence.',
        budget=dict(pilot_seconds=120, total_seconds=600, CPU_threads=1, memory_GB=4,
                    new_EOS_calls=0, new_collision_states=0, automatic_expansion=False),
        forecast='Prior full-radius assembly11s and harmonic solve0.56s. Reuse four constant SDIRK stage factorizations; measure eight pilot steps before production. Time-step solve cost is unmeasured.',
        bindings={str(p):old.old.go.task.digest(p) for p in [Path(__file__), Path(old.__file__),
            old.OUT/'plan.json', old.OUT/'result.json', old.OUT/'all/fine-bank.npz',
            old.native.OUT/'coefficients.npz', old.native.OUT/'sources.npz']}))
    write(OUT/'symbolic.json', symbolic())
    signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    p = Problem()
    setup = time.monotonic()-start
    def approximate(seconds):
        _, energy = p.model.heat.faces(seconds/p.model.heat.geometry.tc)
        return float(max(abs(p.TE@energy)))
    period = 2*np.pi/p.omega*p.model.heat.geometry.tc
    horizon = brentq(lambda t:approximate(t)-.02, 0, period, xtol=1e-15)
    d = p.model.heat.d
    sources = np.load(old.native.OUT/'sources.npz')
    info = dict(classification='Counterexample candidate', coordinate_horizon_seconds=horizon,
        orbital_period_seconds=period, horizon_over_orbit=horizon/period,
        no_motion_endpoint_delta_lnT=approximate(horizon),
        maximum_direct_native_reaction_delta_lnT=float(max(abs(sources['rest_to_internal_heating']/d['thermo'][:, 3]))*horizon),
        maximum_frozen_composition_increment=float(max(abs(sources['dxdt']).ravel())*horizon),
        source_scope='Direct native source tangents only; not a bound after nonlinear fluid/metric/chemical evolution.')
    write(OUT/'input.json', info)
    row = p.evolve(horizon, 8, 'pilot')
    forecast = 1.5*(2*setup+row['seconds']*176/8+20)
    write(OUT/'pilot.json', dict(classification='Counterexample candidate', row=row,
        setup_seconds=setup, seconds=time.monotonic()-start, forecast_seconds=forecast,
        assumption='Pilot cost scaled to176 production steps, two assemblies,20s and50percent margin; fixed factorization overhead makes this conservative but p2 remains unmeasured.'))
    print('BACKGROUND PILOT', info, forecast, flush=True)


def run():
    plan = json.loads((OUT/'plan.json').read_text())
    pilot = json.loads((OUT/'pilot.json').read_text())
    assert not (OUT/'result.json').exists() and pilot['forecast_seconds'] < 600
    for p, digest in plan['bindings'].items():
        assert old.old.go.task.digest(Path(p)) == digest, p
    signal.alarm(int(600-pilot['seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    horizon = json.loads((OUT/'input.json').read_text())['coordinate_horizon_seconds']
    p = Problem()
    rows = [p.evolve(horizon, n, 'p4-'+str(n)) for n in [16, 32, 64]]
    arrays = [np.load(OUT/f'p4-{n}.npz')['temperature'] for n in [16, 32, 64]]
    norm = float(np.max(abs(arrays[-1])))
    first = float(np.max(abs(arrays[0]-arrays[1][::2]))/norm)
    last = float(np.max(abs(arrays[1]-arrays[2][::2]))/norm)
    comparison = dict(previous=first, last=last, order=float(np.log2(first/last)))
    time_pass = last < .01 and comparison['order'] > 1.5
    spatial_pass = False
    if time_pass:
        del p
        other = Problem(2).evolve(horizon, 64, 'p2-64')
        data = np.load(OUT/'p2-64.npz')['temperature']
        comparison['spatial'] = float(np.max(abs(data-arrays[-1]))/norm)
        spatial_pass = comparison['spatial'] < .01
    fine = rows[-1]['history']
    t = np.array([r['seconds'] for r in fine])
    value = np.array([r['maximum_delta_lnT'] for r in fine])
    crossing = float(np.interp(.01, value, t)) if max(value) >= .01 else None
    result = dict(classification='Counterexample candidate', passed=time_pass and spatial_pass,
        temperature_comparison=comparison, linear_interpolated_one_percent_crossing_seconds=crossing,
        final=fine[-1], seconds=time.monotonic()-start,
        total_compute_seconds=time.monotonic()-start+pilot['seconds'],
        full_physical_background=False, full_four_component_convergence=False,
        physical_atmosphere=False, full_goal_complete=False)
    write(OUT/'result.json', result)
    print('BACKGROUND RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'run'])
    globals()[parser.parse_args().action]()
