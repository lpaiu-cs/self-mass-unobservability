"""Propagate saved background drift through actual radial transport to charge.

Counterexample candidate: finite transport-input update only. GR geometry,
native EOS tangent and microscopic electron conduction remain fixed. This is
not a rebuilt hydrostatic star, an atmosphere, or an orbit-long background.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags
import def_orbital_radiative_transport as old
import def_thermal_background_drift as drift

OUT = old.OUT.parent/'def-background-transport-charge'
write = old.write
solve = FunctionType(old.direct_solve.__code__, dict(old.namespace, OUT=OUT))


def symbolic():
    import sympy as s
    ad, rhoL, displacement, slopeT, slopeR, heat = s.symbols('ad rhoL x t r heat')
    temperature = ad*rhoL-displacement*slopeT-heat
    density = temperature/ad+displacement*(slopeT/ad-slopeR)+heat/ad
    assert s.simplify(density-(rhoL-displacement*slopeR)) == 0
    return dict(classification='Proven', passed=True,
        identity='Eulerian density is the Lagrangian baryon/metric perturbation minus displacement times the background density gradient; heat changes temperature through the native energy equation, not by an extra density source.',
        scope='Same linearized conservative GR/native EOS reconstruction; no nonlinear EOS error bound.')


class Problem(old.Problem):
    def __init__(self, degree=4):
        super().__init__(degree)
        self.initial_conductance = self.conductance.copy()
        self.initial_lam = self.lam.copy()
        self.initial_theta = self.theta.copy()
        self.initial_Gq = self.Gq.copy()
        self.initial_GEraw = self.GEraw.copy()
        self.initial_gL = self.gL.copy()
        # Thermal arrays are outer-to-inner; heat-energy edge indices run
        # inner-to-outer. They label the same faces but are not interchangeable.
        self.native_faces = np.load(old.OUT/'all/fine-bank.npz')['faces']
        n = len(self.theta)
        faces = self.native_faces
        assert np.array_equal(self.active, n-faces)
        self.diff = coo_matrix((np.tile([-1., 1.], len(self.active)),
            (np.repeat(np.arange(len(self.active)), 2),
             np.c_[n-faces, n-1-faces].ravel())),
             shape=(len(self.active), n)).tocsr()

    def state(self, steps):
        m = self.model
        d = m.heat.d
        data = np.load(drift.OUT/f'p{m.degree}-{steps}.npz')
        assert np.array_equal(data['active_faces'], self.active)
        assert np.array_equal(data['grid'], m.grid)
        r = m.original.native
        slots = np.minimum(np.arange(len(r)), len(r)-2)
        slopeT = np.diff(d['lnT'][::-1])[slots]/np.diff(r)[slots]
        slopeR = np.diff(np.log(d['raw'][::-1, 0]))[slots]/np.diff(r)[slots]
        ad = self.point['adiabatic_T_rho']
        Rq = diags(1/ad)@self.Tq+diags(r*(slopeT/ad-slopeR))@self.V[0]
        _, _, J = old.old.go.task.fem.source_points(m.heat, self.point)
        b = 1-2*self.point['m']/r
        RE = -diags(1/(r*b))@J
        full = np.zeros(len(m.heat.edges), np.longdouble)
        full[self.active] = data['E']
        temp = (self.Tq@data['q']+self.TE@full).real
        density = (Rq@data['q']+RE@full).real
        assert max(abs(temp-data['temperature'][-1])) < 1e-15
        assert max(abs(temp)) < .05 and max(abs(density)) < .05
        return np.asarray(temp[::-1], float), np.asarray(density[::-1], float)

    def update(self, data):
        d = self.model.heat.d
        faces = self.native_faces
        Tshift, Rshift = data['delta_lnT'], data['delta_lnrho']
        opacity = data['opacity']
        ratio = np.exp(3*Tshift-Rshift)*d['opacity'][:, 0]/opacity
        collision = np.exp(Rshift)*opacity/d['opacity'][:, 0]
        self.conductance = self.initial_conductance.copy()
        self.conductance[:, 0] *= np.sqrt(ratio[faces-1]*ratio[faces])
        self.lam = self.initial_lam.copy()
        self.lam[:, 0] *= np.sqrt(collision[faces-1]*collision[faces])
        delta_theta = self.initial_theta*np.expm1(Tshift[::-1])
        self.theta = self.initial_theta+delta_theta
        self.Gq = self.initial_Gq+self.diff@diags(delta_theta)@self.Tq
        self.GEraw = self.initial_GEraw+(self.diff@diags(delta_theta)@self.TE)[:, self.active]
        self.gL = self.initial_gL+self.diff@(delta_theta*self.TL)
        self.energy_scale = 1/np.maximum(abs(self.GEraw.diagonal()), np.longdouble('1e-100'))
        assert np.all(self.conductance >= 0) and np.all(self.lam > 0)
        assert np.array_equal(self.conductance[:, 1:], self.initial_conductance[:, 1:])


def opacity(model, d, temperature, density, ids):
    return np.array([old.native.two.opacity_parts(model,
        (d['lnd'][i]+density[i], d['lnT'][i]+temperature[i], d['X'][i]))[0] for i in ids])


def inputs(p, steps, model):
    label = f'p{p.model.degree}-{steps}'
    assert not (OUT/(label+'-input.npz')).exists()
    temperature, density = p.state(steps)
    d = p.model.heat.d
    kappa = opacity(model, d, temperature, density, range(len(temperature)))
    np.savez_compressed(OUT/(label+'-input.npz'), delta_lnT=temperature,
                        delta_lnrho=density, opacity=kappa)
    return dict(delta_lnT=temperature, delta_lnrho=density, opacity=kappa)


def prepare():
    assert not OUT.exists()
    OUT.mkdir()
    paths = [Path(__file__), Path(old.__file__), Path(drift.__file__),
             old.OUT/'result.json', old.OUT/'all/fine-bank.npz',
             old.native.OUT/'coefficients.npz', drift.OUT/'audit.json']
    paths += [drift.OUT/(label+'.npz') for label in ['p4-32', 'p4-64', 'p2-64']]
    write(OUT/'plan.json', dict(classification='Counterexample candidate', checkpoint='93b43539f',
        claim='Carry the actually evolved temperature/density drift into radius-dependent Rosseland transport, solve the same orbital heat/GR system, and measure its external scalar-charge change.',
        decision='Resolve the CHANGE of outgoing charge against its own time/space contrasts, not merely convergence of the larger baseline. A small change excludes dominance only of this transport-input update; unresolved or large change does not trigger automatic refinement.',
        changed='At saved endpoints, evaluate the same opacity table at rho*exp(delta_lnrho),T*exp(delta_lnT), retained composition. Update Krad, radiative pole rate and redshifted-temperature factors in all internal-face heat maps. Hold native EOS tangent, GR geometry/operator, electron conduction and exterior wave propagation fixed.',
        scope='Finite transport sensitivity along the saved linear background trajectory. Not full background re-equilibration, updated EOS/geometry, instantaneous coefficient derivative feedback, atmosphere, or an orbit-long physical transfer function. The short-time boundary scalar itself is never used as radiation.',
        cases=['p4-64 harmonics1,2,3', 'p4-32 harmonics1,2,3', 'p2-64 harmonics1,2,3'],
        gates=dict(opacity_baseline_replay=1e-10, linear=1e-9, heat=1e-9,
            balance=2e-13, adjoint=1e-5, change_time=.01, change_space=.02),
        budget=dict(pilot_seconds=120, total_seconds=600, CPU_threads=1, memory_GB=4,
            new_EOS_calls=0, new_time_steps=0, automatic_expansion=False),
        forecast='Prior full GR assembly11s, harmonic solve below1s. Measure32 baseline/updated opacity evaluations before all three native updates, and one joint solve before production.',
        bindings={str(p):old.old.go.task.digest(p) for p in paths}))
    write(OUT/'symbolic.json', symbolic())
    signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    p = Problem()
    setup = time.monotonic()-start
    temperature, density = p.state(64)
    model = old.native.two.radiative.tables.Opacity()
    d = p.model.heat.d
    ids = np.unique(np.linspace(0, len(temperature)-1, 32).astype(int))
    tick = time.monotonic()
    zeros = np.zeros(len(temperature))
    base = opacity(model, d, zeros, zeros, ids)
    shifted = opacity(model, d, temperature, density, ids)
    elapsed = time.monotonic()-tick
    replay = float(max(abs(base/d['opacity'][ids, 0]-1)))
    assert replay < 1e-10
    forecast_inputs = elapsed*3*len(temperature)/(2*len(ids))
    forecast = 1.5*(3*setup+forecast_inputs+9+30)
    write(OUT/'opacity-pilot.json', dict(classification='Counterexample candidate',
        seconds=elapsed, baseline_relative=replay, input_forecast_seconds=forecast_inputs,
        total_pre_solve_forecast_seconds=forecast, sampled_cells=ids.tolist(),
        sampled_relative_opacity_change=(shifted/base-1).tolist()))
    assert forecast < 600, 'Input pilot exceeds the registered total budget'
    data = inputs(p, 64, model)
    p.update(data)
    row = solve(p, 1, 'p4-64')
    measured = time.monotonic()-start
    forecast = max(forecast, measured+1.5*(2*setup+2*forecast_inputs/3+8*row['seconds']+30))
    write(OUT/'pilot.json', dict(classification='Counterexample candidate', row=row,
        seconds=measured, setup_seconds=setup, forecast_seconds=forecast,
        assumption='Two remaining assemblies and measured opacity/solve rates plus30s and50percent margin; p2 has not yet been measured.'))
    print('BACKGROUND CHARGE PILOT', measured, forecast, flush=True)


def repair():
    assert (OUT/'first-failure.json').exists() and not (OUT/'pilot.json').exists()
    plan = json.loads((OUT/'first-plan.json').read_text())
    plan['repair'] = 'The no-op support test exposed reversed native face IDs versus ascending heat-energy edge IDs. Preserve the failure and original source; use the original bank face order for native coefficient and temperature maps, and active edge order only for E columns. Add a bitwise no-op control. No physical input, gate, memory or resolution change.'
    plan['bindings'][str(Path(__file__))] = old.old.go.task.digest(Path(__file__))
    for name in ['first-source.py', 'first-plan.json', 'first-failure.json', 'sparsity.json']:
        plan['bindings'][str(OUT/name)] = old.old.go.task.digest(OUT/name)
    write(OUT/'plan.json', plan)
    signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    p = Problem()
    setup = time.monotonic()-start
    zeros = np.zeros(len(p.theta))
    p.update(dict(delta_lnT=zeros, delta_lnrho=zeros, opacity=p.model.heat.d['opacity'][:, 0]))
    assert (p.Gq-p.initial_Gq).nnz == 0
    assert (p.GEraw-p.initial_GEraw).nnz == 0
    assert np.array_equal(p.gL, p.initial_gL)
    assert np.array_equal(p.conductance, p.initial_conductance)
    assert np.array_equal(p.lam, p.initial_lam)
    write(OUT/'noop-control.json', dict(classification='Counterexample candidate', passed=True,
        maps_and_coefficients_bitwise=True, native_faces_begin=p.native_faces[:3].tolist(),
        active_edges_begin=p.active[:3].tolist()))
    p.update(np.load(OUT/'p4-64-input.npz'))
    row = solve(p, 1, 'p4-64')
    spent = time.monotonic()-start
    prior = 120+json.loads((OUT/'sparsity.json').read_text())['seconds']
    opacity_seconds = json.loads((OUT/'opacity-pilot.json').read_text())['input_forecast_seconds']
    forecast = prior+spent+1.5*(2*setup+opacity_seconds+8*row['seconds']+30)
    write(OUT/'pilot.json', dict(classification='Counterexample candidate', row=row,
        seconds=prior+spent, setup_seconds=setup, forecast_seconds=forecast,
        failed_pilot_charged_seconds=120, diagnosis_seconds=prior-120,
        successful_pilot_seconds=spent,
        assumption='Failed pilot charged at its120s cap; measured diagnosis and repaired pilot. Two assemblies, measured remaining opacity/solve rates,30s and50percent margin.'))
    print('REPAIRED BACKGROUND CHARGE PILOT', spent, forecast, flush=True)


def run():
    plan = json.loads((OUT/'plan.json').read_text())
    pilot = json.loads((OUT/'pilot.json').read_text())
    for path, sha in plan['bindings'].items():
        assert old.old.go.task.digest(Path(path)) == sha, path
    assert not (OUT/'result.json').exists() and pilot['forecast_seconds'] < 600
    signal.alarm(int(600-pilot['seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    model = old.native.two.radiative.tables.Opacity()
    cases = {}
    for degree in [4, 2]:
        p = Problem(degree)
        for steps in ([64, 32] if degree == 4 else [64]):
            label = f'p{degree}-{steps}'
            data = np.load(OUT/(label+'-input.npz')) if label == 'p4-64' else inputs(p, steps, model)
            p.update(data)
            cases[label] = [pilot['row'] if label == 'p4-64' and n == 1 else solve(p, n, label)
                            for n in [1, 2, 3]]
        del p
    prior = json.loads((old.OUT/'result.json').read_text())['cases']
    pair = lambda v: [float(v.real), float(v.imag)]
    changes = {}
    for label, rows in cases.items():
        base = prior['all-p'+label[1]]
        changes[label] = np.array([complex(*a['radiative_charge_gain'])-
            complex(*b['radiative_charge_gain']) for a, b in zip(rows, base)])
    comparisons = []
    for i in range(3):
        ref = changes['p4-64'][i]
        row = cases['p4-64'][i]
        comparisons.append(dict(harmonic=i+1, charge_gain_change=pair(ref),
            delta_alpha_radiative_over_phi0=abs(ref)*row['actual_drive_amplitude']/.001,
            change_relative_to_original_thermal=abs(ref)/abs(complex(*prior['all-p4'][i]['radiative_charge_gain'])),
            time_relative_to_change=float(abs(changes['p4-32'][i]-ref)/max(abs(ref), 1e-100)),
            space_relative_to_change=float(abs(changes['p2-64'][i]-ref)/max(abs(ref), 1e-100))))
    passed = all(r['time_relative_to_change'] < .01 and r['space_relative_to_change'] < .02 for r in comparisons)
    value = dict(classification='Counterexample candidate', passed=passed, cases=cases,
        comparisons=comparisons, seconds=time.monotonic()-start,
        total_compute_seconds=time.monotonic()-start+pilot['seconds'],
        finite_transport_update_solved=True, full_background_feedback=False,
        physical_atmosphere=False, full_goal_complete=False)
    write(OUT/'result.json', value)
    print('BACKGROUND CHARGE RESULT', json.dumps(comparisons), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'repair', 'run'])
    globals()[parser.parse_args().action]()
