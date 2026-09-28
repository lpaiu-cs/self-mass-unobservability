"""Correct the common EOS heat map and rerun the actual heat/GR problems."""
from pathlib import Path
from types import FunctionType
import argparse
import json
import resource
import signal
import time
import numpy as np
import def_native_temperature_closure as closure
import def_orbital_radiative_transport as orbit
import def_thermal_background_drift as drift

OUT = orbit.OUT.parent/'def-native-temperature-reconnect'
write = orbit.write
orbital_solve = FunctionType(orbit.direct_solve.__code__, dict(orbit.namespace, OUT=OUT/'orbital'))
time_evolve = FunctionType(drift.Problem.evolve.__code__, dict(vars(drift), OUT=OUT/'background'))


class Orbital(orbit.Problem):
    def __init__(self, degree=4):
        super().__init__(degree)
        self.closure = closure.correct(self)


class Background(drift.Problem):
    def __init__(self, degree=4):
        super().__init__(degree)
        self.closure = closure.correct(self)


def native_control():
    d = np.load(orbit.native.OUT/'coefficients.npz')
    eos = orbit.native.h.molecular.model.EOS()
    ids = [0, 132, 600, 1175, 2300, 4000, 5734]
    rows = []
    for i in ids:
        lr, lt, x = float(d['lnd'][i]), float(d['lnT'][i]), d['X'][i]
        raw = eos(2, lr, lt, x)
        cv = raw[10]
        cp = d['thermo'][i, 5]
        estimates = []
        for step in [2e-5, 1e-5]:
            minus = eos(2, lr, lt-step, x)
            plus = eos(2, lr, lt+step, x)
            measured = (np.longdouble(plus[2])-np.longdouble(minus[2]))/(2*step)
            estimates.append(float(measured))
        error = abs(estimates[-1]/cv-1)
        assert error < 2e-5, (i, error)
        assert abs(cv/d['thermo'][i, 3]-1) < 1e-7
        rows.append(dict(cell=i, cvT=float(cv), cpT=float(cp), cp_over_cv=float(cp/cv),
            constant_density_native_heat_derivative=estimates, finite_relative_error=float(error),
            old_fixed_density_temperature_relative_bias=float(cv/cp-1)))
    value = dict(classification='Counterexample candidate', passed=True, native_calls=len(ids)*5,
        rows=rows, model_error_confirmed=True,
        scope='Native EOS finite-difference positive control at fixed rho,X; not continuum thermodynamic certification.')
    write(OUT/'native-control.json', value)
    return value


def prepare():
    assert not OUT.exists()
    OUT.mkdir()
    (OUT/'orbital').mkdir()
    (OUT/'background').mkdir()
    paths = [Path(__file__), Path(closure.__file__), Path(orbit.__file__), Path(drift.__file__),
        Path(orbit.old.thermal.__file__), orbit.native.OUT/'coefficients.npz',
        orbit.OUT/'all/fine-bank.npz', orbit.OUT/'result.json', drift.OUT/'plan.json',
        drift.OUT/'result.json']
    write(OUT/'plan.json', dict(classification='Counterexample candidate', checkpoint='a3f16fa71',
        claim='Repair the native EOS heat-temperature closure in the actual orbital and background GR feedback equations, then rerun the previously approved finite domains and original gates.',
        cause='The physical-density temperature map took thermo[:,5]=cp*T where thermo[:,3]=cv*T is required. The fixed-pressure source rho_ref correctly keeps cp. Historical files are frozen, and their affected physical temperature/feedback conclusions are superseded, not their raw numerical records.',
        repair='Shared closure rebuilds TE, GE and GEraw/energy scale; both time and orbital consumers use it. No output-only correction, EOS/table replacement, source-energy change or acceptance relaxation.',
        cases=['orbital p4/p2 harmonics1,2,3', 'same0.496490906ms background p4-16/32/64 and p2-64'],
        gates=dict(native_heat_derivative=2e-5, EOS_chain=1e-7, GR=1e-9, heat=1e-9,
            balance=2e-13, adjoint=1e-5, orbital_space=.02,
            temperature_time=.01, temperature_order=1.5, temperature_space=.01),
        background_scope='Only the original temperature convergence target; velocity/scalar diagnostics retained. Same zero-current state, horizon, frozen coefficients/composition and omitted native outer flux. Physical atmosphere and orbit-long background remain unproved.',
        budget=dict(pilot_seconds=120, total_seconds=600, CPU_threads=1, memory_GB=4,
            native_control_calls=35, new_native_states_for_evolution=0,
            automatic_expansion=False),
        forecast='Recent full assembly about11s, joint harmonic solve0.55s and64-step path about2s. Native controls and corrected stiffness must first be measured; no orbit integration.',
        bindings={str(p):orbit.old.go.task.digest(p) for p in paths}))
    write(OUT/'symbolic.json', closure.symbolic())
    signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    native_control()
    native_seconds = time.monotonic()-start
    tick = time.monotonic()
    p = Orbital()
    setup = time.monotonic()-tick
    row = orbital_solve(p, 1, 'p4')
    seconds = time.monotonic()-start
    forecast = seconds+1.5*(4*setup+5*row['seconds']+25+30)
    write(OUT/'pilot.json', dict(classification='Counterexample candidate', row=row,
        closure=p.closure, seconds=seconds, setup_seconds=setup,
        native_control_seconds=native_seconds, forecast_seconds=forecast,
        assumption='Four remaining assemblies, five harmonic solves,25s for unchanged finite time paths and30s audit allowance,50percent margin. Time-path rate after repair not yet measured.'))
    print('NATIVE TEMPERATURE PILOT', seconds, forecast, p.closure, flush=True)


def run():
    plan = json.loads((OUT/'plan.json').read_text())
    pilot = json.loads((OUT/'pilot.json').read_text())
    assert not (OUT/'result.json').exists() and pilot['forecast_seconds'] < 600
    for path, sha in plan['bindings'].items():
        assert orbit.old.go.task.digest(Path(path)) == sha, path
    signal.alarm(int(600-pilot['seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    cases = {}
    for degree in [4, 2]:
        p = Orbital(degree)
        cases[f'p{degree}'] = [pilot['row'] if degree == 4 and n == 1
            else orbital_solve(p, n, f'p{degree}') for n in [1, 2, 3]]
        del p
    spatial = [abs(complex(*a['radiative_charge_gain'])/complex(*b['radiative_charge_gain'])-1)
        for a, b in zip(cases['p2'], cases['p4'])]
    write(OUT/'orbital-result.json', dict(classification='Counterexample candidate',
        passed=max(spatial)<.02, cases=cases, spatial=spatial))
    horizon = json.loads((drift.OUT/'input.json').read_text())['coordinate_horizon_seconds']
    p = Background()
    rows = [time_evolve(p, horizon, n, f'p4-{n}') for n in [16, 32, 64]]
    del p
    arrays = [np.load(OUT/f'background/p4-{n}.npz')['temperature'] for n in [16, 32, 64]]
    norm = np.max(abs(arrays[-1]))
    previous = float(np.max(abs(arrays[0]-arrays[1][::2]))/norm)
    last = float(np.max(abs(arrays[1]-arrays[2][::2]))/norm)
    comparison = dict(previous=previous, last=last, order=float(np.log2(previous/last)))
    time_passed = last < .01 and comparison['order'] > 1.5
    if time_passed:
        p = Background(2)
        time_evolve(p, horizon, 64, 'p2-64')
        other = np.load(OUT/'background/p2-64.npz')['temperature']
        comparison['spatial'] = float(np.max(abs(other-arrays[-1]))/norm)
    fine = rows[-1]['history']
    ts = np.array([r['seconds'] for r in fine])
    values = np.array([r['maximum_delta_lnT'] for r in fine])
    crossing = float(np.interp(.01, values, ts)) if max(values) >= .01 else None
    value = dict(classification='Counterexample candidate',
        passed=max(spatial)<.02 and time_passed and comparison.get('spatial', 1)<.01,
        orbital_spatial=spatial, temperature_comparison=comparison,
        background_final=fine[-1], linear_interpolated_one_percent_crossing_seconds=crossing,
        seconds=time.monotonic()-start, total_compute_seconds=pilot['seconds']+time.monotonic()-start,
        native_temperature_closure_corrected=True, full_four_component_convergence=False,
        full_physical_background=False, physical_atmosphere=False, full_goal_complete=False)
    write(OUT/'result.json', value)
    print('NATIVE TEMPERATURE RESULT', json.dumps(value), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'run'])
    globals()[parser.parse_args().action]()
