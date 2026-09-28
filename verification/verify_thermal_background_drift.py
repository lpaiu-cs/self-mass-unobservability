"""Independent endpoint, comparison and boundary-balance audit of Phase87."""
from pathlib import Path
import hashlib
import json
import resource
import signal
import time
import numpy as np
import sympy as s
import def_thermal_background_drift as run


def main():
    assert not (run.OUT/'audit.json').exists()
    signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    for path, sha in json.loads((run.OUT/'plan.json').read_text())['bindings'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha, path
    result = json.loads((run.OUT/'result.json').read_text())
    assert result['passed'] and not result['full_four_component_convergence']
    series = {n: json.loads((run.OUT/f'p4-{n}.json').read_text()) for n in [16, 32, 64]}
    spatial = json.loads((run.OUT/'p2-64.json').read_text())
    diagnostics = {}
    for field in ['velocity_mass_RMS_m_s', 'scalar_mass_RMS', 'temperature_mass_RMS']:
        arrays = [np.array([r[field] for r in series[n]['history']]) for n in [16, 32, 64]]
        norm = max(abs(arrays[-1]))
        previous = max(abs(arrays[0]-arrays[1][::2]))/norm
        last = max(abs(arrays[1]-arrays[2][::2]))/norm
        other = np.array([r[field] for r in spatial['history']])
        diagnostics[field] = dict(previous=float(previous), last=float(last),
            order=float(np.log2(previous/last)), spatial=float(max(abs(other-arrays[-1]))/norm))
    p = run.Problem()
    d = p.model.heat.d
    stored = np.load(run.OUT/'p4-64.npz')
    state = tuple(stored[k] for k in ['q', 'w', 'E', 'currents'])
    horizon = result['final']['seconds']
    _, temp, _ = p.read(state, horizon)
    score = float(max(abs(temp-stored['temperature'][-1]))/max(abs(temp)))
    assert score < 1e-12
    full = np.zeros(len(p.model.heat.edges), np.longdouble)
    full[p.active] = stored['E']
    debit = -np.diff(full)
    balance = float(abs(sum(debit))/max(abs(full)))
    assert balance < 2e-13
    crossings = {}
    for n, data in series.items():
        t = np.array([r['seconds'] for r in data['history']])
        v = np.array([r['maximum_delta_lnT'] for r in data['history']])
        j = int(np.flatnonzero(v >= .01)[0])
        crossings[n] = dict(below_seconds=float(t[j-1]), above_seconds=float(t[j]),
            below_delta_lnT=float(v[j-1]), above_delta_lnT=float(v[j]),
            interpolated_seconds=float(np.interp(.01, v, t)))
    # Necessary energy balance for a motionless fixed outer shell in the
    # DECLARED LTE diffusion model; no statement about a true thin atmosphere.
    face = int(np.flatnonzero(p.model.heat.face_ids == len(d['faces_cm'])-2)[0])
    inward = float(np.sum(p.model.heat.amplitude[face]))
    sources = np.load(run.old.native.OUT/'sources.npz')
    nuclear = float(d['dm'][0]*d['A'][0]**2*d['N'][0]**2*sources['rest_to_internal_heating'][0])
    needed = inward+nuclear
    assert needed < 0
    Lin, Lout, source = s.symbols('Lin Lout source', real=True)
    energy_rate = Lin-Lout+source
    assert s.simplify(energy_rate.subs(Lout, Lin+source)) == 0
    boundary = dict(classification='Counterexample candidate',
        internal_face_outward_luminosity_erg_s=inward,
        direct_native_nuclear_thermal_power_erg_s=nuclear,
        steady_required_outer_outward_luminosity_erg_s=needed,
        outer_cell_baryon_mass_fraction=float(d['dm'][0]/sum(d['dm'])),
        outer_T_K=float(np.exp(d['lnT'][0])), neighbor_T_K=float(np.exp(d['lnT'][1])),
        conditional_no_outgoing_only_steady_state=True,
        scope='For this fixed, motionless outer shell, this LTE flux law and retained sources, Udot=Lin-Lout+Q. Since Lin+Q<0, no Lout>=0 preserves its energy. Thin-layer transport and the actual native-exterior connection are not certified; this is not a physical-star cooling bound.')
    value = dict(classification='Counterexample candidate', passed=True,
        temperature_endpoint_reconstruction=score, energy_balance=balance,
        crossings=crossings, component_diagnostics=diagnostics, boundary=boundary,
        maximum_GR_residual=max(x['maximum_linear_residual'] for x in list(series.values())+[spatial]),
        maximum_heat_residual=max(x['maximum_heat_residual'] for x in list(series.values())+[spatial]),
        symbolic_boundary_identity=dict(classification='Proven', passed=True,
            identity='Udot=Lin-Lout+Q requires Lout=Lin+Q for a stationary motionless shell. Lin+Q<0 excludes every outgoing-only boundary in the declared model.'),
        seconds=time.monotonic()-start,
        total_compute_seconds=time.monotonic()-start+result['total_compute_seconds'],
        full_physical_background=False, full_four_component_convergence=False,
        full_goal_complete=False)
    assert value['total_compute_seconds'] < 600
    run.write(run.OUT/'audit.json', value)
    print(json.dumps(value), flush=True)


if __name__ == '__main__':
    main()
