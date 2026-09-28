"""Audit native EOS closure through the solved density, pressure and heat."""
from pathlib import Path
import hashlib
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
import def_native_temperature_reconnect as run
import def_orbital_charge_audit as project


def main():
    assert not (run.OUT/'audit.json').exists()
    signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    for name, sha in json.loads((run.OUT/'plan.json').read_text())['bindings'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == sha, name
    result = json.loads((run.OUT/'result.json').read_text())
    orbital = json.loads((run.OUT/'orbital-result.json').read_text())
    assert result['passed'] and orbital['passed']
    p = run.Orbital()
    m = p.model
    d = m.heat.d
    point = p.point
    r = point['r']
    b = 1-2*point['m']/r
    alpha = -4*point['phi']
    A4 = np.exp(-8*point['phi']**2)
    dlz = r*r*point['v']**2/2-4*np.pi*r*r*A4*point['p']/b-point['m']/(r*b)
    Rq = -diags(r)@p.D[0]-diags(3+dlz+3*alpha*r*point['v'])@p.V[0]
    Rq -= diags(r*point['v']+3*alpha)@p.V[1]
    source = run.orbit.old.go.task.fem.source_points(m.heat, point)
    RE = -diags(1/(r*b))@source[2]
    gr = run.orbit.old.go.task.fem.base.task.h.gr
    geo = gr.G*.1*m.heat.geometry.R**2/gr.C**4
    rho = d['raw'][::-1, 0]
    cp = d['thermo'][::-1, 5]
    ad = d['thermo'][::-1, 4]
    gamma = d['raw'][::-1, 4]
    rt = d['raw'][::-1, 8]
    slots = np.minimum(np.arange(len(r)), len(r)-2)
    slopeT = np.diff(d['lnT'][::-1])[slots]/np.diff(r)[slots]
    stored = np.load(run.OUT/'background/p4-64.npz')
    full = np.zeros(len(m.heat.edges), np.longdouble)
    full[p.active] = stored['E']
    density = (Rq@stored['q']+RE@full).real
    heat = -(source[1]@full)/(rho*geo)
    rho_ref = rt*heat/cp
    pressure = gamma*(density-rho_ref)
    # Pressure-coordinate thermodynamics uses cp and subtracts the entropy
    # density source; no call to the corrected density-coordinate heat map.
    lagrange_T = ad/gamma*pressure+heat/cp
    expected_T = lagrange_T-r*slopeT*(p.V[0]@stored['q']).real
    map_T = (p.Tq@stored['q']+p.TE@full).real
    score = float(max(abs(map_T-expected_T))/max(abs(expected_T)))
    assert score < 1e-7
    assert max(abs(map_T-stored['temperature'][-1])) < 1e-15
    maxima = dict(GR=0., heat=0., adjoint=0.)
    for n in [1, 2, 3]:
        a = np.load(run.OUT/f'orbital/p4-{n}.npz')
        qa, q, E = [a[k] for k in ['qa', 'dq', 'E']]
        z = np.clongdouble(-1j*n*p.omega)
        B = p.original_force-z*z*(p.Mx@p.H)
        defect = p.Kx@q+z*z*(p.Mx@q)-B@E
        er = float(max(abs(defect)/(p.absK@abs(q)+abs(z*z)*(p.absM@abs(q))+abs(B)@abs(E)+1e-100)))
        gain = np.sum(p.conductance*m.heat.geometry.tc*p.lam/(z*(z+p.lam)), axis=1)
        heat_error = float(max(abs(E-gain*(p.Gq@(qa+q)+p.gL+p.GEraw@E)))/max(abs(E)))
        direct = (p.kl+z*z*p.ml)@q-p.fL@E
        reciprocal = -(B.T@qa+p.fL)@E
        adjoint = float(abs(direct-reciprocal)/abs(direct))
        assert er < 1e-9 and heat_error < 1e-9 and adjoint < 1e-5
        for name, value in [('GR', er), ('heat', heat_error), ('adjoint', adjoint)]:
            maxima[name] = max(maxima[name], value)
    previous = json.loads((run.orbit.OUT/'result.json').read_text())['cases']['all-p4']
    changes = [abs(complex(*new['radiative_charge_gain'])/complex(*old['radiative_charge_gain'])-1)
        for old, new in zip(previous, orbital['cases']['p4'])]
    base = json.loads((run.orbit.old.orbit.OUT.parent/'def-full-orbital-exterior/comparator-result.json').read_text())
    project.mp.mp.dps = 70
    mp = project.mp
    combined = [base['cases']['p4-tight'][0]]
    amplitudes = [mp.mpf(0)]
    for a, row in zip(base['cases']['p4-tight'][1:], orbital['cases']['p4']):
        combined.append(dict(charge=[mp.nstr(mp.mpf(str(x))+mp.mpf(str(y)), 65)
            for x, y in zip(a['charge'], row['radiative_charge_gain'])]))
        amplitudes.append(mp.mpf(str(row['actual_drive_amplitude'])))
    fits = [project.projection(combined, degree, amplitudes, mp.mpf('.001')) for degree in [4, 5]]
    assert fits[-1]['residual_norm'] < 1e-70
    diagnostics = {}
    series = {n:json.loads((run.OUT/f'background/p4-{n}.json').read_text())['history'] for n in [16,32,64]}
    other = json.loads((run.OUT/'background/p2-64.json').read_text())['history']
    for field in ['velocity_mass_RMS_m_s', 'scalar_mass_RMS', 'temperature_mass_RMS']:
        arrays = [np.array([r[field] for r in series[n]]) for n in [16,32,64]]
        norm = max(abs(arrays[-1]))
        first = max(abs(arrays[0]-arrays[1][::2]))/norm
        last = max(abs(arrays[1]-arrays[2][::2]))/norm
        spatial = max(abs(np.array([r[field] for r in other])-arrays[-1]))/norm
        diagnostics[field] = dict(last=float(last), order=float(np.log2(first/last)), spatial=float(spatial))
    times = np.array([r['seconds'] for r in series[64]])
    temps = np.array([r['maximum_delta_lnT'] for r in series[64]])
    j = int(np.flatnonzero(temps>=.01)[0])
    value = dict(classification='Counterexample candidate', passed=True,
        pressure_coordinate_temperature_reconstruction_relative=score,
        maximum_original_orbital_residuals=maxima,
        orbital_complex_relative_change_from_old=changes,
        summed_thermal_charge_amplitudes=sum(r['delta_alpha_radiative_over_phi0'] for r in orbital['cases']['p4']),
        degree4_full_residual=fits[0]['residual_norm'],
        previous_degree4_control_change=base['degree4_control_difference'],
        background_component_diagnostics=diagnostics,
        one_percent_crossing_bracket_seconds=[float(times[j-1]),float(times[j])],
        seconds=time.monotonic()-start,
        total_compute_seconds=result['total_compute_seconds']+time.monotonic()-start,
        superseded='Density-coordinate cp heat map and physical temperature/feedback interpretation in Phase81,85-88. Preserve raw results. Phase84 adiabatic comparison and the Phase86 same-state spectral coefficient comparison are unaffected by this correction.',
        full_four_component_convergence=False, full_physical_background=False,
        physical_atmosphere=False, full_goal_complete=False)
    assert value['total_compute_seconds'] < 600
    run.write(run.OUT/'audit.json', value)
    print(json.dumps(value), flush=True)


if __name__ == '__main__':
    main()
