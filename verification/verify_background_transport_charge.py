"""Original-equation and output-change audit for the transport snapshot update."""
from pathlib import Path
import hashlib
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
import def_background_transport_charge as run
import def_orbital_charge_audit as projection


def main():
    assert not (run.OUT/'audit.json').exists()
    signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    plan = json.loads((run.OUT/'plan.json').read_text())
    for name, sha in plan['bindings'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == sha, name
    result = json.loads((run.OUT/'result.json').read_text())
    p = run.Problem()
    m = p.model
    d = m.heat.d
    state = np.load(run.drift.OUT/'p4-64.npz')
    data = np.load(run.OUT/'p4-64-input.npz')
    r = m.original.native
    point = p.point
    b = 1-2*point['m']/r
    alpha = -4*point['phi']
    A4 = np.exp(-8*point['phi']**2)
    dlz = r*r*point['v']**2/2-4*np.pi*r*r*A4*point['p']/b-point['m']/(r*b)
    slots = np.minimum(np.arange(len(r)), len(r)-2)
    slope = np.diff(np.log(d['raw'][::-1, 0]))[slots]/np.diff(r)[slots]
    # Direct baryon/metric expression, independently of the temperature map.
    Rq = -diags(r)@p.D[0]-diags(3+dlz+3*alpha*r*point['v']+r*slope)@p.V[0]
    Rq -= diags(r*point['v']+3*alpha)@p.V[1]
    _, _, J = run.old.old.go.task.fem.source_points(m.heat, point)
    full = np.zeros(len(m.heat.edges), np.longdouble)
    full[p.active] = state['E']
    density = (Rq@state['q']-diags(1/(r*b))@J@full).real[::-1]
    density_error = float(max(abs(density-data['delta_lnrho']))/max(max(abs(density)), 1e-100))
    assert density_error < 1e-10
    p.update(data)
    maxima = dict(GR=0., heat=0., adjoint=0.)
    for n in [1, 2, 3]:
        saved = np.load(run.OUT/f'p4-64-{n}.npz')
        qa, dq, E = [saved[k] for k in ['qa', 'dq', 'E']]
        z = np.clongdouble(-1j*n*p.omega)
        B = p.original_force-z*z*(p.Mx@p.H)
        defect = p.Kx@dq+z*z*(p.Mx@dq)-B@E
        error = float(np.max(abs(defect)/(p.absK@abs(dq)+abs(z*z)*(p.absM@abs(dq))+abs(B)@abs(E)+1e-100)))
        gain = np.sum(p.conductance*m.heat.geometry.tc*p.lam/(z*(z+p.lam)), axis=1)
        theta = p.Gq@(qa+dq)+p.gL+p.GEraw@E
        heat_error = float(max(abs(E-gain*theta))/max(abs(E)))
        direct = (p.kl+z*z*p.ml)@dq-p.fL@E
        reciprocal = -(B.T@qa+p.fL)@E
        adjoint_error = float(abs(direct-reciprocal)/abs(direct))
        for name, value in [('GR', error), ('heat', heat_error), ('adjoint', adjoint_error)]:
            maxima[name] = max(maxima[name], value)
        assert error < 1e-9 and heat_error < 1e-9 and adjoint_error < 1e-5
    projection.mp.mp.dps = 70
    mp = projection.mp
    prior = json.loads((run.old.OUT/'result.json').read_text())['cases']['all-p4']
    baseline = json.loads((run.old.old.orbit.OUT.parent/'def-full-orbital-exterior/comparator-result.json').read_text())
    adiabatic = baseline['cases']['p4-tight']
    combined = [adiabatic[0]]
    change = [dict(charge=[0, 0])]
    amps = [mp.mpf(0)]
    for a, old, new in zip(adiabatic[1:], prior, result['cases']['p4-64']):
        combined.append(dict(charge=[mp.nstr(mp.mpf(str(x))+mp.mpf(str(y)), 65)
                                     for x, y in zip(a['charge'], new['radiative_charge_gain'])]))
        change.append(dict(charge=[mp.nstr(mp.mpf(str(x))-mp.mpf(str(y)), 65)
                                   for x, y in zip(new['radiative_charge_gain'], old['radiative_charge_gain'])]))
        amps.append(mp.mpf(str(new['actual_drive_amplitude'])))
    whole = projection.projection(combined, 4, amps, mp.mpf('.001'))
    projected_change = projection.projection(change, 4, amps, mp.mpf('.001'))
    fifth = projection.projection(combined, 5, amps, mp.mpf('.001'))
    assert fifth['residual_norm'] < 1e-70
    temp = data['delta_lnT']
    rho = data['delta_lnrho']
    Kratio = np.exp(3*temp-rho)*d['opacity'][:, 0]/data['opacity']
    value = dict(classification='Counterexample candidate', audit_passed=True,
        change_numerically_resolved=result['passed'], density_reconstruction_relative=density_error,
        maximum_original_residuals=maxima,
        maximum_delta_lnT=float(max(abs(temp))), maximum_delta_lnrho=float(max(abs(rho))),
        maximum_relative_radiative_K_change=float(max(abs(Kratio-1))),
        outer_relative_radiative_K_change=float(Kratio[0]-1),
        summed_charge_change=sum(r['delta_alpha_radiative_over_phi0'] for r in result['comparisons']),
        degree4_full_residual=whole['residual_norm'], degree4_update=projected_change['residual_norm'],
        previous_degree4_control_change=baseline['degree4_control_difference'],
        seconds=time.monotonic()-start,
        total_compute_seconds=result['total_compute_seconds']+time.monotonic()-start,
        full_background_feedback=False, physical_atmosphere=False, full_goal_complete=False,
        scope='Actual finite update of the radiative transport input and temperature factor only. Original fixed GR/EOS tangent and omitted atmosphere remain. Finite contrasts are not a rigorous error envelope.')
    assert value['total_compute_seconds'] < 600
    run.write(run.OUT/'audit.json', value)
    print(json.dumps(value), flush=True)


if __name__ == '__main__':
    main()
