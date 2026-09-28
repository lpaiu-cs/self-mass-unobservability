"""Compare the solved grey channel and its actual outer spectral transport."""
from pathlib import Path
import argparse
import inspect
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_orbital_radiative_transport as run


def prepare():
    result = json.loads((run.OUT/'result.json').read_text())
    assert result['passed'] and not (run.OUT/'audit-plan.json').exists()
    remaining = 600-result['total_compute_seconds']
    assert remaining > 60
    run.write(run.OUT/'audit-plan.json', dict(classification='Counterexample candidate',
        claim='Project the already solved grey radiative addition into the same declared derivative comparator, and compare its outer Rosseland conductivity with the zero-frequency first angular moment of the saved actual spectral collision operator.',
        decision='If microscopic and native grey conductivities differ, retain their difference as a model discrepancy; do not renormalize or extend one outer state throughout the star. Neither a small grey charge nor a local conductivity match certifies the missing atmosphere or stationary background.',
        formula='K_spectral=c^2/3 * sqrt(C_gamma)^T C_l1^{-1} sqrt(C_gamma). Retain full absorption and finite frequency redistribution in the l=1 collision operator. No Rosseland-to-Planck substitution.',
        budget=dict(hard_seconds=min(90, int(remaining)), CPU_threads=1, memory_GB=4,
                    new_EOS_calls=0, new_frequency_cells=0),
        gates=dict(spectral_solve_residual=1e-10, positive_conductivity=True),
        bindings={str(p):run.old.go.task.digest(p) for p in [Path(__file__),
            run.OUT/'result.json', run.OUT/'plan.json']}))


def project(result):
    import def_orbital_charge_audit as check
    check.mp.mp.dps = 70
    prior = json.loads((run.old.orbit.OUT.parent/'def-full-orbital-exterior/comparator-result.json').read_text())
    baseline = prior['cases']['p4-tight']
    bycase = {}
    residuals = {}
    for label, rows in result['cases'].items():
        combined = [baseline[0]]
        for a, b in zip(baseline[1:], rows):
            values = [check.mp.mpf(str(x))+check.mp.mpf(str(y))
                      for x, y in zip(a['charge'], b['radiative_charge_gain'])]
            combined.append(dict(charge=[check.mp.nstr(x, 60) for x in values]))
        amplitudes = [check.mp.mpf(0)]+[check.mp.mpf(str(r['actual_drive_amplitude'])) for r in rows]
        fits = [check.projection(combined, degree, amplitudes, check.mp.mpf('.001'))
                for degree in [2, 3, 4, 5]]
        residuals[label] = np.array(fits[2]['residual'])
        bycase[label] = dict(fits=fits, summed_thermal_amplitudes=sum(r['delta_alpha_radiative_over_phi0'] for r in rows))
    ref = residuals['all-p4']
    return dict(classification='Counterexample candidate', cases=bycase,
        degree4_thermal_spatial_change=float(np.linalg.norm(ref-residuals['all-p2'])),
        degree4_opaque_component_change=float(np.linalg.norm(ref-residuals['opaque-p4'])),
        previous_adiabatic_control_change=prior['degree4_control_difference'],
        full_goal_complete=False,
        scope='Same six charge-coefficient norm, not an observational likelihood. Native-grey and opaque-component contrasts are not physical error enclosures.')


def spectral():
    import def_photon_radial_gr as radial
    ph = radial.ph
    b, d, _, _ = radial.geometry()
    b = dict(b)
    th = np.load(ph.OUT/'thermo.npz')
    b.update(Cm=th['Cf'], velocity_order=256, q_points=257)
    source = inspect.getsource(ph.photon_operator).replace('Operator(b,8,1,24)', 'Operator(b,2,0,24)')
    ns = dict(vars(ph))
    exec(compile(source, 'spectral-transport-source.py', 'exec'), ns)
    op = ns['photon_operator'](b)
    matrix = (diags(op.L.diagonal().reshape(len(op.u), 2)[:, 1])+op.gain[1]).tocsc()
    root = np.sqrt(op.Ci)
    factor = splu(matrix)
    answer = factor.solve(root)
    residual = np.linalg.norm(matrix@answer-root)/np.linalg.norm(root)
    conductivity = float(b['c'])**2/3*(root@answer)
    assert residual < 1e-10 and conductivity > 0
    np.savez_compressed(run.OUT/'spectral-transport.npz', photon_heat_capacity=op.Ci,
                        moment_response=answer, frequency=op.u)
    return dict(classification='Counterexample candidate', frequency_cells=len(op.u),
        temperature_K=float(b['T']), density_g_cm3=float(b['material_rho']),
        K_spectral_cgs=float(conductivity), K_native_Rosseland_cgs=float(d['K'][0, 0]),
        spectral_over_native=float(conductivity/d['K'][0, 0]), residual=float(residual),
        matrix_nonzeros=matrix.nnz, factor_nonzeros=factor.L.nnz+factor.U.nnz,
        LTE_capacity_relative=float(sum(op.Ci)/(4*float(b['arad'])*float(b['T'])**3)-1),
        scope='One saved actual outer EOS state and retained spectral atomic/RPA model. No arbitrary multiplicative correction, full radial atomic population matching or absorption-thermalization certificate.')


def audit():
    plan = json.loads((run.OUT/'audit-plan.json').read_text())
    assert not (run.OUT/'audit.json').exists()
    for p, h in plan['bindings'].items():
        assert run.old.go.task.digest(Path(p)) == h, p
    signal.alarm(plan['budget']['hard_seconds'])
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    result = json.loads((run.OUT/'result.json').read_text())
    projection = project(result)
    run.write(run.OUT/'combined-projection.json', projection)
    local = spectral()
    seconds = time.monotonic()-start
    run.write(run.OUT/'audit.json', dict(classification='Counterexample candidate', passed=True,
        spectral_transport=local, seconds=seconds,
        total_compute_seconds_charged=seconds+result['total_compute_seconds'],
        full_goal_complete=False))
    print('SPECTRAL TRANSPORT', json.dumps(local), flush=True)
    print('PROJECTED', json.dumps(projection), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'audit'])
    globals()[parser.parse_args().action]()
