"""Full native-radius Rosseland transport coupled to the orbital GR response.

Counterexample candidate: one-temperature, frozen-coefficient transport only.
Rosseland is used for flux transport, never as a Planck absorption mean. The
optically thin layers and omitted atmosphere prevent physical certification.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
import def_orbital_conductive_feedback as old
import def_free_surface_thermal as native

OUT = old.OUT.parent / 'def-orbital-radiative-transport'
ORIGINAL_BANK = old.go.task.BANK
write = old.write


def inputs():
    d = np.load(native.OUT / 'coefficients.npz')
    b = np.load(ORIGINAL_BANK / 'fine-bank.npz')
    faces = np.arange(1, len(d['lnT']))
    T = np.exp(d['lnT'])
    theta = d['A'] * d['N'] * T
    ell = 1 / (d['raw'][:, 0] * d['opacity'][:, 0])
    distance = -np.diff(d['radius_cm'])
    proper = np.sqrt(d['A'][:-1] * d['A'][1:] * d['metric'][:-1] * d['metric'][1:])
    H = distance * proper / abs(np.diff(np.log(theta)))
    kn = np.sqrt(ell[:-1] * ell[1:]) / H
    mode = np.zeros((len(faces), 13))
    rates = np.ones_like(mode)
    mode[:, 0] = np.sqrt(d['K'][:-1, 0] * d['K'][1:, 0]) / 1e5
    # This is a declared grey transport moment time, not an absorption rate.
    rates[:, 0] = 2.99792458e10 / np.sqrt(ell[:-1] * ell[1:])
    mode[b['faces'] - 1, 1:] = b['mode_K_SI']
    rates[b['faces'] - 1, 1:] = b['poles_proper_s']
    luminosity = d['luminosity'][:, 0]
    free = (4 * np.pi * d['faces_cm'][1:-1]**2 *
            d['N'][:-1] * d['N'][1:] * (d['A'][:-1] * d['A'][1:])**2 *
            4 * 5.670400e-5 * (T[:-1] * T[1:])**2)
    return d, b, faces, mode, rates, dict(
        knudsen=kn, mean_free_path_cm=ell, proper_temperature_scale_cm=H,
        background_L_over_cE=abs(luminosity) / free)


class Problem(old.Problem):
    def __init__(self, degree=4, label='all'):
        # Reuse the original conservative source, scalar lift and joint solve.
        # Only the independently saved transport bank is changed for assembly.
        previous = old.go.task.BANK
        old.go.task.BANK = OUT / label
        try:
            super().__init__(degree)
        finally:
            old.go.task.BANK = previous
        self.label = label

    def solve(self, n, label):
        previous = old.OUT
        old.OUT = OUT
        try:
            row = super().solve(n, label)
        finally:
            old.OUT = previous
        row.update(transport_bank=self.label,
                   full_native_internal_radiative_transport=self.label == 'all',
                   grey_LTE_diffusion_candidate=True,
                   full_radial_spectral_photons=False,
                   physical_atmosphere=False)
        # Independent reciprocal readout of the actual stored correction.
        saved = np.load(OUT / f'{label}-{n}.npz')
        z = np.clongdouble(-1j * n * self.omega)
        B = self.load_active - self.Kx @ self.H - z*z*(self.Mx @ self.H)
        direct = (self.kl + z*z*self.ml) @ saved['dq'] - self.fL @ saved['E']
        adjoint = -(B.T @ saved['qa'] + self.fL) @ saved['E']
        row['adjoint_readout_relative'] = float(abs(direct-adjoint)/max(abs(direct), 1e-100))
        assert row['adjoint_readout_relative'] < 1e-5, row['adjoint_readout_relative']
        write(OUT / f'{label}-{n}.json', row)
        return row


def prepare():
    assert not OUT.exists()
    OUT.mkdir()
    write(OUT / 'plan.json', dict(
        classification='Counterexample candidate', checkpoint='c8db24fa9',
        claim='Solve full native internal-radius grey LTE radiative transport plus the already computed microscopic conduction, actual orbital scalar drive and full GR temperature/energy/momentum feedback.',
        decision='Compare its outgoing charge with the known microscopic conduction and adiabatic comparator. A large change motivates missing photon physics; a small result cannot bound spectral/atmosphere errors. Compare all radiative faces with Kn<=0.01 components to locate sensitivity to the diffusion premise. Neither bank is a physical error envelope.',
        model='Use saved radiative Rosseland K=16 sigma T^3/(3 rho kappa_R) and grey flux-relaxation rate c rho kappa_R. This single-temperature moment model does not introduce a Planck absorption coefficient or duplicate LTE photon energy. Frozen geometry, coefficients and background; all native INTERNAL faces, omitted outer atmosphere flux. Retain existing microscopic electron conduction on its computed faces.',
        approximation='Nonuniform background is affine, not certified stationary. Kn is mean free path over proper Tolman-temperature scale. Small Kn alone does not establish absorption thermalization. Coefficient interpolation and native opacity physics remain conditional.',
        cases=['all-p4 harmonics1,2,3', 'all-p2 harmonics1,2,3', 'opaque-p4 harmonics1,2,3'],
        gates=dict(linear=1e-9, heat_law=1e-9, internal_balance=2e-13,
                   adjoint=1e-5, spatial=.02),
        budget=dict(pilot_seconds=120, total_seconds=600, CPU_threads=1,
                    memory_GB=4, new_EOS_calls=0, new_stellar_time_steps=0,
                    automatic_expansion=False),
        forecast='Prior joint setup10.3s, solve0.46s at52029dofs. Extra1722face unknowns and stiff radiative coupling are unmeasured; use one pilot before production. Do not automatically expand upon failure.',
        references=['https://arxiv.org/abs/astro-ph/0601635'],
        bindings={str(p): old.go.task.digest(p) for p in [Path(__file__), Path(old.__file__),
                  Path(old.thermal.__file__), native.OUT/'coefficients.npz',
                  ORIGINAL_BANK/'fine-bank.npz', old.OUT/'combined-projection.json']}))
    signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    d, b, faces, mode, rates, domain = inputs()
    for label in ['all', 'opaque']:
        path = OUT / label
        path.mkdir()
        selected = mode.copy()
        if label == 'opaque':
            selected[domain['knudsen'] > .01, 0] = 0
        use = selected.sum(1) > 0
        np.savez_compressed(path/'fine-bank.npz', faces=faces[use],
                            mode_K_SI=selected[use], poles_proper_s=rates[use],
                            unknown_faces=np.r_[0, faces[~use]])
    np.savez_compressed(OUT/'transport-domain.npz', faces=faces, **domain)
    write(OUT/'domain.json', dict(classification='Counterexample candidate',
        radiative_internal_faces=len(faces), microscopic_conductive_faces=len(b['faces']),
        maximum_Kn=float(max(domain['knudsen'])),
        faces_Kn_above_001=(faces[domain['knudsen'] > .01]).tolist(),
        maximum_background_L_over_cE=float(max(domain['background_L_over_cE'])),
        maximum_transport_time_seconds=float(max(domain['mean_free_path_cm'])/2.99792458e10),
        physical_diffusion_or_LTE_certificate=False))
    write(OUT/'symbolic.json', old.thermal.control())
    p = Problem()
    setup = time.monotonic()-start
    row = p.solve(1, 'all-p4')
    elapsed = time.monotonic()-start
    forecast = 1.5*(3*setup+9*row['seconds']+20)
    write(OUT/'pilot.json', dict(classification='Counterexample candidate', row=row,
        setup_seconds=setup, seconds=elapsed, forecast_seconds=forecast,
        maximum_omega_transport_time=float(p.omega/p.model.heat.geometry.tc *
            max(domain['mean_free_path_cm']/2.99792458e10/(d['A']*d['N']))),
        assumption='Three assemblies and nine solves at measured pilot rate with20s output allowance and50percent margin; p2/opaque rates are estimates.'))
    print('RADIATIVE PILOT', elapsed, forecast, flush=True)


def run():
    plan = json.loads((OUT/'plan.json').read_text())
    pilot = json.loads((OUT/'pilot.json').read_text())
    assert not (OUT/'result.json').exists() and pilot['forecast_seconds'] < 600
    for p, sha in plan['bindings'].items():
        assert old.go.task.digest(Path(p)) == sha, p
    signal.alarm(int(600-pilot['seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
    start = time.monotonic()
    cases = {}
    for label, degree, bank in [('all-p4', 4, 'all'), ('all-p2', 2, 'all'), ('opaque-p4', 4, 'opaque')]:
        p = Problem(degree, bank)
        cases[label] = [pilot['row'] if label == 'all-p4' and n == 1
                       else p.solve(n, label) for n in [1, 2, 3]]
        del p
        if label == 'all-p2':
            spatial = [abs(complex(*a['radiative_charge_gain'])/complex(*b['radiative_charge_gain'])-1)
                       for a, b in zip(cases[label], cases['all-p4'])]
            if max(spatial) >= .02:
                break
    contrast = []
    for i in range(3):
        ref = complex(*cases['all-p4'][i]['radiative_charge_gain'])
        row = dict(harmonic=i+1, spatial=spatial[i])
        if 'opaque-p4' in cases:
            core = complex(*cases['opaque-p4'][i]['radiative_charge_gain'])
            row['opaque_to_all_relative_difference'] = abs(core/ref-1)
        contrast.append(row)
    write(OUT/'result.json', dict(classification='Counterexample candidate',
        passed=max(spatial)<.02, comparisons=contrast, cases=cases,
        seconds=time.monotonic()-start,
        total_compute_seconds=time.monotonic()-start+pilot['seconds'],
        full_native_internal_grey_radiation_coupled=True,
        full_spectral_radial_photons=False, stationary_background_certified=False,
        physical_atmosphere=False, full_goal_complete=False))
    print('RADIATIVE RESULT', json.dumps(contrast), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'run'])
    globals()[parser.parse_args().action]()
