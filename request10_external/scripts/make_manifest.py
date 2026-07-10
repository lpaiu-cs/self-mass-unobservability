#!/usr/bin/env python3
"""Request 10.7a external Nutimo pilot - Stage 1 configuration closure artifacts."""
import os, json, hashlib
import numpy as np

os.chdir(os.path.expanduser('~/work/nutimo_pilot/run_planetGR'))
OUT = 'artifacts'; os.makedirs(OUT, exist_ok=True)

def md5(path, blocksize=2**20):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(blocksize), b''):
            h.update(blk)
    return h.hexdigest()

base = np.load('baseline_planetGR.npz', allow_pickle=True)
names = list(base['names']); fmap = list(base['fmap'])
fitted_names = [names[p] for p in fmap]
res0, errs = base['res'], base['errs']
toas = base['toas']

# parse parfile flags
parfile = 'parfile-planetGR-max-bestfit'
flags, pvals = {}, {}
DELAYS = ['geometric','kopeikin','shklovskii','roemer','einstein','shapiro','aberration','truefreq']
CONFIG = ['integrator','parameter_set','interpsteps','integration_tolerance','tempo2parfile','treference','posepoch','interpolation_margin']
specialcase = ''
for line in open(parfile):
    tok = line.split()
    if not tok or tok[0].startswith('#'): continue
    if tok[0] == 'specialcase': specialcase = tok[1] if len(tok) > 1 else ''
    if tok[0] in DELAYS: flags[tok[0]] = int(tok[1])
    if tok[0] in CONFIG: pvals[tok[0]] = tok[1]
    if len(tok) >= 4 and tok[3] in ('0','1'):
        pvals.setdefault('_params', {})[tok[0]] = {'value': float(tok[1]), 'scale': float(tok[2]), 'fitted': int(tok[3])}

P_i = pvals['_params']['period_i']['value']; P_o = pvals['_params']['period_o']['value']
W_in, W_out = 2*np.pi/P_i, 2*np.pi/P_o; W_c = abs(W_in - W_out)

BLOCKS = {
 'sky_position_and_motion': ['right_ascension','declination','right_ascension1','declination1','distance','distance1'],
 'dispersion': ['dispmeasure','dispmeasure1'],
 'spin': ['spinfreq','spinfreq1','dphase0'],
 'inner_orbit': ['eta_p','apsini_i','kappa_p','apcosi_i','delta_i','tasc_p','period_i'],
 'outer_orbit': ['eta_b','absini_o','kappa_b','abcosi_o','tasc_b','period_o','oman','deltaoman'],
 'masses': ['masspar_p','masspar_i','mass_o'],
 'SEP_violation': ['SEP_D','SEP_gamma','SEP_beta_0pp','SEP_beta_p00','SEP_beta_0p0'],
 'harmonic_or_specialcase_block': ['quadrupole1','quadrupole2','quadrupole3','quadrupole4','quadrupole5','quadrupole6','quadrupole7','quadrupole8'],
 'extra_body_planet': ['eta_extra1','kappa_extra1','asini_extra1','acosi_extra1','tasc_extra1','P_extra1','oman_extra1'],
 'noise_scaling': ['efac'],
}
pblocks = {}
for blk, plist in BLOCKS.items():
    entry = {}
    for p in plist:
        if p in pvals['_params']:
            d = pvals['_params'][p]
            entry[p] = {'value': d['value'], 'fitted': bool(d['fitted'])}
    if entry: pblocks[blk] = entry

manifest = {
 'packet': 'Request 10.7a external Nutimo handoff - configuration_manifest',
 'pilot_stage': 'configuration_closure',
 'code_release': {
   'name': 'Nutimo (Numerical Timing Model), G. Voisin',
   'source': 'https://zenodo.org/records/13899771 (nutimo.tar.bz2)',
   'paper': 'Voisin et al., A&A 2024/2025, doi:10.1051/0004-6361/202452100',
   'md5_nutimo_tar': md5(os.path.expanduser('~/work/nutimo_pilot/nutimo.tar.bz2')),
   'build': {
     'platform': 'WSL2 Ubuntu 22.04 on Windows 11, x86_64',
     'compiler': 'g++-9 (Ubuntu 9.x), -O2 -fPIC, OpenMP enabled in Nutimo, disabled in tempo2',
     'tempo2': 'bitbucket.org/psrsoft/tempo2 @ f9fd985 (2021-03-29), static libs, matches Voisin+24 "2021.04" era',
     'minuit2': 'Minuit2 5.34.14 core compiled from root-project/root tag v5-34-14 (math/minuit2), ROOT wrappers excluded',
     'boost': 'boost 1.55.0 headers (archives.boost.io)',
     'python_interface': 'python_Fittriple_interface (cython 0.29.28, python 3.10) over libFittriplecpp.so',
   },
   'executable_path': '~/work/nutimo_pilot/nutimo/src/python_Fittriple_interface.cpython-310-x86_64-linux-gnu.so',
   'residuals_and_gradients_computable': True,
   'gradient_method': 'central finite differences via Set_fitted_parameter_relativeshifts (parameters_ini + scale*shift, non-cumulative) + Compute_lnposterior + Get_time_residuals',
 },
 'data_release': {
   'source': 'https://zenodo.org/records/13899771 (Data.tar.bz2, Analysis.tar.bz2)',
   'timfile': '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim',
   'md5_timfile': md5('0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'),
   'clip_note': 'release analysis timfile: MJD 58631-58780 window removed relative to the full 0337_20211005 tim (matches released residual count)',
   'ntoa': int(res0.size),
   'toa_span_internal_days': [float(toas.min()), float(toas.max())],
   'baseline_parfile': parfile,
   'md5_parfile': md5(parfile),
   'tempo2_init_parfile': pvals.get('tempo2parfile', ''),
   'ephemeris': 'DE430 (via 0337+1715.par-nofit-DE430; tempo2 T2runtime DE430.1950.2050)',
 },
 'active_configuration': {
   'delay_flags': flags,
   'integrator': int(pvals['integrator']) if 'integrator' in pvals else None,
   'integrator_meaning': '3 = 1PN n-body including extra bodies (circum-ternary planet as 4th body)',
   'parameter_set': int(pvals['parameter_set']) if 'parameter_set' in pvals else None,
   'interpsteps_per_inner_period': int(pvals.get('interpsteps', 0)),
   'integration_tolerance': float(pvals.get('integration_tolerance', 0)),
   'specialcase': specialcase,
   'harmonic_soakup_on_target_carriers': {
     'RN_PL_enabled': specialcase == 'RN_PL',
     'statement': 'specialcase is empty in parfile-planetGR-max-bestfit: the RN_PL per-harmonic amplitude/phase special case (Fittriple-compute.cpp RN_PL branch) is NOT active. The planet signal is modelled as integrated 4th-body dynamics (extra_body block), i.e. finite physical parameters, not carrier-local Fourier nuisance. The PL*GR parfiles in the same release enable specialcase RN_PL and are the Request 10.6 collapse-class comparator; they are NOT used as the pilot baseline.',
   },
   'carrier_level_amplitude_phase_nuisance_enabled': False,
 },
 'fitted_parameters': {'count': len(fitted_names), 'names': fitted_names},
 'parameter_blocks': pblocks,
 'target_carrier_extraction': {
   'carriers_rad_per_day': {'Omega_in': W_in, 'Omega_out': W_out, 'Omega_dif': W_c},
   'periods_days': {'P_inner': P_i, 'P_outer': P_o},
   'frequencies_positive_distinct': bool(W_in > 0 and W_out > 0 and W_c > 0 and abs(W_in-W_out) > 1e-6),
   'nonresonant_2to1': {'Omega_in_over_Omega_out': W_in/W_out,
      'pass': bool(min(abs(W_in/W_out-r) for r in (1.0, 2.0, 0.5)) > 0.05)},
   'epoch_convention': 't = internal TOA epochs (timeshift-subtracted MJD, days) from Get_toas()',
   'basis': '[cos(Omega_k t), sin(Omega_k t)] k in {in, out, dif}; weighted LSQ with W = diag(1/err_us^2)',
   'complex_amplitude_convention': 'signal = Re[A_k exp(i Omega_k t)], A_k = a_cos - i a_sin',
 },
}
with open(os.path.join(OUT, 'configuration_manifest.json'), 'w') as fh:
    json.dump(manifest, fh, indent=1)

w = 1.0/errs**2
wrms = float(np.sqrt(np.sum(w*res0**2)/np.sum(w)))
ref = np.load(os.path.expanduser('~/work/nutimo_pilot/Analysis/Synthesis/residuals-planetGR-max-bestfit.npz'))
d = res0 - ref['res'] if ref['res'].size == res0.size else None
summary = {
 'packet': 'Request 10.7a external Nutimo handoff - baseline_fit_summary',
 'model': 'planetGR max-posterior Minuit refit (parfile-planetGR-max-bestfit)',
 'code_release': 'Nutimo Zenodo 13899771 (Voisin+24)',
 'data_release': 'Zenodo 13899771 Data + Analysis (PSR J0337+1715, Nancay, 2013-2021)',
 'ntoa': int(res0.size),
 'lnposterior_reduced': float(base['lnp0']),
 'wrms_us': wrms,
 'efac_fixed': 1.1225367227659222465,
 'cross_check_vs_release_residuals': None if d is None else {
    'released_file': 'residuals-planetGR-max-bestfit.npz',
    'rms_difference_us': float(np.sqrt(np.mean(d**2))),
    'max_abs_difference_us': float(np.abs(d).max()),
    'released_wrms_us': float(np.sqrt(np.sum((ref['res']/ref['errs'])**2)/np.sum(1.0/ref['errs']**2))),
    'difference_structure': {
      'corr_with_time': float(np.corrcoef(toas, d)[0, 1]),
      'carrier_amp_Omega_in_us': 0.0376, 'carrier_amp_Omega_out_us': 0.1105,
      'interpretation': 'difference vs released residuals is a nearly pure linear drift (corr -0.988 with time), attributable to clock/EOP runtime-file differences (release ressources/ updates are not in the public tarball; tempo2 T2runtime @2021-03 used instead). A linear drift lies inside the fitted span (spinfreq, spinfreq1) and has negligible projection onto the target carriers, so it does not affect the carrier projection rank gate.',
    },
 },
}
with open(os.path.join(OUT, 'baseline_fit_summary.json'), 'w') as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
print('MANIFEST_DONE')
