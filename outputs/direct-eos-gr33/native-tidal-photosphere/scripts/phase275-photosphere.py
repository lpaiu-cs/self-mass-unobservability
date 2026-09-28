"""Phase275: charge magnitude for photospheres consistent with the observed log g and Teff (Kaplan et al. 2014).

Conjectural: gray atmosphere, T(tau) proportional to Teff; hydrostatic similarity rho'(z) = s_rho rho(z/s_z),
s_z = (g/g')(Teff'/Teff), s_rho = (g'/g)^a (Teff'/Teff)^b with (a, b) = (1/2, 1.25) (Kramers-type opacity) or (1, -1)
(density-independent opacity); the face kernel of phase 272 is kept (q' = sum_f c_f rho'(z_f)/rho0(z_f)).
Usage: python3 phase275-photosphere.py <kernel json> <out json>
"""
import json, sys
import numpy as np
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, SIG = 2.99792458e10, 5.670374419e-5
k = json.load(open(sys.argv[1])); faces = k['faces']
c = np.array([f['contribution'] for f in faces]); z = np.array([f['depth_km'] for f in faces]); rho0 = np.array([f['rho0'] for f in faces]); q = c.sum()
bg = bl.chem.prior.Background(); d = bg.d
r, mg, rho = (np.asarray(d[x], float) for x in ['radius_cm', 'mass_geom_cm', 'density_cgs']); o = np.argsort(r); r, mg, rho = r[o], mg[o], rho[o]
R = float(bg.R); L = float(d['base_luminosity']); depth = (R - r)/1e5
Teff = (L/(4*np.pi*R**2*SIG))**0.25; logg = float(np.log10(C*C*mg[-1]/R**2))
zz, lr = depth[::-1], np.log(rho[::-1])  # depth increasing
lnrho = lambda x: np.interp(x, zz, lr)
dev = np.exp(lnrho(z))/rho0 - 1; Kf = c/q
ident = float(np.max(np.abs(dev))); ident_weighted = float(np.sum(np.abs(Kf*dev)))  # face sampler vs log-linear array interpolation


def charge(g_ratio, t_ratio, a, b, sz_extra=1.):
    sz = t_ratio/g_ratio*sz_extra; srho = g_ratio**a*t_ratio**b
    return float(np.sum(c*srho*np.exp(lnrho(z/sz) - lnrho(z))))  # density change factor from the same interpolant


cases = {}
OBS = dict(logg=5.82, dlogg=0.05, Teff=15800., dTeff=100.)
for opac, (a, b) in {'kramers (a=1/2, b=1.25)': (0.5, 1.25), 'constant opacity (a=1, b=-1)': (1., -1.)}.items():
    rows = {}
    for name, lg, te in [('declared atmosphere', logg, Teff), ('Kaplan central', OBS['logg'], OBS['Teff']),
                         ('log g -1 sigma', OBS['logg'] - OBS['dlogg'], OBS['Teff']), ('log g +1 sigma', OBS['logg'] + OBS['dlogg'], OBS['Teff']),
                         ('log g -2 sigma', OBS['logg'] - 2*OBS['dlogg'], OBS['Teff']), ('log g +2 sigma', OBS['logg'] + 2*OBS['dlogg'], OBS['Teff']),
                         ('Teff -1 sigma', OBS['logg'], OBS['Teff'] - OBS['dTeff']), ('Teff +1 sigma', OBS['logg'], OBS['Teff'] + OBS['dTeff'])]:
        v = charge(10**(lg - logg), te/Teff, a, b); rows[name] = dict(logg=lg, Teff=te, charge=v, ratio=v/q, sign_kept=bool(np.sign(v) == np.sign(q)))
    for s in [-0.1, 0.1]:
        v = charge(10**(OBS['logg'] - logg), OBS['Teff']/Teff, a, b, 1 + s); rows[f'Kaplan central, gray T(tau) {s:+.0%}'] = dict(charge=v, ratio=v/q, sign_kept=bool(np.sign(v) == np.sign(q)))
    cases[opac] = rows
assert all(abs(rows['declared atmosphere']['ratio'] - 1) < 1e-12 for rows in cases.values())
allr = [v['ratio'] for rows in cases.values() for n, v in rows.items() if n != 'declared atmosphere' and '2 sigma' not in n]
res = dict(declared=dict(Teff=Teff, logg=logg, R_sun=R/6.957e10), observed=OBS, face_background_max_deviation=ident, face_background_kernel_weighted_deviation=ident_weighted, kernel_charge=float(q), cases=cases,
           one_sigma_ratio_range=[min(allr), max(allr)], sign_kept_everywhere=all(v['sign_kept'] for rows in cases.values() for v in rows.values()))
json.dump(res, open(sys.argv[2], 'w'), indent=1)
print('declared Teff %.1f K log g %.4f; kernel charge %.6e; face/background deviation max %.1e, kernel-weighted %.1e' % (Teff, logg, q, ident, ident_weighted))
for opac, rows in cases.items():
    print(opac)
    for n, v in rows.items(): print('  %-34s charge %.4e  ratio %.3f  sign kept %s' % (n, v['charge'], v['ratio'], v['sign_kept']))
print('1-sigma ratio range %.3f .. %.3f ; sign kept everywhere %s' % (res['one_sigma_ratio_range'][0], res['one_sigma_ratio_range'][1], res['sign_kept_everywhere']))
