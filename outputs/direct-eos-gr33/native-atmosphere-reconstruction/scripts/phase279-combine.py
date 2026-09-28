"""Phase279 part 3: combine the spherical table-opacity gray atmospheres (part 1) with the non-gray LTE density correction (part 2)
through the phase-272 face kernel, and record the photospheric flux constancy of the non-gray models.

Conjectural: physical density at depth z below the P_gas = 1 cut: rho_case(z) = rho_tab,case(z) * r_NG,case(z), where rho_tab is the
spherical gray envelope with the declared EOS/tables (validated on the declared background) and r_NG = rho_nongray/rho_gray of the
plane-parallel LTE H+He model at the same (Teff, log g). Charge: q = sum_f c_f rho_case(z_f)/rho_tab,declared(z_f).
Usage: python phase279-combine.py <kernel json> <gray-sph json> <gray-sph npz> <nongray json> <nongray npz prefix> <out json>
"""
import json, sys
import numpy as np
kern, gj, gz, nj, npre, outp = sys.argv[1:7]
k = json.load(open(kern)); c = np.array([f['contribution'] for f in k['faces']]); zf = np.array([f['depth_km'] for f in k['faces']])*1e5; q = c.sum()
G = json.load(open(gj)); A = np.load(gz); N = json.load(open(nj))
lr = lambda z, zz, rr: np.interp(z, zz, np.log(rr))
key = lambda name: name.replace(' ', '_')
dec_tab = lr(zf, A[f'{key("declared")}__z'], A[f'{key("declared")}__rho'])
sys.path.insert(0, '.')
NG = {}; exec(open('phase279-nongray.py', encoding='utf-8').read().split('# self-checks')[0], NG)  # separate namespace (it defines c, k, h)
out = dict(kernel_charge=float(q), gray_validation=G['validation'], cases={})
for name in G['cases']:
    tab = lr(zf, A[f'{key(name)}__z'], A[f'{key(name)}__rho'])
    Z = np.load(f'{npre}-{key(name)}.npz')
    rng = lr(zf, Z['nongray_z'], Z['nongray_rho']) - lr(zf, Z['gray_z'], Z['gray_rho'])
    q_gray = float(np.sum(c*np.exp(tab - dec_tab))); q_ng = float(np.sum(c*np.exp(tab - dec_tab + rng)))
    # photospheric flux constancy of the stored non-gray model (tau_R <= 100)
    m, T, rho, tauR = Z['nongray_m'], Z['nongray_T'], Z['nongray_rho'], Z['nongray_tauR']
    ab, sc = NG['opacity'](rho, T); J, Hh, B = NG['solve_J'](T, ab, sc, m)
    Hn = np.r_[J[0:1]*0.5, 0.5*(Hh[1:] + Hh[:-1]), Hh[-1:]]; H0 = NG['SIG']*G['cases'][name]['Teff']**4/(4*np.pi)
    ferr = float(np.max(np.abs(np.sum(Hn*NG['w_nu'], 1)[tauR <= 100]/H0 - 1)))
    out['cases'][name] = dict(Teff=G['cases'][name]['Teff'], logg=G['cases'][name]['logg'], gray_charge=q_gray, gray_ratio=q_gray/q,
                              nongray_charge=q_ng, nongray_ratio=q_ng/q, nongray_factor=q_ng/q_gray, sign_kept=bool(np.sign(q_ng) == np.sign(q) and np.sign(q_gray) == np.sign(q)),
                              nongray_flux_error_tau_le_100=ferr, T0_over_Teff_nongray=N['cases'][name]['T0_over_Teff_nongray'])
    print('%-16s gray x%.3f  non-gray x%.3f (factor %.3f)  flux err(tau<=100) %.1e' % (name, q_gray/q, q_ng/q, q_ng/q_gray, ferr), flush=True)
g1 = [v['gray_ratio'] for n_, v in out['cases'].items() if n_ != 'declared' and '2 sigma' not in n_]
n1 = [v['nongray_ratio'] for n_, v in out['cases'].items() if n_ != 'declared' and '2 sigma' not in n_]
out['one_sigma'] = dict(gray=[min(g1), max(g1)], nongray=[min(n1), max(n1)])
out['two_sigma'] = dict(gray=[out['cases']['log g -2 sigma']['gray_ratio'], out['cases']['log g +2 sigma']['gray_ratio']],
                        nongray=[out['cases']['log g -2 sigma']['nongray_ratio'], out['cases']['log g +2 sigma']['nongray_ratio']])
out['sign_kept_everywhere'] = all(v['sign_kept'] for v in out['cases'].values())
json.dump(out, open(outp, 'w'), indent=1)
print(json.dumps({k_: out[k_] for k_ in ['one_sigma', 'two_sigma', 'sign_kept_everywhere']}))
