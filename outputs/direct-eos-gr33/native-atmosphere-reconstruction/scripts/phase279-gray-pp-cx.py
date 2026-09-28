"""Phase279 part 1: plane-parallel gray Eddington atmospheres with the declared native EOS and Rosseland tables.

Conjectural: T^4 = (3/4) Teff^4 (tau + 2/3), dP_tot/dtau = cx g/kappa_R (gravitating density e/c^2 = cx rho, cx = 1.00702), P_gas = P_tot - a T^4/3, dz = dtau/(kappa rho); the depth origin
is the declared envelope's outer cut P_gas = 1 dyn/cm^2 (tau ~ 8e-6), not tau = 0,
the same EOS call (total pressure), opacity parts and composition as the declared envelope (def_native_radiative_envelope.Envelope).
Validation: the declared parameters must reproduce the declared background rho(z) in the charge layers (240-465 km) within 2%.
Charges: face kernel of phase 272 times the density ratio rho_case(z)/rho_declared(z) of the same integrator.
Usage: python3 phase279-gray.py <kernel json> <out json>
"""
import json, sys, time
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, 'verification')
import def_native_radiative_envelope as env
import gr_two_carrier_evolution as two
import def_native_boundary_layer as bl
C = 2.99792458e10
t0 = time.monotonic(); E = env.Envelope(); X = E.X; ARAD = E.arad; SIG = ARAD*C/4
calls = [0]


def rho_kap(Ptot, T):
    calls[0] += 1
    a = E.eos(1, float(np.log(Ptot)), float(np.log(T)), X); rho = float(a[0])
    return rho, float(two.opacity_parts(E.opacity, (np.log(rho), np.log(T), X))[0])


def atmosphere(Teff, logg, tau_min=1e-9, tau_max=300., n=700):
    g = E.cx*10**logg; F = SIG*Teff**4  # the declared hydrostatics uses the rest energy cx rho as gravitating density
    T_of = lambda tau: (0.75*Teff**4*(tau + 2/3))**0.25
    T0 = T_of(tau_min); Pr0 = ARAD*T0**4/3; Pg = 1e-6
    for _ in range(60):  # P_gas(tau_min) = (g - kappa F/c) tau_min/kappa at the top
        rho, kap = rho_kap(Pg + Pr0, T0); new = max((g - kap*F/C)*tau_min/kap, 1e-30)
        if abs(new/Pg - 1) < 1e-12: break
        Pg = new
    def rhs(s, y):
        tau = np.exp(s); T = T_of(tau); rho, kap = rho_kap(y[0], T)
        return [tau*g/kap, tau/(kap*rho)]
    s = np.linspace(np.log(tau_min), np.log(tau_max), n)
    sol = solve_ivp(rhs, (s[0], s[-1]), [Pg + Pr0, 0.], t_eval=s, rtol=1e-9, atol=[1e-30, 1e-3], method='LSODA')
    assert sol.success, sol.message
    tau = np.exp(s); T = T_of(tau); P = sol.y[0]; z = sol.y[1]
    rho = np.array([rho_kap(p, t)[0] for p, t in zip(P, T)])
    Pg = P - ARAD*T**4/3; assert Pg[0] < 1 < Pg[-1]
    z0 = float(np.interp(0., np.log(Pg), z)); tau_cut = float(np.exp(np.interp(0., np.log(Pg), np.log(tau))))
    return dict(tau=tau, T=T, P=P, z=z - z0, rho=rho, g=g, Teff=Teff, tau_cut=tau_cut)


kern = json.load(open(sys.argv[1])); faces = kern['faces']
c = np.array([f['contribution'] for f in faces]); zf = np.array([f['depth_km'] for f in faces])*1e5; q = c.sum()
bg = bl.chem.prior.Background(); d = bg.d
r, rhob, mg = (np.asarray(d[k], float) for k in ['radius_cm', 'density_cgs', 'mass_geom_cm']); o = np.argsort(r); r, rhob, mg = r[o], rhob[o], mg[o]
R = float(bg.R); L = float(d['base_luminosity'])
Teff_d = (L/(4*np.pi*R**2*SIG))**0.25; logg_d = float(np.log10(C*C*mg[-1]/R**2))
dec = atmosphere(Teff_d, logg_d)
lnrho_pp = lambda A, z: np.interp(z, A['z'], np.log(A['rho']))
depth_bg = (R - r)[::-1]; lnrho_bg = np.log(rhob)[::-1]
zz = np.linspace(0, 624e5, 625); lay = (zz >= 240e5) & (zz <= 465e5)
dev = np.exp(lnrho_pp(dec, zz) - np.interp(zz, depth_bg, lnrho_bg)) - 1
res = dict(declared=dict(Teff=Teff_d, logg=logg_d), eos_calls=0,
           validation=dict(max_abs_rel_charge_layers=float(np.max(np.abs(dev[lay]))), max_abs_rel_0_624km=float(np.max(np.abs(dev))),
                           kernel_weighted=float(np.sum(np.abs(c/q)*np.abs(np.exp(lnrho_pp(dec, zf) - np.interp(zf, depth_bg, lnrho_bg)) - 1))),
                           rows=[dict(depth_km=float(z/1e5), rel=float(v)) for z, v in zip(zz[::48], dev[::48])]))
print('declared Teff %.2f log g %.5f ; validation charge layers max %.3e, kernel-weighted %.3e' % (Teff_d, logg_d, res['validation']['max_abs_rel_charge_layers'], res['validation']['kernel_weighted']), flush=True)
OBS = dict(logg=5.82, dlogg=0.05, Teff=15800., dTeff=100.)
cases = [('declared', Teff_d, logg_d), ('Kaplan central', OBS['Teff'], OBS['logg']),
         ('log g -1 sigma', OBS['Teff'], OBS['logg'] - OBS['dlogg']), ('log g +1 sigma', OBS['Teff'], OBS['logg'] + OBS['dlogg']),
         ('log g -2 sigma', OBS['Teff'], OBS['logg'] - 2*OBS['dlogg']), ('log g +2 sigma', OBS['Teff'], OBS['logg'] + 2*OBS['dlogg']),
         ('Teff -1 sigma', OBS['Teff'] - OBS['dTeff'], OBS['logg']), ('Teff +1 sigma', OBS['Teff'] + OBS['dTeff'], OBS['logg'])]
res['cases'] = {}; atm = {}
for name, te, lg in cases:
    A = dec if name == 'declared' else atmosphere(te, lg); atm[name] = A
    ratio = np.exp(lnrho_pp(A, zf) - lnrho_pp(dec, zf)); qc = float(np.sum(c*ratio))
    tau_at = lambda zk: float(np.exp(np.interp(zk, A['z'], np.log(A['tau']))))
    res['cases'][name] = dict(Teff=te, logg=lg, charge=qc, ratio=qc/q, sign_kept=bool(np.sign(qc) == np.sign(q)), tau_cut=A['tau_cut'],
                              tau_at_363km=tau_at(363e5), rho_at_363km=float(np.exp(lnrho_pp(A, 363e5))), depth_tau1_km=float(np.interp(0., np.log(A['tau']), A['z'])/1e5))
    print('%-16s Teff %.0f log g %.3f  charge %.4e ratio %.3f  tau(363km) %.3f  depth(tau=1) %.1f km' % (name, te, lg, qc, qc/q, res['cases'][name]['tau_at_363km'], res['cases'][name]['depth_tau1_km']), flush=True)
allr = [v['ratio'] for n_, v in res['cases'].items() if n_ != 'declared' and '2 sigma' not in n_]
res['one_sigma_ratio_range'] = [min(allr), max(allr)]; res['two_sigma_ratio_range'] = [res['cases']['log g -2 sigma']['ratio'], res['cases']['log g +2 sigma']['ratio']]
res['kernel_charge'] = float(q); res['eos_calls'] = calls[0]; res['seconds'] = time.monotonic() - t0
np.savez_compressed(sys.argv[2].replace('.json', '.npz'), **{f'{n_.replace(" ", "_")}__{k}': v for n_, A in atm.items() for k, v in A.items() if isinstance(v, np.ndarray)})
json.dump(res, open(sys.argv[2], 'w'), indent=1)
print('1-sigma ratio range %.3f .. %.3f ; 2-sigma %.3f .. %.3f ; %d EOS calls, %.1f s' % (*res['one_sigma_ratio_range'], *res['two_sigma_ratio_range'], calls[0], res['seconds']))
