"""Phase287: thermal-relaxation strength of the white dwarf's structural monopole response at orbital frequencies
(declared background only; reuses the loader and the static solver of phase274-tides.py).

Conjectural (per-layer relaxation picture):
 - a layer at radius r relaxes thermally on the overlying thermal time tau_th(r) = int_r^R c_p T dm / L with c_p T = 2.5 P_gas/rho
   (ideal monatomic gas; ionization energy left out, so tau_th is underestimated by up to ~2.5 and more layers count as fast);
 - its fully relaxed limit is bracketed by the isothermal exponent Gamma_T = P_gas/P (ideal gas plus radiation) instead of Gamma1.
For each cut tau_cut, the layers with tau_th < tau_cut take Gamma_T and the static monopole response of phase 274 is solved again:
dS(tau_cut) = dq/eps(relaxed above the cut) - dq/eps(adiabatic). The increment of dS between neighbouring cuts is the relaxation
strength of the shell between them. At drive frequency omega a Debye shell of strength s and time tau adds the quadrature
|s| omega tau/(1 + omega^2 tau^2) <= |s|/2 whatever its sign, so Q(omega) = sum |ds| f(omega tau) estimates the lag and
TV/2 = sum |ds|/2 bounds it within the scanned shells. The scan stops before the modified operator loses positive definiteness.
Results are relative to the adiabatic dq/eps (the same ratio holds for S_struct).
Usage: python3 phase287-relax.py <out json>
"""
import json, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.integrate import cumulative_trapezoid, trapezoid
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, G, DAY, ARAD = 2.99792458e10, 6.67430e-8, 86400., 7.565723e-15
bg = bl.chem.prior.Background(); d = bg.d
keys = ['radius_cm', 'mass_geom_cm', 'pressure_cgs', 'density_cgs', 'gamma1', 'temperature_K', 'phi', 'phi_prime_cm', 'lapse']
r, mg, P, rho, G1, T, phi, dphi, lapse = (np.asarray(d[k], float) for k in keys)
order = np.argsort(r); r, mg, P, rho, G1, T, phi, dphi, lapse = (v[order] for v in (r, mg, P, rho, G1, T, phi, dphi, lapse))
keep = np.r_[True, np.diff(r) > 0]; r, mg, P, rho, G1, T, phi, dphi, lapse = (v[keep] for v in (r, mg, P, rho, G1, T, phi, dphi, lapse))
R = float(bg.R); M = C*C*mg[-1]/G; g = C*C*mg/np.maximum(r, 1.)**2; L = float(d['base_luminosity'])


def matrix(G1v):
    h = np.diff(r); rm = (r[:-1] + r[1:])/2
    A = np.interp(rm, r, G1v*P*r**4); Q = np.interp(rm, r, r**3*np.gradient((3*G1v - 4)*P, r)); Fm = np.interp(rm, r, rho*r**3*g)
    n = len(r); main = np.zeros(n); off = np.zeros(n - 1); F = np.zeros(n)
    main[:-1] += A/h - Q*h/3; main[1:] += A/h - Q*h/3; off += -A/h - Q*h/6
    F[:-1] += -Fm*h/2; F[1:] += -Fm*h/2  # force -eps g per unit mass, eps = 1 (phase 274)
    return main, off, F


def positive_definite(main, off):  # LDL^T pivots of the symmetric tridiagonal matrix
    p = main[0]
    if p <= 0: return False
    for i in range(1, len(main)):
        p = main[i] - off[i - 1]**2/p
        if p <= 0: return False
    return True


def solve(G1v):
    main, off, F = matrix(G1v)
    return spsolve(sp.diags([off, main, off], [-1, 0, 1], format='csc'), F), positive_definite(main, off)


w = -4*phi*lapse; wp = np.gradient(w, r); dm = 4*np.pi*r**2*rho
dq = lambda zeta: float(trapezoid(dm*r*zeta*wp, r)/M)
z0, pd0 = solve(G1); dq0 = dq(z0); assert pd0
assert abs(dq0/-1.840953484516629e-07 - 1) < 1e-9, dq0  # phase 274 stored value

Pg = P - ARAD*T**4/3; GT = Pg/P
th = cumulative_trapezoid((2.5*Pg*4*np.pi*r**2)[::-1], r[::-1], initial=0)[::-1]*-1/L  # overlying thermal time, 0 at the surface
mfrac_above = (M - C*C*mg/G)/M
n_in, n_out = 2*np.pi/(1.62939901*DAY), 2*np.pi/(327.25512703*DAY)
depth_km = lambda tau: float(np.interp(tau, th[::-1], ((R - r)/1e5)[::-1]))
mass_above = lambda tau: float(np.interp(tau, th[::-1], mfrac_above[::-1]))

cuts = np.logspace(-3, 13, 161); rows = []
for tc in cuts:
    m = th < tc
    if not m.any(): continue
    z, pd = solve(np.where(m, GT, G1))
    if not pd: break
    rows.append(dict(tau_cut=float(tc), points=int(m.sum()), depth_km=depth_km(tc), mass_fraction=mass_above(tc), dS_rel=(dq(z) - dq0)/abs(dq0)))
tau = np.array([x['tau_cut'] for x in rows]); dS = np.array([x['dS_rel'] for x in rows])
inc = np.diff(np.r_[0., dS]); tmid = np.sqrt(np.r_[tau[0]/10, tau[:-1]]*tau)  # shell between neighbouring cuts
f = lambda wt: wt/(1 + wt*wt)
out = dict(R=R, M_g=M, L=L, dq_eps_adiabatic=dq0, gamma_T_range=[float(GT.min()), float(GT.max())],
           tau_th_center_s=float(th[1]), cuts_scanned=len(rows), last_cut_s=float(tau[-1]), stopped_by_positivity=len(rows) < len(cuts),
           last_cut_depth_km=rows[-1]['depth_km'], last_cut_mass_fraction=rows[-1]['mass_fraction'], dS_rel_at_last_cut=float(dS[-1]),
           total_variation_rel=float(np.abs(inc).sum()))
for name, om in (('n_in', n_in), ('n_out', n_out), ('2 n_in', 2*n_in)):
    q = np.abs(inc)*f(om*tmid)
    out[name] = dict(omega=om, one_over_omega_s=1/om, depth_km_at_tau_1_over_omega=depth_km(1/om), mass_fraction_above=mass_above(1/om),
                     debye_quadrature_rel=float(q.sum()), shells_within_1e4_of_1_over_omega_TV_half=float(np.abs(inc)[(om*tmid > 1e-4) & (om*tmid < 1e4)].sum()/2),
                     deeper_than_scan_suppression=float(1/(om*tau[-1])))
out['profile'] = [dict(tau_cut=x['tau_cut'], depth_km=x['depth_km'], mass_fraction=x['mass_fraction'], dS_rel=x['dS_rel']) for x in rows[::10]] + [rows[-1]]
json.dump(out, open(sys.argv[1], 'w'), indent=1); print(json.dumps({k: v for k, v in out.items() if k != 'profile'}, indent=1))
