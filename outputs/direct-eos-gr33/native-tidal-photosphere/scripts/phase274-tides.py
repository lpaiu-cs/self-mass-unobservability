"""Phase274: tidal (l>=2) channel and dissipation of the free-fall charge at orbital timescales (declared background only).

Conjectural (Newtonian weak field, compactness 4e-6, linear scalar coupling A = exp(-2 phi^2), alpha = -4 phi, beta = -4):
 1. apsidal constant k2 from Clairaut-Radau  r eta' = -eta^2 + eta + 6 - 6 (rho/rhobar)(eta + 1), eta(0) = 0,
    k2_aps = (3 - eta_R)/(2 (2 + eta_R)); the fluid Love number is 2 k2_aps.
 2. N^2 = g (dlnP/dr / Gamma1 - dlnrho/dr); asymptotic g-mode spacing DeltaPi_l = 2 pi^2/(sqrt(l(l+1)) int N/r dr).
 3. optically thick envelope (gray tau >= 1, no energy sources, L_r = L): K = L/(4 pi r^2 |dT/dr| rho c_p), c_p = 2.5 P_gas/(rho T),
    P_gas = P - a T^4/3; one-way damping depth of an l=2 gravity wave tau_w = int K k_r^3/omega dr, k_r = sqrt(6) N/(omega r)
    (lower bound: the optically thin atmosphere, where diffusion fails, and the core are left out).
 4. thermal time of the overlying envelope tau_th(r) = int_r^R 2.5 P_gas 4 pi r^2 dr / L (ionization energy left out, <~ x2.5).
 5. static structural monopole response: adiabatic radial operator of phase 273 at omega = 0 with the extra force -eps g per unit
    mass, dq/eps = (1/M) int xi w' dm, w = alpha(phi0) * lapse. Self-check: uniform Gamma1 = 5/3 gives zeta = -eps exactly.
 6. J0337 magnitudes (masses, semimajor axes, periods: repository parameter_set=6; eccentricities: Ransom et al. 2014).
Usage: python3 phase274-tides.py <out json>
"""
import json, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.integrate import solve_ivp, cumulative_trapezoid, trapezoid
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, G, MSUN_CM, DAY, SIG, ARAD = 2.99792458e10, 6.67430e-8, 1.476625e5, 86400., 5.670374419e-5, 7.565723e-15
bg = bl.chem.prior.Background(); d = bg.d
keys = ['radius_cm', 'mass_geom_cm', 'pressure_cgs', 'density_cgs', 'gamma1', 'temperature_K', 'phi', 'phi_prime_cm', 'lapse']
r, mg, P, rho, G1, T, phi, dphi, lapse = (np.asarray(d[k], float) for k in keys)
order = np.argsort(r); r, mg, P, rho, G1, T, phi, dphi, lapse = (v[order] for v in (r, mg, P, rho, G1, T, phi, dphi, lapse))
keep = np.r_[True, np.diff(r) > 0]; r, mg, P, rho, G1, T, phi, dphi, lapse = (v[keep] for v in (r, mg, P, rho, G1, T, phi, dphi, lapse))
R = float(bg.R); M = C*C*mg[-1]/G; g = C*C*mg/np.maximum(r, 1.)**2; L = float(d['base_luminosity'])
rj = float(np.asarray(d['radius_cm'], float)[int(d['junction_index']) - 1]); env = r >= rj
res = dict(points=int(len(r)), R=R, M_g=M, M_sun=M/1.98841e33, L=L, junction_radius=rj, junction_depth_km=(R - rj)/1e5)

# 1. apsidal constant
rhobar = 3*mg*C*C/(G*4*np.pi*np.maximum(r, 1.)**3); ratio = np.where(r > 0, rho/rhobar, 1.)
f = lambda lr, e: [-e[0]**2 + e[0] + 6 - 6*np.interp(np.exp(lr), r, ratio)*(e[0] + 1)]
sol = solve_ivp(f, [np.log(r[1]), np.log(R)], [0.], rtol=1e-10, atol=1e-12, max_step=0.01)
eta = float(sol.y[0, -1]); k2 = (3 - eta)/(2*(2 + eta))
poly1 = lambda x: x*x*np.sin(x)/(3*(np.sin(x) - x*np.cos(x)))  # n=1 polytrope rho/rhobar at xi
e1 = solve_ivp(lambda lx, e: [-e[0]**2 + e[0] + 6 - 6*poly1(np.exp(lx))*(e[0] + 1)], [np.log(1e-4), np.log(np.pi)], [0.], rtol=1e-11, atol=1e-13).y[0, -1]
k2_poly1 = (3 - e1)/(2*(2 + e1)); assert abs(k2_poly1 - (15/np.pi**2 - 1)/2) < 1e-6, k2_poly1  # known n=1 apsidal constant
res['love'] = dict(n1_polytrope_check=float(k2_poly1), eta_R=eta, k2_apsidal=k2, love_number=2*k2, central_to_mean_density=float(rho[0]/(3*M/(4*np.pi*R**3))))

# 2. buoyancy and g-mode spacing
N2 = g*(np.gradient(np.log(P), r)/G1 - np.gradient(np.log(rho), r)); N = np.sqrt(np.clip(N2, 0, None))
In = float(trapezoid(np.where(r > 0, N/np.maximum(r, 1.), 0.), r))
dPi = {l: 2*np.pi**2/(np.sqrt(l*(l + 1))*In) for l in (1, 2)}
neg = N2 < 0
res['buoyancy'] = dict(int_N_over_r=In, DeltaPi_l1_s=dPi[1], DeltaPi_l2_s=dPi[2], N_max=float(N.max()), r_N_max_depth_km=float((R - r[np.argmax(N)])/1e5),
                       convective_fraction_of_radius=float(trapezoid(neg.astype(float), r)/R),
                       convective_depth_ranges_km=[[float((R - r[i])/1e5)] for i in np.flatnonzero(np.diff(neg.astype(int)) != 0)][:12])

# 3. traveling-wave test and 4. thermal time (envelope only)
Teff = (L/(4*np.pi*R**2*SIG))**0.25; tau_g = np.clip(4/3*(T/Teff)**4 - 2/3, 0, None); Pg = P - ARAD*T**4/3; thick = env & (tau_g >= 1)
dT = np.abs(np.gradient(T, r)); cp = 2.5*Pg/(rho*T); K = L/(4*np.pi*r**2*np.maximum(dT, 1e-300)*rho*cp)
th = np.zeros_like(r); th[env] = cumulative_trapezoid((2.5*Pg*4*np.pi*r**2)[env][::-1], r[env][::-1], initial=0)[::-1]*-1/L
depth = (R - r)/1e5
n_in, n_out = 2*np.pi/(1.62939901*DAY), 2*np.pi/(327.25512703*DAY)
drives = {'2 n_in (non-rotating tide)': 2*n_in, 'n_in (eccentric)': n_in, 'n_out': n_out, '1 h': 2*np.pi/3600., '10 min': 2*np.pi/600.}
res['drives'] = {}
for name, w in drives.items():
    kr = np.sqrt(6)*N/(w*np.maximum(r, 1.)); m = thick & (N > w)
    tau = float(trapezoid(np.where(m, K*kr**3/w, 0.), r)); Pw = 2*np.pi/w
    res['drives'][name] = dict(omega=w, period_s=Pw, gmode_order_l2=Pw/dPi[2], relative_spacing=dPi[2]/Pw, tau_wave_thick_envelope=tau,
                               depth_tau_th_equals_1_over_omega_km=float(np.interp(1/w, th[env][::-1], depth[env][::-1])) if th[env].max() > 1/w else None)
lay = env & (depth >= 240) & (depth <= 465)
res['thermal'] = dict(tau_th_charge_layers_s=[float(th[lay].min()), float(th[lay].max())], tau_th_base_s=float(th[env].max()),
                      Teff=float(Teff), thick_depth_km=[float(depth[thick].min()), float(depth[thick].max())], K_thick=[float(K[thick].min()), float(K[thick].max())],
                      N_thick=[float(N[thick].min()), float(N[thick].max())], Prad_over_P_charge_layers=[float((1 - Pg/P)[lay].min()), float((1 - Pg/P)[lay].max())])

# 5. static structural monopole response
def static(G1v):
    h = np.diff(r); rm = (r[:-1] + r[1:])/2
    A = np.interp(rm, r, G1v*P*r**4); Q = np.interp(rm, r, r**3*np.gradient((3*G1v - 4)*P, r)); Fm = np.interp(rm, r, rho*r**3*g)
    n = len(r); main = np.zeros(n); off = np.zeros(n - 1); F = np.zeros(n)
    main[:-1] += A/h - Q*h/3; main[1:] += A/h - Q*h/3; off += -A/h - Q*h/6
    F[:-1] += -Fm*h/2; F[1:] += -Fm*h/2  # force -eps g per unit mass, eps = 1
    return spsolve(sp.diags([off, main, off], [-1, 0, 1], format='csc'), F)
z53 = static(np.full_like(G1, 5/3)); dmw = 4*np.pi*r**2*rho
check = float(trapezoid(dmw*np.abs(z53 + 1), r)/trapezoid(dmw, r)); check_max = float(np.max(np.abs(z53[(R - r) > 1e5] + 1)))
assert check < 1e-3, check  # uniform Gamma1 = 5/3: zeta = -eps/(3 Gamma1 - 4) = -1 (up to the discrete hydrostatic balance)
zeta = static(G1); w = -4*phi*lapse; wp = np.gradient(w, r); dm = 4*np.pi*r**2*rho
dq_eps = float(trapezoid(dm*r*zeta*wp, r)/M)
phi_inf, alpha0, beta = 1e-3, -4e-3, -4.
eps_per_dphi = 4*abs(alpha0) + 2*abs(alpha0*beta)
kappa = abs(dq_eps)*eps_per_dphi
res['structural'] = dict(uniform_gamma_check_mass_weighted=check, uniform_gamma_check_max_below_1km=check_max, zeta_center=float(zeta[1]), zeta_surface=float(zeta[-1]), zeta_charge_layers=[float(zeta[lay].min()), float(zeta[lay].max())],
                         dq_per_eps=dq_eps, eps_per_dphi_max=eps_per_dphi, kappa_struct_max=kappa, over_direct_beta=kappa/abs(beta),
                         scaling='dq/eps ~ phi_inf, eps/dphi ~ phi_inf: kappa ~ phi_inf^2 (declared phi_inf = 1e-3)')

# 6. J0337
mp, mc, mo = 1.4378144085, 0.1975363853, 0.4101027069; a_in, a_out = 4.7761915344e11, 1.7648750789e13; e_in, e_out = 6.92e-4, 0.0354
q = mp/mc; x = R/a_in
phi_p = mp*MSUN_CM/a_in; phi_o = mo*MSUN_CM/a_out
mod = {'pulsar, inner eccentricity': phi_p*e_in, 'outer, outer eccentricity': phi_o*e_out, 'outer, inner-orbit motion': phi_o*a_in*mp/(mp + mc)/a_out}
tide = dict(tidal_parameter_pulsar=q*x**3, tidal_parameter_outer=(mo/mc)*(R/a_out)**3, fractional_tidal_force=6*k2*q*x**5,
            apsidal_tide_over_n=15*k2*q*x**5, apsidal_gr_over_n=3*(mp + mc)*MSUN_CM/(a_in*(1 - e_in**2)))
tide['apsidal_tide_over_gr'] = tide['apsidal_tide_over_n']/tide['apsidal_gr_over_n']
wmax = max(mod.values())
res['j0337'] = dict(dphi_static_per_unit_charge=dict(pulsar=phi_p, outer=phi_o), dphi_modulation_per_unit_charge=mod, tides=tide,
                    lagged_monopole_delta_bound=kappa*wmax, paper_b_lag_limit=1.7e-9, margin=1.7e-9/(kappa*wmax),
                    second_order_monopole_periodic=abs(dq_eps)*(q*x**3)**2*6*e_in,
                    nonstatic_inphase_tide=(2*n_in*float(np.sqrt(R**3/(C*C*mg[-1]))))**2)
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(json.dumps(res, indent=1))
