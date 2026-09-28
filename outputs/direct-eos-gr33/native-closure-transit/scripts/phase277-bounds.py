"""Phase277: bounds for the self-consistent GR fixed point, ADM/infinity normalization, nonlinearity and phi_inf dependence.
Conjectural inputs: declared background sampler, the 4x readout grid, stored phase 252/260/270 measurements (quoted below).
Usage: python3 phase277-bounds.py <out json>
"""
import json, sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, G = 2.99792458e10, 6.67430e-8
d = dict(np.load('readout268-quad64-work/gr/field-source-64.npz'))
T = float(d['t'][-1]); D = float(d['drive_duration']); eta = float(d['drive_amplitude']); Rd = float(d['drive_radius'])
Mcm = float(d['M_cm']); Kcm = float(d['K_cm']); edges = np.asarray(d['edges'], float); radius = np.asarray(d['radius'], float)
bg = bl.chem.prior.Background(); b = bg.d
r, rho, P, G1, mg, phi = (np.asarray(b[k], float) for k in ['radius_cm', 'density_cgs', 'pressure_cgs', 'gamma1', 'mass_geom_cm', 'phi'])
o = np.argsort(r); r, rho, P, G1, mg, phi = r[o], rho[o], P[o], G1[o], mg[o], phi[o]
R = float(bg.R); phi_inf = 1e-3; alpha0 = -4*phi_inf; alphaA = Kcm/Mcm
inner = float(edges.min()); rho_max = float(np.interp(inner, r, rho))
L = 2*np.pi*G*rho_max*T**2  # radial Newtonian self-gravity feedback over the evolved interval
F1 = quad(lambda p: (4*p*(1 - p))**4, 0, 1)[0]
xi_max = 4*C*eta*phi.max()*D*F1*Rd/inner  # frozen free-fall displacement (x' ~ 1)
H = P/(rho*C*C*mg/np.maximum(r, 1.)**2); H_min = float(np.min(H[(r >= inner) & (r <= R)]))
cs = np.sqrt(G1*P/rho); cs_min = float(np.min(cs[(r >= inner) & (r <= R)]))
res = dict(T=T, D=D, eta=eta, R_d=Rd, readout_inner_edge=inner, depth_reached_km=(R - inner)/1e5,
    fixed_point=dict(rho_max=rho_max, L_bound=L, relative_residual_bound=L/(1 - L)),
    adm=dict(alpha0=alpha0, alphaA=alphaA, M_over_r=Mcm/Rd, exterior_cross_energy_relative=abs(alpha0*alphaA)*Mcm/Rd,
             mass_normalization_measured={'phase252_common_period': 3.30957e-5, 'phase260_full_period_1x': (2.483473532861 - 2.483417865013)/2.483417865013}),
    infinity=dict(schwarzschild_tail_bound=Mcm/Rd, phase252_exterior_scalar_relative=6.041864583481e-66/1.368965112963e-51),
    nonlinear=dict(dphi_over_phi0=eta*Rd/inner/phi.min(), displacement_over_H=xi_max/H_min, velocity_over_cs=xi_max*C/inner/cs_min, gr_second_order=eta, F1=F1, xi_max=xi_max),
    phi_inf=dict(state_share=0.9999999900051401, geometry_share=9.994859754558131e-09, crossover_if_geometry_phi0=phi_inf*np.sqrt(9.994859754558131e-09/0.9999999900051401)))
q_rich = -2.334965719648403e-51; mn = res['adm']['mass_normalization_measured']['phase260_full_period_1x']
res['continuum'] = dict(compact=q_rich, mass_normalized=q_rich*(1 + mn), per_eta_phi2=q_rich*(1 + mn)/(eta*phi_inf**2))
res['largest_correction'] = max(L, res['adm']['exterior_cross_energy_relative'], res['infinity']['schwarzschild_tail_bound'], res['nonlinear']['dphi_over_phi0'], mn)
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(json.dumps(res, indent=1))
