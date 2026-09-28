"""Phase273 parts 2-4: restoring timescales of the declared background and the long-wavelength force ratio.

Conjectural (Newtonian linear adiabatic radial oscillations; compactness 4e-6, GR corrections neglected).
 - omega_0: lowest radial adiabatic eigenfrequency of the whole background star from
       d/dr[Gamma1 P r^4 zeta'] + r^3 d/dr[(3 Gamma1 - 4) P] zeta + omega^2 rho r^4 zeta = 0,  zeta = delta r / r,
   as a linear finite-element generalized eigenproblem K z = omega^2 M z with natural boundary conditions.
 - omega_ac = c_s/(2H) in the charge layers (depth 138-482 km), c_s^2 = Gamma1 P/rho, H = P/(rho g), g = c^2 m_geom/r^2.
 - Long-wavelength monopole force: for a drive uniform over the star the scalar force per unit mass is
   c^2 (a^2/B^2) |d alpha(phi0)/dr| dphi = 4 c^2 |phi0'| dphi; the pulse force is 4 c^2 phi0 |d dphi/dr|, with
   |d dphi/dr| <= dphi_max max|f'|/(c D max f) for the chain's pulse. Their ratio at equal amplitude is reported.
 - Orbital drives: (omega/omega_0)^2 and (omega/omega_ac)^2 for several periods.
Usage: python3 phase273-modes.py <out json>
"""
import json, sys
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C = 2.99792458e10
bg = bl.chem.prior.Background(); d = bg.d
r = np.asarray(d['radius_cm'], float); mg = np.asarray(d['mass_geom_cm'], float); P = np.asarray(d['pressure_cgs'], float)
rho = np.asarray(d['density_cgs'], float); G1 = np.asarray(d['gamma1'], float); phi = np.asarray(d['phi'], float); dphi = np.asarray(d['phi_prime_cm'], float)
order = np.argsort(r); r, mg, P, rho, G1, phi, dphi = (v[order] for v in (r, mg, P, rho, G1, phi, dphi))
keep = np.r_[True, np.diff(r) > 0]; r, mg, P, rho, G1, phi, dphi = (v[keep] for v in (r, mg, P, rho, G1, phi, dphi))
R = float(bg.R); g = C*C*mg/np.maximum(r, 1.)**2
res = dict(points=int(len(r)), r_min=float(r[0]), r_max=float(r[-1]), R=R, M_geom_cm=float(mg[-1]), central_density=float(rho[0]),
           gamma1_range=[float(G1.min()), float(G1.max())])
# charge layers
depth = (R - r)/1e5; layer = (depth >= 138) & (depth <= 482)
cs = np.sqrt(G1*P/rho); H = P/(rho*g); wac = cs/(2*H)
res['charge_layers'] = dict(c_s_cm_s=[float(cs[layer].min()), float(cs[layer].max())], H_km=[float(H[layer].min()/1e5), float(H[layer].max()/1e5)],
                            omega_ac=[float(wac[layer].min()), float(wac[layer].max())], period_ac_s=[float(2*np.pi/wac[layer].max()), float(2*np.pi/wac[layer].min())],
                            sound_crossing_17km_s=[float(1.7e6/cs[layer].max()), float(1.7e6/cs[layer].min())])
# long-wavelength force ratio at equal amplitude
p = np.linspace(0, 1, 200001); f = (4*p*(1 - p))**4; fp = np.gradient(f, p); D = 0.0017172155589643512
grad_scale = np.max(np.abs(fp))/(C*D*np.max(f))  # max |d dphi/dr| / dphi_max for the chain's pulse (1/cm)
ratio = np.abs(dphi[layer])/(np.abs(phi[layer])*grad_scale)
res['long_wavelength_force_ratio'] = dict(phi0=[float(phi[layer].min()), float(phi[layer].max())], phi0_prime_per_cm=[float(np.abs(dphi[layer]).min()), float(np.abs(dphi[layer]).max())],
                                         pulse_gradient_scale_per_cm=float(grad_scale), ratio=[float(ratio.min()), float(ratio.max())])
# radial adiabatic eigenproblem (linear FEM on the background grid)
h = np.diff(r); rm = (r[:-1] + r[1:])/2
A = np.interp(rm, r, G1*P*r**4)                       # stiffness coefficient
Q = np.interp(rm, r, r**3*np.gradient((3*G1 - 4)*P, r))  # potential term coefficient
Mc = np.interp(rm, r, rho*r**4)                       # mass coefficient
n = len(r); K = np.zeros((n, n)); Mm = np.zeros((n, n))
for e in range(n - 1):
    ke = A[e]/h[e]*np.array([[1, -1], [-1, 1]]); qe = -Q[e]*h[e]/6*np.array([[2, 1], [1, 2]]); me = Mc[e]*h[e]/6*np.array([[2, 1], [1, 2]])
    K[e:e+2, e:e+2] += ke + qe; Mm[e:e+2, e:e+2] += me
w2 = eigh(K, Mm, subset_by_index=[0, 3], eigvals_only=True)
res['radial_modes'] = dict(omega_squared=[float(v) for v in w2], omega=[float(np.sqrt(v)) if v > 0 else None for v in w2],
                           period_s=[float(2*np.pi/np.sqrt(v)) if v > 0 else None for v in w2], dynamical_time_s=float(np.sqrt(R**3/(C*C*mg[-1]))))
w0 = float(np.sqrt(w2[0])) if w2[0] > 0 else None
orbits = {'10 min': 600., '1 h': 3600., '1 d': 86400., '1.629 d (J0337 inner)': 1.629*86400.}
res['orbital_drives'] = {k: dict(omega=2*np.pi/Pd, over_omega0_squared=(2*np.pi/Pd/w0)**2 if w0 else None,
                                 over_omega_ac_squared=(2*np.pi/Pd/float(wac[layer].min()))**2) for k, Pd in orbits.items()}
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(json.dumps(res, indent=1))
