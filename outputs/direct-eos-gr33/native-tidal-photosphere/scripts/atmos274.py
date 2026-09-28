"""Phase274/275 helper: where the charge layers sit in the declared atmosphere (column mass, gray optical depth), N^2 sign in the envelope."""
import json, sys
import numpy as np
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, G, SIG = 2.99792458e10, 6.67430e-8, 5.670374419e-5
bg = bl.chem.prior.Background(); d = bg.d
r, mg, P, rho, G1, T = (np.asarray(d[k], float) for k in ['radius_cm', 'mass_geom_cm', 'pressure_cgs', 'density_cgs', 'gamma1', 'temperature_K'])
o = np.argsort(r); r, mg, P, rho, G1, T = (v[o] for v in (r, mg, P, rho, G1, T)); k = np.r_[True, np.diff(r) > 0]; r, mg, P, rho, G1, T = (v[k] for v in (r, mg, P, rho, G1, T))
R = float(bg.R); L = float(d['base_luminosity']); g = C*C*mg/np.maximum(r, 1.)**2; Teff = (L/(4*np.pi*R**2*SIG))**0.25
depth = (R - r)/1e5; tau = np.clip(4/3*(T/Teff)**4 - 2/3, 0, None); col = P/g
N2 = g*(np.gradient(np.log(P), r)/G1 - np.gradient(np.log(rho), r)); env = depth <= 624.5
out = dict(Teff=Teff, logg_surface=float(np.log10(g[-1])), R_sun=R/6.957e10, T_surface=float(T[-1]), T_over_Teff_surface=float(T[-1]/Teff),
           envelope_N2_negative_depths_km=[float(x) for x in depth[env & (N2 < 0)]][:20], rows=[])
for z in [0, 50, 100, 150, 200, 240, 300, 363, 378, 420, 465, 515, 560, 624]:
    i = int(np.argmin(np.abs(depth - z)))
    out['rows'].append(dict(depth_km=float(depth[i]), rho=float(rho[i]), P=float(P[i]), T=float(T[i]), column_g_cm2=float(col[i]), gray_tau=float(tau[i]),
                            H_km=float(P[i]/(rho[i]*g[i])/1e5), Gamma1=float(G1[i])))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print('Teff %.1f logg %.4f R %.5f Rsun  T_s/Teff %.4f' % (Teff, out['logg_surface'], out['R_sun'], out['T_over_Teff_surface']))
print('envelope N2<0 depths (first 20):', out['envelope_N2_negative_depths_km'])
for row in out['rows']: print('  depth %7.1f km  rho %.3e  P %.3e  T %8.0f  m %.3e g/cm2  tau %.3e  H %.1f km  G1 %.3f' % tuple(row[k] for k in ['depth_km', 'rho', 'P', 'T', 'column_g_cm2', 'gray_tau', 'H_km', 'Gamma1']))
