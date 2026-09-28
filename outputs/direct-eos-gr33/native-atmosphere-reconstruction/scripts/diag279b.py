"""Phase279 diagnostic: declared plane-parallel gray atmosphere versus the declared background at matched depth and matched pressure."""
import sys, time, json
import numpy as np
sys.path.insert(0, 'verification')
exec(open('.phase279-gray.py').read().split("kern = json.load")[0])  # reuse EOS/opacity/atmosphere definitions only
bg = bl.chem.prior.Background(); d = bg.d
r, rhob, mg, Tb, Pb = (np.asarray(d[k], float) for k in ['radius_cm', 'density_cgs', 'mass_geom_cm', 'temperature_K', 'pressure_cgs']); o = np.argsort(r)
r, rhob, mg, Tb, Pb = r[o], rhob[o], mg[o], Tb[o], Pb[o]
R = float(bg.R); L = float(d['base_luminosity']); Teff = (L/(4*np.pi*R**2*SIG))**0.25; logg = float(np.log10(C*C*mg[-1]/R**2))
A = atmosphere(Teff, logg, n=400)
zb = (R - r)[::-1]; iz = lambda arr, z: np.interp(z, zb, arr[::-1])
print('tau_cut %.3e' % A['tau_cut'])
for zk in [0, 50, 100, 200, 240, 300, 363, 420, 466, 515, 600]:
    z = zk*1e5
    Tm, Pm, rm_ = (float(np.interp(z, A['z'], v)) for v in (A['T'], A['P'], A['rho']))
    print('z %4d km  T pp/bg %.5f  P pp/bg %.5f  rho pp/bg %.5f' % (zk, Tm/iz(Tb, z), Pm/iz(Pb, z), rm_/iz(rhob, z)))
for Pk in [1e2, 1e3, 1e4, 1e5, 5e5]:  # matched total pressure: depth and T
    zp = float(np.interp(np.log(Pk), np.log(A['P']), A['z'])); zbk = float(np.interp(np.log(Pk), np.log(Pb[::-1]), zb))
    print('P %.0e: depth pp %.2f km bg %.2f km ; T pp %.1f bg %.1f' % (Pk, zp/1e5, zbk/1e5, float(np.interp(np.log(Pk), np.log(A['P']), A['T'])), float(np.interp(np.log(Pk), np.log(Pb[::-1]), Tb[::-1]))))
