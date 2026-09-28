cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
timeout 300 python3 - <<'PY'
import sys, numpy as np
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C = 2.99792458e10
bg = bl.chem.prior.Background(); d = bg.d
r, rho, P, mg, T, lapse, e = (np.asarray(d[k], float) for k in ['radius_cm', 'density_cgs', 'pressure_cgs', 'mass_geom_cm', 'temperature_K', 'lapse', 'energy_cgs'])
o = np.argsort(r); r, rho, P, mg, T, lapse, e = (v[o] for v in (r, rho, P, mg, T, lapse, e)); R = float(bg.R)
g = C*C*mg/r**2; geff = -np.gradient(P, r)/rho; arad = 7.565733250280002e-15
Prad = arad*T**4/3; Pg = P - Prad; geff_gas = -np.gradient(Pg, r)/rho
for zk in [150, 250, 300, 363, 420, 466, 515, 600]:
    i = int(np.argmin(abs((R - r)/1e5 - zk)))
    print('depth %4.0f km  -dP/dr/rho / g = %.5f   -dPgas/dr/rho / g = %.5f   (e/(rho c^2)) = %.6f  lapse %.8f' % ((R - r[i])/1e5, geff[i]/g[i], geff_gas[i]/g[i], e[i]/(rho[i]*C*C), lapse[i]))
PY
