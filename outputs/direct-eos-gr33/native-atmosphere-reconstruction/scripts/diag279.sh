cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
timeout 900 python3 - <<'PY'
import sys, numpy as np
sys.path.insert(0, 'verification')
import def_native_radiative_envelope as env, def_native_boundary_layer as bl
E = env.Envelope(); two = env.prior.two
print('prior.two module:', two.__name__)
bg = bl.chem.prior.Background(); d = bg.d
r, rho, T, P = (np.asarray(d[k], float) for k in ['radius_cm', 'density_cgs', 'temperature_K', 'pressure_cgs']); o = np.argsort(r); r, rho, T, P = r[o], rho[o], T[o], P[o]
R = float(bg.R); L = float(d['base_luminosity']); SIG = E.arad*2.99792458e10/4; Teff = (L/(4*np.pi*R**2*SIG))**0.25
env_ = (R - r) <= 700e5; rr, rh, TT, PP = r[env_][::-1], rho[env_][::-1], T[env_][::-1], P[env_][::-1]  # from the surface inward
kap = np.array([two.opacity_parts(E.opacity, (np.log(a), np.log(b), E.X))[0] for a, b in zip(rh, TT)])
dz = -np.diff(rr); tau_k = np.r_[0, np.cumsum(0.5*(kap[1:]*rh[1:] + kap[:-1]*rh[:-1])*dz)]
tau_T = np.clip(4/3*(TT/Teff)**4 - 2/3, 0, None)
rho_eos = np.array([E.eos(1, float(np.log(p)), float(np.log(t)), E.X)[0] for p, t in zip(PP, TT)])
for zk in [100, 200, 240, 300, 363, 400, 466, 515, 600]:
    i = int(np.argmin(abs((R - rr)/1e5 - zk)))
    print('depth %5.0f km  T %7.0f  rho %.3e  kappa %.3e  tau_kappa %.3e  tau_T %.3e  rho_EOS(P,T)/rho %.4f' % ((R - rr[i])/1e5, TT[i], rh[i], kap[i], tau_k[i], tau_T[i], rho_eos[i]/rh[i]))
PY
