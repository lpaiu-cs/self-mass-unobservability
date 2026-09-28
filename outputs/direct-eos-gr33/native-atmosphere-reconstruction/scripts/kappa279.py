"""Phase279 check: Rosseland mean of the simple H+He continuum opacity versus the table opacity along the declared gray atmosphere."""
import json, sys
import numpy as np
NG = {}; exec(open('phase279-nongray.py', encoding='utf-8').read().split('# self-checks')[0], NG)
A = np.load('phase279-gray-sph.npz'); z, tau, rho, T = (A[f'declared__{k}'] for k in ['z', 'tau', 'rho', 'T'])
kap_tab = np.gradient(tau, z)/rho  # table Rosseland opacity recovered from dtau/dz = kappa rho
ab, sc = NG['opacity'](rho, T); kR = NG['rosseland'](ab, sc, T)
rows = []
for zk in [100, 200, 240, 300, 363, 420, 466, 515, 600]:
    i = int(np.argmin(abs(z/1e5 - zk))); rows.append(dict(depth_km=float(z[i]/1e5), T=float(T[i]), rho=float(rho[i]), kappa_table=float(kap_tab[i]), kappa_simple_R=float(kR[i]), ratio=float(kR[i]/kap_tab[i])))
    print('depth %5.0f km  T %7.0f  rho %.3e  kappa table %.3e  simple Rosseland %.3e  ratio %.3f' % tuple(rows[-1][k] for k in ['depth_km', 'T', 'rho', 'kappa_table', 'kappa_simple_R', 'ratio']))
lay = (z >= 240e5) & (z <= 465e5); r = kR[lay]/kap_tab[lay]
json.dump(dict(rows=rows, charge_layer_ratio_range=[float(r.min()), float(r.max())]), open('phase279-kappa.json', 'w'), indent=1)
print('charge layers ratio range %.3f .. %.3f' % (r.min(), r.max()))
