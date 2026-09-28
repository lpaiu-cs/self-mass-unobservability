"""Phase271: free-fall (pressureless) reproduction of the baryon redistribution and of the endpoint compact charge.

Conjectural model, derived from the declared theory only (not from the chain's matter equations):
  Einstein-frame static background ds^2 = -a^2 c^2 dt^2 + B^2 dr^2 + r^2 dOmega^2, coupling A(phi) = exp(-2 phi^2),
  alpha = dlnA/dphi = -4 phi, background scalar phi0(r). Dust follows geodesics of A^2 g; to first order in the
  incident field and in the velocity, the radial displacement obeys
      d2xi/dt2 = -c^2 (a^2/B^2) d/dr[alpha(phi0) dphi],   dphi = eta (R_d/r) f((t + x(r)/c)/D),
  f(p) = (4p(1-p))^4 on 0<p<1 (the chain's C3 pulse), x the characteristic coordinate of the drive (Born and response
  fields are 1e-10 and 1e-21 of the incident field and are neglected; pressure is neglected, P/(rho c^2) ~ 2e-9).
  With F1 = int_0^p f, F2 = int_0^p F1 the double time integral is exact:
      xi = 4 c^2 (a^2/B^2) eta R_d D^2 [ (phi0'/r - phi0/r^2) F2(p) + phi0 x' F1(p) / (r c D) ].
  Baryon mass through a face: Phi = 4 pi r^2 B rho0 xi; cell perturbation dM_k = Phi_k - Phi_{k+1}.
Background values (rho0, phi0, a, B) at the faces come from the declared background sampler; x at the faces is
interpolated from the drive's cell values. Only interior cells (boundary-layer geometry) are modelled; atmosphere
cells keep the chain's values.
Modes:
  compare <readout folder> <label> : dM(t, cell) vs the chain's baryon_g at the stored stage times.
  charge  <readout folder> <label> : endpoint charge from the free-fall dM through the unchanged setup/propagator
                                     (cubic Hermite state polynomials; every other source part masked), next to the
                                     chain's baryon-only charge.
"""
import json, shutil, sys, time
from pathlib import Path
from math import comb
import numpy as np
sys.path.insert(0, 'verification')
C = 2.99792458e10
mode, out, label = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
PHASE = sys.argv[4] if len(sys.argv) > 4 else 'x'  # pulse arrival: 'x' = drive_x/c, 'delay' = the source's delay array
z = dict(np.load(out/'gr/source-64.npz')); t = z['t']; edges = z['edges']; rc = z['radius']
eta, Rd, D = float(z['drive_amplitude']), float(z['drive_radius']), float(z['drive_duration'])
import def_native_boundary_layer as bl
ncell = len(np.load('outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz')['r'])  # interior (boundary-layer) cells
bg = bl.chem.prior.Background(); fe = edges[:ncell + 1].copy()
zf = bg.sample(fe); rho, phi, a, B = zf['rho'], zf['phi'], zf['a'], zf['B']
xc = z['drive_x'][:ncell]; xf = np.interp(fe, rc[:ncell], xc)  # linear in r; x is smooth (x' = B/a to 1e-5)
xf[0] = xc[0] - (rc[0] - fe[0])*(xc[1] - xc[0])/(rc[1] - rc[0]); xf[-1] = xc[-1] + (fe[-1] - rc[ncell - 1])*(xc[-1] - xc[-2])/(rc[ncell - 1] - rc[ncell - 2])
dx = np.gradient(xf, fe); dphi0 = np.gradient(phi, fe)
if PHASE == 'delay':  # the source's delay differs from drive_x/c by a constant; use it for the arrival time
    offset = z['delay'][:ncell] - xc/C; assert np.ptp(offset) < 1e-12, np.ptp(offset); xf = xf + C*float(np.mean(offset))
    print('arrival from the delay array: constant offset %.6e s' % float(np.mean(offset)))
# pulse integrals: f = sum_k w_k p^k (k=4..8); F1 = sum w_k p^(k+1)/(k+1); F2 = sum w_k p^(k+2)/((k+1)(k+2))
W = {4: 256., 5: -1024., 6: 1536., 7: -1024., 8: 256.}
f = lambda p: sum(w*p**k for k, w in W.items())
F1 = lambda p: sum(w*p**(k + 1)/(k + 1) for k, w in W.items())
F2 = lambda p: sum(w*p**(k + 2)/((k + 1)*(k + 2)) for k, w in W.items())
I0, J0 = F1(1.), F2(1.)
def pulse_integrals(p):
    p = np.asarray(p, float); inside = (p > 0) & (p < 1); after = p >= 1
    f0 = np.where(inside, f(np.clip(p, 0, 1)), 0.)
    g1 = np.where(inside, F1(np.clip(p, 0, 1)), np.where(after, I0, 0.))
    g2 = np.where(inside, F2(np.clip(p, 0, 1)), np.where(after, J0 + I0*(p - 1), 0.))
    return f0, g1, g2
pref = 4*C*C*(a*a/(B*B))*eta*Rd*D*D
flux_w = 4*np.pi*fe*fe*B*rho
def dM_at(tt):  # dM (g) and d(dM)/dt per interior cell at times tt
    p = (np.asarray(tt, float)[:, None] + xf[None]/C)/D
    f0, g1, g2 = pulse_integrals(p)
    xi = pref*((dphi0/fe - phi/fe**2)*g2 + phi*dx*g1/(fe*C*D))
    xidot = pref*((dphi0/fe - phi/fe**2)*g1 + phi*dx*f0/(fe*C*D))/D
    Phi, Phidot = flux_w*xi, flux_w*xidot
    return Phi[:, :-1] - Phi[:, 1:], Phidot[:, :-1] - Phidot[:, 1:]
report = out/f'freefall-{mode}-{label}.json'  # one record per mode
if mode == 'compare':
    m, _ = dM_at(t); chain = z['baryon_g'][:, :ncell]
    scale = float(np.sum(chain[-1]*m[-1])/np.sum(m[-1]*m[-1]))  # least-squares scale at T (1 if units agree)
    rows = []
    for c in range(ncell):
        if abs(chain[-1, c]) > 1e-3*np.max(abs(chain[-1])):
            rows.append(dict(cell=c, chain_T=float(chain[-1, c]), freefall_T=float(m[-1, c]), relative_T=float(m[-1, c]/chain[-1, c] - 1),
                             max_relative_over_time=float(np.max(abs(m[:, c] - chain[:, c]))/np.max(abs(chain[:, c])))))
    res = dict(label=label, cells=ncell, scale_fit=scale, rows=rows)
    report.write_text(json.dumps(res, indent=1) + '\n')
    print('least-squares scale (free-fall -> chain) at T: %.9f' % scale)
    for r in rows: print('cell %2d  chain %+.6e  free-fall %+.6e  rel(T) %+.3e  max rel over time %.3e' % (r['cell'], r['chain_T'], r['freefall_T'], r['relative_T'], r['max_relative_over_time']))
elif mode == 'charge':
    import read_full_captured_history as cap, extend_retarded_history as ext
    rch = cap.prior; KEYS = rch.KEYS
    y, yd = dM_at(t); h = np.diff(t)[:, None]
    y0, y1, d0, d1 = y[:-1], y[1:], yd[:-1]*h, yd[1:]*h  # normalized cubic Hermite per stage interval
    co = np.stack([y0, d0, 3*(y1 - y0) - 2*d0 - d1, 2*(y0 - y1) + d0 + d1])
    raw = dict(z); d = dict(np.load(out/'gr/field-source-64.npz')); M = float(d['M_cm']); T = float(d['t'][-1])
    def masked(model):
        rr = dict(raw)
        for key in KEYS:
            if key != 'baryon_g': rr['state_coeff_' + key] = np.zeros_like(raw['state_coeff_' + key])
            if key not in KEYS[-2:]: rr['geometry_coeff_' + key] = np.zeros_like(raw['geometry_coeff_' + key])
        if model:
            sc = raw['state_coeff_baryon_g'].copy(); sc[:, :, :ncell] = co; rr['state_coeff_baryon_g'] = sc
        return rr
    cap.bind(rch.endpoint.initialize, OUT=out)(); setup = ext.previous.prior.Response.setup; rows = {}
    for name, model in [('chain_baryon_only', False), ('freefall_interior', True)]:
        folder = out/f'freefall-in-{label}-{name}'; (folder/'gr').mkdir(parents=True)
        np.savez(folder/'gr/source-64.npz', **masked(model))
        s0 = time.monotonic(); mm = ext.Response(); ext.bind(setup, INPUT=folder)(mm, d, 8)
        free = ext.at(mm, mm.source, np.array([T])); rows[name] = -float(free[0][-1, -1])/M
        print(name, rows[name], 'seconds %.1f' % (time.monotonic() - s0), flush=True)
        del mm, free; shutil.rmtree(folder)
    rows['relative'] = rows['freefall_interior']/rows['chain_baryon_only'] - 1
    report.write_text(json.dumps(dict(label=label, cells=ncell, **rows), indent=1) + '\n'); print(json.dumps(rows))
else: raise ValueError(mode)
