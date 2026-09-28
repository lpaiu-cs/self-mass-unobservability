"""Phase272: exact face-density kernel of the free-fall endpoint charge (4x grid) and background-structure families.

Conjectural (model of phase 271). The free-fall displacement is independent of rho0, so the endpoint charge is linear
in the face densities: q = sum_f c_f rho0_f (face f adds +Phi_f to cell f and -Phi_f to cell f-1). Each face's
contribution c_f rho0_f is read through the unchanged setup/propagator with only that face's flux in the baryon state
(cubic Hermite), every other source part masked and atmosphere cells zero. The full interior free-fall charge is read
the same way and must equal the sum of the face contributions.
Modes:
  faces   <readout> <label> <first face> <last face>   contributions of a face range (inclusive) -> kernel-<label>-<a>-<b>.json
  full    <readout> <label>                             interior free-fall charge (all faces) -> kernel-<label>-full.json
  analyze <readout> <label>                             kernel, sign margin, smooth structure families -> kernel-<label>.json
"""
import glob, json, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
C = 2.99792458e10
mode, out, label = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
z = dict(np.load(out/'gr/source-64.npz')); t = z['t']; edges = z['edges']; rc = z['radius']
eta, Rd, D = float(z['drive_amplitude']), float(z['drive_radius']), float(z['drive_duration'])
ncell = len(np.load('outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz')['r'])
fe = edges[:ncell + 1].copy()
if mode != 'analyze':
    import def_native_boundary_layer as bl
    bg = bl.chem.prior.Background(); zf = bg.sample(fe); rho, phi, a, B = zf['rho'], zf['phi'], zf['a'], zf['B']
    xc = z['drive_x'][:ncell]; xf = np.interp(fe, rc[:ncell], xc)
    xf[0] = xc[0] - (rc[0] - fe[0])*(xc[1] - xc[0])/(rc[1] - rc[0]); xf[-1] = xc[-1] + (fe[-1] - rc[ncell - 1])*(xc[-1] - xc[-2])/(rc[ncell - 1] - rc[ncell - 2])
    dx = np.gradient(xf, fe); dphi0 = np.gradient(phi, fe)
    W = {4: 256., 5: -1024., 6: 1536., 7: -1024., 8: 256.}
    f = lambda p: sum(w*p**k for k, w in W.items())
    F1 = lambda p: sum(w*p**(k + 1)/(k + 1) for k, w in W.items())
    F2 = lambda p: sum(w*p**(k + 2)/((k + 1)*(k + 2)) for k, w in W.items())
    I0, J0 = F1(1.), F2(1.)
    p = (t[:, None] + xf[None]/C)/D; inside = (p > 0) & (p < 1); after = p >= 1; pc = np.clip(p, 0, 1)
    f0 = np.where(inside, f(pc), 0.); g1 = np.where(inside, F1(pc), np.where(after, I0, 0.)); g2 = np.where(inside, F2(pc), np.where(after, J0 + I0*(p - 1), 0.))
    pref = 4*C*C*(a*a/(B*B))*eta*Rd*D*D
    xi = pref*((dphi0/fe - phi/fe**2)*g2 + phi*dx*g1/(fe*C*D)); xidot = pref*((dphi0/fe - phi/fe**2)*g1 + phi*dx*f0/(fe*C*D))/D
    w = 4*np.pi*fe*fe*B*rho; Phi, Phidot = w*xi, w*xidot
    import read_full_captured_history as cap, extend_retarded_history as ext
    rch = cap.prior; KEYS = rch.KEYS; raw = dict(z); d = dict(np.load(out/'gr/field-source-64.npz')); M = float(d['M_cm']); T = float(d['t'][-1])
    cap.bind(rch.endpoint.initialize, OUT=out)(); setup = ext.previous.prior.Response.setup
    h = np.diff(t)[:, None]
    def charge(y, yd, tag):
        y0, y1, d0, d1 = y[:-1], y[1:], yd[:-1]*h, yd[1:]*h
        co = np.stack([y0, d0, 3*(y1 - y0) - 2*d0 - d1, 2*(y0 - y1) + d0 + d1])
        rr = dict(raw)
        for key in KEYS:
            if key != 'baryon_g': rr['state_coeff_' + key] = np.zeros_like(raw['state_coeff_' + key])
            if key not in KEYS[-2:]: rr['geometry_coeff_' + key] = np.zeros_like(raw['geometry_coeff_' + key])
        sc = np.zeros_like(raw['state_coeff_baryon_g']); sc[:, :, :ncell] = co; rr['state_coeff_baryon_g'] = sc
        folder = out/f'kernel-in-{label}-{tag}'; (folder/'gr').mkdir(parents=True); np.savez(folder/'gr/source-64.npz', **rr)
        m = ext.Response(); ext.bind(setup, INPUT=folder)(m, d, 8); free = ext.at(m, m.source, np.array([T])); shutil.rmtree(folder)
        return -float(free[0][-1, -1])/M
    def face_pattern(fidx):
        y = np.zeros((len(t), ncell)); yd = np.zeros_like(y)
        if fidx < ncell: y[:, fidx] += Phi[:, fidx]; yd[:, fidx] += Phidot[:, fidx]
        if fidx > 0: y[:, fidx - 1] -= Phi[:, fidx]; yd[:, fidx - 1] -= Phidot[:, fidx]
        return y, yd
if mode == 'faces':
    a0, b0 = int(sys.argv[4]), int(sys.argv[5]); rows = []
    for fidx in range(a0, b0 + 1):
        s0 = time.monotonic(); q = charge(*face_pattern(fidx), f'face{fidx}')
        rows.append(dict(face=fidx, r=float(fe[fidx]), depth_km=float((fe[-1] - fe[fidx])/1e5), rho0=float(rho[fidx]), contribution=q, seconds=time.monotonic() - s0))
        print(fidx, q, flush=True)
    (out/f'kernel-{label}-{a0}-{b0}.json').write_text(json.dumps(rows, indent=1) + '\n')
elif mode == 'full':
    y = Phi[:, :ncell] - Phi[:, 1:ncell + 1]; yd = Phidot[:, :ncell] - Phidot[:, 1:ncell + 1]
    q = charge(y, yd, 'full'); (out/f'kernel-{label}-full.json').write_text(json.dumps(dict(interior_freefall_charge=q), indent=1) + '\n'); print('full', q)
elif mode == 'analyze':
    rows = sorted((r for p_ in glob.glob(str(out/f'kernel-{label}-*-*.json')) for r in json.load(open(p_))), key=lambda r: r['face'])
    assert [r['face'] for r in rows] == list(range(ncell + 1)), [r['face'] for r in rows]
    full = json.load(open(out/f'kernel-{label}-full.json'))['interior_freefall_charge']
    c = np.array([r['contribution'] for r in rows]); depth = np.array([r['depth_km'] for r in rows]); rho0 = np.array([r['rho0'] for r in rows])
    q = c.sum(); K = c/q
    margin = 1/np.abs(K).sum()
    centroid = float(np.sum(np.abs(c)*depth)/np.sum(np.abs(c)))
    lnrho = np.log(rho0); dlnrho = np.gradient(lnrho, -depth*1e5)  # d ln rho / dr (r increases outward = -depth)
    fam = {}
    for dr_km in [-100, -30, -10, 10, 30, 100]:  # profile moved outward by dr (rho'(r) = rho(r - dr))
        rho_new = np.exp(np.interp(-(depth + dr_km), -depth, lnrho))  # rho at depth + dr_km (profile shifted outward)
        fam[f'shift {dr_km:+d} km'] = float(np.sum(c*rho_new/rho0))
    iref = int(np.argmin(abs(depth - centroid)))
    for eps in [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]:  # scale height H -> H(1+eps) about the centroid depth
        rho_new = np.exp(lnrho[iref] + (lnrho - lnrho[iref])/(1 + eps))
        fam[f'scale height {eps:+.0%}'] = float(np.sum(c*rho_new/rho0))
    for s in [-0.5, 0.5]: fam[f'uniform {s:+.0%}'] = float(q*(1 + s))
    res = dict(label=label, interior_freefall_charge=full, sum_of_face_contributions=float(q), additivity_relative=float(q/full - 1),
               sign_margin_Linf=float(margin), positive_share=float(c[c > 0].sum()/q), negative_share=float(c[c < 0].sum()/q), centroid_depth_km=centroid,
               families={k: dict(charge=v, relative=v/q - 1, sign_kept=bool(np.sign(v) == np.sign(q))) for k, v in fam.items()},
               faces=[dict(face=int(r['face']), depth_km=r['depth_km'], rho0=r['rho0'], contribution=r['contribution'], K=float(k_)) for r, k_ in zip(rows, K)])
    (out/f'kernel-{label}.json').write_text(json.dumps(res, indent=1) + '\n')
    print('interior free-fall charge %.6e, sum of faces %.6e (additivity %.1e)' % (full, q, q/full - 1))
    print('sign margin (L-inf relative density change needed to flip the sign): %.4f' % margin, ' centroid depth %.1f km' % centroid)
    print('positive/negative shares %.3f / %.3f' % (res['positive_share'], res['negative_share']))
    for k, v in res['families'].items(): print('%-22s q %.6e  rel %+.4f  sign kept %s' % (k, v['charge'], v['relative'], v['sign_kept']))
    for r in res['faces']:
        if abs(r['K']) > 0.01: print('face %2d depth %6.1f km  K %+.4f' % (r['face'], r['depth_km'], r['K']))
else: raise ValueError(mode)
