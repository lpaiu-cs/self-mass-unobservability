"""Phase279 part 2: non-gray LTE correction of the photospheric density at fixed geometric depth.

Conjectural. Plane-parallel, LTE, H+He composition (X=0.8887, Y=0.1073; metals and lines omitted), ideal Saha EOS.
Opacity per gram: H I bound-free n<=12 (hydrogenic, g_bf=1, stimulated emission), free-free of H II, He II, He III (Kramers, g_ff=1),
He I ground bound-free (threshold 24.587 eV, sigma ~ nu^-2), He II bound-free n<=4 (hydrogenic Z=2), Thomson scattering.
Transfer: Eddington closure per frequency, (1/3) J'' = eps (J - B), top H = J/2, bottom diffusion; coherent electron scattering.
Temperature: Unsold-Lucy correction  dB = (kJ/kB) J - B + (kJ/kB)[3 int_0^m chiH dH dm + 2 dH(0)],  dH = H0 - H.
Hydrostatics: dP_gas/dm = cx g - g_rad, g_rad = (4 pi/c) int chi_nu H_nu dnu; depth dz = dm/rho, origin at P_gas = 1 dyn/cm^2.
The gray reference uses the Rosseland mean of the same opacity with T^4 = (3/4) Teff^4 (tau_R + 2/3); only the ratio
rho_nongray(z)/rho_gray(z) is used (applied to the table-opacity gray atmospheres of part 1).
Self-checks: int B_nu dnu = sigma T^4/pi on the grid; the gray solver with a constant opacity reproduces the Eddington relation.
Usage: python3 phase279-nongray.py <out json>
"""
import json, sys, time
import numpy as np
h, k, c, me, mH, sigT = 6.62607015e-27, 1.380649e-16, 2.99792458e10, 9.1093837e-28, 1.6735575e-24, 6.6524587e-25
SIG = 5.670374419e-5; eV = 1.602176634e-12; RY = 13.605693 * eV
X, Y, CX = 0.8887014795256516, 0.10730973083975379, 1.0070231612220188
t0 = time.monotonic()
# frequency grid with both sides of every edge
edges = [RY/n**2/h for n in range(1, 13)] + [24.587*eV/h] + [4*RY/n**2/h for n in range(1, 5)]
nu = np.geomspace(1e13, 2.5e16, 520)
nu = np.unique(np.r_[nu, [e*(1 - 1e-4) for e in edges], [e*(1 + 1e-4) for e in edges]])
w_nu = np.gradient(nu)  # trapezoid-like weights


def planck(T):  # (n_depth, n_freq)
    x = np.minimum(h*nu[None, :]/(k*T[:, None]), 700)
    return 2*h*nu[None, :]**3/c**2/np.expm1(x)


def saha(rho, T):
    """Ionization equilibrium for H and He at given rho, T: returns n_e, n_HI levels (n<=12), n_HII, n_HeI, n_HeII(levels n<=4), n_HeIII."""
    nH, nHe = rho*X/mH, rho*Y/(4*mH)
    lam = (2*np.pi*me*k*T/h**2)**1.5
    n = np.arange(1, 13); En = RY*(1 - 1/n**2)
    UH = np.sum(2*n[None, :]**2*np.exp(-En[None, :]/(k*T[:, None])), axis=1)
    SH = lam*(1/UH)*np.exp(-RY/(k*T))  # nHII ne / nHI  (U_HII = 1, g factor 2 cancels with electron spin: 2 UHII/UHI)
    SH = 2*SH
    SHe1 = lam*2*(2/1)*np.exp(-24.587*eV/(k*T))  # nHeII ne / nHeI (U_HeII = 2, U_HeI = 1)
    SHe2 = lam*2*(1/2)*np.exp(-54.418*eV/(k*T))  # nHeIII ne / nHeII
    lo, hi = np.full_like(T, 1e-30), nH + 2*nHe
    for _ in range(200):  # bisection on log n_e
        ne = np.sqrt(lo*hi)
        xH = SH/(SH + ne); a1, a2 = SHe1/ne, SHe1*SHe2/ne**2; den = 1 + a1 + a2
        f = nH*xH + nHe*(a1 + 2*a2)/den - ne
        lo, hi = np.where(f > 0, ne, lo), np.where(f > 0, hi, ne)
    ne = np.sqrt(lo*hi); xH = SH/(SH + ne); a1, a2 = SHe1/ne, SHe1*SHe2/ne**2; den = 1 + a1 + a2
    nHI, nHII = nH*(1 - xH), nH*xH; nHeI, nHeII, nHeIII = nHe/den, nHe*a1/den, nHe*a2/den
    levH = nHI[:, None]*2*n[None, :]**2*np.exp(-En[None, :]/(k*T[:, None]))/UH[:, None]
    m4 = np.arange(1, 5); E4 = 4*RY*(1 - 1/m4**2); U4 = np.sum(2*m4[None, :]**2*np.exp(-E4[None, :]/(k*T[:, None])), axis=1)
    levHe2 = nHeII[:, None]*2*m4[None, :]**2*np.exp(-E4[None, :]/(k*T[:, None]))/U4[:, None]
    return dict(ne=ne, levH=levH, nHII=nHII, nHeI=nHeI, nHeII=nHeII, nHeIII=nHeIII, levHe2=levHe2, n_tot=nH + nHe + ne)


def rho_from_Pg(Pg, T):  # ideal gas with ionization: P = (n_nuclei + n_e) k T
    rho = Pg*mH/(k*T)
    for _ in range(60):
        s = saha(rho, T); new = rho*Pg/(s['n_tot']*k*T)
        if np.max(np.abs(new/rho - 1)) < 1e-12: break
        rho = new
    return new


def opacity(rho, T):
    """absorption and scattering per gram, (n_depth, n_freq)."""
    s = saha(rho, T); stim = -np.expm1(-np.minimum(h*nu[None, :]/(k*T[:, None]), 700))
    ab = np.zeros((len(T), len(nu)))
    for n in range(1, 13):
        nun = RY/n**2/h; sig = np.where(nu >= nun, 7.906e-18*n*(nun/nu)**3, 0.)
        ab += s['levH'][:, n - 1:n]*sig[None, :]
    for m_ in range(1, 5):
        nun = 4*RY/m_**2/h; sig = np.where(nu >= nun, 7.906e-18*m_/4*(nun/nu)**3, 0.)  # hydrogenic, Z=2
        ab += s['levHe2'][:, m_ - 1:m_]*sig[None, :]
    nu1 = 24.587*eV/h; ab += s['nHeI'][:, None]*np.where(nu >= nu1, 7.4e-18*(nu1/nu)**2, 0.)[None, :]
    ab *= stim
    ions = s['nHII'] + s['nHeII'] + 4*s['nHeIII']
    ab += 3.692e8*(s['ne']*ions)[:, None]/np.sqrt(T)[:, None]/nu[None, :]**3*stim
    return ab/rho[:, None], (s['ne']*sigT/rho)[:, None]*np.ones((1, len(nu)))


def rosseland(ab, sc, T):
    x = np.minimum(h*nu[None, :]/(k*T[:, None]), 700); dB = x*np.exp(x)/np.expm1(x)**2*nu[None, :]**3  # dB/dT shape
    return np.sum(dB*w_nu, 1)/np.sum(dB/(ab + sc)*w_nu, 1)


def solve_J(T, ab, sc, m):
    """Eddington transfer per frequency on the column-mass grid; returns J (nd, nf) and H at half points (nd-1, nf)."""
    chi = ab + sc; eps = ab/chi; B = planck(T)
    dtau = 0.5*(chi[1:] + chi[:-1])*np.diff(m)[:, None]; nd = len(m)
    a = np.zeros_like(chi); b = np.zeros_like(chi); cc = np.zeros_like(chi); d = np.zeros_like(chi)
    Dl, Du = dtau[:-1], dtau[1:]; Dm = 0.5*(Dl + Du)
    a[1:-1] = 1/(3*Dl*Dm); cc[1:-1] = 1/(3*Du*Dm); b[1:-1] = -(a[1:-1] + cc[1:-1]) - eps[1:-1]; d[1:-1] = -eps[1:-1]*B[1:-1]
    b[0] = -1/(3*dtau[0]) - 0.5 - dtau[0]/2*eps[0]; cc[0] = 1/(3*dtau[0]); d[0] = -dtau[0]/2*eps[0]*B[0]
    a[-1] = -1.; b[-1] = 1.; d[-1] = B[-1] - B[-2]
    for i in range(1, nd):  # Thomas algorithm, vectorized over frequency
        w = a[i]/b[i - 1]; b[i] = b[i] - w*cc[i - 1]; d[i] = d[i] - w*d[i - 1]
    J = np.zeros_like(chi); J[-1] = d[-1]/b[-1]
    for i in range(nd - 2, -1, -1): J[i] = (d[i] - cc[i]*J[i + 1])/b[i]
    Hh = (J[1:] - J[:-1])/(3*dtau)
    return J, Hh, B


DAMP, CLIP = 0.5, 0.02
def atmosphere(Teff, logg, gray, m=None, iters=1500, T_init=None, verbose=False):
    g = CX*10**logg; H0 = SIG*Teff**4/(4*np.pi)
    m = np.geomspace(1e-9, 60., 200) if m is None else m
    # start: the converged gray model (non-gray) or an isothermal guess (gray)
    T = (0.75*Teff**4*(np.full_like(m, 1.) + 2/3))**0.25 if T_init is None else T_init.copy()
    Pg = g*m; rho = rho_from_Pg(Pg, T)
    hist = []
    for it in range(iters):
        ab, sc = opacity(rho, T); kR = rosseland(ab, sc, T)
        tauR = np.r_[0., np.cumsum(0.5*(kR[1:] + kR[:-1])*np.diff(m))] + kR[0]*m[0]
        if gray:
            Tn = (0.75*Teff**4*(tauR + 2/3))**0.25; grad = kR*4*np.pi*H0/c; flux_err = 0.
        else:
            J, Hh, B = solve_J(T, ab, sc, m)
            Hn = np.r_[J[0:1]*0.5, 0.5*(Hh[1:] + Hh[:-1]), Hh[-1:]]  # node flux per frequency
            Htot = np.sum(Hn*w_nu, 1); dH = H0 - Htot
            kB = np.sum(ab*B*w_nu, 1)/np.sum(B*w_nu, 1); kJ = np.sum(ab*J*w_nu, 1)/np.sum(J*w_nu, 1)
            chiH = np.sum((ab + sc)*np.abs(Hn)*w_nu, 1)/np.sum(np.abs(Hn)*w_nu, 1)
            Jt, Bt = np.sum(J*w_nu, 1), np.sum(B*w_nu, 1)
            integ = np.r_[0., np.cumsum(0.5*(chiH[1:]*dH[1:] + chiH[:-1]*dH[:-1])*np.diff(m))] + chiH[0]*dH[0]*m[0]
            dB = kJ/kB*Jt - Bt + kJ/kB*(3*integ + 2*dH[0])
            Bnew = np.maximum(Bt + dB, 0.3*Bt)
            Tn = (np.pi*Bnew/SIG)**0.25; Tn = T*(1 + DAMP*np.clip(Tn/T - 1, -CLIP, CLIP))
            grad = 4*np.pi/c*np.sum((ab + sc)*Hn*w_nu, 1); flux_err = float(np.max(np.abs(dH[:-5]/H0)))  # bottom diffusion nodes excluded
        # hydrostatics in column mass
        geff = g - grad; assert np.all(geff > 0), 'super-Eddington'
        Pg_n = np.r_[0., np.cumsum(0.5*(geff[1:] + geff[:-1])*np.diff(m))] + geff[0]*m[0]
        dT = float(np.max(np.abs(Tn/T - 1))); T = Tn if not gray else Tn
        rho = rho_from_Pg(Pg_n, T); Pg = Pg_n
        hist.append(dict(it=it, dT=dT, flux_err=flux_err))
        if verbose and it % 200 == 0: print('   it %4d dT %.2e flux %.2e T0 %.1f' % (it, dT, flux_err, T[0]), flush=True)
        if gray and dT < 1e-10 and it > 5: break
        if not gray and it > 50 and dT < 1e-7 and flux_err < 1e-4: break
    z = np.r_[0., np.cumsum(0.5*(1/rho[1:] + 1/rho[:-1])*np.diff(m))]
    z0 = float(np.interp(0., np.log(Pg), z))
    return dict(m=m, T=T, Pg=Pg, rho=rho, z=z - z0, tauR=tauR, kR=kR, iterations=len(hist), final=hist[-1], Teff=Teff, logg=logg)


# self-checks
Tt = np.array([8e3, 16e3, 4e4]); Bint = np.sum(planck(Tt)*w_nu, 1); chk_planck = float(np.max(np.abs(Bint/(SIG*Tt**4/np.pi) - 1)))
assert chk_planck < 5e-3, chk_planck
res = dict(planck_grid_check=chk_planck, n_freq=len(nu), cases={})
lay = np.linspace(0., 700e5, 701)
for name, te, lg in [('declared', 15898.27, 5.73962), ('Kaplan central', 15800., 5.82), ('log g -1 sigma', 15800., 5.77), ('log g +1 sigma', 15800., 5.87),
                     ('log g -2 sigma', 15800., 5.72), ('log g +2 sigma', 15800., 5.92), ('Teff -1 sigma', 15700., 5.82), ('Teff +1 sigma', 15900., 5.82)]:
    G = atmosphere(te, lg, gray=True); N = atmosphere(te, lg, gray=False, T_init=G['T'])
    lr = lambda A, z: np.interp(z, A['z'], np.log(A['rho']))
    ratio = np.exp(lr(N, lay) - lr(G, lay))
    res['cases'][name] = dict(Teff=te, logg=lg, gray_iterations=G['iterations'], nongray_iterations=N['iterations'], nongray_final=N['final'],
        T0_over_Teff_gray=float(G['T'][0]/te), T0_over_Teff_nongray=float(N['T'][0]/te),
        depth_tau1_gray_km=float(np.interp(0., np.log(G['tauR']), G['z'])/1e5), depth_tau1_nongray_km=float(np.interp(0., np.log(N['tauR']), N['z'])/1e5),
        rho_ratio_z=dict(z_km=(lay/1e5).tolist(), ratio=ratio.tolist()),
        T_ratio_at_tau=[dict(tauR=tq, ratio=float(np.interp(np.log(tq), np.log(N['tauR']), N['T'])/np.interp(np.log(tq), np.log(G['tauR']), G['T']))) for tq in [1e-3, 1e-2, 0.1, 1., 10.]])
    np.savez_compressed(sys.argv[1].replace('.json', f'-{name.replace(" ", "_")}.npz'), **{f'gray_{kk}': v for kk, v in G.items() if isinstance(v, np.ndarray)}, **{f'nongray_{kk}': v for kk, v in N.items() if isinstance(v, np.ndarray)})
    print('%-15s gray it %d, nongray it %d (flux %.1e, dT %.1e); T0/Teff gray %.4f nongray %.4f; tau=1 depth gray %.1f nongray %.1f km; rho ratio 240/363/465 km %.3f %.3f %.3f; %.0fs' % (
        name, G['iterations'], N['iterations'], N['final']['flux_err'], N['final']['dT'], res['cases'][name]['T0_over_Teff_gray'], res['cases'][name]['T0_over_Teff_nongray'],
        res['cases'][name]['depth_tau1_gray_km'], res['cases'][name]['depth_tau1_nongray_km'], ratio[240], ratio[363], ratio[465], time.monotonic() - t0), flush=True)
json.dump(res, open(sys.argv[1], 'w'), indent=1)
