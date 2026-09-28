"""Phase278: whole-star transit and long-time relaxation of the compact charge (linear adiabatic radial model).

Conjectural model (declared background only, Newtonian self-gravity, compactness 4e-6):
  M zeta'' + K zeta = F(t),  zeta = xi/r,  phase-273 finite elements (consistent mass, natural boundaries),
  F = int rho r^3 f v dr,  f = -c^2 (a^2/B^2) d/dr[alpha(phi0) dphi] = 4 c^2 (a^2/B^2) d/dr[phi0 dphi],
  dphi = eta R_d [s(p_in) - s(p_out)]/r, p_in = (t + x/c)/D, p_out = (t - (x - 2 x0)/c)/D, s(p) = (4p(1-p))^4 on (0,1),
  x = -int_r^{R_d} (B/a) dr (optical coordinate), x0 = x(0): the regular spherical pulse reflected through the centre.
Transit (0-1.5 s): average-acceleration Newmark; element mass perturbations dM_e = Phi_l - Phi_r, Phi = 4 pi r^2 B rho0 xi,
and their running time integrals I_e(t). Relaxation (t >= 1 s): full eigen-decomposition of (K, M), analytic free modes.
Readout (flat exterior, exact monopole): q(u) = (1/M) sum_e w_e (c/2r_e) [I_e(u + r_e/c) - I_e(u - r_e/c)], w = alpha(phi0) a,
u = t_R - R_d/c (t_R: readout time at R_d, the chain's readout radius). Free modes: q(u) = sum_k Q_k c_k(u),
Q_k = (1/M) sum_e w_e dM_{e,k} sinc(omega_k r_e/c). Static-window limit: (1/M) sum_e w_e dM_e.
Usage: python3 phase278-transit.py <out json>
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.linalg import eigh
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
C, G = 2.99792458e10, 6.67430e-8
t0 = time.monotonic(); log = lambda *a: print('%7.1fs' % (time.monotonic() - t0), *a, flush=True)
bg = bl.chem.prior.Background(); d = bg.d
keys = ['radius_cm', 'mass_geom_cm', 'pressure_cgs', 'density_cgs', 'gamma1', 'phi', 'phi_prime_cm', 'lapse']
r, mg, P, rho, G1, phi, dphi, a = (np.asarray(d[k], float) for k in keys)
o = np.argsort(r); r, mg, P, rho, G1, phi, dphi, a = (v[o] for v in (r, mg, P, rho, G1, phi, dphi, a))
k = np.r_[True, np.diff(r) > 0]; r, mg, P, rho, G1, phi, dphi, a = (v[k] for v in (r, mg, P, rho, G1, phi, dphi, a))
R = float(bg.R); n = len(r)
fs = np.load('readout268-quad64-work/gr/field-source-64.npz')
Rd, D, eta, T, Mcm = float(fs['drive_radius']), float(fs['drive_duration']), float(fs['drive_amplitude']), float(fs['t'][-1]), float(fs['M_cm'])
Mg = C*C*Mcm/G
B = np.where(r > 0, 1/np.sqrt(np.clip(1 - 2*mg/np.maximum(r, 1.), 1e-12, None)), 1.)
xp = B/a
x = np.empty(n); x[-1] = -(Rd - R)*xp[-1]
x[:-1] = x[-1] - np.cumsum(((xp[1:] + xp[:-1])/2*np.diff(r))[::-1])[::-1]
x0 = float(x[0])
# finite elements (phase 273 operator, tridiagonal)
h = np.diff(r); rm = (r[:-1] + r[1:])/2
A = np.interp(rm, r, G1*P*r**4); Q = np.interp(rm, r, r**3*np.gradient((3*G1 - 4)*P, r)); Mc = np.interp(rm, r, rho*r**4)
Kd = np.zeros(n); Ko = np.zeros(n - 1); Md = np.zeros(n); Mo = np.zeros(n - 1)
Kd[:-1] += A/h - Q*h/3; Kd[1:] += A/h - Q*h/3; Ko += -A/h - Q*h/6
Md[:-1] += Mc*h/3; Md[1:] += Mc*h/3; Mo += Mc*h/6
Ks = sp.diags([Ko, Kd, Ko], [-1, 0, 1], format='csc'); Ms = sp.diags([Mo, Md, Mo], [-1, 0, 1], format='csc')
# force ingredients at element midpoints
ip = lambda v: np.interp(rm, r, v)
phm, dphm, am, Bm, xm, rhom, xpm = ip(phi), ip(dphi), ip(a), ip(B), ip(x), ip(rho), ip(xp)
wload = rhom*rm**3*h/2  # consistent (midpoint) load weight per node


def s(p): return np.where((p > 0) & (p < 1), (4*p*(1 - p))**4, 0.)
def sp_(p): return np.where((p > 0) & (p < 1), 4*(4*p*(1 - p))**3*(4 - 8*p), 0.)


def force(t):
    pin = (t + xm/C)/D; pout = (t - (xm - 2*x0)/C)/D
    u = s(pin) - s(pout); du = (sp_(pin) + sp_(pout))*xpm/(C*D)
    dph = eta*Rd*u/rm; ddph = eta*Rd*(du/rm - u/rm**2)
    f = 4*C*C*(am**2/Bm**2)*(dphm*dph + phm*ddph)
    F = np.zeros(n); fe = f*wload; F[:-1] += fe; F[1:] += fe
    return F


Phi_w = 4*np.pi*r**2*B*rho*r  # Phi = Phi_w * zeta
rme = np.r_[rm, R]  # element centres plus the pseudo-element outside R that receives the surface mass flux
wE = -4*np.r_[np.interp(rm, r, phi)*am, phi[-1]*a[-1]]  # w = alpha(phi0) a
def dM(z): Ph = Phi_w*z; return np.r_[Ph[:-1] - Ph[1:], Ph[-1]]  # mass conserving: sum = Phi(0) = 0
rmax_c = float(rm.max()/C)

# ---- transit: Newmark ----
schedule = [(0.5, 1e-5), (1.5, 1e-4)]
zeta = np.zeros(n); vel = np.zeros(n); acc = np.zeros(n); tt = 0.; I = np.zeros(n); dMold = np.zeros(n); static_hist = [0.]
store_t = [0.]; store_I = [I.copy()]; ffree = None; mid_state = None
for t_end, dt in schedule:
    Keff = splu((Ks + (4/dt**2)*Ms).tocsc())
    nsteps = int(round((t_end - tt)/dt))
    for i in range(nsteps):
        tn = tt + dt
        Fn = force(tn) if tn < (2*abs(x0))/C + 2*D else np.zeros(n)
        rhs = Fn + Ms @ ((4/dt**2)*zeta + (4/dt)*vel + acc)
        zn = Keff.solve(rhs); an = (4/dt**2)*(zn - zeta) - (4/dt)*vel - acc; vel = vel + dt/2*(acc + an); acc = an; zeta = zn; tt = tn
        dMn = dM(zeta); I += dt/2*(dMold + dMn); dMold = dMn
        if tt < 0.01 or (tt < 0.8 and (i + 1) % max(1, int(round(1e-4/dt))) == 0) or (i + 1) % max(1, int(round(5e-4/dt))) == 0:
            store_t.append(tt); store_I.append(I.copy()); static_hist.append(float(np.sum(wE*dMn))/Mg)
        if abs(tt - T) < dt/2: ffree = dict(t=tt, dM=dMn.copy(), zeta=zeta.copy())
        if abs(tt - 1.0) < dt/2: mid_state = (tt, zeta.copy(), vel.copy())
    log('transit to %.3f s with dt %.0e done' % (tt, dt))
store_t = np.array(store_t); store_I = np.array(store_I); static_hist = np.array(static_hist); zeta1, vel1 = zeta.copy(), vel.copy()
log('stored I', store_I.shape)


def I_at(times):  # I_e at per-element times (vector over elements), linear interpolation in the stored table
    times = np.asarray(times); j = np.clip(np.searchsorted(store_t, times) - 1, 0, len(store_t) - 2)
    w1 = np.clip((times - store_t[j])/(store_t[j + 1] - store_t[j]), 0, 1); cols = np.arange(len(times))
    out = store_I[j, cols]*(1 - w1) + store_I[j + 1, cols]*w1
    return np.where(times <= 0, 0., out)


def q_ret(tR):
    u = tR - Rd/C
    return float(np.sum(wE*(C/(2*rme))*(I_at(u + rme/C) - I_at(u - rme/C)))/Mg)


from numpy.polynomial import polynomial as Pn
sc = Pn.polypow([0, 4, -4], 4); F1c = Pn.polyint(sc); F2c = Pn.polyint(F1c)
def Fk(c, p): return np.where(p <= 0, 0., np.where(p < 1, Pn.polyval(np.clip(p, 0, 1), c), Pn.polyval(1., c) + (Pn.polyval(1., Pn.polyder(c)) if c is F2c else 0.)*(p - 1)))
def xi_ff(t):
    p = (t + x/C)/D; rr = np.maximum(r, 1.)
    return 4*C*C*(a**2/B**2)*eta*Rd*D*D*((dphi/rr - phi/rr**2)*Fk(F2c, p) + phi*xp*Fk(F1c, p)/(rr*C*D))
res = dict(R=R, R_d=Rd, D=D, eta=eta, T=T, M_g=Mg, x0=x0, center_arrival_s=-x0/C, exit_s=-2*x0/C + D, nodes=n)
# validation at T and the static-window value
qT = q_ret(T); lay = (R - r >= 240e5) & (R - r <= 465e5); xf = xi_ff(ffree['t']); xz = ffree['zeta']*r
res['validation'] = dict(freefall_displacement_rel_Linf_charge_layers=float(np.max(np.abs(xz[lay] - xf[lay]))/np.max(np.abs(xf[lay]))), q_T=qT, chain_continuum_mass_normalized=-2.3350180598198245e-51,
                                        relative=qT/-2.3350180598198245e-51 - 1, static_window_T=float(np.sum(wE*ffree['dM'])/Mg) if ffree else None)
log('validation q(T) %.6e  rel %+.3e  static %.3e' % (qT, res['validation']['relative'], res['validation']['static_window_T']))
tR_grid = np.r_[np.linspace(0, 0.02, 81)[1:], np.arange(0.02, 1.5, 2e-4)[1:]]
qR = np.array([q_ret(t) for t in tR_grid])
sel = np.r_[np.arange(0, 80), np.arange(80, len(tR_grid), 10)]
res['transit'] = dict(tR=tR_grid[sel].tolist(), q=qR[sel].tolist(), static_t=store_t[::20].tolist(), static_q=static_hist[::20].tolist(),
                      static_max_abs=float(np.max(np.abs(static_hist))), q_max_abs=float(np.max(np.abs(qR))), t_max_abs=float(tR_grid[np.argmax(np.abs(qR))]),
                      q_after_exit_window=float(qR[np.searchsorted(tR_grid, 0.75)]), q_min=float(qR.min()), q_max=float(qR.max()),
                      sign_changes=int(np.sum(np.diff(np.sign(qR[qR != 0])) != 0)))
log('transit max |q| %.3e at %.3f s, sign changes %d' % (res['transit']['q_max_abs'], res['transit']['t_max_abs'], res['transit']['sign_changes']))

# ---- relaxation: modes ----
Kdense = Ks.toarray(); Mdense = Ms.toarray()
w2, V = eigh(Kdense, Mdense); del Kdense, Mdense
log('eigh done, omega^2 range %.3e .. %.3e' % (w2[0], w2[-1]))
assert w2[0] > 0
om = np.sqrt(w2); t1 = 1.0
# state at t1 from the Newmark run: rerun cheaply is not needed; use the stored end state if t1 == end, else project the end state
t_state = tt
ak = V.T @ (Ms @ zeta1); bk = (V.T @ (Ms @ vel1))/om
PV = Phi_w[:, None]*V; dMk = np.vstack([PV[:-1] - PV[1:], PV[-1:]]); del PV
Qk = ((wE[:, None]*dMk)*np.sinc((om[None, :]*rme[:, None]/C)/np.pi)).sum(0)/Mg
Qs = (wE[:, None]*dMk).sum(0)/Mg
def q_modal(tR, Qv=Qk, a_=None, b_=None, ts=None):
    a_ = ak if a_ is None else a_; b_ = bk if b_ is None else b_; ts = t_state if ts is None else ts
    u = np.asarray(tR) - Rd/C; out = np.empty(len(u))
    for i in range(0, len(u), 400):
        ph = om[None, :]*(u[i:i + 400, None] - ts); out[i:i + 400] = (np.cos(ph)*a_ + np.sin(ph)*b_) @ Qv
    return out
def q_mean(tA, tB):  # exact time average of the free-mode charge over readout times [tA, tB]
    ua, ub = tA - Rd/C - t_state, tB - Rd/C - t_state
    return float(np.sum(Qk*(ak*(np.sin(om*ub) - np.sin(om*ua)) - bk*(np.cos(om*ub) - np.cos(om*ua)))/om)/(tB - tA))
if mid_state is not None:  # modal propagation 1.0 s -> end versus Newmark at the end (time-discretization check)
    tm, zm, vm = mid_state; am_ = V.T @ (Ms @ zm); bm_ = (V.T @ (Ms @ vm))/om
    z_prop = V @ (am_*np.cos(om*(t_state - tm)) + bm_*np.sin(om*(t_state - tm)))
    res['relaxation_check'] = dict(static_window_newmark=float(np.sum(wE*dM(zeta1))/Mg), static_window_modal=float(np.sum(wE*dM(z_prop))/Mg),
                                   zeta_rel_L2=float(np.linalg.norm(z_prop - zeta1)/np.linalg.norm(zeta1)))
over = np.array([1.5])  # consistency at the end of the transit table (window fully after t_state needs u - r/c >= t_state)
res['relaxation'] = dict(t_state=t_state, omega_lowest=om[:6].tolist(), period_lowest=(2*np.pi/om[:6]).tolist(), energy_modal=float(0.5*np.sum(om**2*(ak**2 + bk**2))))
lateR = np.r_[np.linspace(t_state + Rd/C + rmax_c + 0.01, 60, 3000), np.geomspace(60, 1e4, 3000)[1:]]
qL = q_modal(lateR); qLs = q_modal(lateR, Qs)
amp = np.abs(Qk)*np.sqrt(ak**2 + bk**2); top = np.argsort(amp)[::-1][:8]
res['relaxation'].update(tR_first=float(lateR[0]), q_first=float(qL[0]), q_rms=float(np.sqrt(np.mean(qL**2))), q_max_abs=float(np.max(np.abs(qL))),
    mean_9000_10000s=q_mean(9000., 1e4), mean_first_to_1e4=q_mean(float(lateR[0]), 1e4), amplitude_sum=float(amp.sum()),
    dominant_modes=[dict(k=int(i), period_s=float(2*np.pi/om[i]), amplitude=float(amp[i]), Q=float(Qk[i])) for i in top],
    static_window_rms=float(np.sqrt(np.mean(qLs**2))), samples=dict(tR=lateR[::60].tolist(), q=qL[::60].tolist()))
log('relaxation first q %.3e, rms %.3e, max %.3e, mean(9000-1e4 s) %.3e' % (qL[0], res['relaxation']['q_rms'], res['relaxation']['q_max_abs'], res['relaxation']['mean_9000_10000s']))
json.dump(res, open(sys.argv[1], 'w'), indent=1)
log('done')
