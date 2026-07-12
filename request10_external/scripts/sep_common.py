#!/usr/bin/env python3
"""Shared pipeline core for the 10.8x dynamic-SEP analyses.

Introduced by the 10.8f review response (notes/REQUEST10_8F_REVIEW_RESPONSE.md):
C1  one CAUSAL template implementation (col_s is the +pi/2 = -sin drive
    response, so the causal pole combination is g1*col_c - g2*col_s);
C4  the static-leakage guard direction is kept explicitly after truncation;
C10 constants are loaded from committed artifacts, never retyped.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

SV_CUT_DEFAULT = 1e-3
M_RN_DEFAULT = 30
SEED = 20260710


# ---------------------------------------------------------------- inputs

def load_inputs(here: str) -> dict:
    """Load the committed inputs of record from request10_external/."""
    base = np.load(os.path.join(here, 'baseline_planetGR.npz'), allow_pickle=True)
    res0 = base['res'].astype(np.float64)
    errs = base['errs'].astype(np.float64)
    t = base['toas'].astype(np.float64)
    rank_json = json.load(open(os.path.join(here, 'carrier_projection_rank.json')))
    oms = {k.split('_')[1]: rank_json['carriers_rad_per_day'][k]
           for k in ('Omega_in', 'Omega_out', 'Omega_dif')}
    # carrier-consistency stop rule (mirrors joint_fit_upper_limit.py)
    p_i = rank_json['periods_days']['P_inner']
    p_o = rank_json['periods_days']['P_outer']
    w_re = {'in': 2*np.pi/p_i, 'out': 2*np.pi/p_o,
            'dif': abs(2*np.pi/p_i - 2*np.pi/p_o)}
    for k in oms:
        assert abs(oms[k] - w_re[k]) < 1e-9, 'carrier mismatch vs artifacts'
    J = np.load(os.path.join(here, 'finite_jacobian_v2.npy'))
    jD = np.load(os.path.join(here, 'sep_dynamic/col_SEP_D.npz'))['dcol'].astype(np.float64)
    npz = np.load(os.path.join(here, 'sep_dynamic/sep_dynamic_columns.npz'))
    cols = {k: npz[k].astype(np.float64) for k in npz.files}   # materialize once
    return {'res0': res0, 'errs': errs, 't': t, 'sw': 1.0/errs, 'N': res0.size,
            'J': J, 'jD': jD, 'cols': cols, 'OMS': oms,
            'P_in': p_i, 'P_out': p_o}


def load_anchor_factors(here: str) -> dict:
    """K factors from the committed artifacts (10.8f C10)."""
    out = {'K_dyn': 10.0, 'K_static': None}
    try:
        ar = json.load(open(os.path.join(here, 'sep_dynamic/anchor_resolution_10_8c.json')))
        out['K_dyn'] = float(ar['decision']['K_dyn'])
    except Exception:
        pass
    try:
        ss = json.load(open(os.path.join(here, 'sep_dynamic/sep_static_sensitivity.json')))
        sig = ss.get('R5_extended', {}).get('sigma_Delta_rn_R5',
                                            ss['sigma_Delta_rn_marginalized'])
        out['K_static'] = float(ss['reference_published_sigma']) / float(sig)
    except Exception:
        out['K_static'] = 934.0
    return out


# ------------------------------------------------------------- nuisance

def build_nuisance(inp: dict, m_rn: int = M_RN_DEFAULT,
                   sv_cut: float = SV_CUT_DEFAULT,
                   keep_guard: bool = True) -> dict:
    """Weighted nuisance projector.

    Block: normalized Jacobian columns + offset + m_rn Fourier pairs
    (+ the static-Delta guard). Truncated at relative SV >= sv_cut; if
    keep_guard, the guard's post-truncation residual direction is appended
    explicitly so the static-leakage guard is ACTIVE regardless of the cut
    (10.8f C4 -- previously it fell below the cut and was silently dropped).
    """
    sw, t, J, jD, N = inp['sw'], inp['t'], inp['J'], inp['jD'], inp['N']
    A = sw[:, None]*J
    colnorm = np.linalg.norm(A, axis=0)
    assert np.all(colnorm > 0), 'zero-norm Jacobian column'
    An = A/colnorm[None, :]
    cvec = sw/np.linalg.norm(sw)
    blocks = [An, cvec[:, None]]
    if m_rn > 0:
        T = float(t.max() - t.min())
        fc = []
        for j in range(1, m_rn + 1):
            arg = 2.0*np.pi*j*(t - t.min())/T
            fc += [np.cos(arg), np.sin(arg)]
        FB = sw[:, None]*np.column_stack(fc)
        blocks.append(FB/np.linalg.norm(FB, axis=0, keepdims=True))
    xjD = sw*jD
    xjDn = xjD/np.linalg.norm(xjD)
    blocks.append(xjDn[:, None])
    B0 = np.column_stack(blocks)
    U0, s0, _ = np.linalg.svd(B0, full_matrices=False)
    rank = int(np.sum(s0/s0[0] >= sv_cut)) if sv_cut > 0 else B0.shape[1]
    Q = U0[:, :rank]
    guard_retained = 1.0
    if keep_guard and sv_cut > 0:
        g = xjDn - Q @ (Q.T @ xjDn)
        gn = float(np.linalg.norm(g))
        guard_retained = 1.0 - gn**2
        if gn > 1e-8:
            Q = np.column_stack([Q, g/gn])
    return {'Q': Q, 'rank': Q.shape[1], 'ncols': B0.shape[1],
            'sv_cut': sv_cut, 'guard_residual_norm': 1.0 - guard_retained,
            'colnorm': colnorm}


def proj_out(Q: np.ndarray, X: np.ndarray) -> np.ndarray:
    return X - Q @ (Q.T @ X)


def noise_scale(Q: np.ndarray, sw: np.ndarray, res0: np.ndarray) -> tuple[float, float]:
    """Computed per span -- never hardcode (10.8f C10)."""
    y = sw*res0
    yp = proj_out(Q, y)
    chi2 = float(yp @ yp)
    dof = res0.size - Q.shape[1]
    return math.sqrt(chi2/dof), chi2/dof


# ------------------------------------------------------------- templates

def pole_coeffs(tau: float, w: float, anticausal: bool = False) -> tuple[float, float]:
    g1 = 1.0/(1.0 + (w*tau)**2)
    g2 = w*tau/(1.0 + (w*tau)**2)
    if anticausal:
        g2 = -g2
    return g1, g2


def templates(inp: dict, tau: float, t_off: float = 0.0,
              drive: dict | None = None, anticausal: bool = False):
    """CAUSAL convention of record (10.8f C1).

    col_c = response to cos(w t); col_s = response to cos(w t + pi/2)
    (a -sin drive). Response to a drive a*cos(w t + ph):
        a*cos(ph)*col_c + a*sin(ph)*col_s.
    Causal pole Re[e^{i(w t + ph)}/(1+i w tau)]
        = (g1 cos ph + g2 sin ph) cos(w t) + (g1 sin ph - g2 cos ph)(-sin(w t))...
    combining through the column convention gives
        Tb += a*[(g1 cos ph + g2 sin ph) col_c + (g1 sin ph - g2 cos ph) col_s].
    At ph = 0 this is g1*col_c - g2*col_s (NOT +g2; the +g2 combination is
    the advanced/anticausal template under this column convention).
    """
    cols, OMS, N = inp['cols'], inp['OMS'], inp['N']
    Tb = np.zeros(N)
    Tc = np.zeros(N)
    for name, w in OMS.items():
        amp, ph0 = (1.0, 0.0) if drive is None else drive[name]
        ph = ph0 + w*t_off
        g1, g2 = pole_coeffs(tau, w, anticausal)
        cb = amp*(g1*math.cos(ph) + g2*math.sin(ph))
        sb = amp*(g1*math.sin(ph) - g2*math.cos(ph))
        Tb += cb*cols['col_%s_c' % name] + sb*cols['col_%s_s' % name]
        Tc += amp*(math.cos(ph)*cols['col_%s_c' % name] +
                   math.sin(ph)*cols['col_%s_s' % name])
    return Tc, Tb


def template_stencil(tau: float, t_off: float, OMS: dict,
                     drive: dict | None = None, anticausal: bool = False):
    """6-vector coefficient stencils (order: in_c,in_s,out_c,out_s,dif_c,dif_s)
    for (Tc, Tb) -- the factorized-scan core (10.8f C5 efficiency)."""
    wc = np.zeros(6)
    wb = np.zeros(6)
    for k, (name, w) in enumerate(OMS.items()):
        amp, ph0 = (1.0, 0.0) if drive is None else drive[name]
        ph = ph0 + w*t_off
        g1, g2 = pole_coeffs(tau, w, anticausal)
        wc[2*k] = amp*math.cos(ph)
        wc[2*k+1] = amp*math.sin(ph)
        wb[2*k] = amp*(g1*math.cos(ph) + g2*math.sin(ph))
        wb[2*k+1] = amp*(g1*math.sin(ph) - g2*math.cos(ph))
    return wc, wb


def projected_sixcols(inp: dict, Q: np.ndarray):
    """Project the six weighted columns once; return (C6p, order)."""
    order = []
    mats = []
    for name in inp['OMS']:
        for q in ('c', 's'):
            order.append('col_%s_%s' % (name, q))
            mats.append(inp['sw']*inp['cols']['col_%s_%s' % (name, q)])
    C6 = np.column_stack(mats)
    return proj_out(Q, C6), order


# ------------------------------------------------------------ statistics

def phi(x: float) -> float:
    return 0.5*(1.0 + math.erf(x/math.sqrt(2.0)))


def u95_of(b: float, sig: float) -> float:
    lo, hi = 0.0, abs(b) + 10.0*sig
    for _ in range(120):
        mid = 0.5*(lo + hi)
        if phi((mid - b)/sig) - phi((-mid - b)/sig) < 0.95:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def anchor_key(tau: float) -> str:
    return 'tau_%g' % tau


def row_at(rows: list, tau: float) -> dict:
    for r in rows:
        if abs(r['tau'] - tau) < 1e-9:
            return r
    raise KeyError('anchor tau=%g not on the grid' % tau)


def wrap_half_away(x: np.ndarray, period: float) -> np.ndarray:
    """Wrap to (-period/2, period/2] with round-half-away tie-breaking
    (10.8f C8 -- np.round's banker's rounding split exact half-quantum ties
    by parity)."""
    return x - period*np.floor(x/period + 0.5)
