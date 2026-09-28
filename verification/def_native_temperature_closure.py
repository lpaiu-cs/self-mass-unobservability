"""Native EOS temperature closure in physical-density GR coordinates.

Frozen historical producers retain their hashes. New calculations call this
shared correction after their original assembly; all affected heat maps are
rebuilt together rather than fixing only a reported temperature.
"""
import numpy as np
from scipy.sparse import coo_matrix, diags


def correct(problem):
    p = problem
    m = p.model
    d = m.heat.d
    rho = d['raw'][::-1, 0]
    cvT = d['thermo'][::-1, 3]
    cpT = d['thermo'][::-1, 5]
    assert np.all(cvT > 0) and np.all(cpT > 0)
    # Import lazily: this function accepts both the time and orbital producers.
    import def_gr_temperature_feedback as prior
    gr = prior.go.task.fem.base.task.h.gr
    geo = gr.G*.1*m.heat.geometry.R**2/gr.C**4
    _, loss, J = prior.go.task.fem.source_points(m.heat, p.point)
    r = p.point['r']
    b = 1-2*p.point['m']/r
    ad = p.point['adiabatic_T_rho']
    p.TE = (-diags(ad/(r*b))@J-diags(1/(rho*geo*cvT))@loss).astype(np.clongdouble)
    active = m.heat.face_ids
    n = len(r)
    faces = n-active
    difference = coo_matrix((np.tile([-1., 1.], len(active)),
        (np.repeat(np.arange(len(active)), 2), np.c_[n-faces, n-1-faces].ravel())),
        shape=(len(active), n)).tocsr()
    p.GE = (difference@diags(p.theta)@(p.TE-p.Tq@m.H)).astype(np.clongdouble)
    if hasattr(p, 'GEraw'):
        dtype = p.GEraw.dtype
        value = (difference@diags(p.theta)@p.TE)[:, active]
        p.GEraw = value.real.astype(dtype) if np.issubdtype(dtype, np.floating) else value.astype(dtype)
        p.energy_scale = 1/np.maximum(abs(p.GEraw.diagonal()), np.longdouble('1e-100'))
    return dict(minimum_cp_over_cv=float(min(cpT/cvT)), maximum_cp_over_cv=float(max(cpT/cvT)),
        maximum_chain_error=float(max(abs(cpT/cvT-(1-ad*d['raw'][::-1, 8])))),
        corrected_maps=['TE', 'GE', 'GEraw when present', 'energy_scale when present'])


def symbolic():
    import sympy as s
    rp, rt, up, ut, pressure_over_rho = s.symbols('rp rt up ut P_rho', nonzero=True)
    heat, density = s.symbols('Q density')
    cv = ut-up*rt/rp
    cp = ut-pressure_over_rho*rt
    ad = (pressure_over_rho-up/rp)/cv
    assert s.simplify(cp-cv*(1-ad*rt)) == 0
    rho_ref = rt*heat/cp
    temperature = ad*(density-rho_ref)+heat/cp
    assert s.simplify(temperature-(ad*density+heat/cv)) == 0
    assert s.simplify(temperature.subs(density, 0)-heat/cv) == 0
    return dict(classification='Proven', passed=True,
        identity='Delta lnT = ad*Delta lnrho + Q/(cv*T). Equivalently ad*(Delta lnrho-rho_ref)+Q/(cp*T), rho_ref=rho_T_at_P*Q/(cp*T). cp/cv=1-ad*rho_T_at_P.',
        preserved_cp_use='The fixed-pressure entropy-density source rho_ref still uses cp, as before. Only the temperature map written in total-density coordinates needs cv.',
        scope='Thermodynamic differential identity at fixed composition; no physical EOS error certificate.')
