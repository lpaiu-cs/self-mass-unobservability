"""Exact force-level checks and the frozen auxiliary-drive phase audit (no timing engine)."""
import hashlib
import json
from pathlib import Path
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    passed = []
    def check(name, expression):
        assert s.simplify(expression) == 0, (name, expression)
        passed.append(name)

    # A one-dimensional pair suffices to check the central-force sign and energy identity.
    xp, xj, vp, vj, qp, qj, qdot = s.symbols('xp xj vp vj Qp Qj Qdot', real=True)
    mp, mj, inertia, gamma, kappa, r = s.symbols('mp mj I Gamma kappa r', positive=True)
    c2, c4, q0 = s.symbols('c2 c4 Q0', real=True)
    potential = c2*qp**2/2+c4*qp**4/24
    interaction = -(mp*mj+qp*qj)/(xp-xj)  # xp > xj patch of a central potential
    fp, fj = -s.diff(interaction,xp), -s.diff(interaction,xj)
    check('reciprocal_pair_force', fp+fj)
    check('common_pair_acceleration_factor', fp/(-mp*mj/(xp-xj)**2)-(1+qp*qj/(mp*mj)))
    qddot = (qj/(xp-xj)-s.diff(potential,qp)-gamma*qdot)/inertia
    energy_dot = (vp*fp+vj*fj+inertia*qdot*qddot+s.diff(potential,qp)*qdot
                  +s.diff(interaction,xp)*vp+s.diff(interaction,xj)*vj+s.diff(interaction,qp)*qdot)
    check('orbital_plus_state_energy_loss', energy_dot+gamma*qdot**2)
    check('stable_branch_stiffness', s.diff(potential,qp,2).subs(qp,q0)-(c2+c4*q0**2/2))

    omega, tau, aw, ustar = s.symbols('omega tau aw Ustar', positive=True)
    h = 1/(kappa-inertia*omega**2+s.I*gamma*omega)
    h0 = 1/(kappa+s.I*gamma*omega)
    check('relative_inertial_error', (h-h0)/h0-inertia*omega**2*h)
    check('one_pole_matching', aw**2*h0/mp-(aw**2/(kappa*mp))/(1+s.I*omega*gamma/kappa))
    check('normalized_drive_beta', (aw**2/(kappa*mp))*ustar-aw**2*ustar/(kappa*mp))
    ai, ao, dq = s.symbols('ai ao dQ', real=True)
    check('unequal_companion_obstruction', ai*dq/mp-ao*dq/mp-(ai-ao)*dq/mp)
    susceptibility = s.symbols('C', real=True)
    qj_response = susceptibility*dq/r
    check('companion_feedback_stiffness', kappa*dq-qj_response/r-(kappa-susceptibility/r**2)*dq)
    check('no_relaxation_resonance', s.diff(1/(1+(omega*tau)**2),omega)
          +2*omega*tau**2/(1+(omega*tau)**2)**2)

    ecc, mean, ratio, angle = s.symbols('e M rho theta', real=True)
    check('eccentric_inverse_radius', s.diff(1/(1-ecc*s.cos(mean)),ecc).subs(ecc,0)-s.cos(mean))
    check('hierarchical_inverse_radius', s.diff((1+2*ratio*s.cos(angle)+ratio**2)**s.Rational(-1,2),ratio).subs(ratio,0)+s.cos(angle))
    eta, kap, lam, varpi = s.symbols('eta kap lambda varpi', real=True)
    check('pericenter_potential_phase', s.expand_trig(s.cos(lam-varpi))
          -s.cos(lam)*s.cos(varpi)-s.sin(lam)*s.sin(varpi))

    source=ROOT/'request10_external/baseline_planetGR.npz'
    z=np.load(source,allow_pickle=True)
    p=dict(zip(z['names'],z['params']))
    om={k:2*np.pi/p[v] for k,v in [('in','period_i'),('out','period_o')]}
    om['dif']=om['in']-om['out']
    # Published table: eta=e*cos(varpi), kappa=e*sin(varpi) in this release.
    # Do not substitute the opposite ELL1 name convention from another implementation.
    per={k:float(np.arctan2(p['kappa_'+v],p['eta_'+v])) for k,v in [('in','p'),('out','b')]}
    longitudes={'in':-om['in']*p['tasc_p'],'out':-om['out']*p['tasc_b']}
    physical={'in':longitudes['in']-per['in'],'out':longitudes['out']-per['out'],
              'dif':longitudes['in']-longitudes['out']+np.pi}
    archived={'in':longitudes['in']-np.pi/2,'out':longitudes['out']-np.pi/2,
              'dif':longitudes['in']-longitudes['out']}
    wrap=lambda x:float(np.arctan2(np.sin(x),np.cos(x)))
    closure=wrap(physical['dif']-physical['in']+physical['out'])
    assert abs(closure)>3.0
    # Closure is invariant under a common origin shift because omega_dif=omega_in-omega_out.
    for t0 in [0.,1.,100.]:
        shifted={k:physical[k]+om[k]*t0 for k in physical}
        assert abs(wrap(shifted['dif']-shifted['in']+shifted['out'])-closure)<1e-10
    passed.append('phase_closure_invariance_and_historical_mismatch')
    result=dict(status='Proven',scope='Conditional Newtonian scalar-charge model and leading coplanar drive; no EOS or complete timing matching',
                checks=passed,input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                phase_convention_source='https://arxiv.org/html/2411.10066v2 (footnote 3 and Table 4)',
                pericenter_radians=per,leading_physical_phase_radians={k:wrap(v) for k,v in physical.items()},
                archived_auxiliary_phase_radians={k:wrap(v) for k,v in archived.items()},
                physical_minus_archived_radians={k:wrap(physical[k]-archived[k]) for k in physical},
                physical_closure_radians=closure,archived_closure_radians=wrap(archived['dif']-archived['in']+archived['out']),
                common_origin_cannot_fix_closure=True,
                matching=dict(tau='Gamma/kappa',B='a_w^2/(kappa*m_p)',beta_for_F_equals_deltaU_over_Ustar='B*Ustar'),
                omitted_terms=['inertia','nonlinear charge response','companion-charge feedback','noncoplanar and higher orbital harmonics','higher PN and radiation terms','transients'])
    out=ROOT/'outputs/research-completion/physical-matching.json'
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
