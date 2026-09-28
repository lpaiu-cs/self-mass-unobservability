"""Moving-surface adapter for the verified dynamic vacuum boundary.

This is a scalar matching condition, not the missing material pressure or
atmosphere evolution law. No new EOS or radial integration is performed.
"""
from pathlib import Path
import json
import sympy as sp
import def_radiative_exterior as exterior


def residual(delta_phi_lagrangian, delta_Rphi_prime_lagrangian, xi_over_R,
             q, R2phi_second, impedance, h, incoming):
    return (delta_Rphi_prime_lagrangian-impedance*delta_phi_lagrangian
            -incoming*h.conjugate()*(impedance.conjugate()-impedance)
            -xi_over_R*(R2phi_second-impedance*q))


def induced_outgoing(delta_phi_lagrangian_contrast, xi_over_R_contrast, q, h):
    """Same incoming field/background: avoid subtracting two large wave amplitudes."""
    return (delta_phi_lagrangian_contrast-xi_over_R_contrast*q)/h


def check():
    phi,prime,xi,q,second,Z,rhs=sp.symbols('phi prime xi q second Z rhs')
    lagphi=phi+xi*q;lagprime=prime+xi*second
    assert sp.expand(lagprime-Z*lagphi-rhs-xi*(second-Z*q)-(prime-Z*phi-rhs))==0
    path=exterior.OUT/'result.json';result=json.loads(path.read_text());assert result['passed']
    errors=[]
    for row in result['rows']:
        if row['case']!='stellar' or not row['omega_R_over_c']:continue
        solution=row['solves'][-1];h=complex(*solution['h']);Z=complex(*solution['impedance'])
        mu,q=row['mu'],row['scalar_flux'];second=-(2+2*mu/(1-2*mu))*q
        incoming=.13-.21j;outgoing=.08+.19j;xi=.003-.004j
        phi=outgoing*h+incoming*h.conjugate()
        prime=outgoing*h*Z+incoming*h.conjugate()*Z.conjugate()
        err=abs(residual(phi+xi*q,prime+xi*second,xi,q,second,Z,h,incoming))
        contrast=.01+.02j
        err=max(err,abs(induced_outgoing(contrast*h+xi*q,xi,q,h)-contrast))
        assert err<1e-14;errors.append(err)
    return dict(classification='Proven',passed=True,symbolic_identity=True,
        scope='Linear moving-surface change of variables and same-incoming outgoing contrast; not material free-surface dynamics.',
        numerical_control_errors=errors,source_sha256=exterior.e.digest(Path(__file__)),
        dynamic_exterior_result_sha256=exterior.e.digest(path))


if __name__=='__main__':
    target=exterior.OUT/'moving-surface.json';assert not target.exists()
    value=check();exterior.e.write(target,value);print(json.dumps(value))
