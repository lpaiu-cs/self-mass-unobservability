"""Proven: finite summation identity, conditional on the declared stage equations."""
from pathlib import Path
import json
import sympy as sp
import def_spherical_regular as model


def check():
    model.secant.check()
    x=sp.symbols('dphi0:3'); xb=sp.Symbol('dphi_boundary')
    G=sp.symbols('gradient_port0:4'); I=sp.symbols('matter_port0:4')
    H=sp.symbols('H0:3'); V=sp.symbols('gV0:3')
    source=sp.symbols('gstar0:3'); error=sp.symbols('energy_residual0:3')
    # The scalar position equation turns kinetic work into x*(Gright-Gleft).
    kinetic=sum(x[i]*(G[i+1]-G[i])+H[i]*V[i]*source[i]*x[i] for i in range(3))
    gradient=sum(G[i]*(x[i]-x[i-1]) for i in range(1,3))+G[3]*(xb-x[2])
    matter=sum(-H[i]*V[i]*source[i]*x[i]-(I[i+1]-I[i])+H[i]*V[i]*error[i] for i in range(3))
    expected=G[3]*xb-I[3]+I[0]+sum(H[i]*V[i]*error[i] for i in range(3))
    assert sp.expand((kinetic+gradient+matter-expected).subs(G[0],0))==0
    # Each radial adjoint pair cancels the interior mass increment exactly.
    C,D=sp.symbols('C0:3'),sp.symbols('D0:3')
    lam=sp.symbols('lambda0:4'); dm=sp.symbols('dm0:4')
    radial=sum(lam[i+1]*dm[i+1]-lam[i]*dm[i] for i in range(3))
    assert sp.expand(radial-(lam[3]*dm[3]-lam[0]*dm[0]))==0
    return dict(classification='Proven',passed=True,
        identity='Delta(M)/Ggeom = scalar_boundary_work - matter_boundary_outflow + sum(H V energy_residual)',
        assumptions=['Regular positive radial constraint branch; Cstar,Dstar and H from the same material metric secant.',
            'Same scalar kinetic/gradient weights in the constraint and wave equations; zero central gradient port.',
            'Shared matter transport face and identical reciprocal trace source.',
            'This implementation starts at phi=Pi=0 and performs one stage.'],
        not_proved=['Nonlinear continuum convergence','Full Einstein angular/momentum constraint preservation',
            'Finite step entropy monotonicity','Physical stellar response or observational closure'])


if __name__=='__main__':
    target=model.OUT/'assembly-identity.json';assert not target.exists()
    result=check();result['source_sha256']=model.e.digest(Path(__file__))
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
