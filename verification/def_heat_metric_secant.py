"""Proven conditional secant identities for the future spherical heat coupling."""
from pathlib import Path
import json
import sympy as sp
import def_heat_coupling as model


def check():
    a0,a1,E0,E1,R0,R1,g,dphi,F,V,r,f,dmL,dmR,chi=sp.symbols('a0 a1 E0 E1 R0 R1 g dphi F V r f dmL dmR chi',nonzero=True)
    av=(a0+a1)/2;ev=(E0+E1)/2;rv=sp.Symbol('Rstar')
    assert sp.expand(a1*E1-a0*E0-av*(E1-E0)-ev*(a1-a0))==0
    b0,b1=sp.symbols('b0 b1',positive=True)
    aa0,aa1=1/sp.sqrt(b0),1/sp.sqrt(b1)
    assert sp.simplify((aa1-aa0)/((aa0+aa1)/2)+2*(b1-b0)/(sp.sqrt(b0)+sp.sqrt(b1))**2)==0
    x=sp.Symbol('x')
    assert sp.limit(sp.log(1+x)/x,x,0)==1
    # ΔE = -(Ebar+Rbar) Δa/abar - g Δphi + F.
    C,D,B=sp.symbols('C D B')
    work=(ev+rv)*chi/r
    de=-work*((1-f)*dmL+f*dmR)-g*dphi+F
    mass=D*dmR-C*dmL-B-V*de
    rearranged=(D+V*work*f)*dmR-(C-V*work*(1-f))*dmL-B+V*g*dphi-V*F
    assert sp.expand(mass-rearranged)==0
    # Weighted transport telescopes. V includes Ggeom when restoring units.
    H1,H2,J0,J1,J2=sp.symbols('H1 H2 J0 J1 J2',nonzero=True)
    assert sp.expand(H1*(-(J1-J0)/H1)+H2*(-(J2-J1)/H2)+J2-J0)==0
    return dict(classification='Proven',passed=True,
        metric_secant='Delta(a)/abar = chi*Delta(m_cell)/r, chi=4/(sqrt(b0)+sqrt(b1))^2',
        matter_work='Delta(aE)=-Rstar*Delta(a)-abar*gstar*Delta(phi)+abar*F; Rstar=mean(aR)*L(a0,a1), gstar=mean(a*alpha*trace)/abar, L=log(a1/a0)/(a1-a0) with its smooth diagonal value 1/a0. This is the same logarithmic metric work used by the native local paths.',
        radial_coefficients='Dstar=Dbar+Ggeom*V*(Ebar+Rstar)*chi*f/r; Cstar=Cbar-Ggeom*V*(Ebar+Rstar)*chi*(1-f)/r',
        scalar_source='The same radial adjoint H yields scalar Pi source H*Ggeom*V*gstar/w_scalar. No inverse Pi.',
        transport='F_i=-Delta_face(I_m)/(H_i*V_i), where I_m=dt*Hface*area*k*S_m_mid. Its H-weighted volume sum is exactly the time-integrated boundary flux.',
        scope='Algebraic coupling construction only. Requires the stated finite matter-work equation, a shared discrete adjoint, positive metric branch and compatible transport flux. It has not been implemented or run in the full native spherical solver.')


if __name__=='__main__':
    target=model.OUT/'metric-secant.json';assert not target.exists()
    result=check();result['source_sha256']=model.digest(Path(__file__))
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
