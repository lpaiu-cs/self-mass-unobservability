"""Canonical energy form of the original continuous fluid/scalar equations.

q=(xi/r, delta_phi), p=(Delta_p/p, delta_Phi). The physical pressure,
not the entropy-lifted auxiliary pressure, is used when transporting heat.
"""
from functools import lru_cache
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as s
import def_gr_matrix_propagation as prior

OUT=prior.OUT.parent/'def-gr-canonical';write=prior.write


@lru_cache(maxsize=1)
def symbolic():
    r,m,p,e,phi,v,G,A,C,ep,Gp=s.symbols('r m p e phi v Gamma A4 C e_prime Gamma_prime',nonzero=True)
    z,eta,f,V,rr,loss,J=s.symbols('z eta f V rho_ref loss J')
    beta=s.Integer(-4);alpha=beta*phi;b=1-2*m/r;w=e+p
    mp=4*s.pi*r*r*A*e+r*r*b*v*v/2
    g=m/(r*r*b)+4*s.pi*r*A*p/b+r*v*v/2+alpha*v
    F=4*s.pi*A/b*(alpha*(e-3*p)+r*v*(e-p))-2*(r-m)*v/(r*r*b)
    Cp=C*(2*m/(r*r*b)+4*s.pi*r*A*(p-e)/b)
    variables=[r,m,p,e,phi,v,G,A,C]
    rates=[1,mp,-w*g,ep,v,F,Gp,4*alpha*A*v,Cp]
    derivative=lambda expr:sum(s.diff(expr,a)*da for a,da in zip(variables,rates))
    clean=lambda M:M.applyfunc(s.cancel)
    xi=r*z;Dm=r*r*b*v*f-(4*s.pi*r*r*A*p+r*r*b*v*v/2)*xi+J
    dl=(Dm/r-m*xi/r**2)/b
    xp=-eta/G-rr-2*z-dl-3*alpha*f
    inc=[xi,Dm,p*eta,w*(eta/G+rr)-loss,f,V,4*alpha*A*f]
    args=[r,m,p,e,phi,v,A]
    Dg=sum(s.diff(g,x)*dx for x,dx in zip(args,inc))
    DF=sum(s.diff(F,x)*dx for x,dx in zip(args,inc))
    eta_prime=-g*(w*(eta/G+rr)-loss+p*eta+w*xp)/p-w*Dg/p+w*g*eta/p
    rhs=s.Matrix([(xp-z)/r,eta_prime,V+xp*v,DF+xp*F])
    H=rhs.jacobian([z,eta,f,V]);source=rhs.jacobian([rr,loss,J])
    T=s.Matrix([[1,0,0,0],[0,0,1,0],[r*v,1,0,0],[r*F,0,0,1]])
    Tp=T.applyfunc(derivative)
    H=clean(T.inv()*(H*T-Tp));source=clean(T.inv()*source)
    Aq=H[:2,:2];Bq=H[:2,2:];Cq=H[2:,:2];Dq=H[2:,2:]
    assert clean(Bq-s.diag(-1/(G*r),1))==s.zeros(2)
    R=s.diag(-4*s.pi*r**3*A*p*C/b,r*r*C)
    sc=s.cancel(-R[0,0]*Dq[0,1]);S=s.Matrix([[0,sc],[-sc,0]])
    Acan=clean(Aq-Bq*R.inv()*S/2);Bcan=clean(Bq*R.inv())
    X=clean(R.applyfunc(derivative)+R*Dq+S*Bq/2)
    Dcan=clean(X*R.inv())
    Ccan=clean(R*Cq+S.applyfunc(derivative)/2+S*Aq/2-X*R.inv()*S/2)
    residual_d=clean(Dcan+Acan.T);residual_c=clean(Ccan-Ccan.T)
    assert residual_d==s.zeros(2),str(residual_d)
    assert residual_c==s.zeros(2),str(residual_c)
    W=clean(R*s.diag(-w*r/(p*C*C),1/(C*C)))
    assert clean(W-s.diag(4*s.pi*r**4*A*w/(b*C),r*r/C))==s.zeros(2)
    gsource=source[:2,:];hsource=clean(R*source[2:,:]+S*gsource/2)
    assert clean(gsource-s.Matrix([[-1/r,0,-1/(r*r*b)],[0,0,0]]))==s.zeros(2,3)
    expr=[*Acan,*Bcan.inv(),*Ccan,*W,*gsource,*hsource,*R,*S]
    fn=s.lambdify([r,m,p,e,phi,v,G,A,C,ep,Gp],expr,'numpy',cse=True)
    text=dict(classification='Proven',passed=True,
        variables='q=(xi/r, Eulerian delta_phi); p=(physical Delta_p/p, Eulerian delta_Phi)',
        pressure='eta_ad=physical_Delta_p/p+Gamma*rho_ref. Move this source out of the derivative equation before the weak form.',
        canonical_momentum='P=R*p+S*q/2; R=diag(-4*pi*r^3*A4*p*C/b,r^2*C), C=N*sqrt(b). S is antisymmetric.',
        equations='q_prime=Acan*q+Bcan*P+g; P_prime=Ccan*q-Acan.T*P+W*q_tt+h; Bcan and Ccan symmetric.',
        kinetic='W=diag(4*pi*r^4*A4*(e+p)/(b*C),r^2/C), positive in matter and positive in the scalar vacuum channel.',
        energy='V=integral[(q_prime-Acan*q).T*Bcan^-1*(q_prime-Acan*q)+q.T*Ccan*q]/2; kinetic=integral q_t.T*W*q_t/2.',
        load='Weak thermal force is integral[(test_prime-Acan*test).T*Bcan^-1*g-test.T*h], with physical natural pressure boundary.',
        heat_momentum='With H=(H_zeta,0), u=q+H gives M*u_tt+K*u=f+K*H. Recover pressure from physical q and g, then add Gamma*rho_ref only for the old eta_ad readout.',
        scope='Continuous equations at beta=-4 with the declared mechanical GR background identities. No claim that the old nodal midpoint pencil is self-adjoint, or that a new spatial/time discretization has passed.',
        source_fields=['rho_ref','loss','J'],expressions=[str(x) for x in expr])
    return fn,text


def main():
    assert not OUT.exists();OUT.mkdir();signal.alarm(120);start=time.monotonic()
    write(OUT/'derivation-plan.json',dict(classification='Conjectural',checkpoint='29fe9acd',
        claim='Derive and check the complete canonical energy form, physical pressure and thermal weak source of the existing continuous GR equations before modifying spatial inertia.',
        budget=dict(symbolic_seconds=120,new_evolution_paths=0,new_EOS_calls=0),
        source=prior.prior.prior.digest(Path(__file__))))
    _,result=symbolic();result['seconds']=time.monotonic()-start
    write(OUT/'symbolic.json',result);signal.alarm(0);print('CANONICAL',result['seconds'],result['passed'],flush=True)


if __name__=='__main__':main()
