"""An equivalent canonical momentum with finite pressure-free-surface coefficients."""
from functools import lru_cache
import sympy as s
import def_gr_canonical as old


@lru_cache(maxsize=1)
def symbolic():
    _,prior=old.symbolic()
    r,m,p,e,phi,v,G,A,C,ep,Gp=s.symbols('r m p e phi v Gamma A4 C e_prime Gamma_prime',nonzero=True)
    symbols=[r,m,p,e,phi,v,G,A,C,ep,Gp]
    locals={str(x):x for x in symbols};locals['pi']=s.pi
    values=[s.sympify(x,locals=locals) for x in prior['expressions']]
    matrix=lambda i,n,m:s.Matrix(n,m,values[i:i+n*m])
    Ao=matrix(0,2,2);Bi=matrix(4,2,2);Co=matrix(8,2,2);W=matrix(12,2,2)
    gs=matrix(16,2,3);hs=matrix(22,2,3);R=matrix(28,2,2);S=matrix(32,2,2)
    U=s.Matrix([[0,-S[0,1]/2],[-S[0,1]/2,0]])
    b=1-2*m/r;alpha=-4*phi;w=e+p
    mp=4*s.pi*r*r*A*e+r*r*b*v*v/2
    g=m/(r*r*b)+4*s.pi*r*A*p/b+r*v*v/2+alpha*v
    F=4*s.pi*A/b*(alpha*(e-3*p)+r*v*(e-p))-2*(r-m)*v/(r*r*b)
    Cp=C*(2*m/(r*r*b)+4*s.pi*r*A*(p-e)/b)
    variables=[r,m,p,e,phi,v,G,A,C];rates=[1,mp,-w*g,ep,v,F,Gp,4*alpha*A*v,Cp]
    derivative=lambda x:sum(s.diff(x,a)*da for a,da in zip(variables,rates))
    clean=lambda M:M.applyfunc(s.cancel)
    An=clean(Ao-Bi.inv()*U)
    Cn=clean(Co+U.applyfunc(derivative)+U*Ao+Ao.T*U-U*Bi.inv()*U)
    hn=clean(hs+U*gs);T=clean(S/2+U)
    assert Cn==Cn.T and U==U.T and T[0,0]==T[0,1]==T[1,1]==0
    # No reciprocal pressure remains in the differential energy coefficients.
    assert all(s.simplify(s.denom(x).subs(p,0))!=0 for x in [*An,*Cn,*hn])
    x1,x2,dx1,dx2=s.symbols('x1 x2 dx1 dx2');q=s.Matrix([x1,x2]);dq=s.Matrix([dx1,dx2])
    vo=(dq-Ao*q).T*Bi*(dq-Ao*q)+q.T*Co*q
    vn=(dq-An*q).T*Bi*(dq-An*q)+q.T*Cn*q
    boundary=2*dq.T*U*q+q.T*U.applyfunc(derivative)*q
    assert clean(vn-vo-boundary)==s.zeros(1)
    expr=[*An,*Bi,*Cn,*W,*gs,*hn,*R,*(2*T)]
    fn=s.lambdify(symbols,expr,'numpy',cse=True)
    return fn,dict(classification='Proven',passed=True,
        change='P_new=P_old+U*q, U symmetric, U12=U21=-S12/2. Thus P_new,fluid=R11*physical_pressure, with no subtractive scalar-pressure cancellation.',
        energy_identity='new_energy_density-old_energy_density=(q.T*U*q)_r/2. The boundary term vanishes at regular center, zero-density surface and in vacuum.',
        pressure_regular=True,scope='Same continuous GR and source equations under the declared background identities; discrete quadrature/background interpolation may break exact gauge equivalence and are checked on actual paths.',
        storage='Last four expressions are2*T, not an antisymmetric S; reconstruction always subtracts this matrix times q/2.',expressions=[str(x) for x in expr])
