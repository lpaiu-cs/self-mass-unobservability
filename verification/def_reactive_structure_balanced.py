"""Conservative mass defect avoids cancellation of reversible work terms."""
from pathlib import Path
from types import FunctionType
import inspect
import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu
import def_reactive_structure as old

OUT=old.OUT/'balanced'


def mass_identity():
    r,m,p,e,phi,v,beta=sp.symbols('r m p e phi v beta',real=True)
    xi,f,V,J,eta,ga,rr,loss=sp.symbols('xi f V J eta ga rr loss')
    b=1-2*m/r;A4=sp.exp(2*beta*phi**2);alpha=beta*phi;w=e+p
    M=4*sp.pi*r*r*A4*e+r*r*b*v*v/2
    g=m/(r*r*b)+4*sp.pi*r*A4*p/b+r*v*v/2+alpha*v
    F=4*sp.pi*A4/b*(alpha*(e-3*p)+r*v*(e-p))-2*(r-m)*v/(r*r*b)
    H=r*r*b*v*f-(4*sp.pi*r*r*A4*p+r*r*b*v*v/2)*xi
    Dm=H+J;dl=(Dm/r-m*xi/r**2)/b
    xp=-eta/ga-rr-2*xi/r-dl-3*alpha*f
    inc=[xi,Dm,p*eta,w*(eta/ga+rr)-loss,f,V]
    DM=sum(sp.diff(M,t)*u for t,u in zip([r,m,p,e,phi,v],inc))
    Hp=sum(sp.diff(H,t)*u for t,u in zip([r,m,p,phi,v,xi,f],[1,M,-w*g,v,F,xp,V+xp*v]))
    target=-(4*sp.pi*r*A4*w/b+r*v*v)*J-4*sp.pi*r*r*A4*loss
    assert sp.simplify(DM+xp*M-Hp-target)==0
    return dict(classification='Proven',passed=True,
        definition='J=Delta m-[r^2*b*Phi*f-(4*pi*r^2*A^4*p+r^2*b*Phi^2/2)*xi]',
        equation='J_prime=-(4*pi*r*A^4*(e+p)/b+r*Phi^2)*J-4*pi*r^2*A^4*loss',
        integrating_factor='N*a; (N*a*J)_prime=-4*pi*r^2*A^4*N*a*loss',
        ADM_first_law='At fixed phi_infinity: dot M_geom=-(G/c^4)*1e-7*sum(dm*A^2*N^2*epsilon_nu) for dm in grams and epsilon_nu in erg/g/s.',
        scope='Linear quasi-static first law. Nuclear rest energy becomes internal energy; only specified escaped neutrino energy is removed at the initial zero-flux state.')


def coefficients(bg,rho_ref,loss,fn):
    r,m,p,e,phi,v,N,ga=[bg[k] for k in ['r','m','p','e','phi','v','N','gamma']]
    b=1-2*m/r;A4=np.exp(-8*phi*phi);alpha=-4*phi;w=e+p
    val=fn(r,m,p,e,phi,v,-4.);M,g,F=val[:3];dg,dF=val[9:15],val[15:21]
    basis=np.zeros((5,6,len(r)));basis[:,:5,:]=np.eye(5)[:,:,None]
    z,eta,f,V,j=basis;xi=r*z
    rr=np.zeros_like(z);rr[-1]=rho_ref
    le=np.zeros_like(z);le[-1]=loss
    H=r*r*b*v*f-(4*np.pi*r*r*A4*p+r*r*b*v*v/2)*xi
    Dm=H+r*j;dl=(Dm/r-m*xi/r**2)/b
    xp=-eta/ga-rr-2*z-dl-3*alpha*f
    inc=[xi,Dm,p*eta,w*(eta/ga+rr)-le,f,V]
    DG=lambda der:sum(q*u for q,u in zip(der,inc))
    bracket=g*(2*z+dl+3*alpha*f+eta)-DG(dg)
    bracket[-1]+=np.divide(g*loss,w,out=np.zeros_like(w),where=w!=0)
    jp=-(4*np.pi*r*A4*w/b+r*v*v+1/r)*j-4*np.pi*r*A4*le
    with np.errstate(divide='ignore',invalid='ignore'):
        result=np.array([(xp-z)/r,w/p*bracket-g*eta,V+xp*v,DG(dF)+xp*F,jp])
    return np.moveaxis(result,-1,0),bracket,M,F


ns=dict(vars(old),OUT=OUT,__file__=__file__,old=old,coefficients=coefficients,splu=splu,mass_identity=mass_identity)
source=inspect.getsource(old.solve)
before="    y=spsolve(matrix.multiply((1/scale)[:,None]).tocsc(),rhs/scale).reshape(n+1,5)"
after="""    lu=splu(matrix.multiply((1/scale)[:,None]).tocsc())
    answer=lu.solve(rhs/scale)
    extended=matrix.astype(np.longdouble)
    for _ in range(2):
        defect=rhs.astype(np.longdouble)-extended@answer.astype(np.longdouble)
        answer+=lu.solve(np.asarray(defect/scale,float))
    y=answer.reshape(n+1,5)"""
assert source.count(before)==1;source=source.replace(before,after)
# Convert the conservative defect to Delta m/r for exterior rows/readout only.
anchor="    ev,jac=ext\n"
replacement="""    ev,jac=ext
    b=1-2*saved['m'][-1]/R
    convert=np.array([-R*R*b*v*v/2,0,R*b*v,0,1.])
"""
source=source.replace(anchor,replacement)
source=source.replace("    bcs=[", "    phiinfty+=jac[2,0]*(convert-np.array([0,0,0,0,1.]))\n    bcs=[")
source=source.replace("    z,eta,f,V,dd=y[-1];", "    z,eta,f,V,_=y[-1];dd=float(convert@y[-1]);")
exec(compile(source,__file__,'exec'),ns)
main=inspect.getsource(old.main)
main=main.replace("paths=[Path(__file__),", "paths=[Path(__file__),Path(old.__file__),old.OUT/'result.json',")
main=main.replace("symbolic=proof,", "symbolic=proof,mass_identity=mass_identity(),")
main=main.replace("expected=-power*h.gr.G*.1/h.gr.C**4", "expected=-power*h.gr.G*1e-7/h.gr.C**4")
main=main.replace("claim='Connect", "preserved_failure='Original independent-mass solve lost precision through reversible-work cancellation. Original first-law reporting also used the energy-density conversion 0.1 for luminosity, whose correct conversion is 1e-7. Preserve its failed results.',\n        claim='Connect")
exec(compile(main,__file__,'exec'),ns)


if __name__=='__main__':ns['main']()
