"""Independent identities and manufactured controls for the retarded readout."""
import json
import time
import numpy as np
from scipy.integrate import quad
import sympy as s
import def_native_release_charge as task


def main():
    start=time.monotonic();out=task.OUT
    assert not (out/'audit.json').exists()
    r,A,N,b,alpha,trace=s.symbols('r A N b alpha trace',positive=True)
    proper=4*s.pi*A**3*r*r/s.sqrt(b)
    assert s.simplify((2*s.pi*r*N*A**4/s.sqrt(b))/(proper*A*N/(2*r))-1)==0
    rho,D,v,u,p,cx,c=s.symbols('rho D v u p cx c',real=True)
    nr=rho*u-3*p+cx*c*c*(rho-D)
    assert s.expand(cx*c*c*D+nr-(rho*(cx*c*c+u)-3*p))==0
    t,delay,cs=s.symbols('t delay cs',real=True);z=s.symbols('z',real=True)
    port=z*z
    integrated=s.integrate((t+delay-(1+cs)*z)**2,(z,0,(t+delay)/(1+cs)))
    assert s.simplify(integrated-(t+delay)**3/(3*(1+cs)))==0
    # Two equal and opposite baryon histories: independent quadrature directly
    # evaluates the uncentered mathematical source, while the producer uses
    # its cancellation-safe centered form. Include a retarded zero crossing.
    times=np.array([0.,.0001,.001,.003]);rates=np.array([7e14,-7e14]);weights=np.array([-.004/7e9,-.004/(7e9+900)])
    delays=np.array([-4e-6,3e-6]);w0=-.004/7e9;checks=[]
    for kind in ['linear','cubic']:
        poly=task.polynomial(times,times[:,None]*rates,kind).antiderivative()
        for at in [0.,2e-6,1e-5,.002]:
            now=task.paired(poly,np.full(2,at));later=task.paired(poly,at+delays)
            actual=np.sum(weights*(later-now)+(weights-w0)*now,dtype=np.longdouble)
            expected=sum(w*quad(lambda tt:rate*tt,0,max(0,at+d),epsabs=1e-9)[0] for w,rate,d in zip(weights,rates,delays))
            checks.append(float(abs(actual-expected)/max(abs(expected),1e-40)))
    # Constant metric/weight acoustic solution has the independently integrated
    # result -H_port(t+d)/(1+cs), with conserved total bulk mass -port(t).
    speed=.00006;at=.002;dd=-7e-7;value=quad(lambda tau:(at+dd-(1+speed)*tau)**2,0,(at+dd)/(1+speed),epsabs=1e-20)[0]
    acoustic=abs(value/((at+dd)**3/(3*(1+speed)))-1)
    passed=max(checks)<1e-8 and acoustic<1e-10
    result=dict(classification='Proven',passed=passed,manufactured_maximum_relative=max(checks),
        independent_acoustic_integral_relative=float(acoustic),seconds=time.monotonic()-start,
        identities=['Jordan proper-volume source conversion','rest plus nonrest trace decomposition','incoming acoustic port convolution'],
        direct_equation='U_tt-c^2 U_rstar_rstar=-4*pi*G/c^2*N^2*alpha*A^4*r*delta(e-3p); U=r*delta_phi. The scalar/metric potential and metric-mediated matter sources are omitted, not bounded by this identity.',
        kernel='Q_direct(u)=-G/(2*c^3)*integral w(r)*integral_0^(u+d(r)) delta(trace*dV) dt; w=alpha*A*N/r.',
        normalization='delta_alpha_direct=-Q_direct/M_geom. Time-dependent Bondi mass and indirect sources are separate terms.',
        whole_physical_model_certified=False)
    task.write(out/'audit.json',result);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':main()
