"""Identical null variations driven by continuous cumulative photon emission.

Moving source-ray pressure terms are exact total derivatives. Retain their
launch/end terms instead of asking an adaptive ODE to resolve every source jet.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import propagate_smooth_characteristics as base

b=base.b;p=base.p;C=p.C;G=p.G;LD=p.LD
OLD=base.OUT;base.OUT=b.OUT=p.OUT=OUT=Path('native-primitive-photon263-work')
read,write,sha=base.read,base.write,base.sha


class Metric(b.Metric):
    def at(self,t,r):
        t,r=np.broadcast_arrays(t,r);shape=t.shape;t=t.ravel();r=r.ravel();d=self.d
        mass,N,B,a,Phi=d.bg.metric(r/d.model.m.R);mass*=d.model.m.R;Phi/=d.model.m.R
        U=self.wave(t,d.optical(r))[0];tau=d.optical(r)/C;zero=np.zeros(len(t),int)
        integral=self.scalar.convolution(t,tau,zero,self.qcoeff[:,0])
        nr=mass/(r*r*B)+r*Phi*Phi/2;lr=r*Phi*Phi/2-mass/(r*r*B)
        F0=np.zeros(len(t),LD);F1=F0.copy();Fr=F0.copy();kernel=F0.copy()
        for k,mu0 in enumerate(self.mu):
            v=np.clip((r-d.r0)/(self.reach[k]-d.r0),0,1)
            age=d.T*self.inverse[k].sol(v)[0];age=np.where(r>self.reach[k],2*d.T,age)
            bins=np.full(len(t),self.bins[k],int);F=self.aw[k]*self.emission.primitive(np.maximum(t-age,0),bins)
            direction=np.sqrt(np.maximum(1-(1-mu0*mu0)*(N*d.r0/(d.z0['lapse'][0]*r))**2,0))
            F0+=F;F1+=direction*F;Fr+=(1-direction*direction)/direction*(1/r-nr)*F
            kernel+=self.aw[k]*self.emission.convolution(t,age,bins,self.coeff[:,k])
        residual=np.interp(t,self.clock,self.residual);fac=LD(G)/LD(C)**4;K=1/(r*a)
        J=residual-fac*F0;scalar=Phi*U
        lam=scalar+J*K;nu=scalar-integral-J*K-fac*kernel
        values=dict(nu=nu,lam=lam,zeta=nu-lam,N=N,b=B,a=a,nr=nr,lr=lr,scalar=scalar,
            residual=residual,F0=fac*F0,F1=fac*F1,Fr=fac*Fr,K=K)
        return {k:np.asarray(v,float).reshape(shape) for k,v in values.items()}


class Photons(base.Photons):
    def integrate(self,age,te,owner,scale):
        n=len(te);r0,mu0=self.ray(np.zeros(n),owner);launch=self.metric.at(te,r0)
        def endpoint(mu,z):return mu*z['scalar']-mu*mu*z['K']*z['residual']-mu*z['K']*z['F1']
        initial=np.r_[np.zeros(n),-endpoint(mu0,launch)/scale]
        def rhs(s,y):
            r,mu=self.ray(age*s,owner);z=self.metric.at(te+age*s,r)
            v=C*z['a']*mu;md=(1-mu*mu)*C*z['a']*(1/r-z['nr'])
            h=y[n:]+endpoint(mu,z)/scale
            delay=-z['zeta']/scale-(1-mu*mu)*h/(mu*mu)
            qdot=z['residual']*z['K']*(-v*(1+mu*mu)/(r*z['b'])+2*mu*md)
            qdot+=z['K']*(v/(r*z['b'])*(z['F0']-mu*z['F1'])+md*z['F1']+mu*v*z['Fr'])
            qdot+=z['scalar']*(-md+v/r*(1+mu-(1-mu)/z['b']))
            return np.r_[age*delay,age*qdot/scale]
        sol=solve_ivp(rhs,[0.,1.],initial,method='DOP853',rtol=2e-8,atol=2e-11)
        assert sol.success,sol.message
        r,mu=self.ray(age,owner);z=self.metric.at(te+age,r);h=sol.y[n:,-1]+endpoint(mu,z)/scale
        work=(z['nu']-launch['nu'])/scale+h
        sol.y=np.r_[sol.y[:n,-1],h,work][:,None]
        return sol


install_previous=base.install
def install():
    install_previous();Photons.propagate=base.Photons.propagate;base.Photons=Photons;p.Metric=Metric


def symbolic():
    u,w,v,md,wd,K,r,B,F,L,e,R=sp.symbols('u w v md wd K r B F L e R',nonzero=True)
    # kdot pressure term = -e*K*u*(u+w)*(d/dt F), Fdot=(1-u/w)*L.
    assert sp.factor(e*K*L*(-u*w+u**3/w)+e*K*u*(u+w)*(1-u/w)*L)==0
    J=R-e*F;Adot=-K*v/(r*B)*u*(u+w)+K*md*(2*u+w)+K*u*wd
    old=J*K*(-v*(1+u*u)/(r*B)+2*u*md)+e*F*Adot
    new=R*K*(-v*(1+u*u)/(r*B)+2*u*md)+e*K*F*(v*(1-u*w)/(r*B)+w*md+u*wd)
    assert sp.simplify(old-new)==0
    result=base.symbolic();result.update(primitive_transform='q=h-mu*Phi*U+mu^2*K*C_infinity+G/c^4*K*mu*sum weight*mu_source*F(t-travel). Keep its nonzero launch and endpoint terms.',
        crossing='Fdot=(1-mu/mu_source)*L. The pressure coefficient contains this exact factor; no division by a nearly tangent crossing is needed.',
        work='Reconstructed work is an identity, not an independent energy-integral check. Actual direct-jet comparison remains required.')
    # Use the independently proved outgoing-scalar identity too.
    import propagate_outgoing_characteristics as outgoing
    result['outgoing_identity']=outgoing.symbolic()['outgoing_identity']
    # Import changes module output globals; restore this execution's owner.
    base.OUT=b.OUT=p.OUT=OUT
    return result


def prepare():
    base.prepare();plan=read(OUT/'plan.json')
    plan['claim']='Complete the identical returned-metric photon variations using cumulative emission and scalar values after exact integration by parts of all source time jets and moving-ray pressure terms.'
    plan['bindings'][str(Path(__file__))]=sha(__file__)
    plan['bindings']['verification/propagate_outgoing_characteristics.py']=sha('verification/propagate_outgoing_characteristics.py')
    plan['transform']='The pressure numerator cancels the ray-crossing derivative factor. No small-denominator crossing division, source smoothing, new clock or relaxed ODE tolerance.'
    write(OUT/'plan.json',plan);write(OUT/'symbolic.json',symbolic())


def check():
    install();p.initialize();old=b.Metric(8);new=Metric(8);d=new.d
    t=np.linspace(.001,.999,128)*d.T;r=d.inverse(.23*C*t)
    a=old.at(t,r);z=new.at(t,r)
    errors={k:float(np.max(abs(z[k]-a[k]))/max(np.max(abs(a[k])),1e-290)) for k in ['nu','lam','zeta','N','b','a','nr','lr']}
    assert max(errors.values())<1e-10,errors
    # Independent evaluation of the full transformed rhs identity at actual
    # metric points, including residual-dot and outgoing scalar derivative.
    # Time differentiation is a control only; production uses value primitives.
    mu=np.linspace(.02,.97,len(t));v=C*z['a']*mu;md=(1-mu*mu)*C*z['a']*(1/r-z['nr'])
    dt=1e-7*d.T
    def shift(tt,rr,uu):
        q=new.at(tt,rr);return -uu*q['scalar']+uu*uu*q['K']*q['residual']+uu*q['K']*q['F1']
    derivative=(shift(t+dt,r+v*dt,mu+md*dt)-shift(t-dt,r-v*dt,mu-md*dt))/(2*dt)
    direct=-v*a['nup']-mu*mu*a['lt']+derivative
    qdot=z['residual']*z['K']*(-v*(1+mu*mu)/(r*z['b'])+2*mu*md)
    qdot+=z['K']*(v/(r*z['b'])*(z['F0']-mu*z['F1'])+md*z['F1']+mu*v*z['Fr'])
    qdot+=z['scalar']*(-md+v/r*(1+mu-(1-mu)/z['b']))
    rel=float(np.max(abs(qdot-direct))/max(np.max(abs(qdot)),1e-290))
    assert rel<.002,rel
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,metric_values_relative=errors,
        independent_actual_transformed_rhs_relative=rel,actual_direct_trajectory_comparison_pending=True))


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(1800 if action=='pilot' else 600)
    try:
        if action!='prepare':
            for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
        if action in ['prepare','check']:globals()[action]()
        else:
            install();base.install=lambda:None;getattr(base,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
