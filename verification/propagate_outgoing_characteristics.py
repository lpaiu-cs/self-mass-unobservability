"""Also remove the declared outgoing scalar's time jet, exactly."""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import propagate_smooth_characteristics as base

p=base.p; b=base.b; C=p.C
OLD=base.OUT; base.OUT=b.OUT=p.OUT=OUT=Path('native-outgoing-photon263-work')
read,write,sha=base.read,base.write,base.sha


class Photons(base.Photons):
    def pieces(self,t,r,z):
        d=self.d;U,Ut,_=self.metric.wave(t,d.optical(r))
        Phi=d.bg.metric(r/d.model.m.R)[4]/d.model.m.R
        scalar=Phi*U; lamp=z['lam']-scalar
        nr=z['nup']-Phi*(-Ut/(C*z['a'])+U*(1/z['b']-1)/r)
        lr=z['lamp']-Phi*(-Ut/(C*z['a'])-U*(1+1/z['b'])/r)
        return scalar,lamp,nr,lr

    def integrate(self,age,te,owner,scale):
        n=len(te);r0,mu0=self.ray(np.zeros(n),owner);launch=self.metric.at(te,r0)
        scalar,lamp,_,_=self.pieces(te,r0,launch)
        initial=np.r_[np.zeros(n),(-mu0*scalar+mu0*mu0*lamp)/scale]
        def rhs(s,y):
            r,mu=self.ray(age*s,owner);z=self.metric.at(te+age*s,r)
            scalar,lamp,nr,lr=self.pieces(te+age*s,r,z)
            md=(1-mu*mu)*C*z['a']*(1/r-z['nr']);v=C*z['a']*mu
            h=y[n:]+(mu*scalar-mu*mu*lamp)/scale
            delay=-z['zeta']/scale-(1-mu*mu)*h/(mu*mu)
            qdot=-v*nr+2*mu*md*lamp+mu*mu*v*lr
            qdot+=scalar*(-md+v/r*(1+mu-(1-mu)/z['b']))
            return np.r_[age*delay,age*qdot/scale]
        sol=solve_ivp(rhs,[0.,1.],initial,method='DOP853',rtol=2e-8,atol=2e-11)
        assert sol.success,sol.message
        r,mu=self.ray(age,owner);z=self.metric.at(te+age,r);scalar,lamp,_,_=self.pieces(te+age,r,z)
        h=sol.y[n:,-1]+(mu*scalar-mu*mu*lamp)/scale
        work=(z['nu']-launch['nu'])/scale+h
        sol.y=np.r_[sol.y[:n,-1],h,work][:,None]
        return sol


install_previous=base.install
def install():
    install_previous(); Photons.propagate=base.Photons.propagate
    base.Photons=Photons


def symbolic():
    mu,md,v,Phi,U,Ut,r,B=sp.symbols('mu md v Phi U Ut r B',nonzero=True)
    # v=c*a*mu and along an outgoing characteristic dU/dt=(1-mu)*Ut.
    old=mu*(1-mu*mu)*Phi*Ut-v*Phi*U*(1/B-1)/r+2*mu*md*Phi*U-mu*mu*v*Phi*U*(1+1/B)/r
    removed=mu*(1+mu)*Phi*(1-mu)*Ut+(1+2*mu)*md*Phi*U-mu*(1+mu)*Phi*(1+1/B)*v*U/r
    new=Phi*U*(-md+v/r*(1+mu-(1-mu)/B))
    assert sp.simplify(old-removed-new)==0
    result=base.symbolic();result.update(outgoing_transform='q=h-mu*Phi*U+mu^2*lambda_mass. The outgoing U(t-x/c) derivative cancels exactly; initial q=-mu0*Phi0*U_launch+mu0^2*lambda_mass_launch.',outgoing_identity=True)
    return result


def prepare():
    base.prepare();plan=read(OUT/'plan.json')
    plan['claim']='Complete the identical returned-metric photon characteristics after exactly removing both mass time jets and outgoing scalar time jets from the integrated variable.'
    plan['additional_transform']='q=h-mu*Phi*U+mu^2*lambda_mass; valid for the SAME declared outgoing scalar continuation only. Keep all initial/endpoint terms and reconstruct h before every arrival-delay evaluation.'
    plan['bindings'][str(Path(__file__))]=sha(__file__)
    write(OUT/'plan.json',plan);write(OUT/'symbolic.json',symbolic())


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(1800 if action=='pilot' else 600)
    try:
        if action!='prepare':
            for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
        if action=='prepare':prepare()
        else:
            install();base.install=lambda:None
            getattr(base,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
