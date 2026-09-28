"""Conservative nonrest Killing energy for the native radial gas release.

Counterexample candidate: closes the discrete gas/photon/inner-port energy
ledger on the saved metric. Native equilibrium chemistry is still a premise.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline
import def_native_radial_thermo as prior

task=prior.task
old=task.old
C=task.C
OUT=prior.OUT.parent/'def-native-energy-release'
write=task.write


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0a71a62d1',
        claim='Evolve conserved baryon, radial momentum and rest-subtracted Killing energy, with native EOS recovery. Close gas, escaped dilute inventory, inner transfer and opposite photon work before imposing GR mass constraints.',
        equations='K=a*(E-cX*D)+(a-a_surface)*cX*D, F_K=(K+a*p)*v; divergence uses the same spherical lapse/area flux and proper volume as baryons. The only gas energy source is the paired elastic photon work. No measured residual is added as heat.',
        EOS='Reuse the147 density states on all3 old entropy planes; add entropy planes[-0.02,0.03,0.1,0.3,0.7,1.5,3]. Native isentropes and entropy derivatives of pressure/u/T, no dilute extrapolation. Positive entropy range allows conservative numerical dissipation rather than forcing entropy advection.',
        native_budget=dict(calls=5000,seconds=180,temperature_K=[100,1e7]),
        evolution_budget=dict(seconds=240,CPU_threads=1,memory_GB=2),
        paths=['224-cell full-horizon pilot','896-cell saved-horizon run','1792-cell saved-horizon run'],
        same_horizon_seconds=.0034344311179287023,domain_m=[-200,1200],CFL=.35,
        gates=dict(native_EOS_relative=.002,primitive_energy_relative=1e-8,baryon_ledger_relative=1e-10,energy_ledger_over_response=1e-8,outside_mass_refinement=.02,trace_refinement=.02),
        stop='Stop on EOS/root/positivity, independent-control or budget failure. No automatic entropy/density extension, grid refinement, horizon enlargement or tolerance relaxation. Save failed states.',
        unchanged_scope='Frozen metric, LTE gas chemistry and optically thin instantaneous elastic scattering remain assumptions; full GR and physical chemical/radiative closure remain required.',
        bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(prior.__file__),Path(task.__file__),prior.OUT/'eos.npz',prior.OUT/'cells-1792.npz',prior.OUT/'flow-audit.json']}))


def bank():
    assert not (OUT/'eos.npz').exists();start=time.monotonic();signal.alarm(180)
    saved=np.load(prior.OUT/'eos.npz');x=saved['x'];sigma=np.sort(np.r_[saved['sigma'],[-.02,.03,.1,.3,.7,1.5,3.]])
    fan=task.prior.Fan(call_cap=5000,reuse=True);raw=np.zeros((len(sigma),len(x),21));temp=np.zeros((len(sigma),len(x)));done=np.zeros(temp.shape,bool)
    for j,z in enumerate(sigma):
        match=np.where(saved['sigma']==z)[0]
        if len(match):raw[j]=saved['raw'][match[0]];temp[j]=saved['T'][match[0]];done[j]=True
    def root(xx,z,guess):
        for _ in range(24):
            assert np.log(100)<guess<np.log(1e7),'Native temperature domain'
            a=fan.call(np.log(float(saved['rho0']))+xx,guess);err=(a[3]-float(saved['s0'])-z*float(saved['sunit']))*np.exp(guess)/a[10]
            if abs(err)<2e-12:return a,np.exp(guess)
            guess-=np.clip(err,-.3,.3)
        raise AssertionError('Native entropy root did not converge')
    count=0
    for j,z in enumerate(sigma):
        for i in range(len(x)-1,-1,-1):
            if done[j,i]:continue
            known=np.where(done[:,i])[0];k=known[np.argmin(abs(sigma[known]-z))]
            guess=np.log(temp[k,i])+np.clip((z-sigma[k])*float(saved['sunit'])*temp[k,i]/raw[k,i,10],-.8,.8)
            if i+1<len(x) and done[j,i+1]:guess=np.log(temp[j,i+1])+(raw[j,i+1,1]/raw[j,i+1,0]-raw[j,i+1,9])/raw[j,i+1,10]*(x[i]-x[i+1])
            raw[j,i],temp[j,i]=root(float(x[i]),float(z),float(guess));done[j,i]=True;count+=1
            if count%16==0:
                np.savez_compressed(OUT/'eos-progress.npz',x=x,sigma=sigma,raw=raw,T=temp,done=done,native_calls=fan.calls)
            if count==16:
                elapsed=time.monotonic()-start;forecast=elapsed*(7*len(x)/16)*1.25
                write(OUT/'native-budget.json',dict(measured16_states_seconds=elapsed,forecast_seconds=forecast,assumption='Same-density-step native roots across7 new entropy planes,25percent margin; high entropy costs not measured.'))
                assert forecast<175,'Native bank forecast exceeds budget'
    np.savez_compressed(OUT/'eos.npz',x=x,sigma=sigma,raw=raw,T=temp,rho0=saved['rho0'],cx=saved['cx'],s0=saved['s0'],sunit=saved['sunit'])
    eos=EOS();checks=[]
    for xx,z in [(-.037,-.009),(-.7,.015),(-3.3,.05),(-10.3,.2),(-.7,.5),(-3.3,1.),(-10.3,2.),(-17.7,.5),(-.14,.06),(-.86,.2),(-1.07,1.1),(-5.2,2.7)]:
        p,u,gamma,T,kap=eos(np.array([np.exp(xx)]),np.array([z]));a,tt=root(xx,z,float(np.log(T[0])))
        errors=[float(abs(p[0]*eos.rho0*C*C/a[1]-1)),float(abs(u[0]*C*C/a[2]-1)),float(abs(gamma[0]/a[4]-1)),float(abs(T[0]/tt-1))]
        checks.append(dict(log_density=xx,sigma=z,relative=errors))
    passed=max(max(r['relative']) for r in checks)<.002
    row=dict(classification='Counterexample candidate',passed=passed,controls=checks,native_calls=fan.calls,new_states=count,reused_states=441,seconds=time.monotonic()-start,
        minimum_T=float(temp.min()),maximum_T=float(temp.max()),maximum_relative=max(max(r['relative']) for r in checks))
    write(OUT/'eos.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert passed,'Native conservative EOS controls failed'


class EOS:
    def __init__(self):
        self.d=np.load(OUT/'eos.npz');d=self.d;raw=d['raw'];self.x=d['x'];self.sigma=d['sigma'];self.rho0=float(d['rho0']);self.cx=float(d['cx']);self.sunit=float(d['sunit']);self.floor=np.exp(self.x[0]);self.top=np.exp(self.x[-1])
        ad=(raw[:,:,1]/raw[:,:,0]-raw[:,:,9])/raw[:,:,10];ds=self.sunit*d['T']/raw[:,:,10]
        self.fun=[CubicHermiteSpline(self.x,np.log(raw[:,:,1]/(self.rho0*C*C)),raw[:,:,4],axis=1),
                  CubicHermiteSpline(self.x,raw[:,:,2]/C**2,raw[:,:,1]/raw[:,:,0]/C**2,axis=1),
                  CubicHermiteSpline(self.x,np.log(d['T']),ad,axis=1)]
        self.ds=[CubicSpline(self.x,raw[:,:,6]*ds,axis=1),CubicSpline(self.x,self.sunit*d['T']/C**2,axis=1),CubicSpline(self.x,ds,axis=1)]
        self.kap=CubicSpline(self.x,raw[:,:,13]/raw[:,:,0],axis=1)

    def __call__(self,rho,sigma):
        active=rho>=self.floor
        assert np.max(rho)<=self.top*(1+1e-9),'High density outside native table'
        assert np.all((sigma[active]>=self.sigma[0]-1e-9)&(sigma[active]<=self.sigma[-1]+1e-9)),'Entropy outside native table'
        xx=np.log(np.maximum(rho,self.floor));z=np.clip(sigma,self.sigma[0],self.sigma[-1]);i=np.clip(np.searchsorted(self.sigma,z)-1,0,len(self.sigma)-2);j=np.arange(len(rho));h=self.sigma[i+1]-self.sigma[i];f=(z-self.sigma[i])/h
        h00=2*f**3-3*f*f+1;h10=f**3-2*f*f+f;h01=-2*f**3+3*f*f;h11=f**3-f*f
        def blend(v,derivative):return h00*v[i,j]+h01*v[i+1,j]+h*(h10*derivative[i,j]+h11*derivative[i+1,j])
        lp,u,lt=[blend(fn(xx),der(xx)) for fn,der in zip(self.fun,self.ds)]
        gamma=blend(self.fun[0](xx,1),self.ds[0](xx,1));kk=self.kap(xx);kap=(1-f)*kk[i,j]+f*kk[i+1,j]
        p=np.exp(lp)*active;u*=active;T=np.exp(lt)
        assert np.all(gamma[active]>1) and np.all(u[active]>0),'Thermodynamic stability'
        return p,u,gamma,T,np.maximum(kap,0)*6.6524587321e-25/1.66053906660e-24*active


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','bank']);globals()[p.parse_args().action]()
