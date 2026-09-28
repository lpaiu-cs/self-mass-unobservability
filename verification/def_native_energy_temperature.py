"""Reuse native states in the density/temperature coordinates of energy inversion."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
import def_native_energy_release as old

OUT=old.OUT
C=old.C


class EOS:
    def __init__(self,path=None):
        self.d=np.load(path or OUT/'temperature-columns.npz');d=self.d;raw=d['raw'];self.x=d['x'];self.rho0=float(d['rho0']);self.cx=float(d['cx']);self.floor=np.exp(self.x[0]);self.top=np.exp(self.x[-1]);self.sunit=float(d['sunit'])
        self.curves=[];self.low=[];self.high=[]
        for j in range(len(self.x)):
            section=slice(d['offsets'][j],d['offsets'][j+1]);a=raw[section];lt=d['logT'][section];T=np.exp(lt);self.low.append(lt[0]);self.high.append(lt[-1])
            values=[np.log(a[:,1]/(self.rho0*C*C)),a[:,2]/C**2,a[:,3]]
            thermal=[a[:,6],a[:,10]/C**2,a[:,10]/T]
            radial=[a[:,5],a[:,9]/C**2,-a[:,1]/a[:,0]*a[:,6]/T]
            self.curves.append(([CubicHermiteSpline(lt,v,q,extrapolate=False) for v,q in zip(values,thermal)],
                                [PchipInterpolator(lt,v,extrapolate=False) for v in radial],PchipInterpolator(lt,a[:,13]/a[:,0],extrapolate=False)))
        self.low=np.array(self.low);self.high=np.array(self.high)
        # Batch PPoly evaluation across ragged native temperature columns.
        # The32-unit separation exceeds the entire registered lnT range.
        self.starts=np.concatenate([f[0].x[:-1] for f,g,k in self.curves])
        self.search=np.concatenate([f[0].x[:-1]+32*j for j,(f,g,k) in enumerate(self.curves)])
        self.coeff=np.concatenate([np.stack([q.c for q in f+g+[k]]) for f,g,k in self.curves],axis=2)

    def limits(self,rho):
        x=np.log(np.maximum(rho,self.floor));i=np.clip(np.searchsorted(self.x,x)-1,0,len(self.x)-2)
        return np.maximum(self.low[i],self.low[i+1]),np.minimum(self.high[i],self.high[i+1])

    def evaluate(self,rho,lt):
        active=rho>=self.floor;assert max(rho)<=self.top*(1+1e-9)
        x=np.log(np.maximum(rho,self.floor));i=np.clip(np.searchsorted(self.x,x)-1,0,len(self.x)-2);lo,hi=self.limits(rho)
        assert np.all((lt[active]>=lo[active]-1e-12)&(lt[active]<=hi[active]+1e-12)),'Temperature outside native column overlap'
        lt=np.clip(lt,lo,hi);z=(x-self.x[i])/(self.x[i+1]-self.x[i]);h=self.x[i+1]-self.x[i];v=np.zeros((2,3,len(rho)));vr=v.copy();vt=v.copy();vrt=v.copy();kap=np.zeros((2,len(rho)))
        for side in range(2):
            at=np.searchsorted(self.search,lt+32*(i+side),side='right')-1;dt=lt-self.starts[at];a,b,c,d=self.coeff[:,:,at].transpose(1,0,2)
            values=((a*dt+b)*dt+c)*dt+d;derivatives=(3*a*dt+2*b)*dt+c
            v[side]=values[:3];vt[side]=derivatives[:3];vr[side]=values[3:6];vrt[side]=derivatives[3:6];kap[side]=values[6]
        h00=2*z**3-3*z*z+1;h01=-2*z**3+3*z*z;h10=z**3-2*z*z+z;h11=z**3-z*z
        value=h00*v[0]+h01*v[1]+h*(h10*vr[0]+h11*vr[1])
        thermal=h00*vt[0]+h01*vt[1]+h*(h10*vrt[0]+h11*vrt[1])
        radial=((6*z*z-6*z)*v[0]+(-6*z*z+6*z)*v[1])/h+(3*z*z-4*z+1)*vr[0]+(3*z*z-2*z)*vr[1]
        p=np.exp(value[0])*active;u=value[1]*active;cvT=thermal[1];chiT=thermal[0];chir=radial[0]
        gamma=chir+p/np.maximum(rho,self.floor)*chiT**2/cvT
        assert np.all(cvT[active]>0) and np.all(gamma[active]>1),'Native temperature interpolation stability'
        kk=np.maximum(0,(1-z)*kap[0]+z*kap[1])*active*6.6524587321e-25/1.66053906660e-24
        return p,u,gamma,np.exp(lt),kk,cvT,value[2]

    def __call__(self,rho,lt):return self.evaluate(rho,lt)[:5]


def audit():
    assert not (OUT/'temperature-audit.json').exists();start=time.monotonic();signal.alarm(10)
    old.write(OUT/'temperature-plan.json',dict(classification='Counterexample candidate',
        prior_entropy_interpolation_failed=True,actual_first_step_sigma=.7989337086070264,
        decision='Preserve the failed entropy interpolation and first temperature-column overlap failure. Use the separately registered native column repair, with the same rho,T Hermite derivatives and unchanged12 original independent controls.',
        new_native_states=0,control_call_cap=80,seconds=10,gate=.002,
        bindings={str(p.relative_to(old.old.ROOT)):old.old.photons.digest(p) for p in [Path(__file__),OUT/'temperature-columns.npz',OUT/'column-result.json',OUT/'eos.json',OUT/'first-step-native.json']}))
    eos=EOS();seed=old.EOS();fan=old.task.prior.Fan(call_cap=80,reuse=True);checks=[]
    for row in json.loads((OUT/'eos.json').read_text())['controls']:
        xx,z=row['log_density'],row['sigma'];rho=np.array([np.exp(xx)]);_,_,_,T,_=seed(rho,np.array([z]));lt=float(np.log(T[0]));target=float(eos.d['s0'])+z*eos.sunit
        for _ in range(12):
            a=fan.call(np.log(eos.rho0)+xx,lt);error=(a[3]-target)*np.exp(lt)/a[10]
            if abs(error)<2e-12:break
            lt-=np.clip(error,-.3,.3)
        else:raise AssertionError('Native independent control root')
        p,u,gamma,T,_=eos(rho,np.array([lt]));errors=[float(abs(p[0]*eos.rho0*C*C/a[1]-1)),float(abs(u[0]*C*C/a[2]-1)),float(abs(gamma[0]/a[4]-1))]
        checks.append(dict(log_density=xx,sigma=z,temperature_K=float(np.exp(lt)),relative=errors))
    passed=max(max(r['relative']) for r in checks)<.002
    result=dict(classification='Counterexample candidate',passed=passed,controls=checks,maximum_relative=max(max(r['relative']) for r in checks),native_calls=fan.calls,seconds=time.monotonic()-start)
    old.write(OUT/'temperature-audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':audit()
