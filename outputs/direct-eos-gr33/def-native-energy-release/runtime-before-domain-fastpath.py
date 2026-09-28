"""Query the native EOS only where the actual conservative trajectory needs it."""
import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
import def_native_energy_temperature as parent

OUT=parent.OUT


class EOS(parent.EOS):
    def __init__(self):
        path=OUT/'runtime-columns.npz';super().__init__(path if path.exists() else None)
        self.prior_calls=int(self.d['runtime_native_calls']) if 'runtime_native_calls' in self.d else 0
        self.fan=parent.old.task.prior.Fan(call_cap=3000-self.prior_calls,reuse=True)
        self.cache=[]
        for j in range(len(self.x)):
            sl=slice(self.d['offsets'][j],self.d['offsets'][j+1]);self.cache.append(dict(zip(self.d['logT'][sl],self.d['raw'][sl])))

    def limits(self,rho):return np.full_like(rho,np.log(100.)),np.full_like(rho,np.log(1e7))

    def save(self):
        ts=[];rows=[];offsets=[0]
        for cache in self.cache:
            keys=sorted(cache);ts.extend(keys);rows.extend([cache[t] for t in keys]);offsets.append(len(ts))
        path=OUT/'runtime-columns.tmp.npz'
        np.savez_compressed(path,x=self.x,logT=ts,raw=rows,offsets=offsets,rho0=self.rho0,cx=self.cx,sunit=self.sunit,s0=self.d['s0'],runtime_native_calls=self.prior_calls+self.fan.calls)
        path.replace(OUT/'runtime-columns.npz')

    def extend(self,j,lo,hi):
        cache=self.cache[j];oldlo=min(cache);oldhi=max(cache);newlo=max(np.log(100.)+1e-10,lo-.04);newhi=min(np.log(1e7)-1e-10,hi+.04)
        def call(t):
            if t not in cache:cache[t]=self.fan.call(np.log(self.rho0)+self.x[j],float(t))
            return cache[t]
        def refine(l,r,depth):
            a,b=call(l),call(r);mid=(l+r)/2;n=call(mid)
            lp=CubicHermiteSpline([l,r],np.log([a[1],b[1]]),[a[6],b[6]])
            u=CubicHermiteSpline([l,r],[a[2],b[2]],[a[10],b[10]])
            p=float(np.exp(lp(mid)));cv=float(u(mid,1));cr=(a[5]+b[5])/2;gamma=cr+p/n[0]*float(lp(mid,1))**2/cv
            err=max(abs(p/n[1]-1),abs(float(u(mid))/n[2]-1),abs(cv/n[10]-1),abs(cr/n[5]-1),abs(gamma/n[4]-1))
            if err>.0005:
                assert depth<7,'Runtime native column depth'
                refine(l,mid,depth+1);refine(mid,r,depth+1)
        if lo<oldlo:refine(newlo,oldlo,0)
        if hi>oldhi:refine(oldhi,newhi,0)
        ts=sorted(cache);a=np.array([cache[t] for t in ts]);T=np.exp(ts);C=parent.C
        values=[np.log(a[:,1]/(self.rho0*C*C)),a[:,2]/C**2,a[:,3]]
        thermal=[a[:,6],a[:,10]/C**2,a[:,10]/T];radial=[a[:,5],a[:,9]/C**2,-a[:,1]/a[:,0]*a[:,6]/T]
        self.curves[j]=([CubicHermiteSpline(ts,v,q,extrapolate=False) for v,q in zip(values,thermal)],
                        [PchipInterpolator(ts,v,extrapolate=False) for v in radial],PchipInterpolator(ts,a[:,13]/a[:,0],extrapolate=False))
        self.low[j]=ts[0];self.high[j]=ts[-1]

    def evaluate(self,rho,lt):
        active=rho>=self.floor;x=np.log(np.maximum(rho,self.floor));i=np.clip(np.searchsorted(self.x,x)-1,0,len(self.x)-2)
        changed=False
        try:
            for side in range(2):
                for j in np.unique((i+side)[active]):
                    mask=active&(i+side==j);lo=float(min(lt[mask]));hi=float(max(lt[mask]))
                    if lo<self.low[j] or hi>self.high[j]:self.extend(j,lo,hi);changed=True
        except Exception:
            self.save();raise
        if changed:
            self.starts=np.concatenate([f[0].x[:-1] for f,g,k in self.curves])
            self.search=np.concatenate([f[0].x[:-1]+32*j for j,(f,g,k) in enumerate(self.curves)])
            self.coeff=np.concatenate([np.stack([q.c for q in f+g+[k]]) for f,g,k in self.curves],axis=2)
            self.save()
        lo,hi=parent.EOS.limits(self,rho);lt=np.where(active,lt,np.clip(lt,lo,hi))
        return super().evaluate(rho,lt)
