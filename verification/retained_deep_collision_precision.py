"""Same saved deep cubic EOS/rates in arbitrary precision; no new native states."""
import mpmath as mp
import numpy as np
from def_native_hydrogen_exchange import C as SPEED, K


def number(x):
    if isinstance(x,(float,np.floating)):
        p,q=x.as_integer_ratio();return mp.mpf(p)/q
    return mp.mpf(x)


def coefficients(state,k,z,factor,dps,pressure_only=False):
    m=state.m;model=m.model;b=model.bulk;e=b.eos;base=e.base;row=m.material.point(k)
    with mp.workdps(dps):
        n=m.nb;C=number(SPEED)
        eps=number(1e-26)*number(factor);cx=number(model.cx)
        h=[number(v) for v in row['h']];running=mp.mpf(0)
        for j in range(n):running+=number(z[0,j]);h[j+1]-=eps*running
        Q=[[number(row['Q'][i,j])+eps*number(z[i,j])*int(row['active'][j]) for j in range(n)] for i in range(4)]
        result=np.empty((2,n) if pressure_only else (2,n,m.q,m.nf),object);max_root=mp.mpf(0)
        def polynomial(co,j,t):
            index=int(np.clip(np.searchsorted(co.x,float(t),side='right')-1,0,len(co.x)-2))
            x=t-number(co.x[index]);c=co.c[:,index,j]
            values=[];der=[]
            for idx in np.ndindex(c.shape[1:]):
                p=[number(c[(i,)+idx]) for i in range(4)]
                values.append(((p[0]*x+p[1])*x+p[2])*x+p[3]);der.append((3*p[0]*x+2*p[1])*x+p[2])
            return np.array(values,object).reshape(c.shape[1:]),np.array(der,object).reshape(c.shape[1:])
        fraction=number(e.f)
        for j in range(n):
            mass=number(model.mass0[j])-h[j+1]+h[j]
            rho=mass/number(b.volume[j]);a=number(b.d['a'][j])
            beta=Q[1][j]/(cx*mass*C*C)
            eta=Q[3][j]/(mass*number(b.d['thermo'][j,4])*number(b.d['y0'][j]))-1
            u=(Q[2][j]-(a-number(model.m.a0))*cx*C*C*mass)/(a*mass)-cx*(beta*C)**2/2
            x=(h[j]-h[j+1])/number(model.mass0[j]);x=number(e.density_shift[j])+(1+number(e.density_shift[j]))*x
            xi=sum(number(v)*q for v,q in zip(model.mech.xi[j],h))
            w=(1+eta-number(base.ratios[0]))/number(np.diff(base.ratios)[0])
            t=number(row['theta'][j])+number(e.theta0[j]);assert 0<=w<=1
            ur=number(e.saved['ur'][0,j])*(1-fraction)+number(e.saved['ur'][1,j])*fraction
            for _ in range(12):
                T=number(base.d['T'][j])*mp.exp(t);q,qt=polynomial(base.uc,j,t)
                v=[number(base.u0[j,i])+number(base.R[j,i])*T*mp.mpf('1.5')+q[i] for i in range(2)]
                dv=[number(base.R[j,i])*T*mp.mpf('1.5')+qt[i] for i in range(2)]
                residual=v[0]*(1-w)+v[1]*w+ur*x-number(e.inventory[1,j])*xi-u
                step=residual/(dv[0]*(1-w)+dv[1]*w);t-=step
                if abs(step)<mp.power(10,-dps+12):break
            else:raise AssertionError('High precision conserved thermal root')
            max_root=max(max_root,abs(step));assert number(base.temperature_bounds[0])<=t<=number(base.temperature_bounds[1])
            if pressure_only:
                lp,_=polynomial(base.pc,j,t);T=number(base.d['T'][j])*mp.exp(t)
                pp=[mp.exp(lp[i])*number(base.d['rho'][j])*T for i in range(2)]
                pr=number(e.saved['pr'][0,j])*(1-fraction)+number(e.saved['pr'][1,j])*fraction
                gas=(pp[0]*(1-w)+pp[1]*w+pr*x-number(e.inventory[0,j])*xi)*number(b.volume[j])
                result[:,j]=[gas,gas+cx*mass*(beta*C)**2];continue
            lp,_=polynomial(base.rc,j,t);y=number(base.d['y0'][j])*(1+eta);assert 0<y<1
            T=number(base.d['T'][j])*mp.exp(t);rate=[]
            for channel in range(2):
                rr=[]
                for f in range(m.nf):
                    boltz=number(base.d['Einf'][f])/(number(base.d['a'][j])*number(K)*T) if channel else 0
                    p=[mp.exp(lp[i,channel,f]-boltz) if base.mask[j,i,channel,f] else mp.mpf(0) for i in range(2)]
                    dr=number(e.spectral['density'][0,j,channel,f])*(1-fraction)+number(e.spectral['density'][1,j,channel,f])*fraction
                    df=number(e.spectral['frequency'][0,j,channel,f])*(1-fraction)+number(e.spectral['frequency'][1,j,channel,f])*fraction
                    correction=(1+dr*x-number(e.rinventory[channel,j,f])*xi)*(1+df*number(e.frequency_shift[j]))
                    assert correction>0
                    rr.append((p[0]*(1-w)+p[1]*w)*(y if channel==0 else 1-y)*correction)
                rate.append(rr)
            scale=a*C*rho*number(b.d['thermo'][j,4])
            for angular,mu in enumerate(m.mu):
                for f in range(m.nf):
                    values=[]
                    for channel in range(2):
                        df=number(e.spectral['frequency'][0,j,channel,f])*(1-fraction)+number(e.spectral['frequency'][1,j,channel,f])*fraction
                        values.append(scale*rate[channel][f]*(1-beta*number(mu)*(1+df)))
                    result[0,j,angular,f]=values[1];result[1,j,angular,f]=values[0]-values[1]
        return result,float(max_root)


def difference(state,k,z,factor,dps,pressure_only=False):
    with mp.workdps(dps):
        a,ra=coefficients(state,k,np.zeros_like(z),0.,dps,pressure_only);b,rb=coefficients(state,k,z,factor,dps,pressure_only)
        delta=np.array([np.longdouble(str(v)) for v in (b-a).ravel()]).reshape(a.shape)
        center=np.array([np.longdouble(str(v)) for v in a.ravel()]).reshape(a.shape)
        gross=np.array([np.longdouble(str(abs(x)+abs(y))) for x,y in zip(a.ravel(),b.ravel())]).reshape(a.shape)
        return delta,center,gross*float(16*mp.eps),max(ra,rb)


class Pressure:
    """Finite conserved pressure; reuse the existing extended atmospheric owner."""
    def __init__(self,state):self.s=state;self.rows=[]
    def __call__(self,material,k,z,field):
        import time
        started=time.monotonic();s=self.s;m=s.m;nb=m.nb;LD=np.longdouble
        assert not np.any(field),'This source contains no new metric perturbation'
        s.precision(False);m.model.flow.seed=np.asarray(m.model.flow.seed,float);m.material.point(k)
        assert np.array_equal(material.point(k)['Q'],m.material.point(k)['Q'])
        s.precision(True)
        def pressure():
            f=m.model.flow;f.eos.y=m.y;p,u,*_=f.eos.evaluate(m.rho,m.lt);pg=p*LD(f.eos.rho0)*LD(SPEED)**2
            gas=np.r_[m.model.bulk.eos.gas(m.theta,m.eta)[0],pg]*m.material.V;radial=gas.copy()
            radial[:nb]+=2*m.model.kinetic()
            H=m.rho*(LD(f.eos.cx)+u)*LD(f.eos.rho0)*LD(SPEED)**2+pg
            radial[nb:]+=H*m.beta*m.beta/(1-m.beta*m.beta)*m.material.V[nb:]
            return np.array([gas,radial],LD)
        owner=m.coefficients;m.coefficients=pressure
        try:a=s.coefficients(k,np.zeros_like(z),0.);b=s.coefficients(k,z,1.)
        finally:m.coefficients=owner
        delta=b-a;error=16*np.finfo(LD).eps*(abs(a)+abs(b))
        coarse,center,_,_=difference(s,k,z,1.,70,True)
        fine,center,rounding,root=difference(s,k,z,1.,90,True)
        ownership=float(np.max(abs(a[:,:nb]-center))/max(np.max(abs(center)),1.))
        comparison=float(np.max(abs(coarse-fine))/max(np.max(abs(fine)),LD(1e-300)))
        assert ownership<1e-9 and comparison<1e-10 and root<1e-60
        delta[:,:nb]=fine;error[:,:nb]=abs(coarse-fine)+rounding
        changed=np.any(z!=0,axis=0);changed[:nb]|=(m.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:nb])])!=0
        assert np.all(delta[:,~changed]==0);error[:,~changed]=0
        self.rows.append(dict(k=k,steps=material.steps,owner=ownership,precision=comparison,root=root,seconds=time.monotonic()-started))
        self.last=dict(center=a,increment=delta,error=error,deep70=coarse,deep90=fine)
        s.precision(False);return a,delta/LD(1e-26),error/LD(1e-26)
