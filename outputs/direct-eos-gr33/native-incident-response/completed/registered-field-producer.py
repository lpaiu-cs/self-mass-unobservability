"""Counterexample candidate: compact incoming scalar Cauchy data and GR input.

This is an external-input connection test, NOT the wide-binary benchmark.
Use the existing optical Green operator for the first potential return; retain
the higher-return and interpolation budgets separately from physical errors.
"""
from pathlib import Path
import json,resource,signal,sys,time
import numpy as np
import sympy as sp
from numpy.polynomial import legendre as leg
import def_retained_native_acoustic as native
import def_native_characteristic_gr as wave
import def_native_dynamic_lapse as metric

OUT=Path('native-incident-drive155-work');FIELDS=OUT/'fields';METRIC=OUT/'metric'
PHOTON=OUT/'photons';MATERIAL=OUT/'material';GR=OUT/'gr'
read,write,sha=native.read,native.write,native.sha
C=wave.C;G=wave.base.G;LD=np.longdouble;ETA=1e-30
CAPS=dict(fields=120,photon_pilot=90,photon_production=1000,material_pilot=60,material_production=500,readout=180,audit=60)


def pulse(s,derivative=0):
    s=np.asarray(s);x=np.clip(s,0,1)
    # 256*s^4*(1-s)^4, compact C3 incoming packet; peak=1.
    coef=np.r_[np.zeros(4),256*np.array([1.,-4.,6.,-4.,1.])]
    v=np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(coef,derivative))
    return np.where((s>0)&(s<1),v,0.)


def symbolic():
    r,a,R,t,x,D,A,c=sp.symbols('r a R t x D A c',positive=True);g=sp.Function('g')
    U=A*R*g((t+x/c)/D)
    assert sp.simplify(sp.diff(U,t,2)/c**2-sp.diff(U,x,2))==0
    assert sp.simplify(sp.diff(U,t)-c*sp.diff(U,x))==0
    z=sp.symbols('z');p=256*z**4*(1-z)**4
    assert p.subs(z,sp.Rational(1,2))==1
    assert all(sp.diff(p,z,k).subs(z,b)==0 for k in range(4) for b in [0,1])
    return dict(classification='Proven',passed=True,incoming='Uin=eta*r0*g((t+x/c)/D), x(r0)=0; Ut=c*Ux. At t=0 support is outside the star,0<x<cD.',
        potential='U=Uin+L(Uin)+remainder; L=-Gret V. If norm(L)<=etaV<1, norm(remainder)<=etaV^2/(1-etaV)*norm(Uin). This is conditional on that operator bound, not a physical EOS certificate.',
        physical_scope='Initial scalar Cauchy packet on the same frozen nonzero background, not an applied thermal source, static parameter fit or companion-matched orbital input.',
        lapse=metric.symbolic())


def prepare():
    assert not OUT.exists()
    for p in [FIELDS,METRIC/'corrected',PHOTON,MATERIAL,GR]:p.mkdir(parents=True)
    paths=[Path(__file__),Path(wave.__file__),Path(wave.base.__file__),Path(metric.__file__),
           wave.base.flow.INPUT/'balanced-20.npz',native.prior.EV/'source-128.npz',
           Path('retained-static-response154-work/result.json'),Path('retained-static-response154-work/independent.json')]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='efe578325ca021051c1ad1f68856b2337e4e5c29',
        claim='Apply genuinely nonzero incoming scalar Cauchy data, its consistent first-order polar metric, and reciprocal conformal work to the CURRENT retained material/photon response; return the resulting material source to scalar charge.',
        drive=dict(amplitude=ETA,shape='256*s^4*(1-s)^4 on0<s<1',duration='half the existing3.4344311179287023ms interval',initial_support='Outside the actual initial photon edge, in0<x<cD; no initial interior scalar or fluid increments. Exterior scalar lapse can be nonzero initially.',
            companion_matched=False,reason='Connection test for the missing external-input operator. Do not rescale it into the declared wide-binary response or count it as the complete objective.'),
        field='Exact incoming free optical characteristic plus one retarded V return on the actual initial anisotropic GR coefficients. Use the existing characteristic-cut Green quadrature at4/8 radial orders and33 source times; vacuum optical panels extend to2*c*T only for the scalar/lapse calculation, not new fluid cells. Initial exterior data are retained.',
        lapse='delta_m=r^2*b*Phi*f (J=0 for the direct incident component), delta_lambda=r*Phi*f. Integrate the actual canonical gas/photon polar lapse constraint inward from the compact exterior scalar lapse integral, with lapse fixed at infinity. No arbitrary zero lapse at the numerical edge.',
        stages='Evaluate the primary incoming waveform and derivatives at actual SDIRK/SSP stage times. Interpolate only the small first potential return; do not approximate the main pulse with17 saved driving samples.',
        frequency='The unperturbed stationary Eref drift is zero. Use the odd conservative frequency-source derivative,0.5*(F(omega)-F(-omega)), removing the non-differentiable upwind numerical-diffusion term proportional to abs(omega). This is an explicitly declared linear finite-frequency representation, not the derivative of the old upwind scheme. Preserve its exact number/energy/work ledgers.',
        matter='Reuse the existing fine retained trajectory, corrected-EOS collision banks, actual stage-aligned free material equations and compensated increment representation. Fixed deep incoming photon perturbation remains a controlled boundary, not a solved full-core response. Native state Jacobian and nonlinear Einstein closure remain open.',
        budget=CAPS,CPU_threads=1,virtual_GiB=3,new_native_calls=0,new_background_steps=0,
        gates=dict(field_quadrature=.002,field_remainder_over_incident=1e-6,null=1e-12,frequency=1e-12,linear=1e-12,
            time=.02,conservation=1e-8,material_directional=.002,physical_branch=.01,small_state=1e-6,pressure=.002,independent_GR=1e-9),
        dispatch='First measure fields and equal-horizon4/8 photon prefixes. Dispatch only if2x measured remaining step cost including the recorded late-step floor fits1000s. Then measure2-step free-material prefixes with prior full raw-call counts before500s. Reuse all accepted prefixes.',
        stop='Keep amplitude, field family, clocks, physical interval and gates fixed. Stop on a failed gate or forecast; do not automatically refine, add periods, weaken criteria or reinterpret the pulse as an orbit.',
        bindings={str(p):sha(p) for p in paths}))
    write(OUT/'symbolic.json',symbolic())


class Driver:
    def __init__(self,order=8,build=False):
        self.order=order;self.response=m=wave.base.Response();self.bg=m.bg;self.model=m.model;self.geo=m.geo
        self.edges=np.asarray(m.faces);self.r0=float(self.edges[-1]);self.n=len(self.edges)-1
        self.T=float(np.load(native.prior.EV/'source-128.npz')['t'][-1]);self.D=self.T/2
        self.origin=C*self.geo(np.array([self.geo.physical(np.array([self.r0]))[0][0]-self.model.m.RJ]))[1][0]
        rJ=np.r_[self.model.bulk.d['r'],self.model.m.r]
        self.centers=self.geo.metric(rJ-self.model.m.RJ)[3]*self.model.m.R
        self.q=q=wave.base.flow.initial.Quadrature(self.edges,order)
        self.r=q.r.ravel();self.x=self.optical(self.r);self.z=m.coeff(self.r)
        self.a=self.z['lapse']*np.sqrt(self.z['b']);self.zc=m.coeff(self.centers)
        self.xc=self.optical(self.centers);self.z0=m.coeff(np.array([self.r0]))
        self.ex=np.linspace(0,2*C*self.T,17);self.er=self.inverse(self.ex)
        self.eq=wave.base.flow.initial.Quadrature(self.er,order)
        self.rout=self.eq.r.ravel();self.xout=self.optical(self.rout);self.zout=m.coeff(self.rout)
        self.tx=np.r_[self.optical(self.edges[:1]),self.xc,0.,self.ex[1:]]
        self.tr=np.r_[self.edges[0],self.centers,self.r0,self.er[1:]]
        self.times=np.linspace(0,self.T,33)
        if not build:
            z=np.load(FIELDS/f'born-g{order}.npz');self.born={k:z[k] for k in ['U','Ut','Ux']}
            assert np.array_equal(z['x'],self.tx) and np.array_equal(z['t'],self.times)
        self.partial=leg.legvander((self.centers-self.edges[:-1])/q.h-1,order-1)
        self.inv=np.linalg.inv(leg.legvander(q.x,order-1)).astype(LD)
        self.cache={}

    def optical(self,r):
        rJ=self.geo.physical(np.asarray(r))[0]
        return C*self.geo(rJ-self.model.m.RJ)[1]-self.origin

    def inverse(self,x):
        r=self.r0+np.asarray(x).copy()
        for _ in range(4):
            z=self.bg.fields(r);r-=(self.optical(r)-x)*z['lapse']*np.sqrt(z['b'])
        return r

    def wave(self,t,x):
        s=(t+np.asarray(x)/C)/self.D;u=ETA*self.r0*pulse(s);ut=ETA*self.r0/self.D*pulse(s,1)
        if hasattr(self,'born'):
            j=min(max(np.searchsorted(self.times,t,side='right')-1,0),31);w=(t-self.times[j])/(self.times[j+1]-self.times[j])
            out=[np.interp(x,self.tx,(1-w)*self.born[k][j]+w*self.born[k][j+1]) for k in ['U','Ut','Ux']]
            return u+out[0],ut+out[1],ut/C+out[2]
        return u,ut,ut/C

    def born_return(self):
        m=wave.Response.__new__(wave.Response);m.order=self.order;m.t=self.times
        m.xfaces=np.r_[self.optical(self.edges),self.ex[1:]];m.mid=(m.xfaces[:-1]+m.xfaces[1:])/2;m.half=np.diff(m.xfaces)/2
        m.x=np.r_[self.x,self.xout];m.ids=np.repeat(np.arange(len(m.mid)),self.order)
        aout=self.zout['lapse']*np.sqrt(self.zout['b'])
        m.dx=np.r_[(self.q.h[:,None]*self.q.w/self.a.reshape(self.q.r.shape)).ravel(),
                   (self.eq.h[:,None]*self.eq.w/aout.reshape(self.eq.r.shape)).ravel()]
        xx=(m.x.reshape(-1,self.order)-m.mid[:,None])/m.half[:,None]
        m.inverse=np.linalg.inv(leg.legvander(xx,self.order-1));m.gx,m.gw=leg.leggauss((self.order+3)//2)
        m.tx=self.tx;m.distance=abs(m.tx[:,None]-m.x[None,:])/C;m.sign=np.sign(m.tx[:,None]-m.x[None,:])
        V=np.r_[self.z['V'],self.zout['V']];incoming=ETA*self.r0*pulse((self.times[:,None]+m.x/C)/self.D)
        source=-m.dx*V*incoming;start=time.monotonic();values=m.propagate(source)
        self.born=dict(zip(['U','Ut','Ux'],values));eta=float(C*self.T/2*np.sum(m.dx*abs(V),dtype=LD))
        # The estimate is conditional on this finite coefficient representation;
        # it does not enclose errors in the EOS/background or radial continuum.
        remainder=eta*eta/(1-eta);assert eta<.01 and remainder<1e-6
        np.savez_compressed(FIELDS/f'born-g{self.order}.npz',t=self.times,x=self.tx,r=self.tr,**self.born,
            source_radius=np.r_[self.r,self.rout],source_x=m.x,V=V,dx=m.dx,eta=eta)
        return dict(order=self.order,seconds=time.monotonic()-start,potential_norm_estimate=eta,
            higher_returns_over_incident=remainder,maximum_first_return_over_incident=float(np.max(abs(values[0]))/(ETA*self.r0)))

    def variations(self,r,z,x,t):
        U,Ut,Ux=self.wave(t,x);a=z['lapse']*np.sqrt(z['b']);f=U/r;ft=Ut/r;fr=Ux/(a*r)-U/r**2
        lam=r*z['Phi']*f;lt=r*z['Phi']*ft;dm=r*z['b']*lam;volume=3*z['alpha']*f+lam
        P=z['Pg']+z['Pr'];dp=-z['Kg']*volume-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*lam
        nr=(1+8*np.pi*r*r*z['A4']*P)*dm/(r*r*z['b']**2)+4*np.pi*r*z['A4']*(dp+4*z['alpha']*P*f)/z['b']+r*z['Phi']*fr
        return f,ft,fr,lam,lt,nr

    def at(self,t):
        t=float(t)
        if t in self.cache:return self.cache[t]
        f,ft,fr,lam,lt,nr=self.variations(self.r,self.z,self.x,t)
        Uout=self.wave(t,self.xout)[0];U0=self.wave(t,np.array([0.]))[0][0]
        integrand=2*self.zout['Phi']*Uout/(self.rout*self.zout['b'])
        exterior=np.sum(self.eq.h*(integrand.reshape(self.eq.r.shape)@self.eq.w),dtype=LD)
        boundary=self.z0['Phi'][0]*U0-exterior
        integ,faces=self.q.integrate(nr);nu=boundary-(faces[-1]-integ)
        centers=np.sum((nu@self.inv.T)*self.partial,axis=1,dtype=LD).astype(float)
        f,ft,fr,lam,lt,nr=self.variations(self.centers,self.zc,self.xc,t)
        alpha=self.zc['alpha'];u=alpha*f;ut=alpha*ft;up=-4*self.zc['Phi']*f+alpha*fr
        row=dict(delta_phi=f,delta_u=u,delta_u_t=ut,delta_u_prime=up,delta_lambda=lam,
            delta_lambda_rate=lt,delta_nu=centers,delta_nu_prime=nr,
            delta_log_lapse=centers+u,delta_log_radial_length=lam+u,delta_log_areal_radius=u,
            delta_log_volume=lam+3*u,delta_log_speed=centers-lam)
        self.cache={t:row};return row

    def view(self,t,clock):
        """Stage-local arrays consumed by the unchanged coupled owners.

        Values interpolate to the exact queried field; u/lambda slopes equal
        its queried time derivatives. These arrays are never saved as history.
        """
        j=max(0,min(np.searchsorted(clock,t,side='left')-1,len(clock)-2));h=clock[j+1]-clock[j];w=(t-clock[j])/h
        row=self.at(t);g={k:np.broadcast_to(v,(len(clock),len(v))).copy() for k,v in row.items()}
        for key,der in [('delta_u','delta_u_t'),('delta_lambda','delta_lambda_rate')]:
            g[key][j]=row[key]-w*h*row[der];g[key][j+1]=row[key]+(1-w)*h*row[der]
        g['delta_lambda_interval_rate']=np.broadcast_to(row['delta_lambda_rate'],(len(clock)-1,self.n))
        return g


def fields():
    native.prior.initialize();start=time.monotonic();rows=[];saved={}
    for order in [4,8]:
        d=Driver(order,True);rows.append(d.born_return());clock=np.linspace(0,d.T,17)
        records=[d.at(t) for t in clock];values={k:np.array([r[k] for r in records]) for k in records[0]}
        values['delta_lambda_interval_rate']=np.diff(values['delta_lambda'],axis=0)/np.diff(clock)[:,None]
        values['delta_nu_interval_rate']=np.diff(values['delta_nu'],axis=0)/np.diff(clock)[:,None]
        values.update(t=clock,radius_E=d.centers)
        for n in [64,128]:np.savez_compressed(METRIC/f'corrected/metric-{n}-g{order}.npz',**values)
        saved[order]=values
    errors={k:float(np.max(abs(saved[4][k]-saved[8][k]))/max(np.max(abs(saved[8][k])),1e-290))
            for k in ['delta_u','delta_lambda','delta_log_lapse','delta_log_speed','delta_u_prime','delta_nu_prime']}
    # Cauchy interior data are unchanged; exterior pulse sets a small lapse.
    assert np.max(abs(saved[8]['delta_phi'][0]))==0
    result=dict(classification='Counterexample candidate',passed=max(errors.values())<.002,rows=rows,quadrature=errors,
        initial_interior_scalar_zero=True,initial_exterior_scalar_nonzero=True,
        initial_interior_lapse_maximum=float(np.max(abs(saved[8]['delta_nu'][0]))),
        amplitude=ETA,duration_seconds=d.D,horizon_seconds=d.T,seconds=time.monotonic()-start,
        potential_remainder_is_conditional_finite_operator_estimate=True,full_physical_input_error_enclosed=False,
        actual_input_applied_to_material_photons=False,companion_matched=False,full_goal_complete=False)
    write(METRIC/'result.json',result);print(json.dumps(result));assert result['passed'],errors


def symbolic_check():
    assert read(OUT/'prepare-receipt.json')['error']=='AssertionError()'
    write(OUT/'symbolic.json',symbolic())
    write(OUT/'symbolic-repair.json',dict(classification='Proven',passed=True,
        change='Keep c symbolic in the proof. Floating c introduced unequal rounded reciprocal/square constants in an exact SymPy equality; the numerical wave code and all physical parameters are unchanged.',
        original_sha256=sha(OUT/'registered-prepare-producer.py'),producer_sha256=sha(__file__)))


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','fields','symbolic_check'];resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    native.deadline(CAPS[action] if action in CAPS else 30);start=time.monotonic();cpu=time.process_time();error=None
    try:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():
            receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
            write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
                peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
