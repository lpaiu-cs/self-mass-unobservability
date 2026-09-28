"""Counterexample candidate: positive angular photons coupled to native H/heat.

Retain the saved16 radial cells and native spectral tables. Transport angular
occupations, rather than reconstructing incompatible staggered P1 moments.
"""
from pathlib import Path
import json
import signal
import sys
import time
from types import FunctionType
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.interpolate import CubicSpline
import def_native_causal_nonlinear as old

prior=old.prior
C=prior.C
OUT=prior.OUT.parent/'def-native-angular-photons'
write=old.write
sha=prior.sha


class SupportedTable(old.Table):
    __init__=FunctionType(old.Table.__init__.__code__,dict(old.Table.__init__.__globals__,OUT=OUT/'neutral-support'))
    temperature_bounds=(-.06,.06)


class ThermalTable(SupportedTable):
    temperature_bounds=(-.06,.18)

    def __init__(self):
        self.d=d=np.load(prior.OUT/'bank-16-8.npz');v=np.load(OUT/'thermal-support/bank.npz')
        self.raw=v['raw'];self.rates=v['rates'];self.ratios=v['ratios'];self.n=len(d['r'])
        def poly(a):return CubicSpline(v['offsets'],np.moveaxis(a,2,0),extrapolate=False)
        rho=d['rho'][:,None];T=d['T'][:,None];temps=T[:,:,None]*np.exp(v['offsets'])[None,None,:]
        self.R=self.raw[:,:,1,1]/(rho*T);self.u0=self.raw[:,:,1,2]-1.5*self.R*T
        self.uc=poly(self.raw[:,:,:,2]-self.u0[:,:,None]-1.5*self.R[:,:,None]*temps)
        self.pc=poly(np.log(self.raw[:,:,:,1]/(rho[:,:,None]*temps)));self.nec=poly(self.raw[:,:,:,13]/1.66053906660e-24)
        self.mask=self.rates[:,:,1]>0;logs=np.log(np.maximum(self.rates,1e-300))
        logs[:,:,:,1]+=d['Einf'][None,None,None,:]/(d['a'][:,None,None,None]*prior.K*temps[:,:,:,None]);self.rc=poly(logs)

    def polynomial(self,co,t):
        ids=np.clip(np.searchsorted(co.x,t,side='right')-1,0,len(co.x)-2);c=co.c[:,ids,np.arange(self.n)]
        x=(t-co.x[ids]).reshape((self.n,)+(1,)*(c.ndim-2))
        return ((c[0]*x+c[1])*x+c[2])*x+c[3],(3*c[0]*x+2*c[1])*x+c[2]


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='f2be22103',
        claim='Replace the unresolved staggered P1 photon field by nonnegative angular occupations in the same actual native interior, and evolve those occupations together with nonlinear hydrogen and heat.',
        decision='A positive, conservative angular solution with bounded time/angle comparisons may supply the spectrum/flux for the moving atmosphere; no final charge from the old incompatible moment reconstruction.',
        reuse='Same16 native radii,152 frequencies,96 nonlinear EOS/rate states and3.434ms horizon. No repeated native bank, no fluid or full-GR rerun.',
        model='Axisymmetric angular finite volumes at fixed Killing frequency in the saved Jordan metric. Conservative radial and angular upwind streaming; nonlinear reciprocal H bound-free exchange and the normalized angle-dependent Thomson kernel. Fixed density/metric and other local initial ionic inventories remain.',
        initial='Preserve local Planck mean, saved grey flux and K/J=1/3. For flux ratio<=1/2 use two constant hemispheres. Where flux is slightly larger, use a positive quadratic outgoing hemisphere with the same three discrete moments and zero incoming hemisphere. This supplies previously unspecified higher angular moments, not a uniquely inferred observed field.',
        boundary='Vacuum incoming at the outer surface. Inner incoming angular spectrum fixed to its actual saved background; inner outgoing is transported. Export the actual outgoing surface occupation and both Killing energy and photon-number port histories.',
        discretization='Backward Euler positive streaming, followed by conservative nonlinear implicit local exchange. First-order splitting; independent64/128 time and8/16 angle comparisons. Photon positivity is checked without clipping.',
        gates=dict(energy=1e-8,local_exchange=1e-9,positive_photons=True,moment_variance=1e-12,time_trace=.02,time_surface_spectrum=.02,angle_trace=.02,angle_surface_spectrum=.02),
        budget=dict(pilot_angles=8,pilot_steps=32,paths=[[8,64],[8,128],[16,128]],transport_seconds=120,native_audit_calls=80,native_audit_seconds=15,CPU_threads=1,memory_GB=2),
        forecast='Phase109 global nonlinear sparse solve took41.7s. Angular transport here factors only128/256 streaming unknowns once and solves local2variable material equations; actual cost unmeasured. Measure pilot before production; forecast with3x angular scaling margin.',
        stop='No automatic radius,frequency,time,angle or native-support enlargement. Preserve a failed path. Positive angular moments do not certify radial/frequency convergence, physical spectrum, moving gas, all reaction channels or full GR/scalar feedback.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),prior.OUT/'bank-16-8.npz',old.OUT/'bank.npz',old.OUT/'repaired-bank.json',old.OUT/'audit.json']}))
    import sympy as s
    r,a,mu,ap=s.symbols('r a mu ap',real=True);area=r*r/a**2
    divergence=(2*r/a**2-2*r*r*ap/a**3)*mu+s.diff(area*(1/r-ap/a)*(1-mu*mu),mu)
    assert s.simplify(divergence)==0
    u,v=s.symbols('u v',real=True);kernel=s.Rational(1,2)*(1+s.Rational(1,2)*s.legendre(2,u)*s.legendre(2,v))
    assert s.integrate(kernel,(u,-1,1))==1
    assert s.integrate(kernel*u,(u,-1,1))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        streaming='div_(r,mu)[r^2/a^2*(mu,(1/r-a_prime/a)*(1-mu^2))]=0. Conserved occupation measure is r^2*B/a^3 dr dmu.',
        Thomson='Kernel=(1+P2(mu)*P2(mu_prime)/2)/2 is nonnegative and normalized; its first angular outgoing moment is zero. Discrete bin-averaged P2 retains normalization.',
        positivity='Upwind transport has nonnegative off-diagonal rates; backward Euler is positivity preserving for the conservative transport generator. A nonnegative angular density has J>=0,H^2<=J*K and K<=J by Cauchy-Schwarz. These algebraic statements do not prove continuum resolution.'))


class Model:
    def __init__(self,angles=8):
        assert angles>=8 and angles%2==0
        self.eos=ThermalTable() if (OUT/'thermal-support/bank.npz').exists() else SupportedTable() if (OUT/'neutral-support/bank.npz').exists() else old.Table()
        self.d=d=self.eos.d;self.n=n=len(d['r']);self.m=m=len(d['Einf']);self.angles=q=angles
        self.edges_mu=np.linspace(-1,1,q+1);self.mu=(self.edges_mu[:-1]+self.edges_mu[1:])/2;dm=np.diff(self.edges_mu);self.w=dm/2
        self.mu2=self.mu**2+dm**2/12;self.P2=(3*self.mu2-1)/2;self.meanP22=float(self.w@(self.P2**2))
        r,a,B,rho,T=[d[k] for k in ['r','a','B','rho','T']];edges=d['edges'];dx=np.diff(edges)
        self.W=r*r*B/a**3*dx;self.area=edges**2/d['face_a']**2;self.volume=4*np.pi*r*r*B*dx
        self.J0=1/np.expm1(d['Einf'][None,:]/(a*prior.K*T)[:,None]);self.scale=self.J0.max(0)
        bg=prior.Background();L=bg.luminosity(dict(r=r,a=a,B=B,rho=rho,T=T));Erad=(self.J0*d['num']*d['Einf']).sum(1)/a**4
        fraction=L/(4*np.pi*r*r*a*a*C*Erad)
        self.initial=self.J0[:,None,:]*self.shape(fraction)[:,:,None]
        ji=1/np.expm1(d['Einf']/(d['face_a'][0]*prior.K*d['face_T'][0]));ei=(ji*d['num']*d['Einf']).sum()/d['face_a'][0]**4
        fi=d['face_L'][0]/(4*np.pi*edges[0]**2*d['face_a'][0]**2*C*ei)
        self.incoming=ji[None,:]*self.shape(np.array([fi]))[0,:,None]*(self.mu>0)[:,None]
        self.u0=self.eos.gas(np.zeros(n),np.zeros(n))[1];self.cv=d['thermo'][:,0]
        self.en=d['num'][None,:]*d['Einf']/(a**4*rho)[:,None]
        self.num=d['num'][None,:]/(a**3*rho*d['thermo'][:,4]*d['y0'])[:,None]
        self.photon_energy_weight=4*np.pi*self.W[:,None,None]*self.w[None,:,None]*d['num']*d['Einf']
        self.gas_weight=self.volume*rho*a
        self.factor=a*C*rho*d['thermo'][:,4]
        ids=np.arange(n*q).reshape(n,q);rr=[];cc=[];vv=[];boundary=np.zeros((n,q,m))
        def add(i,j,v):rr.append(int(i));cc.append(int(j));vv.append(float(v))
        for k,mu in enumerate(self.mu):
            for face in range(n+1):
                flux=C*self.area[face]*mu
                if face==0:
                    if mu>0:boundary[0,k]+=flux/self.W[0]*self.incoming[k]
                    else:add(ids[0,k],ids[0,k],flux/self.W[0])
                elif face==n:
                    if mu>0:add(ids[-1,k],ids[-1,k],-flux/self.W[-1])
                    # No exterior incoming occupation.
                else:
                    donor=face-1 if mu>0 else face
                    add(ids[face-1,k],ids[donor,k],-flux/self.W[face-1])
                    add(ids[face,k],ids[donor,k],flux/self.W[face])
        # The integrated angular face coefficient uses the exact same radial
        # area difference; its divergence cancels radial focusing on I=const.
        for i in range(n):
            for k in range(1,q):
                flux=C*(self.area[i+1]-self.area[i])/2*(1-self.edges_mu[k]**2)
                assert flux>=0
                add(ids[i,k-1],ids[i,k-1],-flux/(self.W[i]*dm[k-1]))
                add(ids[i,k],ids[i,k-1],flux/(self.W[i]*dm[k]))
        self.A=sparse.coo_matrix((vv,(rr,cc)),shape=(n*q,n*q)).tocsc();self.boundary=boundary
        off=self.A-sparse.diags(self.A.diagonal());assert off.data.min(initial=0)>=0
        self.initial_moment_error=float(max(np.max(abs(self.mean(self.initial)/self.J0-1)),np.max(abs(self.mean(self.initial*self.mu[None,:,None])/self.J0-fraction[:,None])),np.max(abs(self.mean(self.initial*self.mu2[None,:,None])/self.J0-1/3))))
        assert self.initial_moment_error<1e-12

    def mean(self,I):return np.einsum('iqf,q->if',I,self.w)

    def shape(self,f):
        shape=1+2*f[:,None]*np.where(self.mu>0,1.,-1.)[None,:]
        # Match the same J,H,K in a positive outgoing-only polynomial where
        # the retained grey flux slightly exceeds the half-isotropic ratio.
        pos=self.mu>0;z=self.mu[pos];dm=2*self.w[pos]
        p1=2*z-1;p2=6*(z*z+dm*dm/12)-6*z+1
        M=np.array([[dm@(z*p1),dm@(z*p2)],[dm@(self.mu2[pos]*p1),dm@(self.mu2[pos]*p2)]])
        for i in np.flatnonzero(f>.5):
            aa,bb=np.linalg.solve(M,np.array([f[i]-.5,0.]));shape[i]=0.;shape[i,pos]=2*(1+aa*p1+bb*p2)
        assert np.min(shape)>=0 and max(abs(shape@self.w-1))<1e-12
        return shape

    def exchange(self,I,theta,eta,h,exponential=False):
        Jstar=self.mean(I);Kstar=self.mean(I*self.P2[None,:,None]);oldeta=eta.copy();oldtheta=theta.copy();uold=self.eos.gas(theta,eta)[1]
        def equations(t,y,derivatives=True):
            p,u,ut,uy,*_=self.eos.gas(t,y)
            tt,yy=((t+oldtheta)/2,(y+oldeta)/2) if exponential else (t,y)
            ne=self.eos.gas(tt,yy)[6];ab,em,at,et,ay,ey=self.eos.radiation(tt,yy)
            z=h*self.factor[:,None];chi=ab-em
            # Evaluate the small net photon exchange directly; subtracting two
            # almost equal occupations loses the much smaller neutral budget.
            if exponential:
                at,et,ay,ey=[v/2 for v in [at,et,ay,ey]];zz=z*chi;assert zz.min()>-50,'Unbounded stimulated amplification'
                psi,psip=self.phi1(zz);change=z*(em-chi*Jstar)*psi
                J=np.exp(-zz)*Jstar+z*em*psi
                jt=z*(et-(at-et)*Jstar)*psi+z*(em-chi*Jstar)*psip*z*(at-et)
                jy=z*(ey-(ay-ey)*Jstar)*psi+z*(em-chi*Jstar)*psip*z*(ay-ey)
            else:
                den=1+z*chi;assert den.min()>0,'Implicit amplification denominator'
                change=z*(em-chi*Jstar)/den;J=Jstar+change
                jt=z*(et-(at-et)*J)/den;jy=z*(ey-(ay-ey)*J)/den
            resE=(u-uold+(self.en*change).sum(1))/self.cv;resY=y-oldeta-(self.num*change).sum(1)
            if not derivatives:return np.array([resE,resY]),J,ab,em,ne
            a11=(ut+(self.en*jt).sum(1))/self.cv;a12=(uy+(self.en*jy).sum(1))/self.cv
            a21=-(self.num*jt).sum(1);a22=1-(self.num*jy).sum(1)
            return np.array([resE,resY]),(a11,a12,a21,a22)
        for iteration in range(12):
            res,der=equations(theta,eta);err=float(np.max(abs(res)))
            if err<2e-11:break
            a,b,c,d=der;det=a*d-b*c;assert np.all(det>0)
            dt=(-d*res[0]+b*res[1])/det;dy=(c*res[0]-a*res[1])/det
            for cut in range(12):
                tt=theta+dt*2.**(-cut);yy=eta+dy*2.**(-cut)
                lo,hi=getattr(self.eos,'temperature_bounds',(-.06,.06))
                if tt.min()<lo or tt.max()>hi or yy.min()<self.eos.ratios[0]-1 or yy.max()>self.eos.ratios[-1]-1:continue
                try:check=equations(tt,yy,False)[0]
                except AssertionError:continue
                if np.max(abs(check))<err:break
            else:raise AssertionError(('Nonlinear exchange support/line search',err,float(dt.min()),float(dt.max()),float(dy.min()),float(dy.max()),float(theta.min()),float(theta.max()),float(eta.min()),float(eta.max())))
            theta,eta=tt,yy
        else:raise AssertionError('Exchange Newton cap')
        res,J,ab,em,ne=equations(theta,eta,False)
        chi=h*self.factor[:,None]*(ab-em);sc=h*self.d['a']*C*ne*6.6524587321e-25
        if exponential:
            g=np.exp(-sc);g2=np.exp(-sc*(1-.5*self.meanP22))
            scattered=g[:,None,None]*I+(-np.expm1(-sc))[:,None,None]*Jstar[:,None,:]+((g2-g)/self.meanP22)[:,None,None]*self.P2[None,:,None]*Kstar[:,None,:]
            new=np.exp(-chi[:,None,:])*scattered+h*self.factor[:,None,None]*em[:,None,:]*self.phi1(chi)[0][:,None,:]
        else:
            Knew=Kstar/(1+chi+sc[:,None]*(1-.5*self.meanP22))
            new=(I+h*self.factor[:,None,None]*em[:,None,:]+sc[:,None,None]*(J[:,None,:]+.5*self.P2[None,:,None]*Knew[:,None,:]))/(1+chi[:,None,:]+sc[:,None,None])
        assert new.min()>=0,'Negative angular occupation'
        assert np.max(abs(self.mean(new)-J)/(abs(J)+1e-200))<1e-12
        return new,theta,eta,iteration,float(np.max(abs(res)))

    @staticmethod
    def phi1(z):
        small=abs(z)<1e-4;value=np.empty_like(z);der=np.empty_like(z);x=z[small]
        value[small]=1+x*(-.5+x*(1/6+x*(-1/24+x/120)))
        der[small]=-.5+x*(1/3+x*(-1/8+x/30))
        x=z[~small];value[~small]=-np.expm1(-x)/x;der[~small]=(np.exp(-x)*(x+1)-1)/(x*x)
        return value,der

    def exact_stream(self,I,h):
        rate=float(max(-self.A.diagonal()));P=sparse.eye(self.A.shape[0],format='csc')+self.A/rate
        coefficient=np.exp(-rate*h);value=I.copy();result=coefficient*value;port_integral=0.;accumulated=0.
        for k in range(1,64):
            accumulated+=self.port(value)[0]/rate
            value=(P@value.reshape(self.n*self.angles,self.m)).reshape(I.shape)+self.boundary/rate
            coefficient*=rate*h/k;result+=coefficient*value;port_integral+=coefficient*accumulated
            if coefficient*(k+1)<1e-18:break
        else:raise AssertionError('Positive streaming exponential series cap')
        assert result.min()>=0
        return result,port_integral

    def midpoint_exchange(self,I,theta,eta,h):
        # A source substep is refined only when the nonlinear midpoint has no
        # admissible accepted state, never to rescue a failed global comparison.
        for level in range(5):
            v=I.copy();tt=theta.copy();yy=eta.copy();maximum=0.;iterations=0
            try:
                for _ in range(2**level):
                    v,tt,yy,it,err=self.exchange(v,tt,yy,h/2**level,True);maximum=max(maximum,err);iterations=max(iterations,it)
            except AssertionError:continue
            self.maximum_local_substeps=max(self.maximum_local_substeps,2**level)
            return v,tt,yy,iterations,maximum
        raise AssertionError('Registered16 local midpoint-substep cap')

    def port(self,I):
        inner=np.where((self.mu>0)[:,None],self.incoming,I[0]);outer=np.where((self.mu>0)[:,None],I[-1],0.)
        Hin=np.einsum('qf,q->f',inner,self.w*self.mu);Hout=np.einsum('qf,q->f',outer,self.w*self.mu)
        fac=4*np.pi*C*self.area[[0,-1]]
        energy=(fac[0]*Hin-fac[1]*Hout)@(self.d['num']*self.d['Einf'])
        return float(energy),outer,fac[1]*Hout*self.d['num']*self.d['Einf']

    def run(self,steps,label,second_order=False):
        assert not (OUT/(label+'.npz')).exists();start=time.monotonic();h=prior.END/steps
        lu=None if second_order else splu(sparse.eye(self.n*self.angles,format='csc')-h*self.A);I=self.initial.copy();theta=np.zeros(self.n);eta=np.zeros(self.n);self.maximum_local_substeps=1
        E0=float(np.sum(I*self.photon_energy_weight));boundaryE=0.;balance=0.;local=0.;maximum_iterations=0;min_variance=1.;minimum=I.min();failed=None
        p0=self.eos.gas(theta,eta)[0];hist=[];ports=[];snapshots=[]
        def record(t):
            p,u,*_=self.eos.gas(theta,eta);J=self.mean(I);H=self.mean(I*self.mu[None,:,None]);K=self.mean(I*self.mu2[None,:,None])
            variance=np.divide(J*K-H*H,J*J,out=np.zeros_like(J),where=J>1e-140)
            trace=self.d['rho']*(u-self.u0)-3*(p-p0)
            hist.append(dict(t=t,theta=theta.copy(),eta=eta.copy(),J=J,H=H,K=K,trace=trace))
            port=self.port(I);ports.append(port[1]);return float(variance.min()),float((trace*self.volume*self.d['a']).sum())
        record(0.);snapshots.append(I.copy());snapshot_times=[0.]
        try:
            for step in range(steps):
                if second_order:
                    transported,delta=self.exact_stream(I,h/2);boundaryE+=delta
                    I,theta,eta,nit,err=self.midpoint_exchange(transported,theta,eta,h)
                    I,delta=self.exact_stream(I,h/2);boundaryE+=delta
                else:
                    rhs=(I+h*self.boundary).reshape(self.n*self.angles,self.m)
                    transported=lu.solve(rhs).reshape(I.shape);assert transported.min()>=0,'Negative streaming occupation'
                    boundaryE+=h*self.port(transported)[0]
                    I,theta,eta,nit,err=self.exchange(transported,theta,eta,h)
                maximum_iterations=max(maximum_iterations,nit);local=max(local,err)
                minimum=min(minimum,I.min());var,trace=record((step+1)*h);min_variance=min(min_variance,var)
                u=self.eos.gas(theta,eta)[1];energy=float(np.sum((I-self.initial)*self.photon_energy_weight)+(self.gas_weight*(u-self.u0)).sum())
                balance=max(balance,abs(energy-boundaryE))
                if (step+1)%max(1,steps//32)==0:snapshots.append(I.copy());snapshot_times.append((step+1)*h)
        except Exception as exc:failed=repr(exc)
        arrays={k:np.array([x[k] for x in hist]) for k in hist[0]};traces=arrays['trace']@(self.volume*self.d['a'])
        denominator=max(abs(boundaryE),max(abs(traces)),1.)
        np.savez_compressed(OUT/(label+'.npz'),**arrays,trace_energy=traces,snapshots=snapshots,snapshot_times=snapshot_times,
            outgoing_occupation=ports,mu=self.mu,angular_weights=self.w,r=self.d['r'],volume=self.volume,initial=self.initial,final=I)
        result=dict(classification='Counterexample candidate',passed=bool(failed is None and balance/denominator<1e-8 and local<1e-9 and minimum>=0 and min_variance>=-1e-12),
            failure=failed,angles=self.angles,steps=steps,completed_steps=len(hist)-1,seconds=time.monotonic()-start,
            energy_balance_relative=float(balance/denominator),maximum_local_exchange_residual=local,maximum_Newton_iterations=maximum_iterations,
            minimum_occupation=float(minimum),minimum_moment_variance=min_variance,initial_three_moment_relative=self.initial_moment_error,
            maximum_logT_change=float(np.max(abs(arrays['theta']))),maximum_relative_neutral_change=float(np.max(abs(arrays['eta']))),
            trace_endpoint_erg=float(traces[-1]),second_order_splitting=second_order,maximum_local_source_substeps=self.maximum_local_substeps,
            actual_positive_angular_transport=True,spatial_convergence=False,frequency_convergence=False,
            moving_atmosphere_connected=False,full_GR_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def pilot():
    signal.alarm(30);result=Model(8).run(32,'pilot-8-32');forecast=result['seconds']*(2+4+8)*3
    write(OUT/'measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_remaining_seconds=forecast,
        remaining_seconds=120-result['seconds'],eligible=bool(result['passed'] and forecast<120-result['seconds']),
        assumption='Linear step/angle scaling with3x margin;16angle LU/local-exchange cost unmeasured.'))
    signal.alarm(0)


def run(prefix='',second_order=False):
    b=json.loads((OUT/'thermal-measured-budget.json').read_text());assert b['eligible'];signal.alarm(int(b['remaining_seconds']));start=time.monotonic();rows=[]
    for q,n in [(8,64),(8,128),(16,128)]:
        row=Model(q).run(n,f'{prefix}angles-{q}-steps-{n}',second_order);rows.append(row)
        if not row['passed']:break
    result=dict(classification='Counterexample candidate',passed=False,paths=rows,seconds=time.monotonic()-start)
    if len(rows)==3 and all(x['passed'] for x in rows):
        a=np.load(OUT/(prefix+'angles-8-steps-64.npz'));b=np.load(OUT/(prefix+'angles-8-steps-128.npz'));c=np.load(OUT/(prefix+'angles-16-steps-128.npz'));d=np.load(prior.OUT/'bank-16-8.npz')
        errors={}
        def flux(v):return np.einsum('tqf,q->tf',v['outgoing_occupation'],v['angular_weights']*v['mu'])
        for key,x,y,stride in [('time',a,b,2),('angle',b,c,1)]:
            errors[key+'_trace']=float(np.max(abs(x['trace_energy']-y['trace_energy'][::stride]))/np.max(abs(y['trace_energy'])))
            xx,yy=flux(x),flux(y)[::stride];weight=d['num']*d['Einf'];errors[key+'_surface_spectrum']=float(np.max(abs(xx-yy)@weight)/np.max(yy@weight))
        result.update(passed=max(errors.values())<.02,comparisons=errors)
    write(OUT/(prefix+'result.json'),result);signal.alarm(0);print(json.dumps(result),flush=True)


def run_stable():
    write(OUT/'arithmetic-reassessment.json',dict(classification='Counterexample candidate',
        failure='The64step path stopped at an exchange residual2.0038e-11 just above the2e-11 internal Newton threshold. Temperature step was roundoff-sized. Subtracting J_new-J_star loses the tiny neutral-photon balance.',
        repair='Evaluate the algebraically identical deltaJ=z*(emission-net_absorption*J_star)/(1+z*net_absorption) directly. Preserve the failed5step prefix and source; rerun this0.21s incomplete path to retain an exact cumulative boundary ledger. No tolerance,physical support or grid change.',
        source_sha256=sha(__file__)))
    run('stable-')


def second_order():
    assert not (OUT/'second-order-result.json').exists();start=time.monotonic();signal.alarm(30)
    write(OUT/'time-method-reassessment.json',dict(classification='Counterexample candidate',
        failure='All positive angular paths complete, but first-order transport/reaction splitting fails the unchanged2percent time gate: trace11.218percent and surface spectrum8.464percent.8/16angle comparisons pass. Preserve those results.',
        repair='Same64/128 global paths. Symmetric exact-stream half steps and a conservative exponential-midpoint material/photon exchange. Keep actual angular positivity checks; allow at most16 local source substeps only to obtain an admissible nonlinear midpoint. No abundance/temperature clipping or new grid.',
        positivity='Streaming exponential is a positive Poisson series. At frozen midpoint coefficients, bound-free emission and the exact Thomson angular semigroup are nonnegative. Nonlinear material end states and every numerical angular density are also checked.',
        pilot='A new32step pilot measures nonlinear midpoint cost before the same three production paths. Total new run cap90s within the original120s transport allocation after observed previous costs below12s.',source_sha256=sha(__file__)))
    pilot=Model(8).run(32,'second-order-pilot',True);forecast=pilot['seconds']*14*3
    write(OUT/'second-order-budget.json',dict(pilot_seconds=pilot['seconds'],forecast_seconds=forecast,eligible=bool(pilot['passed'] and forecast<90-pilot['seconds'])))
    assert pilot['passed'] and forecast<90-pilot['seconds'];signal.alarm(max(1,int(90-(time.monotonic()-start))));run('second-order-',True)


def repair_support(recover=False):
    target=OUT/'neutral-support';start=time.monotonic();signal.alarm(20)
    if recover:
        assert target.exists() and not (target/'bank.npz').exists()
        write(target/'completion-plan.json',dict(classification='Counterexample candidate',
            root_failure='The shared molecular constraint updated fields only when CURRENT abundance exceeded1e-18. The frozen H2 target3.773e-12 fell to3.406e-20 and was never restored, leaving H+ excess3.331e-12 and failing the unchanged1e-12 full-population gate.',
            repair='Update a molecular field when either current OR target abundance exceeds the declared1e-18 trace threshold. Fix the shared Ions.constrain owner used by both fixed-inventory and reactive native paths. No population tolerance or iteration cap change.',
            resource='The failed48state construction exited before saving its partial bank/counter; exact calls unavailable, bounded by its400call cap. Preserve this operational failure. New completion cap400calls/20s after fixing the shared root cause; save partial state and call counter on failure. Prior existing96state bank remains reused.',
            original_constraint_sha256=sha(OUT/'constraint-producer-before-fix.py'),fixed_constraint_sha256=sha(prior.old.old.cold.old.__file__),source_sha256=sha(__file__)))
    else:
        assert not target.exists();target.mkdir()
        write(target/'plan.json',dict(classification='Counterexample candidate',
        failure='Positive angular pilot stops after5steps because H0 drops below0.25 of its initial value; no radiation positivity or conservation failure occurred before the EOS-support stop.',
        repair='Keep the old upper2.5 chemical plane and all physical radii,T domain,frequencies and gates. Evaluate one lower plane at1e-4 of the local initial neutral H;48 new native states. Actual finite-reaction depletion, not extra spatial refinement or clipped species.',
        new_native_calls=400,new_native_seconds=20,temperature_offsets=[-.06,0.,.06],new_neutral_ratios=[1e-4,2.5],
            stop='No automatic further support extension. If temperature or rates fail, preserve the exact physical support failure and reassess the model/energy source.',source_sha256=sha(__file__)))
    d=np.load(prior.OUT/'bank-16-8.npz');v=dict(np.load(old.OUT/'bank.npz'));native=old.old.Native(cap=400)
    raws=v['raw'].copy();rates=v['rates'].copy()
    done=np.zeros((len(d['r']),3),bool)
    try:
        for j in range(len(d['r'])):
            old.setup(native,d,j)
            for k,t in enumerate(v['offsets']):
                y=native.y0*1e-4;s=native.state(0.,np.log(d['T'][j])+t,y);chi,em=prior.coefficients(native,s,d['Einf']/d['a'][j])
                raws[j,0,k]=s['raw'];rates[j,0,k]=[(chi+em)/y,em/(1-y)];done[j,k]=True
    except Exception as exc:
        np.savez_compressed(target/'partial-bank.npz',raw=raws,rates=rates,done=done,native_calls=native.ion.calls)
        write(target/'completion-failure.json',dict(error=repr(exc),native_calls=native.ion.calls,seconds=time.monotonic()-start,cell=j,temperature_index=k));raise
    np.savez_compressed(target/'bank.npz',raw=raws,rates=rates,offsets=v['offsets'],ratios=np.array([1e-4,2.5]),native_calls=native.ion.calls)
    eos=SupportedTable();checks=[]
    for j in range(len(d['r'])):
        old.setup(native,d,j);theta=.03 if j%2 else -.02;eta=-.95 if j%2 else -.99
        s=native.state(0.,np.log(d['T'][j])+theta,native.y0*(1+eta));chi,em=prior.coefficients(native,s,d['Einf']/d['a'][j]);ab=chi+em
        tt=np.full(len(d['r']),theta);yy=np.full(len(d['r']),eta);p,u,*_=eos.gas(tt,yy);aa,ee,*_=eos.radiation(tt,yy);mask=(ab>0)&(em>0)
        checks.append(dict(cell=j,constitutive=float(max(abs(p[j]/s['raw'][1]-1),abs(u[j]/s['raw'][2]-1))),rate=float(max(np.max(abs(aa[j,mask]/ab[mask]-1)),np.max(abs(ee[j,mask]/em[mask]-1))))))
    result=dict(classification='Counterexample candidate',passed=bool(max(x['constitutive'] for x in checks)<.002 and max(x['rate'] for x in checks)<.002),
        checks=checks,native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=sha(__file__))
    write(target/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def complete_support():repair_support(True)


def repaired_pilot():
    assert json.loads((OUT/'neutral-support/result.json').read_text())['passed'];signal.alarm(30)
    result=Model(8).run(32,'pilot-repaired-8-32');forecast=result['seconds']*(2+4+8)*3
    spent=result['seconds']+json.loads((OUT/'pilot-8-32.json').read_text())['seconds']
    write(OUT/'repaired-measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_remaining_seconds=forecast,
        remaining_seconds=120-spent,eligible=bool(result['passed'] and forecast<120-spent),assumption='Same predeclared3x-margin scaling.'))
    signal.alarm(0)


def diagnose_support():
    signal.alarm(15);d=np.load(prior.OUT/'bank-16-8.npz');native=old.old.Native(cap=80)
    j=int(np.argmin(abs(np.log(d['rho'])+18.520191072923673)));lt=9.74137215461128;rows=[]
    for correction in [False,True]:
        old.setup(native,d,j);y=native.y0*1e-4
        if correction:
            ne=native.base['eos'][13]/1.66053906660e-24
            ratio=np.log1p(d['rho'][j]*native.nH*(native.y0-y)/ne)
            native.prefix['fields'][0,:316]+=native.charges*ratio
        start=native.ion.calls
        try:
            s=native.state(0.,lt,y);row=dict(passed=True,population_error=s['population_error'])
        except AssertionError as exc:
            target=native.target.copy();target[0,:2]=native.epsH*np.array([y,1-y]);last=native.ion.states[-1]['number_fractions'];err=abs(last-target);index=np.unravel_index(err.argmax(),err.shape)
            row=dict(passed=False,error=repr(exc),maximum_population_index=list(map(int,index)),absolute_error=float(err.max()),target=float(target[index]),actual=float(last[index]))
        row.update(electron_seed_correction=correction,calls=native.ion.calls-start);rows.append(row)
    write(OUT/'neutral-support/diagnostic.json',dict(classification='Counterexample candidate',cell=j,rows=rows,native_calls=native.ion.calls,source_sha256=sha(__file__)))
    print(json.dumps(rows),flush=True);signal.alarm(0)


def diagnose_molecules():
    d=np.load(prior.OUT/'bank-16-8.npz');native=old.old.Native(cap=25);j=11
    j=int(np.argmin(abs(np.log(d['rho'])+18.520191072923673)));old.setup(native,d,j)
    try:native.state(0.,9.74137215461128,native.y0*1e-4)
    except AssertionError as exc:error=repr(exc)
    else:error=None
    row=dict(classification='Counterexample candidate',cell=j,error=error,
        initial_molecules=native.base['molecular_H_fractions'].tolist(),last_molecules=native.ion.states[-1]['molecular_H_fractions'].tolist(),
        native_calls=native.ion.calls,source_sha256=sha(__file__))
    write(OUT/'neutral-support/molecular-diagnostic.json',row);print(json.dumps(row),flush=True)


def thermal_support():
    target=OUT/'thermal-support';assert not target.exists();target.mkdir();start=time.monotonic();signal.alarm(30)
    write(target/'plan.json',dict(classification='Counterexample candidate',
        failure='After the shared H2 constraint fix, the same angular path reaches1.717ms; cell10 at378.125km heats from18058.5K to19076.8K and the next nonlinear trial reaches the old+0.06 logT table edge. Positive occupations and energy exchange still hold. Preserve this failure.',
        decision='Complete the physical finite-temperature EOS input for the observed heating rather than clip T, stop the requested horizon, or loosen a numerical gate. The warmer Planck photons are actually streamed from deeper saved layers. Radial numerical diffusion and the partial optical channels are not yet bounded, so this is not evidence of true stellar heating.',
        repair='Reuse every existing native state; append exactly logT offsets0.12 and0.18 at the two existing H planes,64states. Use piecewise cubic interpolation of the slowly varying prefactors, retaining exact ideal-T and photon Boltzmann terms. Independent0.09/0.15 controls, unchanged0.002 gate.',
        budget=dict(native_calls=800,native_seconds=30),stop='No further thermal node/domain addition in this stage. The physical trace remains unaccepted until radial/frequency transport is resolved.',source_sha256=sha(__file__)))
    d=np.load(prior.OUT/'bank-16-8.npz');v=dict(np.load(OUT/'neutral-support/bank.npz'));native=old.old.Native(cap=800)
    raw=np.zeros((len(d['r']),2,5,21));rate=np.zeros((len(d['r']),2,5,2,len(d['Einf'])));raw[:,:,:3]=v['raw'];rate[:,:,:3]=v['rates'];done=np.zeros((len(d['r']),2,2),bool)
    offsets=np.r_[v['offsets'],.12,.18]
    try:
        for j in range(len(d['r'])):
            old.setup(native,d,j)
            for k,factor in enumerate(v['ratios']):
                for it,t in enumerate(offsets[3:]):
                    y=native.y0*factor;s=native.state(0.,np.log(d['T'][j])+t,y);chi,em=prior.coefficients(native,s,d['Einf']/d['a'][j])
                    raw[j,k,it+3]=s['raw'];rate[j,k,it+3]=[(chi+em)/y,em/(1-y)];done[j,k,it]=True
    except Exception as exc:
        np.savez_compressed(target/'partial-bank.npz',raw=raw,rates=rate,done=done,native_calls=native.ion.calls)
        write(target/'failure.json',dict(error=repr(exc),native_calls=native.ion.calls,seconds=time.monotonic()-start));raise
    np.savez_compressed(target/'bank.npz',raw=raw,rates=rate,offsets=offsets,ratios=v['ratios'],native_calls=native.ion.calls)
    eos=ThermalTable();checks=[];controls=[]
    for j in range(len(d['r'])):
        old.setup(native,d,j);theta=.09 if j%2 else .15;eta=-.95 if j%2 else .4
        s=native.state(0.,np.log(d['T'][j])+theta,native.y0*(1+eta));chi,em=prior.coefficients(native,s,d['Einf']/d['a'][j]);ab=chi+em
        tt=np.full(len(d['r']),theta);yy=np.full(len(d['r']),eta);p,u,*_=eos.gas(tt,yy);aa,ee,*_=eos.radiation(tt,yy);mask=(ab>0)&(em>0)
        checks.append(dict(cell=j,constitutive=float(max(abs(p[j]/s['raw'][1]-1),abs(u[j]/s['raw'][2]-1))),rate=float(max(np.max(abs(aa[j,mask]/ab[mask]-1)),np.max(abs(ee[j,mask]/em[mask]-1))))))
        controls.append(dict(theta=theta,eta=eta,raw=s['raw'],absorption=ab,emission=em))
    np.savez_compressed(target/'controls.npz',**{k:np.array([a[k] for a in controls]) for k in controls[0]},native_calls=native.ion.calls)
    result=dict(classification='Counterexample candidate',passed=bool(max(a['constitutive'] for a in checks)<.002 and max(a['rate'] for a in checks)<.002),checks=checks,
        native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=sha(__file__))
    write(target/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def thermal_pilot():
    assert json.loads((OUT/'thermal-support/result.json').read_text())['passed'];signal.alarm(30)
    result=Model(8).run(32,'pilot-thermal-8-32');forecast=result['seconds']*(2+4+8)*3
    spent=result['seconds']+sum(json.loads((OUT/name).read_text())['seconds'] for name in ['pilot-8-32.json','pilot-repaired-8-32.json'])
    write(OUT/'thermal-measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_remaining_seconds=forecast,
        remaining_seconds=120-spent,eligible=bool(result['passed'] and forecast<120-spent),assumption='Same predeclared3x-margin scaling.'))
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
