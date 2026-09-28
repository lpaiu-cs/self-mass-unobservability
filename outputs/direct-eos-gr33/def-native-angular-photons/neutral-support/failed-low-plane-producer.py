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
import def_native_causal_nonlinear as old

prior=old.prior
C=prior.C
OUT=prior.OUT.parent/'def-native-angular-photons'
write=old.write
sha=prior.sha


class SupportedTable(old.Table):
    __init__=FunctionType(old.Table.__init__.__code__,dict(old.Table.__init__.__globals__,OUT=OUT/'neutral-support'))


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
        self.eos=SupportedTable() if (OUT/'neutral-support/bank.npz').exists() else old.Table()
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

    def exchange(self,I,theta,eta,h):
        Jstar=self.mean(I);Kstar=self.mean(I*self.P2[None,:,None]);oldeta=eta.copy();uold=self.eos.gas(theta,eta)[1]
        def equations(t,y,derivatives=True):
            p,u,ut,uy,_,_,ne,_,_=self.eos.gas(t,y);ab,em,at,et,ay,ey=self.eos.radiation(t,y)
            z=h*self.factor[:,None];den=1+z*(ab-em);assert den.min()>0,'Implicit amplification denominator'
            J=(Jstar+z*em)/den;change=J-Jstar
            resE=(u-uold+(self.en*change).sum(1))/self.cv;resY=y-oldeta-(self.num*change).sum(1)
            if not derivatives:return np.array([resE,resY]),J,ab,em,ne
            jt=z*(et-(at-et)*J)/den;jy=z*(ey-(ay-ey)*J)/den
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
                if np.max(abs(tt))>.06 or yy.min()<self.eos.ratios[0]-1 or yy.max()>self.eos.ratios[-1]-1:continue
                try:check=equations(tt,yy,False)[0]
                except AssertionError:continue
                if np.max(abs(check))<err:break
            else:raise AssertionError(('Nonlinear exchange support/line search',err,float(dt.min()),float(dt.max()),float(dy.min()),float(dy.max()),float(theta.min()),float(theta.max()),float(eta.min()),float(eta.max())))
            theta,eta=tt,yy
        else:raise AssertionError('Exchange Newton cap')
        res,J,ab,em,ne=equations(theta,eta,False)
        chi=h*self.factor[:,None]*(ab-em);sc=h*self.d['a']*C*ne*6.6524587321e-25
        Knew=Kstar/(1+chi+sc[:,None]*(1-.5*self.meanP22))
        new=(I+h*self.factor[:,None,None]*em[:,None,:]+sc[:,None,None]*(J[:,None,:]+.5*self.P2[None,:,None]*Knew[:,None,:]))/(1+chi[:,None,:]+sc[:,None,None])
        assert new.min()>=0,'Negative angular occupation'
        assert np.max(abs(self.mean(new)-J)/(abs(J)+1e-200))<1e-12
        return new,theta,eta,iteration,float(np.max(abs(res)))

    def port(self,I):
        inner=np.where((self.mu>0)[:,None],self.incoming,I[0]);outer=np.where((self.mu>0)[:,None],I[-1],0.)
        Hin=np.einsum('qf,q->f',inner,self.w*self.mu);Hout=np.einsum('qf,q->f',outer,self.w*self.mu)
        fac=4*np.pi*C*self.area[[0,-1]]
        energy=(fac[0]*Hin-fac[1]*Hout)@(self.d['num']*self.d['Einf'])
        return float(energy),outer,fac[1]*Hout*self.d['num']*self.d['Einf']

    def run(self,steps,label):
        assert not (OUT/(label+'.npz')).exists();start=time.monotonic();h=prior.END/steps
        lu=splu(sparse.eye(self.n*self.angles,format='csc')-h*self.A);I=self.initial.copy();theta=np.zeros(self.n);eta=np.zeros(self.n)
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
                rhs=(I+h*self.boundary).reshape(self.n*self.angles,self.m)
                transported=lu.solve(rhs).reshape(I.shape);assert transported.min()>=0,'Negative streaming occupation'
                boundaryE+=h*self.port(transported)[0]
                I,theta,eta,nit,err=self.exchange(transported,theta,eta,h);maximum_iterations=max(maximum_iterations,nit);local=max(local,err)
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
            trace_endpoint_erg=float(traces[-1]),actual_positive_angular_transport=True,spatial_convergence=False,frequency_convergence=False,
            moving_atmosphere_connected=False,full_GR_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def pilot():
    signal.alarm(30);result=Model(8).run(32,'pilot-8-32');forecast=result['seconds']*(2+4+8)*3
    write(OUT/'measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_remaining_seconds=forecast,
        remaining_seconds=120-result['seconds'],eligible=bool(result['passed'] and forecast<120-result['seconds']),
        assumption='Linear step/angle scaling with3x margin;16angle LU/local-exchange cost unmeasured.'))
    signal.alarm(0)


def run():
    b=json.loads((OUT/'repaired-measured-budget.json').read_text());assert b['eligible'];signal.alarm(int(b['remaining_seconds']));start=time.monotonic();rows=[]
    for q,n in [(8,64),(8,128),(16,128)]:
        row=Model(q).run(n,f'angles-{q}-steps-{n}');rows.append(row)
        if not row['passed']:break
    result=dict(classification='Counterexample candidate',passed=False,paths=rows,seconds=time.monotonic()-start)
    if len(rows)==3 and all(x['passed'] for x in rows):
        a=np.load(OUT/'angles-8-steps-64.npz');b=np.load(OUT/'angles-8-steps-128.npz');c=np.load(OUT/'angles-16-steps-128.npz');d=np.load(prior.OUT/'bank-16-8.npz')
        errors={}
        def flux(v):return np.einsum('tqf,q->tf',v['outgoing_occupation'],v['angular_weights']*v['mu'])
        for key,x,y,stride in [('time',a,b,2),('angle',b,c,1)]:
            errors[key+'_trace']=float(np.max(abs(x['trace_energy']-y['trace_energy'][::stride]))/np.max(abs(y['trace_energy'])))
            xx,yy=flux(x),flux(y)[::stride];weight=d['num']*d['Einf'];errors[key+'_surface_spectrum']=float(np.max(abs(xx-yy)@weight)/np.max(yy@weight))
        result.update(passed=max(errors.values())<.02,comparisons=errors)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def repair_support():
    target=OUT/'neutral-support';assert not target.exists();target.mkdir();start=time.monotonic();signal.alarm(20)
    write(target/'plan.json',dict(classification='Counterexample candidate',
        failure='Positive angular pilot stops after5steps because H0 drops below0.25 of its initial value; no radiation positivity or conservation failure occurred before the EOS-support stop.',
        repair='Keep the old upper2.5 chemical plane and all physical radii,T domain,frequencies and gates. Evaluate one lower plane at1e-4 of the local initial neutral H;48 new native states. Actual finite-reaction depletion, not extra spatial refinement or clipped species.',
        new_native_calls=400,new_native_seconds=20,temperature_offsets=[-.06,0.,.06],new_neutral_ratios=[1e-4,2.5],
        stop='No automatic further support extension. If temperature or rates fail, preserve the exact physical support failure and reassess the model/energy source.',source_sha256=sha(__file__)))
    d=np.load(prior.OUT/'bank-16-8.npz');v=dict(np.load(old.OUT/'bank.npz'));native=old.old.Native(cap=400)
    raws=v['raw'].copy();rates=v['rates'].copy()
    for j in range(len(d['r'])):
        old.setup(native,d,j)
        for k,t in enumerate(v['offsets']):
            y=native.y0*1e-4;s=native.state(0.,np.log(d['T'][j])+t,y);chi,em=prior.coefficients(native,s,d['Einf']/d['a'][j])
            raws[j,0,k]=s['raw'];rates[j,0,k]=[(chi+em)/y,em/(1-y)]
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


def repaired_pilot():
    assert json.loads((OUT/'neutral-support/result.json').read_text())['passed'];signal.alarm(30)
    result=Model(8).run(32,'pilot-repaired-8-32');forecast=result['seconds']*(2+4+8)*3
    spent=result['seconds']+json.loads((OUT/'pilot-8-32.json').read_text())['seconds']
    write(OUT/'repaired-measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_remaining_seconds=forecast,
        remaining_seconds=120-spent,eligible=bool(result['passed'] and forecast<120-spent),assumption='Same predeclared3x-margin scaling.'))
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
