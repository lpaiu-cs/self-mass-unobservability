"""Counterexample candidate: unsplit native H/heat and angular photon time solve."""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy.linalg.lapack import dgbtrf,dgbtrs
import def_native_angular_photons as prior

OUT=prior.OUT.parent/'def-native-unsplit-photons'
write=prior.write
sha=prior.sha
GAMMA=1-1/np.sqrt(2)


def timeout(*_):raise TimeoutError('Registered wall-time cap')


class Model(prior.Model):
    def __init__(self,angles=8):
        super().__init__(angles)
        n,q=self.n,self.angles;self.N=n*q;self.ids=np.arange(self.N).reshape(n,q)
        self.rows,self.cols=np.indices((q,q));self.S=np.ones((q,1))*self.w[None,:]+.5*self.P2[:,None]*(self.P2*self.w)[None,:]
        self.R=np.eye(q)-self.S
        self.scfactor=self.d['a']*prior.C*6.6524587321e-25
        self.en_scaled=self.en*self.scale;self.num_scaled=self.num*self.scale
        self.stream=self.A.tocoo()
        self.worst_linear_residual=0.

    def evaluate(self,theta,eta,xref,uref,yref,h,jacobian=True):
        n,q=self.n,self.angles;N=self.N
        lo,hi=self.eos.temperature_bounds
        assert np.all((theta>=lo)&(theta<=hi)) and np.all((eta>=self.eos.ratios[0]-1)&(eta<=self.eos.ratios[-1]-1)),'Native support'
        p,u,ut,uy,pt,py,ne,net,ney=self.eos.gas(theta,eta)
        ab,em,at,et,ay,ey=self.eos.radiation(theta,eta);chi=ab-em;ct=at-et;cy=ay-ey
        sc=self.scfactor*ne;sct=self.scfactor*net;scy=self.scfactor*ney
        base=np.zeros((3*q+1,N));A=self.stream
        np.add.at(base,(2*q+A.row-A.col,A.col),-h*A.data)
        base[2*q]+=1.
        for i in range(n):
            base[2*q+self.rows-self.cols,i*q+self.cols]+=h*sc[i]*self.R
        x=np.empty_like(xref);collision=np.empty((n,self.m))
        jet=np.zeros((n,2*n));jey=np.zeros_like(jet)
        rhs=(xref+h*self.boundary/self.scale).reshape(N,self.m)
        for f in range(self.m):
            band=base.copy(order='F');band[2*q]+=np.repeat(h*self.factor*chi[:,f],q)
            lu,piv,info=dgbtrf(band,q,q,overwrite_ab=1);assert info==0,('Photon factorization',info)
            rr=rhs[:,f]+np.repeat(h*self.factor*em[:,f]/self.scale[f],q)
            xx,info=dgbtrs(lu,q,q,rr[:,None],piv);assert info==0
            xx=xx[:,0].reshape(n,q);x[:,:,f]=xx;J=xx@self.w
            collision[:,f]=self.factor*(em[:,f]/self.scale[f]-chi[:,f]*J)
            if not jacobian:continue
            scattering=xx@self.S.T-xx
            bt=self.factor[:,None]*(et[:,f,None]/self.scale[f]-ct[:,f,None]*xx)+sct[:,None]*scattering
            by=self.factor[:,None]*(ey[:,f,None]/self.scale[f]-cy[:,f,None]*xx)+scy[:,None]*scattering
            right=np.zeros((N,2*n));right[self.ids,np.arange(n)[:,None]]=h*bt;right[self.ids,n+np.arange(n)[:,None]]=h*by
            dx,info=dgbtrs(lu,q,q,right,piv);assert info==0
            dJ=np.einsum('iqp,q->ip',dx.reshape(n,q,2*n),self.w)
            dc=-self.factor[:,None]*chi[:,f,None]*dJ
            dc[np.arange(n),np.arange(n)]+=self.factor*(et[:,f]/self.scale[f]-ct[:,f]*J)
            dc[np.arange(n),n+np.arange(n)]+=self.factor*(ey[:,f]/self.scale[f]-cy[:,f]*J)
            jet+=h*self.en_scaled[:,f,None]*dc/self.cv[:,None];jey-=h*self.num_scaled[:,f,None]*dc
        res=np.r_[(u-uref+h*(self.en_scaled*collision).sum(1))/self.cv,eta-yref-h*(self.num_scaled*collision).sum(1)]
        if not jacobian:return res,x,u
        jet[np.arange(n),np.arange(n)]+=ut/self.cv;jet[np.arange(n),n+np.arange(n)]+=uy/self.cv
        jey[np.arange(n),n+np.arange(n)]+=1
        return res,x,u,np.vstack([jet,jey])

    def implicit(self,xref,uref,yref,h,theta,eta):
        theta=theta.copy();eta=eta.copy();n=self.n
        for iteration in range(12):
            res,x,u,jac=self.evaluate(theta,eta,xref,uref,yref,h);err=float(np.max(abs(res)))
            if err<2e-11:break
            delta=np.linalg.solve(jac,-res)
            for cut in range(12):
                tt=theta+delta[:n]*2.**(-cut);yy=eta+delta[n:]*2.**(-cut)
                try:check=self.evaluate(tt,yy,xref,uref,yref,h,False)[0]
                except AssertionError:continue
                if np.max(abs(check))<err:break
            else:raise AssertionError(('Unsplit Newton line search',err,float(theta.min()),float(theta.max()),float(eta.min()),float(eta.max()),float(np.max(abs(delta)))))
            theta,eta=tt,yy
        else:raise AssertionError(('Unsplit Newton cap',err))
        assert x.min()>=0,('Negative unsplit photon',float(x.min()))
        return x,u,theta,eta,iteration,err

    def run(self,steps,label,stop_after=None):
        assert not (OUT/(label+'.npz')).exists();start=time.monotonic();h=prior.prior.END/steps
        x=self.initial/self.scale;theta=np.zeros(self.n);eta=np.zeros(self.n);u=self.u0.copy()
        boundary=0.;balance=0.;local=0.;iterations=0;minimum=0.;failed=None;hist=[];ports=[];snapshots=[];snapshot_times=[]
        p0=self.eos.gas(theta,eta)[0];count=steps if stop_after is None else stop_after
        def record(t):
            I=x*self.scale;p=self.eos.gas(theta,eta)[0];J=self.mean(I);H=self.mean(I*self.mu[None,:,None]);K=self.mean(I*self.mu2[None,:,None])
            trace=self.d['rho']*(u-self.u0)-3*(p-p0)
            hist.append(dict(t=t,theta=theta.copy(),eta=eta.copy(),J=J,H=H,K=K,trace=trace))
            ports.append(self.port(I)[1]);return float((trace*self.volume*self.d['a']).sum())
        record(0.);snapshots.append(self.initial.copy());snapshot_times.append(0.)
        try:
            for step in range(count):
                x1,u1,t1,y1,it1,e1=self.implicit(x,u,eta,GAMMA*h,theta,eta)
                ratio=(1-GAMMA)/GAMMA
                x2,u2,t2,y2,it2,e2=self.implicit(x+ratio*(x1-x),u+ratio*(u1-u),eta+ratio*(y1-eta),GAMMA*h,t1,y1)
                boundary+=h*((1-GAMMA)*self.port(x1*self.scale)[0]+GAMMA*self.port(x2*self.scale)[0])
                x,u,theta,eta=x2,u2,t2,y2;iterations=max(iterations,it1,it2);local=max(local,e1,e2);minimum=min(minimum,float(x.min()))
                record((step+1)*h)
                energy=float(np.sum((x*self.scale-self.initial)*self.photon_energy_weight)+(self.gas_weight*(u-self.u0)).sum())
                balance=max(balance,abs(energy-boundary))
                if (step+1)%max(1,steps//32)==0 or step+1==count:snapshots.append(x*self.scale);snapshot_times.append((step+1)*h)
        except Exception as exc:failed=repr(exc)
        arrays={k:np.array([v[k] for v in hist]) for k in hist[0]};trace=arrays['trace']@(self.volume*self.d['a'])
        denominator=max(abs(boundary),float(max(abs(trace))),1.)
        np.savez_compressed(OUT/(label+'.npz'),**arrays,trace_energy=trace,snapshots=snapshots,snapshot_times=snapshot_times,
            outgoing_occupation=ports,mu=self.mu,angular_weights=self.w,r=self.d['r'],volume=self.volume,initial=self.initial,final=x*self.scale)
        result=dict(classification='Counterexample candidate',passed=bool(failed is None and balance/denominator<1e-8 and local<1e-9 and minimum>=0),
            failure=failed,angles=self.angles,steps=steps,planned_executed_steps=count,completed_steps=len(hist)-1,seconds=time.monotonic()-start,
            energy_balance_relative=float(balance/denominator),maximum_nonlinear_residual=local,maximum_Newton_iterations=iterations,
            minimum_scaled_photon=minimum,maximum_logT_change=float(np.max(abs(arrays['theta']))),maximum_relative_neutral_change=float(np.max(abs(arrays['eta']))),
            trace_endpoint_erg=float(trace[-1]),unsplit_material_photons=True,spatial_convergence=False,frequency_convergence=False,
            moving_atmosphere_connected=False,full_GR_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='fa2c644d0',
        claim='Replace operator splitting by an actual simultaneous nonlinear material and angular-photon time solve; adjudicate the remaining4.889percent full-trace time mismatch.',
        decision='Accept this fixed-background time response only if the existing2percent trace and surface-spectrum gates pass; then resolve radial/frequency error and connect the same photon port to the moving atmosphere.',
        reuse='All existing16 radii,152 Killing frequencies,8/16 angular bins,64/128 time paths,3.434ms horizon,initial angular occupations,thermal native table and boundaries unchanged. No new EOS table or larger grid.',
        equations='At fixed material trial, solve the banded photon equations exactly at each frequency; use their derivatives to form the32variable gas Schur system. Apply a conservative two-stage L-stable SDIRK method with gamma=1-1/sqrt(2). Both stages couple streaming,Thomson scattering,H absorption/emission,neutral inventory and total gas energy.',
        positivity='Higher-order stiff time integration is not unconditionally positive. Check every accepted stage and endpoint without clipping. A negative photon or unsupported native state is a failure, not permission to refine automatically.',
        gates=dict(time_trace=.02,time_surface_spectrum=.02,angle_trace=.02,angle_surface_spectrum=.02,energy=1e-8,nonlinear_residual=1e-9,positive_photons=True,jacobian_directional=1e-5),
        budget=dict(pilot_steps_each=2,pilot_paths=[[8,64],[16,128]],pilot_seconds=30,production_seconds=120,CPU_threads=1,memory_GB=2,native_audit_calls=24,native_audit_seconds=10),
        forecast='Measure both angular path costs. Production prediction uses measured seconds/step; upper bound assumes2x measured cost. Later nonlinear stiffness remains unmeasured.',
        stop='Preserve failures. No automatic256step path,extra radius/frequency,EOS support expansion or weaker gate. A failed prediction or120s cap stops production for reassessment.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(prior.old.__file__),prior.OUT/'thermal-support/bank.npz',prior.prior.OUT/'bank-16-8.npz',prior.OUT/'second-order-result.json']}))
    import sympy as s
    z=s.symbols('z');g=1-1/s.sqrt(2);R=(1+(1-2*g)*z)/(1-g*z)**2
    assert s.simplify(R.subs(z,0))==1 and s.simplify(s.diff(R,z).subs(z,0))==1 and s.simplify(s.diff(R,z,2).subs(z,0))==1
    assert s.limit(R,z,s.oo)==0 and s.simplify(2*g*g-(1-2*g)**2)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Two-stage stability function and second-order consistency for smooth autonomous ODEs. Does not prove stiff order,positivity or physical convergence.',
        stability_function=str(R),A_stability='On z=i*y, denominator modulus squared minus numerator modulus squared equals gamma^4*y^4>=0; the only pole is in the right half-plane.',
        conservation='Apply each stage to gas total internal energy and the photon occupation with identical source and boundary weights; absorption/emission cancels in their sum.'))


def check():
    assert not (OUT/'operator-check.json').exists();signal.signal(signal.SIGALRM,timeout);signal.alarm(15);start=time.monotonic()
    m=Model(8);n=m.n;theta=np.linspace(.003,.012,n);eta=np.linspace(-.2,.1,n);h=GAMMA*prior.prior.END/64
    xref=m.initial/m.scale;res,x,u,jac=m.evaluate(theta,eta,xref,m.u0,np.zeros(n),h)
    ab,em,*_=m.eos.radiation(theta,eta);ne=m.eos.gas(theta,eta)[6];chi=ab-em
    collision=m.factor[:,None,None]*(em[:,None,:]/m.scale-chi[:,None,:]*x)
    scattering=m.scfactor[:,None,None]*ne[:,None,None]*(np.einsum('qk,ikf->iqf',m.S,x)-x)
    dx=(m.A@x.reshape(m.N,m.m)).reshape(x.shape)+m.boundary/m.scale+collision+scattering
    photon=float(np.max(abs(x-xref-h*dx))/max(1.,float(np.max(abs(x)))))
    energy=float(np.sum(dx*m.scale*m.photon_energy_weight)-(m.gas_weight*(m.en_scaled*m.mean(collision)).sum(1)).sum())
    boundary=m.port(x*m.scale)[0];conservation=abs(energy-boundary)/max(abs(boundary),1.)
    rng=np.random.default_rng(111);directions=[]
    for _ in range(2):
        v=rng.normal(size=2*n);v/=np.max(abs(v));eps=1e-6
        plus=m.evaluate(theta+eps*v[:n],eta+eps*v[n:],xref,m.u0,np.zeros(n),h,False)[0]
        minus=m.evaluate(theta-eps*v[:n],eta-eps*v[n:],xref,m.u0,np.zeros(n),h,False)[0]
        fd=(plus-minus)/(2*eps);exact=jac@v;directions.append(float(np.max(abs(fd-exact))/max(1.,float(np.max(abs(exact))))))
    result=dict(classification='Counterexample candidate',passed=bool(photon<1e-11 and conservation<1e-9 and max(directions)<1e-5),
        photon_equation_relative=photon,instantaneous_energy_exchange_relative=conservation,Schur_directional_relative=directions,
        seconds=time.monotonic()-start,source_sha256=sha(__file__))
    write(OUT/'operator-check.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def pilot():
    assert json.loads((OUT/'operator-check.json').read_text())['passed'];signal.signal(signal.SIGALRM,timeout);signal.alarm(30);start=time.monotonic();rows=[]
    for q,n in [(8,64),(16,128)]:
        row=Model(q).run(n,f'pilot-{q}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)<2 or not all(v['passed'] for v in rows) else rows[0]['seconds']/2*192+rows[1]['seconds']/2*128
    write(OUT/'measured-budget.json',dict(rows=rows,forecast_seconds=forecast,upper_seconds=None if forecast is None else forecast*2,
        eligible=bool(forecast is not None and forecast*2<120),seconds=time.monotonic()-start))
    signal.alarm(0)


def run():
    assert json.loads((OUT/'measured-budget.json').read_text())['eligible'];signal.signal(signal.SIGALRM,timeout);signal.alarm(120);start=time.monotonic();rows=[]
    for q,n in [(8,64),(8,128),(16,128)]:
        row=Model(q).run(n,f'angles-{q}-steps-{n}');rows.append(row)
        if not row['passed']:break
    result=dict(classification='Counterexample candidate',passed=False,paths=rows,seconds=time.monotonic()-start)
    if len(rows)==3 and all(x['passed'] for x in rows):
        a=np.load(OUT/'angles-8-steps-64.npz');b=np.load(OUT/'angles-8-steps-128.npz');c=np.load(OUT/'angles-16-steps-128.npz');d=np.load(prior.prior.OUT/'bank-16-8.npz');errors={}
        def flux(v):return np.einsum('tqf,q->tf',v['outgoing_occupation'],v['angular_weights']*v['mu'])
        for key,x,y,stride in [('time',a,b,2),('angle',b,c,1)]:
            errors[key+'_trace']=float(np.max(abs(x['trace_energy']-y['trace_energy'][::stride]))/np.max(abs(y['trace_energy'])))
            xx,yy=flux(x),flux(y)[::stride];weight=d['num']*d['Einf'];errors[key+'_surface_spectrum']=float(np.max(abs(xx-yy)@weight)/np.max(yy@weight))
        result.update(passed=max(errors.values())<.02,comparisons=errors)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
