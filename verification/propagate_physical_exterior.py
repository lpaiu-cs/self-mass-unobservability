"""Actual incident-metric photon worldlines and stress source on the saved EOS.

Counterexample candidate. Preserve separate reference-emission and geometric
responses. This finite-radius Vlasov source is not a solved scalar/GR boundary.
"""
from pathlib import Path
import fcntl, json, os, resource, sys, time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from numpy.polynomial import legendre as leg
import read_returned_exterior as saved
import def_native_incident_drive as incident
import def_native_dynamic_lapse as lapse

OUT=Path('native-dynamic-photon261-work')
ACTUAL=Path('native-short-return259-work')
CHARGE=Path('native-short-return-charge259-work/full/charge')
read,write,sha,LD=saved.read,saved.write,saved.sha,np.longdouble
C,G=incident.C,incident.G
SETTINGS={'fine':(8,8,8),'angular':(4,8,8),'temporal':(8,4,8),'geometry':(8,8,4)}


def initialize():
    out=OUT/'initialization'
    for p in list((CHARGE/'sweep-0').rglob('*.npz'))+[CHARGE/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=out/p.relative_to(CHARGE);dst.parent.mkdir(parents=True,exist_ok=True)
        if not dst.exists():os.link(p,dst)
        else:assert sha(p)==sha(dst)
    for n in ['sweep-1/photons','sweep-1/material','gr']:(out/n).mkdir(parents=True,exist_ok=True)
    with Path('.native-reader-initialization.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        saved.charge.bind(saved.charge.base.endpoint.initialize,OUT=out)()


def symbolic():
    mu,h,dr,q,nu,ul,lam,nr,lr,r,e=sp.symbols('mu h dr q nu ul lam nr lr r e',nonzero=True)
    dmu=(1-mu**2)/mu*(h+q*dr);dH=nu+ul+h
    angular=dH+dr/r-nu-nr*dr-mu*dmu/(1-mu**2)
    assert sp.simplify((angular-ul).subs(q,1/r-nr))==0
    # The Einstein-frame radial energy measure of a packet is H*sqrt(b)/N.
    s=sp.symbols('s');N,B=sp.symbols('N B',positive=True)
    measure=e*sp.exp(s*dH)*sp.sqrt(B)*sp.exp(-s*lam)/(N*sp.exp(s*nu))
    assert sp.simplify(sp.diff(measure,s).subs(s,0)/(e*sp.sqrt(B)/N)-(dH-nu-lam))==0
    rho,pr,pt=sp.symbols('rho pr pt')
    assert sp.expand(rho-rho*mu**2-2*rho*(1-mu**2)/2)==0
    return dict(classification='Proven',passed=True,
        fixed_radius='h=delta_lnH-delta_nu-delta_lnL; dh/dr=-delta_nu_prime-mu*delta_lambda_t/(c*N*sqrt(b)); delta_t_prime=[-zeta-(1-mu^2)*h/mu^2]/(c*N*sqrt(b)*mu).',
        fixed_time='delta_r=-c*N*sqrt(b)*mu*delta_t; delta_mu=(1-mu^2)/mu*[h+(1/r-nu0_prime)*delta_r]; delta_lnH=delta_nu+u_launch+h; delta_lnL=u_launch.',
        source='For any smooth radial test psi, delta integral(4pi*r^2*rho_E*psi dr)=sum E0*sqrt(b)/N*[psi*(delta_lnH-delta_nu-delta_lambda)+delta_r*(psi_prime-(lambda0_prime+nu0_prime)*psi)]. Radial flux and pressure additionally vary mu and mu^2.',
        scope='Linearized null kinematics, angular invariant and packet stress measure. No scalar/mass closure, boundary replay or final charge theorem.')


class Metric:
    def __init__(self,order):
        self.d=d=incident.Driver(order);self.order=order
        self.x,self.w=leg.leggauss(order);self.edges=d.ex
        self.nodes=(self.edges[:-1,None]+np.diff(self.edges)[:,None]*(self.x+1)/2)
        self.weights=np.diff(self.edges)[:,None]*self.w/2
        self.radius=d.inverse(self.nodes.ravel()).reshape(self.nodes.shape)
        _,N,b,a,Phi=d.bg.metric(self.radius.ravel()/d.model.m.R)
        self.coeff=(2*(Phi/d.model.m.R)*a/(self.radius.ravel()*b)).reshape(self.nodes.shape)
        self.scale=abs(incident.ETA*d.r0*float(d.z0['Phi'][0]));assert self.scale>0

    def wave(self,t,x,part='all'):
        """Same saved U values, with exact derivatives of its bilinear Born part."""
        t,x=np.broadcast_arrays(t,x);d=self.d
        s=(t+x/C)/d.D;u=incident.ETA*d.r0*incident.pulse(s)
        ut=incident.ETA*d.r0/d.D*incident.pulse(s,1);ux=ut/C
        if part=='primary':return u,ut,ux
        j=np.clip(np.searchsorted(d.times,t,side='right')-1,0,len(d.times)-2)
        k=np.clip(np.searchsorted(d.tx,x,side='right')-1,0,len(d.tx)-2)
        dt=d.times[j+1]-d.times[j];dx=d.tx[k+1]-d.tx[k]
        f=(t-d.times[j])/dt;g=(x-d.tx[k])/dx;v=d.born['U']
        lo=(1-g)*v[j,k]+g*v[j,k+1];hi=(1-g)*v[j+1,k]+g*v[j+1,k+1]
        ub=(1-f)*lo+f*hi;utb=(hi-lo)/dt
        uxb=((1-f)*(v[j,k+1]-v[j,k])+f*(v[j+1,k+1]-v[j+1,k]))/dx
        return (ub,utb,uxb) if part=='born' else (u+ub,ut+utb,ux+uxb)

    def partial(self,t,lo,hi,part):
        width=np.maximum(hi-lo,0);x=lo[:,None]+width[:,None]*(self.x+1)/2
        r=self.d.inverse(x.ravel());_,N,b,a,Phi=self.d.bg.metric(r/self.d.model.m.R)
        coeff=(2*(Phi/self.d.model.m.R)*a/(r*b)).reshape(x.shape)
        u,ut,_=self.wave(t[:,None],x,part)
        return np.array([np.sum(width[:,None]*self.w/2*coeff*v,axis=1) for v in [u,ut]])

    def integral(self,t,x,part):
        upper=np.maximum(C*(self.d.D-t),0) if part=='primary' else np.full_like(t,self.edges[-1])
        upper=np.minimum(upper,self.edges[-1]);lo=np.minimum(x,upper)
        j=np.clip(np.searchsorted(self.edges,lo,side='right')-1,0,len(self.edges)-2)
        k=np.clip(np.searchsorted(self.edges,upper,side='left')-1,0,len(self.edges)-2)
        u,ut,_=self.wave(t[:,None,None],self.nodes[None,:,:],part)
        included=(self.edges[:-1][None,:]>=lo[:,None])&(self.edges[1:][None,:]<=upper[:,None])
        value=np.array([np.sum(v*self.coeff*self.weights*included[:,:,None],axis=(1,2)) for v in [u,ut]])
        # Full panels at an exact lower knot are already included.
        full_left=(lo==self.edges[j])&(upper>=self.edges[j+1])
        left=np.where(full_left,self.edges[j+1],lo)
        value+=self.partial(t,left,np.minimum(self.edges[j+1],upper),part)
        right=np.maximum(self.edges[k],lo)
        right=np.where((k==j)|(upper==self.edges[k+1]),upper,right)
        value+=self.partial(t,right,upper,part)
        return value

    def at(self,t,r):
        t,r=np.broadcast_arrays(t,r);shape=t.shape;t=t.ravel();r=r.ravel();d=self.d
        x=d.optical(r);x=np.clip(x,0,self.edges[-1])
        mass,N,b,a,Phi=d.bg.metric(r/d.model.m.R);mass*=d.model.m.R;Phi/=d.model.m.R
        u,ut,ux=self.wave(t,x);integ=self.integral(t,x,'primary')+self.integral(t,x,'born')
        lam=Phi*u;nu=lam-integ[0];lt=Phi*ut;nt=lt-integ[1]
        nup=Phi*(ux/a+u*(1/b-1)/r)
        nr=mass/(r*r*b)+r*Phi**2/2;lr=r*Phi**2/2-mass/(r*r*b)
        return {k:v.reshape(shape) for k,v in dict(nu=nu,lam=lam,zeta=-integ[0],lt=lt,nt=nt,nup=nup,N=N,b=b,a=a,nr=nr,lr=lr).items()}


class Photons:
    def __init__(self,angular,geometry):
        self.metric=Metric(geometry);self.d=d=self.metric.d
        self.mu,self.mw,self.bins,self.rays,self.invariant=lapse.Lapse().rays(angular,d.T)
        assert self.invariant<1e-10
        background=incident.native.prior.EV/'coupled-128.npz'
        p=np.load(background);self.clock=np.linspace(0,d.T,17)
        ids=np.array([np.argmin(abs(p['snapshot_t']-t)) for t in self.clock])
        assert np.max(abs(p['snapshot_t'][ids]-self.clock))<1e-18
        outer=p['snapshot_I'][ids].sum(1)[:,-1];b=d.model.bulk
        self.lum=LD(2*np.pi*C)*d.model.area[-1]*(outer[:,b.mu>0]@(b.d['num']*b.d['Einf']))
        self.emission={}
        for label in ['high','low']:
            for n in [64,128]:
                path=saved.full.source.saved(n) if label=='high' else ACTUAL/f'sweep-1/photons/return-{n}.npz'
                p=np.load(path);edges=p['actual_step_edges'];count=2*(len(edges)-1)
                self.emission[label,n]=saved.Emission(edges,p['accepted_angular_luminosity'][:count])

    def ray(self,age,owner):
        age,owner=np.broadcast_arrays(age,owner);flat=age.ravel();ids=owner.ravel().astype(int)
        values=np.empty((2,len(flat)))
        for lo in range(0,len(flat),2048):
            hi=min(lo+2048,len(flat));v=self.rays.sol(flat[lo:hi]/self.d.T).reshape(3,len(self.mu),hi-lo)
            values[:,lo:hi]=v[:2,ids[lo:hi],np.arange(hi-lo)]
        return self.d.r0*values[0].reshape(age.shape),values[1].reshape(age.shape)

    def cohorts(self,now,order,cells=None):
        edges=self.clock[self.clock<=now];assert edges[-1]==now
        gx,gw=leg.leggauss(order);ids=np.arange(len(edges)-1) if cells is None else np.asarray(cells)
        dt=np.diff(edges)[ids];te=(edges[ids,None]+dt[:,None]*(gx+1)/2).ravel()
        tw=(dt[:,None]*gw/2).ravel();lum=np.column_stack([np.interp(te,self.clock,np.asarray(v,float)) for v in self.lum.T])
        owner=np.tile(np.arange(len(self.mu)),len(te));launch=np.repeat(te,len(self.mu))
        energy=(tw[:,None]*lum[:,self.bins]*self.mu*self.mw).ravel().astype(LD)
        return launch,owner,energy

    def reference(self,label,n,now,order):
        emission=self.emission[label,n];gx,gw=leg.leggauss(order)
        lo=emission.edges[:-1];hi=np.minimum(emission.edges[1:],now);keep=hi>lo
        h=hi[keep]-lo[keep];te=(lo[keep,None]+h[:,None]*(gx+1)/2).ravel()
        lum=(emission.a[keep,None,:]+(te.reshape(-1,order)-lo[keep,None])[:,:,None]*emission.b[keep,None,:]).reshape(-1,4)
        energy=((h[:,None]*gw/2).ravel()[:,None]*lum[:,self.bins]*self.mu*self.mw).ravel()
        owner=np.tile(np.arange(len(self.mu)),len(te));r,mu=self.ray(now-np.repeat(te,len(self.mu)),owner)
        _,N,b,_,_=self.d.bg.metric(r/self.d.model.m.R);weight=energy*np.sqrt(b)/N
        values=np.asarray([[np.sum(weight*(self.d.r0/r)**k*v,dtype=LD) for v in [np.ones_like(mu),mu,mu*mu]] for k in [0,1,2]])
        expected=emission.primitives(np.array([now]))[0][0]@(np.arange(1,8,2,dtype=LD)/32)
        port=float(abs(np.sum(energy,dtype=LD)-expected)/max(emission.absolute@(np.arange(1,8,2,dtype=LD)/32),LD('1e-290')))
        assert port<1e-12,(label,n,port)
        return values,dict(component=label,clock=n,port_relative=port,reference_energy_erg=float(expected))

    def propagate(self,now,order,cells=None):
        start=time.monotonic();te,owner,e0=self.cohorts(now,order,cells);age=now-te;n=len(te)
        scale=self.metric.scale
        def rhs(s,y):
            rr,mu=self.ray(age*s,owner);z=self.metric.at(te+age*s,rr)
            # delta_t/scale and h/scale; the background ray is reused exactly.
            h=y[n:2*n];hdot=(-C*z['a']*mu*z['nup']-mu*mu*z['lt'])/scale
            td=(-z['zeta']/scale-(1-mu*mu)*h/(mu*mu))
            work=(z['nt']-mu*mu*z['lt'])/scale
            return np.r_[age*td,age*hdot,age*work]
        sol=solve_ivp(rhs,[0.,1.],np.zeros(3*n),method='DOP853',rtol=2e-8,atol=2e-11)
        assert sol.success,sol.message
        r,mu=self.ray(age,owner);z=self.metric.at(np.full(n,now),r)
        launch=self.metric.at(te,np.full(n,self.d.r0))
        u0=self.d.z0['alpha'][0]*self.metric.wave(te,np.zeros(n))[0]/self.d.r0
        h=scale*sol.y[n:2*n,-1];delay=scale*sol.y[:n,-1]
        dr=-C*z['a']*mu*delay;dm=(1-mu*mu)/mu*(h+(1/r-z['nr'])*dr)
        dH=z['nu']+u0+h;launchH=launch['nu']+u0;work=scale*sol.y[2*n:,-1]
        norm=max(float(np.max(abs(dH))),scale*1e-10)
        energy_identity=float(np.max(abs(dH-launchH-work))/norm)
        invariant=dH+dr/r-z['nu']-z['nr']*dr-mu*dm/(1-mu*mu)-u0
        inv=float(np.max(abs(invariant))/norm)
        values=[]
        for power in [0,1,2]:
            psi=(self.d.r0/r)**power;prime=-power*psi/r
            density=psi*(dH-z['nu']-z['lam'])+dr*(prime-(z['nr']+z['lr'])*psi)
            weight=e0*np.sqrt(z['b'])/z['N']
            values.append([np.sum(weight*v,dtype=LD) for v in [density,mu*density+psi*dm,mu*mu*density+2*psi*mu*dm]])
        result=dict(emission_t=te,owner=owner,background_packet_energy_erg=e0,radius_cm=r,direction=mu,
            delta_radius_cm=dr,delta_direction=dm,delta_log_H=dH,delta_arrival_seconds=delay,
            launch_delta_log_H=launchH,integrated_log_H_work=work,metric_nu=z['nu'],metric_lambda=z['lam'],
            stress_weak_moments_erg=np.asarray(values,LD))
        row=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,packets=n,rhs_evaluations=sol.nfev,
            time=now,energy_identity_relative=energy_identity,angular_invariant_relative=inv,
            physical_launch_energy_increment_erg=float(np.sum(e0*launchH,dtype=LD)),
            instantaneous_packet_energy_increment_erg=float(np.sum(e0*dH,dtype=LD)),
            propagated_work_erg=float(np.sum(e0*work,dtype=LD)),
            maximum_radius_shift_cm=float(np.max(abs(dr))),maximum_direction_change=float(np.max(abs(dm))),
            stress_weak_moments_erg=np.asarray(values,float).tolist())
        assert energy_identity<.002 and inv<1e-10,row
        return result,row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(ACTUAL/'pipeline-status.json')['state']=='completed'
    files=[Path(__file__),Path(incident.__file__),Path(lapse.__file__),ACTUAL/'result.json',
           CHARGE/'result.json',Path('native-short-return-exterior259-work/full/audit.json'),
           incident.native.prior.EV/'coupled-128.npz']
    files += [incident.FIELDS/f'born-g{q}.npz' for q in [4,8]]
    for n in [64,128]:files += [saved.full.source.saved(n),ACTUAL/f'sweep-1/photons/return-{n}.npz']
    write(OUT/'plan.json',dict(classification='Conjectural',
        previous_goal_turn='Progress: corrected complete119/231same-solution mass and conditional charge accepted and archived in Phase260.',
        claim='Propagate the actual background photons through the incident metric used by the accepted coupled high history; retain energy work, angular and radial motion together as an exterior stress distribution for the next GR boundary solve.',
        decision='Replace the frozen physical photon source by its actual time-dependent first variation. Do not infer the new final charge until the matching mass/scalar constraints and returned-metric component are applied.',
        reuse='Same corrected EOS vacuum, original17background emission knots, four angular bins, exact original high/low emission arrays and horizon. No fluid/EOS replay and no new physical grid.',
        metric='Use the identical primary plus saved bilinear Born U values. Differentiate that declared interpolation consistently; compare its jets and boundary values with the saved incident owner. These are finite-representation controls, not uniform derivative certification.',
        representation='Worldline and packet-stress variations, kept separate from large background arrays. Body debit uses launch H; later change in H is explicit metric work, not a fictitious additional launch flux.',
        settings=SETTINGS,gates=dict(energy_identity=.002,angular_invariant=1e-10,quadrature=.002,time=.02,port=1e-12),
        budgets=dict(prepare=120,check=180,pilot=900,each_production=7200,collect=600),virtual_GiB=16,CPU_threads_per_process=1,
        pilot='Measure three selected real emission cells at terminal time. Use2x worst measured per-packet cost for all17readouts before launching bounded production. Preserve completed snapshots on failure.',
        stop='No automatic horizon, physical resolution, or tolerance changes. A failed original gate leaves the physical exterior and final charge unadjudicated.',
        limits='Incident-metric geometric source only; returned-metric geometric source, reciprocal scalar work, background mass remap and same-boundary GR return remain required. No historical diagnostic charge is added.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',symbolic())


def check():
    initialize();m=Metric(8);d=m.d
    ts=np.linspace(0,d.T,17);xs=np.linspace(0,2*C*d.T,19)
    value=max(np.max(abs(m.wave(t,xs)[0]-d.wave(float(t),xs)[0])) for t in ts)/(incident.ETA*d.r0)
    assert value<1e-12,value
    # Independent analytic homogeneous-lapse ray: conformal u has no effect
    # on the null trajectory, while a spatially constant nu changes its clock.
    mu0=np.array([.03,.2,.6,.95]);end=.7
    def rhs(t,y):
        r=np.sqrt(1+2*mu0*t+t*t);mu=(mu0+t)/r;nu=t*(1-t)
        return -nu*np.ones_like(mu)
    sol=solve_ivp(rhs,[0,end],np.zeros(4),method='DOP853',rtol=1e-12,atol=1e-14)
    r=np.sqrt(1+2*mu0*end+end*end);mu=(mu0+end)/r
    delta=-mu*sol.y[:,-1];expected=mu*(end**2/2-end**3/3)
    flat=float(np.max(abs(delta-expected)));assert flat<1e-12
    old=[];new=[];boundary=[]
    for t in ts:
        U,Ut,Ux=d.wave(float(t),np.array([0.]));v=m.wave(t,0.)
        old.append([Ut[0],C*Ux[0]]);new.append([v[1],C*v[2]])
        z=m.at(np.array([t]),np.array([d.r0]))
        ext=np.sum(d.eq.h*((2*d.zout['Phi']*d.wave(float(t),d.xout)[0]/(d.rout*d.zout['b'])).reshape(d.eq.r.shape)@d.eq.w),dtype=LD)
        boundary.append([float(z['nu'][0]),float(d.z0['Phi'][0]*U[0]-ext)])
    old,new,boundary=map(np.asarray,[old,new,boundary])
    jets=np.max(abs(new-old),axis=0)/np.max(abs(old),axis=0)
    face=float(np.max(abs(boundary[:,0]-boundary[:,1]))/max(np.max(abs(boundary[:,1])),m.scale*1e-10))
    assert max(jets)<.002 and face<.002,(jets,face)
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,same_U_value_relative=float(value),
        flat_homogeneous_lapse_error=flat,incident_boundary_value_relative=face,incident_jet_representation_relative=jets.tolist(),
        scalar_metric_consistent_derivatives=True,full_physical_exterior=False))


def pilot():
    initialize();m=Photons(8,8);rows=[]
    for cell in [0,7,15]:
        z,row=m.propagate(m.d.T,8,[cell]);np.savez_compressed(OUT/f'pilot-{cell}.npz',**z);rows.append(row)
        write(OUT/'pilot-progress.json',dict(rows=rows))
    per=max(v['seconds']/v['packets'] for v in rows)
    packets=8*32*sum(range(1,17));upper=2*per*packets+120
    result=dict(classification='Counterexample candidate',rows=rows,forecast_upper_seconds_per_full_path=upper,eligible=upper<7200)
    write(OUT/'pilot.json',result);print(json.dumps(result),flush=True);assert result['eligible'],result


def production(action):
    assert read(OUT/'check.json')['passed'] and read(OUT/'pilot.json')['eligible']
    initialize();a,t,g=SETTINGS[action];m=Photons(a,g);folder=OUT/action;folder.mkdir();rows=[]
    for j,now in enumerate(m.clock[1:],1):
        z,row=m.propagate(float(now),t);ports=[]
        for label in ['high','low']:
            moments,port=m.reference(label,128,float(now),t);z[f'{label}_reference_stress_weak_moments_erg']=moments;ports.append(port)
        z['high_physical_stress_weak_moments_erg']=z['high_reference_stress_weak_moments_erg']+z['stress_weak_moments_erg']
        row.update(reference_ports=ports,high_physical_stress_weak_moments_erg=np.asarray(z['high_physical_stress_weak_moments_erg'],float).tolist(),
            actual_same_solution_emission_applied=True,returned_metric_geometric_response_complete=False)
        np.savez_compressed(folder/f'snapshot-{j}.npz',**z)
        write(folder/f'snapshot-{j}.json',row);rows.append(row)
        write(folder/'progress.json',dict(completed=j,total=16,latest=row));print(json.dumps(dict(action=action,completed=j,seconds=row['seconds'])),flush=True)
    write(folder/'result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,full_physical_exterior=False))


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    if action=='prepare':cap=120
    else:
        plan=read(OUT/'plan.json');cap=plan['budgets'].get(action,7200)
        for p,h in plan['bindings'].items():assert sha(p)==h,p
    incident.native.deadline(cap)
    try:
        if action in SETTINGS:production(action)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
