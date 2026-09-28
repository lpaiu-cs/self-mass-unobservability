"""Counterexample candidate: the emitted radiation metric scatters the input.

Compute the missing exterior mixed term from actual constant-stage photon
emission. Differentiate the moving shell as a distribution before quadrature.
This is a scalar response on the direct radiation metric, not full nonlinear GR.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from numpy.polynomial import legendre as leg
import read_native_incident_infinity as previous
import def_native_dynamic_lapse as metric

OUT=Path('native-exterior-incident159-work')
EV=previous.prior.EV
COARSE=previous.ROOT/'native-retained-tail/supported-temperature/evolution'
C=previous.C;G=previous.G;LD=np.longdouble
read,write,sha=previous.read,previous.write,previous.sha
CAPS=dict(prepare=30,check=45,pilot=60,production=240,bounds=90,audit=30)
SETTINGS=[(128,8,8,8),(128,4,8,8),(128,8,4,8),(128,8,8,4),(64,8,8,8)]


def symbolic():
    x,t,c,r,m,a,b=sp.symbols('x t c r m a b',positive=True)
    Z=sp.Function('Z')(t,x);g=sp.Function('g');U=g(t+x/c)
    first=2*Z*sp.diff(U,x,2)+sp.diff(Z,x)*sp.diff(U,x)+sp.diff(Z,t)*sp.diff(U,t)/c**2
    divergence=sp.diff(Z*sp.diff(U,x),x)+sp.diff(Z*sp.diff(U,x),t)/c
    assert sp.simplify(first-divergence)==0
    derivative=-1/(a*b*r*r)-1/(a*r*r)-2*m/(a*r**3*b)
    assert sp.simplify((derivative+2/(a*b*r*r)).subs(b,1-2*m/r))==0
    # The scalar operator on a constant scalar vanishes in any metric.
    f=sp.symbols('f');assert sp.diff(r*f,r)-(r*f)/r==0
    return dict(classification='Proven',passed=True,
        metric='For a packet e=G E/c^4 outside r, lambda_B=-e/(r a); zeta_B=e*[I(r,rp)+1/(r a)-mu_p^2/(rp a_p)], I=int_r^rp ds/(a b s^2). Outside the packet both direct radiation metric increments vanish. Its interior mass debit is paired with its emitted energy.',
        derivative='The smooth zeta_x=-2e/(b r^2). Its shell term is -e*(1-mu_p^2)/(rp a_p)*delta(x-xp).',
        operator='S=2*zeta*(U_xx-Vgeom*U)+zeta_x*(U_x-a*U/r)+zeta_t*U_t/c^2, Vgeom=a*a_r/r. This is variation of the scalar wave operator on the direct radiation metric.',
        primary='For the exact incoming primary U=g(t+x/c), the first three derivative terms equal (partial_x+partial_t/c)(zeta*U_x). At outgoing infinity their causal-edge terms cancel, leaving -int_0^u zeta(0,t)*U_x(0,t)dt at the inner edge of the exterior domain.',
        reduced='The remaining volume integrand in dr dt per unit e is [2/(b*r^3)-4*zeta_unit*a*m/(b*r^3)]*U. The moving-shell integrand in dt per e is (1-mu_p^2)*U(t,xp)/rp^2. Multiply the sum of boundary, volume and shell integrals by c/2.',
        scope='These are exterior operator and distribution identities. Reciprocal scalar-metric mixed constraints, interior background operator variation and complete physical errors are not certified.')


def prepare(repair=False):
    if repair:
        failed=read(OUT/'prepare-receipt.json');assert failed['error'].startswith('FileNotFoundError')
        assert sha(OUT/'prepare-failed-producer.py')==failed['source_sha256'] and not (OUT/'plan.json').exists()
    else:assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(previous.__file__),Path(previous.incident.__file__),Path(metric.__file__),
           previous.OUT/'final-result.json',previous.OUT/'charge-parts.npz',previous.BACKGROUND,
           previous.incident.FIELDS/'born-g8.npz',previous.wave.base.flow.INPUT/'balanced-20.npz']
    for n in [64,128]:
        folder=COARSE if n==64 else EV
        files.extend([folder/f'accepted-ports-{n}.npz',folder/f'source-{n}.npz'])
    if repair:files.extend([OUT/'prepare-failed-producer.py',OUT/'prepare-receipt.json'])
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='242b5e4fc',
        previous_goal_turn='Progress: Phase158 applied actual coupled emission at fixed-operator null infinity and preserved the physical closure gaps.',
        claim='Apply the actual already-emitted background radiation metric to the declared incoming scalar primary and read its outgoing correction, including the moving shell and exterior-domain boundary. Assess the smaller photon-energy/arrival and Born/potential terms against the resulting charge scale.',
        decision='Does the omitted physical exterior interaction change the sign or size of the Phase158 matter-mediated charge? A changed sign or a failed control requires reconsidering that conditional conclusion, not more fluid replay.',
        reuse='Use the actual saved constant-stage64/128 angular luminosities, fixed four outgoing angular bins, initial exact vacuum, original input eta/duration and17 readout times. No fluid or EOS calls. Quadrature nodes integrate existing emission intervals; they are not new evolution steps.',
        method='Reuse the existing curved outward null rays. Split emission-time quadrature at its existing owner intervals. Split interaction time at the incoming front and the outgoing readout cone; integrate the distributional shell term explicitly. Eliminate derivative cancellation analytically before quadrature.',
        scope='The directly emitted radiation metric acting on the incoming primary is one missing mixed term. It does not include all scalar-metric reciprocal mixed constraints or the evolving interior operator. Do not call the complete physical exterior or final goal solved.',
        settings=SETTINGS,gates=dict(quadrature=.002,background_history=.02,port=1e-12,ray=1e-10,root_fraction=2e-12,optical_fraction=1e-10),
        budgets=CAPS,total_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,new_EOS_calls=0,new_fluid_steps=0,
        pilot='Measure endpoint cohorts in existing emission cells0,31,63 at full8-order settings. Forecast all nine nonzero/cut readouts and five settings with twice the worst per-cell cost plus20s; require the remaining240s cap. Reuse saved endpoint cohort values if the corresponding production path needs them.',
        stop='Keep pulse, intervals, background histories, orders and gates fixed. Stop on a failed forecast or physics/numerical gate. No automatic extra Born order, new clock, longer interval or fluid replay.',
        input_path_repair='The retained64-step history remains in native-retained-tail/supported-temperature/evolution, as its existing infinity owner specifies. Only128 resides in native-retained-completion/evolution. Initial registration stopped before any physical calculation; preserve its receipt/source and bind the actual existing64-step owner.' if repair else None,
        bindings={str(p):sha(p) for p in files}))
    write(OUT/'symbolic.json',symbolic());(OUT/'registered-producer.py').write_bytes(Path(__file__).read_bytes())


class Interaction:
    def __init__(self,angular):
        self.driver=d=previous.incident.Driver(8);self.metric=m=metric.Lapse()
        self.r0=d.r0;self.D=d.D;self.T=d.T;self.M=float(d.bg.z['ADM_mass'])
        self.mu,self.mw,self.bins,self.rays,self.invariant=m.rays(angular,self.D)
        assert self.invariant<1e-10 and self.r0==m.r0
        self.a0=float(d.z0['lapse'][0]*np.sqrt(d.z0['b'][0]));self.root_error=0.
        def rhs(y,_):
            r=self.r0+C*self.D*y;_,N,b,a,_=m.bg.metric(np.array([r/m.model.m.R]))
            return [float(1/a[0]),float(C*self.D*self.r0/(a[0]*b[0]*r*r))]
        self.vacuum=solve_ivp(rhs,[0,1],[0.,0.],rtol=4e-13,atol=2e-15,dense_output=True)
        assert self.vacuum.success

    def fields(self,r):
        shape=np.shape(r);mass,N,b,a,Phi=self.metric.bg.metric(np.asarray(r).ravel()/self.metric.model.m.R)
        return tuple(v.reshape(shape) for v in [mass*self.metric.model.m.R,N,b,a,Phi/self.metric.model.m.R])

    def optical(self,r):
        return C*self.D*self.vacuum.sol((np.asarray(r)-self.r0)/(C*self.D))[0]

    def primitive(self,r):
        return self.vacuum.sol((np.asarray(r)-self.r0)/(C*self.D))[1]/self.r0

    def inverse(self,x):
        r=self.r0+self.a0*np.asarray(x)
        for _ in range(4):r-=(self.optical(r)-x)*self.fields(r)[3]
        return r

    def ray(self,age,owners):
        age,owners=np.broadcast_arrays(age,owners);flat=age.ravel();ids=owners.ravel().astype(int)
        assert np.min(flat,initial=0)>-1e-17 and np.max(flat,initial=0)<self.D+1e-17
        result=np.empty((3,len(flat)))
        for lo in range(0,len(flat),8192):
            hi=min(len(flat),lo+8192)
            v=self.rays.sol(np.clip(flat[lo:hi]/self.D,0,1)).reshape(3,len(self.mu),hi-lo)
            result[:,lo:hi]=v[:,ids[lo:hi],np.arange(hi-lo)]
        return tuple(v.reshape(age.shape) for v in [self.r0*result[0],result[1],result[2]/self.r0])

    def root(self,target,limit,owners,sign):
        target,limit,owners=np.broadcast_arrays(target,limit,owners);rp,_,_=self.ray(limit,owners)
        active=limit+sign*self.optical(rp)/C>target
        age=np.minimum(limit,np.maximum(0,target/(1+sign*self.mu[owners])))
        for _ in range(6):
            rp,mu,_=self.ray(age,owners);residual=age+sign*self.optical(rp)/C-target
            age=np.where(active,np.clip(age-residual/(1+sign*mu),0,limit),limit)
        rp,_,_=self.ray(age,owners);residual=np.where(active,age+sign*self.optical(rp)/C-target,0.)
        error=float(np.max(abs(residual),initial=0)/self.D);self.root_error=max(self.root_error,error)
        assert error<2e-12,error
        return age

    def cohort(self,te,u,temporal,radial):
        te=np.asarray(te)[:,None];owners=np.broadcast_to(np.arange(len(self.mu)),(len(te),len(self.mu)))
        te=np.broadcast_to(te,owners.shape);limit=self.D-te
        tx,tw=leg.leggauss(temporal);rx,rw=leg.leggauss(radial)
        U=lambda t,x:previous.incident.ETA*self.r0*previous.incident.pulse((t+x/C)/self.D)
        Ux=lambda t:previous.incident.ETA*self.r0/(C*self.D)*previous.incident.pulse(t/self.D,1)
        # Boundary of the EXTERIOR domain, not an invented surface drive.
        end=np.maximum(np.minimum(u,self.D)-te,0);age=end[...,None]*(tx+1)/2
        rp,mu,ip=self.ray(age,owners[...,None]);ap=self.fields(rp)[3]
        z0=ip+1/(self.r0*self.a0)-mu*mu/(rp*ap)
        boundary=-np.sum(end[...,None]*tw/2*z0*Ux(te[...,None]+age),axis=-1)
        front=self.root(self.D-te,limit,owners,1)
        causal=self.root(np.maximum(u-te,0),limit,owners,-1) if u<self.D else limit
        end=np.minimum(front,causal);age=end[...,None]*(tx+1)/2
        rp,mu,_=self.ray(age,owners[...,None]);tt=te[...,None]+age;xp=self.optical(rp)
        shell=np.sum(end[...,None]*tw/2*(1-mu*mu)/rp**2*U(tt,xp),axis=-1)
        upper=np.maximum(np.minimum.reduce([limit,(self.D+u)/2-te,causal]),0)
        cuts=np.sort(np.stack([np.zeros_like(upper),np.minimum(front,upper),np.clip(u-te,0,upper),upper],axis=-1),axis=-1)
        lo=cuts[...,:-1];width=np.diff(cuts,axis=-1)
        age=lo[...,None]+width[...,None]*(tx+1)/2
        rp,mu,ip=self.ray(age,owners[...,None,None]);ap=self.fields(rp)[3];tt=te[...,None,None]+age
        rlow=self.inverse(C*np.maximum(tt-u,0));rhigh=np.minimum(rp,self.inverse(C*np.maximum(self.D-tt,0)))
        span=np.maximum(rhigh-rlow,0);rr=rlow[...,None]+span[...,None]*(rx+1)/2
        mass,_,b,a,_=self.fields(rr);zeta=ip[...,None]-self.primitive(rr)+1/(rr*a)-(mu*mu/(rp*ap))[...,None]
        volume_integrand=(2/(b*rr**3)-4*zeta*a*mass/(b*rr**3))*U(tt[...,None],self.optical(rr))
        volume=np.sum(width[...,None]*tw/2*np.sum(span[...,None]*rw/2*volume_integrand,axis=-1),axis=(-1,-2))
        return C/2*np.stack([boundary,volume,shell],axis=-1)


def emission(n):
    folder=COARSE if n==64 else EV
    p=np.load(folder/f'accepted-ports-{n}.npz');L=p['angular_luminosity'];h=float(p['h'])
    source=np.load(folder/f'source-{n}.npz');aw=np.arange(1,8,2,dtype=LD)/32
    cum=np.r_[LD(0),np.cumsum(L.astype(LD)@aw)*LD(h)]
    error=float(np.max(abs(cum-source['outer_cumulative_energy_erg']))/max(abs(cum[-1]),LD(1)))
    assert error<1e-12 and np.min(L)>=0
    return L,h,error


def cells(model,n,u,temporal,radial,selected=None):
    L,h,port=emission(n);gx,gw=leg.leggauss(temporal)
    count=int(round(min(u,model.D)/h));ids=np.arange(count) if selected is None else np.asarray(selected)
    assert np.all(ids<count)
    te=((ids[:,None]+(gx+1)/2)*h).ravel();start=time.monotonic()
    k=model.cohort(te,u,temporal,radial).reshape(len(ids),temporal,len(model.mu),3)
    energy=(h*gw/2)[None,:,None]*L[ids][:,None,model.bins]*model.mu*model.mw
    q=-LD(G)/LD(C)**4/LD(model.M)*np.sum(energy[...,None].astype(LD)*k.astype(LD),axis=(1,2),dtype=LD)
    return q,dict(steps=n,port_relative=port,seconds=time.monotonic()-start,cells=len(ids),root_fraction=model.root_error)


def check():
    previous.prior.initialize();m=Interaction(8);r=np.linspace(m.r0,m.r0+C*m.D,19)
    error=float(np.max(abs(m.optical(r)-m.driver.optical(r)))/(C*m.D));assert error<1e-10,error
    back=float(np.max(abs(m.inverse(m.optical(r))-r))/(C*m.D));assert back<1e-10
    # The ray mass primitive and the independent radial vacuum ODE agree.
    age=np.linspace(0,m.D,17)[:,None];owner=np.arange(len(m.mu))[None,:]
    rp,mu,I=m.ray(age,owner);agreement=float(np.max(abs(I-m.primitive(rp)))*m.r0)
    assert agreement<1e-10
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,ray_invariant=m.invariant,
        optical_fraction=error,inverse_fraction=back,independent_radial_primitive_error=agreement,
        scope='Actual frozen vacuum geometry and rays; no continuous EOS or evolving metric certificate.'))


def pilot():
    previous.prior.initialize();m=Interaction(8);costs=[];saved=[]
    for j in [0,31,63]:
        q,row=cells(m,128,m.D,8,8,[j]);costs.append(row['seconds']);saved.append(q[0])
    # Eight distinct positive readouts untilD; later primary outputs saturate.
    # Conservative all-full-order forecast even for the cheaper control paths.
    upper=2*max(costs)*sum(range(8,65,8))*5+20
    result=dict(classification='Counterexample candidate',point_seconds=costs,upper_remaining_seconds=upper,eligible=upper<CAPS['production'])
    np.savez_compressed(OUT/'pilot.npz',cell_ids=[0,31,63],charge_parts=saved,u=m.D)
    write(OUT/'pilot.json',result);print(json.dumps(result),flush=True);assert result['eligible']


def production():
    assert read(OUT/'pilot.json')['eligible'];previous.prior.initialize();start=time.monotonic();models={};paths=[];rows=[]
    for n,a,t,r in SETTINGS:
        if a not in models:models[a]=Interaction(a)
        m=models[a];clock=np.linspace(0,m.T,17);parts=[np.zeros(3,LD)]
        for u in clock[1:9]:
            reuse=n==128 and a==t==r==8 and u==m.D
            if reuse:
                pilot=np.load(OUT/'pilot.npz');keep=np.setdiff1d(np.arange(64),pilot['cell_ids'])
                q,row=cells(m,n,u,t,r,keep);total=q.sum(axis=0,dtype=LD)+pilot['charge_parts'].sum(axis=0,dtype=LD)
            else:q,row=cells(m,n,u,t,r);total=q.sum(axis=0,dtype=LD)
            parts.append(total);rows.append(dict(row,angular=a,temporal=t,radial=r,retarded_time=float(u)))
        parts=np.asarray(parts+parts[-1:]*8,LD);data=dict(t=clock,charge_parts=parts,charge=parts.sum(axis=1,dtype=LD))
        np.savez_compressed(OUT/f'charge-{n}-a{a}-t{t}-r{r}.npz',**data);paths.append(data)
    fine=paths[0];norm=max(float(np.max(abs(fine['charge']))),1e-290)
    controls={key:float(np.max(abs(d['charge']-fine['charge']))/norm) for key,d in zip(['angular','temporal_quadrature','radial','background_history'],paths[1:])}
    component_norm=np.maximum(np.max(abs(fine['charge_parts']),axis=0),LD('1e-290'))
    component_controls={key:[float(v) for v in np.max(abs(d['charge_parts']-fine['charge_parts']),axis=0)/component_norm]
        for key,d in zip(['angular','temporal_quadrature','radial','background_history'],paths[1:])}
    body=np.load(previous.OUT/'charge-parts.npz')['body'];e0=np.load(previous.BACKGROUND)['epsilon'].astype(LD)
    correction=fine['charge']/(1-e0);combined=body+correction
    np.savez_compressed(OUT/'applied-charge.npz',t=fine['t'],previous_body=body,exterior_interaction=correction,combined=combined,
        exterior_charge_parts=fine['charge_parts'],direct=np.load(previous.OUT/'charge-parts.npz')['direct'])
    passed=max(max([controls[k]]+component_controls[k]) for k in ['angular','temporal_quadrature','radial'])<.002 and max([controls['background_history']]+component_controls['background_history'])<.02
    result=dict(classification='Counterexample candidate',passed=passed,controls=controls,component_controls=component_controls,
        endpoint_exterior_parts=fine['charge_parts'][-1].tolist(),endpoint_exterior=float(correction[-1]),previous_body_endpoint=float(body[-1]),
        combined_endpoint=float(combined[-1]),exterior_over_previous_body=float(correction[-1]/body[-1]),
        actual_background_emission_consumed=True,moving_shell_term_included=True,exterior_domain_boundary_included=True,
        scalar_mixed_constraint_feedback_complete=False,evolving_interior_operator_included=False,
        physical_exterior_closed=False,full_goal_complete=False,seconds=time.monotonic()-start,rows=rows)
    # JSON cannot encode NumPy longdouble; preserve wide values in NPZ.
    result['endpoint_exterior_parts']=[float(v) for v in result['endpoint_exterior_parts']]
    write(OUT/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);assert passed


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','prepare_repair','check','pilot','production']
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    cap=CAPS['prepare']-read(OUT/'prepare-receipt.json')['seconds'] if action=='prepare_repair' else CAPS[action]
    previous.incident.native.deadline(cap);start=time.monotonic();cpu=time.process_time();error=None
    try:
        if not action.startswith('prepare'):
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action=='prepare_repair':prepare(True)
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
