"""Counterexample candidate: evolve the inventory-matched initial projection.

The new metric is installed in every material/photon operator. It remains
frozen during this short experiment: Krr is initial data, not evolved GR.
"""
from pathlib import Path
from types import FunctionType, MethodType
import inspect
import json
import signal
import sys
import textwrap
import time
import numpy as np
from numpy.polynomial import legendre as leg
import def_native_initial_constraints as initial

previous=initial.previous;old=previous.old;C=initial.C;G=initial.G
write=initial.write;sha=initial.sha;green=previous.feedback.prior.green
OUT=initial.OUT.parent/'def-native-projected-evolution'
INPUT=initial.OUT/'finite-volume';LD=np.longdouble


class Background:
    def __init__(self,R):
        self.R=R;self.z=z=np.load(INPUT/'balanced-20.npz');self.edges=z['edges']
        self.rend=float(self.edges[-1]);self.M=float(z['ADM_mass'])/R;self.K=float(z['scalar_K'])/R
        self.q=initial.Quadrature(self.edges,20)
        inverse=np.linalg.inv(leg.legvander(self.q.x,19)).astype(LD)
        self.co={k:z['panel_'+k].astype(LD)@inverse.T for k in ['mass','phi','Phi','lapse','nu_prime']}
        # The exterior starts beyond ALL initial photons, not at the old gas surface.
        from def_native_photon_exterior import vacuum
        mu=float(z['mass_faces'][-1]/self.rend)
        b=1-2*mu;N=float(np.exp(z['nu_faces'][-1]))
        flux=float(z['scalar_flux_faces'][-1]/(N*np.sqrt(b)*self.rend))
        self.vm,self.vk,self.solution=vacuum.background(mu,flux,2e-13)
        self.gx,self.gw=leg.leggauss(16)

    def metric(self,x):
        r=np.asarray(x,float)*self.R;y=r/self.rend
        assert np.min(y)>=1-1e-14,'Vacuum used inside initial material/photon domain'
        mass,c=self.solution.sol(1/y);b=1-2*mass/y;N=c/np.sqrt(b)
        return mass*self.rend/self.R,N,b,c,self.vk/(c*y*y)*self.R/self.rend

    def fields(self,r):
        r=np.asarray(r,float);ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,len(self.q.h)-1)
        x=(np.minimum(r,self.rend)-self.edges[ids])/self.q.h[ids]-1;V=leg.legvander(x,19)
        out={k:np.asarray(np.sum(co[ids]*V,axis=-1,dtype=LD),float) for k,co in self.co.items()}
        outside=r>self.rend
        if outside.any():
            rr=r[outside];mass,N,b,cc,v=self.metric(rr/self.R)
            y=rr/self.rend;nodes=1+(y[:,None]-1)*(self.gx+1)/2
            grad=self.metric(nodes.ravel()*self.rend/self.R)[-1]/self.R
            out['phi'][outside]=float(self.z['phi_faces'][-1])+(rr-self.rend)/2*(grad.reshape(nodes.shape)@self.gw)
            out['mass'][outside]=mass*self.R;out['lapse'][outside]=N;out['Phi'][outside]=v/self.R
            out['nu_prime'][outside]=mass*self.R/(rr*rr*b)+rr*(v/self.R)**2/2
        out['b']=1-2*out['mass']/r
        return out


class Geometry(green.Geometry):
    def __init__(self,m):
        self.m=m;self.bg=m.bg;self.gx,self.gw=leg.leggauss(16)

    def metric(self,x):
        rj=self.m.RJ+np.asarray(x,float);r=rj/self.m.As
        for _ in range(3):
            z=self.bg.fields(r);phi=z['phi'];A=np.exp(-2*phi*phi);alpha=-4*phi
            den=1+alpha*r*z['Phi'];r-=(A*r-rj)/(A*den)
        z=self.bg.fields(r);phi=z['phi'];A=np.exp(-2*phi*phi);alpha=-4*phi
        den=1+alpha*r*z['Phi'];a=A*z['lapse'];B=1/(np.sqrt(z['b'])*den)
        return alpha*a/r,a,B,r/self.m.R,phi

    def physical(self,r):
        z=self.bg.fields(r);A=np.exp(-2*z['phi']**2);alpha=-4*z['phi'];den=1+alpha*r*z['Phi']
        a=A*z['lapse'];B=1/(np.sqrt(z['b'])*den)
        ap=a*(z['nu_prime']+alpha*z['Phi'])/(A*den)
        return A*r,a,B,ap


def remap(I,number,E,ratio,volume_ratio):
    """Positive two-node packet remap; a conservative end-box moment repair.

    No packet is dropped at a finite-frequency endpoint. A tiny positive
    spectral tilt restores the clipped endpoint energy at fixed packet count.
    This declared finite-bin projection is not a spectral continuum proof.
    """
    I=np.asarray(I);n,q,nf=I.shape;packets=(I*number).reshape(n*q,nf)
    dest=E[None,:]*np.repeat(ratio,q)[:,None]
    hi=np.clip(np.searchsorted(E,dest),1,nf-1);lo=hi-1
    fraction=np.clip((dest-E[lo])/(E[hi]-E[lo]),0,1);mapped=np.zeros_like(packets)
    rows=np.arange(n*q)[:,None]
    np.add.at(mapped,(rows,lo),packets*(1-fraction));np.add.at(mapped,(rows,hi),packets*fraction)
    count=packets.sum(1);active=count>0;mean=np.zeros(n*q);target=np.zeros(n*q);variance=np.ones(n*q)
    mean[active]=(mapped@E)[active]/count[active];target[active]=(packets*dest).sum(1)[active]/count[active]
    variance[active]=(mapped*(E[None,:]-mean[:,None])**2).sum(1)[active]/count[active]
    tilt=(target-mean)/variance;factor=1+tilt[:,None]*(E[None,:]-mean[:,None]);assert factor.min()>0
    mapped*=factor
    out=mapped.reshape(I.shape)/number*ratio[:,None,None]**3*volume_ratio[:,None,None]
    errN=np.max(abs(mapped.sum(1)[active]/count[active]-1))
    errE=np.max(abs((mapped@E)[active]/(target[active]*count[active])-1))
    return out,dict(number=float(errN),energy=float(errE),endpoint_tilt=float(np.max(abs(factor-1))))


class Table(previous.Table):
    def __init__(self,base,theta0,density_shift,frequency_shift):
        self.__dict__.update(base.__dict__);self.theta0=theta0;self.density_shift=density_shift;self.frequency_shift=frequency_shift
        lo,hi=self.temperature_bounds;self.temperature_bounds=(lo-theta0,hi-theta0)

    def gas(self,theta,eta):
        x=self.x;self.x=self.density_shift+(1+self.density_shift)*x
        try:return super().gas(theta+self.theta0,eta)
        finally:self.x=x

    def radiation(self,theta,eta):
        x=self.x;self.x=self.density_shift+(1+self.density_shift)*x
        try:v=list(super().radiation(theta+self.theta0,eta))
        finally:self.x=x
        df=self.spectral['frequency'][0]*(1-self.f)+self.spectral['frequency'][1]*self.f
        for k in [0,1]:
            factor=1+df[:,k]*self.frequency_shift[:,None];assert factor.min()>0
            for i in [k,k+2,k+4]:v[i]*=factor
        return v


class Coupled(previous.Coupled):
    def __init__(self):
        super().__init__();b=self.bulk;f=self.flow;m=self.m;n=b.n;zero=np.zeros(n)
        d0=b.d;oldvol=np.r_[b.volume,4*np.pi*m.RJ*m.RJ*m.vol];olda=np.r_[d0['a'],m.a]
        oldr=np.r_[d0['r'],m.r];oldedges=np.r_[d0['edges'][:-1],m.rf];oldf0=dict(self.f0)
        oldincoming=b.incoming.copy();oldinnera=float(self.f0['af'][0])
        rho,v,lt,y=f.primitive(f.initial);f.eos.y=y;_,_,gamma,_,_=f.eos(rho,lt)
        ref=np.load(INPUT/'balanced-initial-state.npz');m.bg=Background(m.R);bg=m.bg;geo=Geometry(m)
        re=ref['radius_E'];ef=bg.edges[-(n+self.n+1):]
        assert len(ef)==n+self.n+1
        radii,a,B,ap=geo.physical(re);faces,af,Bf,_=geo.physical(ef)
        q=initial.Quadrature(ef,20);z=bg.fields(q.r.ravel());A=np.exp(-2*z['phi']**2)
        volume=np.asarray(q.h*(4*np.pi*q.r**2*(A**3/np.sqrt(z['b'])).reshape(q.r.shape)@q.w),float)
        vr=oldvol/volume;dr=np.log(vr)
        anchors=json.loads((INPUT/'audit.json').read_text())['native_anchors']
        theta0=np.log(np.array([x['temperature_K'] for x in anchors])/d0['T'])
        # Keep native interpolation anchors immutable; new evolution theta is
        # relative to its projected T0. Density/frequency corrections use the
        # already audited first-order native derivative bank.
        b.eos=Table(b.eos,theta0,vr[:n]-1,-np.log(a[:n]/olda[:n]))
        b.d=d=dict(d0);d.update(r=radii[:n],edges=faces[:n+1],a=a[:n],B=B[:n],rho=d0['rho']*vr[:n],
            T=d0['T']*np.exp(theta0),phi=np.asarray(ref['phi'][:n],float),face_a=af[:n+1],face_B=Bf[:n+1])
        b.volume=volume[:n];b.W=b.volume/(4*np.pi*d['a']**3);b.area=faces[:n+1]**2/af[:n+1]**2
        b.A=old.transport(b.W,b.area,b.mu,b.edges_mu);b.stream=b.A.tocoo();b.scfactor=d['a']*C*6.6524587321e-25
        b.photon_energy_weight=4*np.pi*b.W[:,None,None]*b.w[None,:,None]*d['num']*d['Einf']
        b.initial,deepremap=remap(b.initial,d['num'],d['Einf'],a[:n]/olda[:n],vr[:n])
        b.incoming,innerremap=remap(oldincoming[None],d['num'],d['Einf'],np.array([af[0]/oldinnera]),np.ones(1));b.incoming=b.incoming[0]
        b.boundary=np.zeros_like(b.initial);pos=b.mu>0
        b.boundary[0,pos]=C*b.area[0]*b.mu[pos,None]/b.W[0]*b.incoming[pos]
        self.initial_I,atmoremap=remap(self.initial_I,d['num'],d['Einf'],a[n:]/olda[n:],vr[n:])
        # Fixed Einstein radii are the Cauchy mesh; map every owner to its new
        # Jordan radii. RJ is a reference area normalization, not a gas wall.
        m.RJ,m.a0,m.As=map(float,[geo.physical(np.array([m.R]))[0][0],geo.physical(np.array([m.R]))[1][0],np.exp(-2*bg.fields(np.array([m.R]))['phi'][0]**2)])
        m.r=radii[n:];m.rf=faces[n:];m.x=m.r-m.RJ;m.xf=m.rf-m.RJ
        m.a=a[n:];m.B=B[n:];m.ap=ap[n:];m.af=af[n:];m.Bf=Bf[n:]
        m.area=(m.rf/m.RJ)**2;m.vol=volume[n:]/(4*np.pi*m.RJ**2);m.dx=float(np.min(np.diff(m.rf)))
        newlt=lt+(gamma-1)*dr[n:];newrho=rho*vr[n:]
        f.initial_temperature=newlt.copy();f.seed=newlt.copy();f.eos.y=y
        f.initial=f.conserved(newrho,v,newlt,y,m.a)[0];f.background_cell=np.array([newrho,np.zeros(self.n),newlt])
        tshift=np.r_[theta0,newlt-lt];face_dr=np.interp(oldedges,oldr,dr);face_t=np.interp(oldedges,oldr,tshift)
        for key in ['background_left','background_right']:
            value=getattr(f,key).copy();value[0]*=np.exp(face_dr[n:]);value[2]+=face_t[n:];setattr(f,key,value)
        self.join0=np.r_[f.background_left[:,0],f.eos.y0];f.join_state=self.join0.copy()
        f.hydro=MethodType(hydro_ns['rhs'],f)
        self.W=volume[n:]/(4*np.pi*m.a**3);self.area=m.rf**2/m.af**2
        self.A=old.transport(self.W,self.area,self.mu,b.edges_mu)
        self.pm=self.A[self.pos][:,self.neg];self.mm=self.A[self.neg][:,self.neg];self.pp=self.A[self.pos][:,self.pos];self.lu={}
        self.energy_weight=4*np.pi*self.W[:,None,None]*self.w[None,:,None]*d['num']*d['Einf']
        self.gas_scale=4*np.pi*m.RJ*m.RJ*f.eos.rho0*C*C
        self.mech.xi=self.mech.xi*(d0['r']**2*d0['B']/(d['r']**2*d['B']))[:,None]/vr[:n,None]
        self.set_material(0.);p0,b.u0,*_=b.eos.gas(zero,zero)
        # Preserve the old physical gas gradient plus the declared change in
        # native pressure. Do not reimpose hydrostatic balance by definition.
        dp=p0-oldf0['p0'];grad_delta=np.gradient(dp,d['r'])
        dg=ap[:n]/a[:n]-oldf0['ap']/oldf0['a']
        oldenthalpy=d0['rho']*(self.cx*C*C+oldf0['u0'])+oldf0['p0']
        deltaenthalpy=(d['rho']-d0['rho'])*(self.cx*C*C+oldf0['u0'])+d['rho']*(b.u0-oldf0['u0'])+dp
        support=oldf0['initial_support']*(a[:n]/B[:n])/(oldf0['a']/oldf0['B'])
        support-=a[:n]/B[:n]*(grad_delta+deltaenthalpy*ap[:n]/a[:n]+oldenthalpy*dg)
        self.f0=dict(oldf0,r=d['r'],edges=d['edges'],volume=b.volume,a=d['a'],B=d['B'],ap=ap[:n],af=af[:n+1],Bf=Bf[:n+1],
            rho=d['rho'],rho_face=oldf0['rho_face']*np.exp(face_dr[:n+1]),p0=p0,u0=b.u0,initial_support=support,pgprime=oldf0['pgprime']+grad_delta)
        self.edge=d['edges'];self.area_gas=4*np.pi*self.edge**2*self.f0['af'];self.face_rho=self.f0['rho_face']
        self.face_K=np.interp(self.edge,d['r'],self.mech.K[0]);self.face_p=np.interp(self.edge,d['r'],p0)
        f.eos.y=np.array([f.eos.y0]);self.face_p[-1]=f.eos(np.array([self.join0[0]]),np.array([self.join0[2]]))[0][0]*f.eos.rho0*C*C
        self.base_momentum_flux=self.area_gas*self.face_p
        self.installation=dict(remap=[deepremap,atmoremap,innerremap],
            mass=float(np.max(abs(d['rho']*b.volume/self.mass0-1))),
            density=float(np.max(abs(np.r_[d['rho'],newrho*f.eos.rho0]-ref['density'])/np.maximum(ref['density'],1e-300))),
            maximum_density_shift=float(np.max(abs(dr))),maximum_temperature_shift=float(np.max(abs(tshift))),
            initial_support_change=float(np.max(abs(support-oldf0['initial_support']))),
            corrected_initial_geometry=True,fixed_metric=True,full_GR_evolution=False)
        assert b.area[-1]==self.area[0]


hydro=previous.prior.base.hydro
hydro=previous.replace(hydro,'C*m.dx/m.vol','C*np.diff(m.rf)/m.vol')
hydro=previous.replace(hydro,'dt=.35*m.dx/np.max(C*m.a/m.B*(abs(v)+cs+1e-100))','dt=.35*np.min(np.diff(m.rf)/(C*m.a/m.B*(abs(v)+cs+1e-100)))')
hydro_ns=dict(previous.prior.base.namespace);exec(compile(hydro,__file__,'exec'),hydro_ns)
runner=previous.prior.base.runner
runner=previous.replace(runner,'p0=b.eos.base.gas(np.zeros(b.n),np.zeros(b.n))[0]',"p0=self.f0['p0']")
run_ns=dict(previous.Coupled.run.__globals__,OUT=OUT);exec(compile(runner,__file__,'exec'),run_ns)
Coupled.run=run_ns['run']


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='f9987fc1a',
        claim='Install the accepted inventory-matched initial GR/scalar projection in all actual material and photon owners, complete the same3.434ms64/128 trajectories and recompute their direct charge.',
        decision='Determine whether the old positive direct response survives a consistent physical initial projection. No old numerical lower bound is inherited.',
        reuse='All19+512 cells,152 frequencies,8 angles, original native tables and first-order density/frequency derivatives. Reuse Phase119 only as a comparison; do not repeat its trajectories.',
        model='Frozen corrected initial metric, exact proper volumes, fixed Einstein mesh mapped to Jordan radius, conserved initial baryons and angular photon number/energy moments, native entropy anchors and conservative frequency remap. Deep theta is relative to new T0. Keep old native gas-gradient reconstruction plus its pressure correction; never reimpose equilibrium. Atmosphere uses the existing native cold EOS and actual pressure/gravity fluxes. Same shared material/photon face.',
        limitations='Initial Krr and scalar balance are retained data but dynamic metric/scalar feedback is not evolved. Local Gamma continuation, first-order deep EOS derivatives, finite spectral endpoint projection, atmosphere Gamma temperature adjustment and inherited spatial gradient reconstruction remain approximations. No final physical charge or observable claim.',
        budget=dict(check_seconds=35,pilot_seconds=40,production_seconds=500,readout_seconds=60,CPU_threads=1,memory_GB=3,paths=[64,128],new_native_bank_calls=0),
        forecast='Phase119 production took333s. Measure two full new steps per clock; resume their prefixes only if1.8x projected remainder plus10s fits500s. Late nonlinear cost unmeasured.',
        gates=dict(inventory=1e-12,photon_moments=1e-12,projected_density=1e-12,native_initial=.002,energy=1e-8,baryon=1e-10,source=1e-9,time_direct=.02,time_total=.02,history_direct=.02,cadence_direct=.02,quadrature_direct=.02),
        stop='Preserve every failure. No automatic third path, finer grid, longer time, table expansion or weaker gate.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(initial.__file__),Path(previous.__file__),INPUT/'balanced-20.npz',INPUT/'balanced-initial-state.npz',INPUT/'audit.json',previous.OUT/'coupled-128.npz']}))
    (OUT/'expanded-hydro.py').write_text(hydro);(OUT/'expanded-run.py').write_text(runner)


def check():
    assert not (OUT/'check.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(35)
    model=Coupled();b=model.bulk;f=model.flow;zero=np.zeros(b.n);ref=np.load(INPUT/'balanced-initial-state.npz')
    row=dict(model.installation);p,u,*_=b.eos.gas(zero,zero)
    anchors=json.loads((INPUT/'audit.json').read_text())['native_anchors']
    row['native_pressure']=float(np.max(abs(p/np.array([x['native_pressure'] for x in anchors])-1)))
    row['native_energy']=float(np.max(abs(u/np.array([x['native_specific_internal_energy'] for x in anchors])-1)))
    state=model.material_state(b.u0,zero);uu,tt,yy=model.recover_material(state,0.,zero)
    row['initial_recovery_logT']=float(np.max(abs(tt)));row['initial_recovery_neutral']=float(np.max(abs(yy)))
    photons=np.concatenate([b.initial,model.initial_I]);a=np.r_[b.d['a'],model.m.a]
    E=np.einsum('iqf,q,f->i',photons,b.w,b.d['num']*b.d['Einf'])/a**4
    P=np.einsum('iqf,q,f->i',photons,b.w*b.mu2,b.d['num']*b.d['Einf'])/a**4
    row['initial_photon_energy']=float(np.max(abs(E/ref['initial_photon_energy']-1)))
    row['initial_photon_pressure']=float(np.max(abs(P/ref['initial_photon_radial_pressure']-1)))
    # Nonzero physical source and paired material/photon ports are evaluated
    # before dispatch, without changing their acceptance based on the result.
    V=f.primitive(f.initial);_,right=f.reconstruct(V,0.);model.join_flux(right[:,0],V[3,0])
    force=model.material_rhs(b.u0,zero,zero)[1]
    row['maximum_initial_force_dyn']=float(np.max(abs(force)));row['initial_join_mass_g_s']=float(model.mflux[0])
    import sympy as s
    N,E0,a0,a1,V0,V1=s.symbols('N E0 a0 a1 V0 V1',positive=True)
    assert s.simplify((N*a1**3*V0/(a0**3*V1))*V1/a1**3-N*V0/a0**3)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Photon packet-number volume/redshift transformation; matched material and photon face fluxes telescope as before. No continuum or dynamic GR theorem.'))
    row.update(classification='Counterexample candidate',seconds=time.monotonic()-start,
        passed=bool(row['mass']<1e-12 and row['density']<1e-12 and max(row['initial_photon_energy'],row['initial_photon_pressure'])<1e-12 and max(row['native_pressure'],row['native_energy'])<.002 and max(max(r['number'],r['energy']) for r in row['remap'])<1e-12))
    write(OUT/'check.json',row);print(json.dumps(row),flush=True);signal.alarm(0);assert row['passed']


def pilot():
    assert json.loads((OUT/'check.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(40);start=time.monotonic();rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'pilot-{steps}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)!=2 or not all(r['passed'] for r in rows) else sum(r['seconds']/2*(r['steps']-2) for r in rows)
    upper=None if forecast is None else 1.8*forecast+10
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,upper_seconds=upper,eligible=bool(upper is not None and upper<500),seconds=time.monotonic()-start));signal.alarm(0)


def production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'production.json').exists()
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():
        assert sha(OUT/'pilot-producer.py' if p==str(Path(__file__)) else p)==h,p
    assert sha(__file__)==json.loads((OUT/'restart-receipt.json').read_text())['source_sha256']
    start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(500);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'coupled-{steps}',restart=f'pilot-{steps}');rows.append(row)
        if not row['passed']:break
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=bool(len(rows)==2 and all(r['passed'] for r in rows)),paths=rows,seconds=time.monotonic()-start,corrected_initial_state_evolved=True,full_GR_feedback=False,final_charge_solved=False));signal.alarm(0)


def receipt():
    write(OUT/'restart-receipt.json',dict(classification='Counterexample candidate',source_sha256=sha(__file__),
        issue='Pre-production inspection found that a restart would evaluate the history-only pressure baseline using the current perturbed density. Both completed two-step pilots began at zero displacement, so their reference was correct.',
        correction='Read the installed immutable f0 pressure baseline when recording histories, including restarts. No RHS, initial state, accepted step or gate changes. Continue the saved pilot prefixes.',
        pilot_producer_sha256=sha(OUT/'pilot-producer.py')))
    (OUT/'expanded-run.py').write_text(runner)


if __name__=='__main__':globals()[sys.argv[1]]()
