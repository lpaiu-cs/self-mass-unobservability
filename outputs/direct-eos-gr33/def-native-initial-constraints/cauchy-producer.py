"""Counterexample candidate: inventory-preserving initial Einstein constraints.

The scalar Cauchy field and its zero normal time derivative are prescribed.
They are not constrained to be stationary. A local isentropic Gamma closure
is explicit; matching these constraints does not certify full native evolution.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from numpy.polynomial import legendre as leg
import def_native_boundary_layer as previous

OUT=previous.OUT.parent/'def-native-initial-constraints'
write=previous.write;sha=previous.sha;C=previous.C;G=previous.feedback.prior.G
LD=np.longdouble


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='537d57c',
        claim='Construct an actual initial mass/lapse solution for the material inventory and photon moments used by Phase119, and supply the missing momentum-constraint extrinsic curvature. Retain compensated differences for the downstream physical evolution.',
        equations='Pi_phi=0; m_prime+r*Phi^2*m=4*pi*r^2*A^4*G*E/c^4+r^2*Phi^2/2. nu_prime=m/(r^2*b)+4*pi*r*A^4*G*Pr/(c^4*b)+r*Phi^2/2. Krr=4*pi*r*A^4*G*F/(c^5*sqrt(b)); m_t=-4*pi*r^2*N*sqrt(b)*A^4*G*F/c^4.',
        sources='Retain the full saved stellar energy/pressure profiles below the causal shell. In that shell use the actual native gas and angular photon moments at all19+512 centers, preserving the saved deep stratification between native anchors. Exterior initial radiation and gas stop at the declared outer numerical face; this is a compact Cauchy field, not an eternal luminosity.',
        inventory='Hold the coordinate baryon measure A^3*rho/sqrt(b) fixed. Thus rho/rho_ref=sqrt(b/b_ref). At every point the gas uses the isentropic Gamma-law anchored to its actual initial energy, pressure and bulk modulus. This is an explicit local constitutive continuation, not replacement by a new globally certified EOS. Native anchor checks are required before any downstream use.',
        scalar='Preserve the stored scalar values and gradients with a C1 Hermite interpolant, continued by the old exact exterior profile. Pi_phi=0 is free Cauchy data; scalar acceleration must be evaluated, not silently set to zero. Flux/extrinsic curvature is evaluated where the actual angular Cauchy field is known; unknown deeper heat flux is not certified.',
        decision='Export corrected initial material/metric variables and evaluate their scalar-stationarity defect. If the correction or that defect invalidates the old tiny-charge interpretation, do not certify or rerun the old subtraction unchanged. Continue with the physical projection/coupling indicated by the measured defect.',
        gates=dict(quadrature_metric=1e-11,inventory=1e-14,picard_density=1e-17,analytic_uniform_star=1e-11,native_anchor=.002),
        budget=dict(seconds=120,orders=[12,20],max_picard=8,CPU_threads=1,memory_GB=3,new_fluid_steps=0,native_calls=100),
        measured_basis='The existing19+512 model construction and moment inspection took9.59s including WSL startup. Constraint integration is vectorized over the saved source knots, with two fixed quadrature orders and at most8 Picard passes. No old fluid history is repeated.',
        stop='Preserve a failed projection or gate. Do not automatically increase mesh, quadrature, horizon or Picard cap. Initial constraint matching alone is not final charge, hydrostatic/scalar stationarity, full EOS certification or dynamic GR.',
        references=['https://arxiv.org/abs/gr-qc/0201064 (RGPS equations233-235)','https://arxiv.org/abs/gr-qc/9707041 (equations2.23-2.26)'],
        bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),previous.OUT/'geometry.npz',previous.OUT/'coupled-128.npz',previous.OUT/'result.json',previous.chem.prior.star.OUT/'background.npz']}))


class Data:
    def __init__(self):
        self.model=model=previous.Coupled();m=model.m;b=model.bulk;self.bg=m.bg;d=self.bg.d
        self.R=self.bg.R;self.rad=d['radius_cm'];self.rend=float(previous.feedback.prior.green.Geometry(m).metric(np.array([m.rf[-1]-m.RJ]))[3][0]*m.R)
        self.phi=CubicHermiteSpline(self.rad,d['phi']-d['phi'][-1],d['phi_prime_cm'],extrapolate=False)
        self.phis=float(d['phi'][-1]);self.geometry=geo=previous.feedback.prior.green.Geometry(m)
        self.r=np.r_[b.d['r'],m.r];_,aa,BB,re,ph=geo.metric(self.r-m.RJ);self.re=re*m.R
        self.begin=float(geo.metric(np.array([b.d['edges'][0]-m.RJ]))[3][0]*m.R)
        self.join=float(geo.metric(np.array([m.rf[0]-m.RJ]))[3][0]*m.R)
        rho,v,lt,y=model.flow.primitive(model.flow.initial)
        assert max(abs(v))==0
        p,u,*_=b.eos.gas(np.zeros(b.n),np.zeros(b.n));f=model.flow;f.eos.y=y
        pp,uu,gg,_,_=f.eos(rho,lt)
        self.rho=np.r_[b.d['rho'],rho*f.eos.rho0]
        self.p=np.r_[p,pp*f.eos.rho0*C*C]
        self.eg=np.r_[b.d['rho']*(model.cx*C*C+u),rho*(f.eos.cx+uu)*f.eos.rho0*C*C]
        # K is the dimensional dP/dlnrho at fixed entropy and inventory.
        gam=np.r_[model.mech.K[0]/p,gg]
        self.gamma=np.where(self.p>0,gam,5/3)
        self.photons=np.concatenate([b.initial,model.initial_I])
        energy=b.d['num']*b.d['Einf']
        self.er=np.einsum('iqf,q,f->i',self.photons,b.w,energy)/aa**4
        self.pr=np.einsum('iqf,q,f->i',self.photons,b.w*b.mu2,energy)/aa**4
        self.flux=C*np.einsum('iqf,q,f->i',self.photons,b.w*b.mu,energy)/aa**4
        env=previous.chem.prior.Background().env
        self.arad=float(env['Prad'][0]*3/env['T'][0]**4)
        self.indices=np.flatnonzero((self.rad>self.begin)&(self.rad<self.join))
        # Use the original stratification times a continuous anchor correction
        # in deep cells. Direct positive reconstruction serves the atmosphere.
        self.deep={};self.atmo={};n=b.n
        z=self.base(self.re[:n])
        for key,values in [('eg',self.eg),('pg',self.p),('rho',self.rho),('er',self.er),('pr',self.pr),('gamma',self.gamma)]:
            self.deep[key]=PchipInterpolator(np.r_[self.begin,self.re[:n],self.join],np.r_[0.,np.log(values[:n]/z[key]),np.log(values[n]/self.base(np.array([self.join]))[key][0])])
            use=(self.rho[n:]>0) if key in ['eg','pg','rho','gamma'] else np.ones(model.n,bool)
            self.atmo[key]=PchipInterpolator(np.r_[self.join,self.re[n:][use],self.rend],np.log(np.r_[values[n],values[n:][use],values[n:][use][-1]]))
        self.flux_interp=PchipInterpolator(self.re,self.flux,extrapolate=True)
        self.cuts=np.unique(np.r_[self.rad,self.begin,self.join,self.re,self.rend])

    def scalar(self,r):
        r=np.asarray(r);inside=r<=self.R;phi=np.empty_like(r);v=np.empty_like(r);vp=np.empty_like(r)
        phi[inside]=self.phis+self.phi(r[inside]);v[inside]=self.phi(r[inside],1);vp[inside]=self.phi(r[inside],2)
        if (~inside).any():
            rr=r[~inside];x=rr/self.R;mass,N,bb,cc,vv=self.bg.metric(x);v[~inside]=vv/self.R
            gx,gw=leg.leggauss(12);nodes=1+(x[:,None]-1)*(gx+1)/2
            phi[~inside]=self.phis+(x-1)/2*(self.bg.metric(nodes.ravel())[-1].reshape(nodes.shape)@gw)
            vp[~inside]=-2*(rr-mass*self.R)/(rr*rr*bb)*v[~inside]
        return phi,v,vp

    def base(self,r):
        d=self.bg.d;r=np.asarray(r)
        def get(k):return np.interp(r,self.rad,d[k])
        er=self.arad*np.exp(np.interp(r,self.rad,np.log(d['temperature_K'])))**4
        return dict(eg=get('energy_cgs')-er,pg=get('pressure_cgs')-er/3,er=er,pr=er/3,
            rho=np.exp(np.interp(r,self.rad,np.log(d['density_cgs']))),gamma=get('gamma1'))

    def source(self,r):
        r=np.asarray(r);z=self.base(r);deep=(r>=self.begin)&(r<self.join);atmo=r>=self.join
        for key in z:
            z[key][deep]*=np.exp(self.deep[key](r[deep]))
            z[key][atmo]=np.exp(self.atmo[key](r[atmo]))
        for key in ['eg','pg','rho']:z[key][r>self.R]=0.
        z['gamma'][r>self.R]=5/3
        return z


class Quadrature:
    def __init__(self,edges,order):
        self.edges=np.asarray(edges);self.x,self.w=leg.leggauss(order);self.h=np.diff(edges)/2
        self.r=(edges[:-1,None]+self.h[:,None]*(self.x+1))
        V=leg.legvander(self.x,order-1);Q=[]
        for j in range(order):
            co=leg.legint(np.eye(order)[j]);Q.append(leg.legval(self.x,co)-leg.legval(-1,co))
        self.Q=np.array(Q).T@np.linalg.inv(V)

    def integrate(self,f):
        f=np.asarray(f,LD).reshape(self.r.shape);full=self.h.astype(LD)*(f@self.w.astype(LD))
        prefix=np.r_[LD(0),np.cumsum(full,dtype=LD)]
        at=prefix[:-1,None]+self.h[:,None].astype(LD)*(f@self.Q.astype(LD).T)
        return at,prefix

    def sample(self,values,r):
        # Barycentric polynomial evaluation at requested points within each
        # Gauss panel. The saved panel knots prevent cross-interface smoothing.
        ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,len(self.h)-1)
        xx=(np.asarray(r)-self.edges[ids])/self.h[ids]-1
        inv=np.linalg.inv(leg.legvander(self.x,len(self.x)-1))
        co=np.asarray(values,LD)@inv.astype(LD).T
        return np.sum(co[ids]*leg.legvander(xx,len(self.x)-1),axis=1,dtype=LD)


def solve(data,order):
    q=Quadrature(data.cuts,order);r=q.r;flat=r.ravel();phi,v,vp=[z.reshape(r.shape).astype(LD) for z in data.scalar(flat)]
    A=np.exp(-2*phi*phi);s={k:z.reshape(r.shape).astype(LD) for k,z in data.source(flat).items()}
    old=data.bg.sample(flat/data.R);mref=old['m'].reshape(r.shape).astype(LD)*LD(data.R);bref=1-2*mref/r
    assert min(s['eg'].flat)>=0 and min(s['pg'].flat)>=0 and min(bref.flat)>0
    integrating,face_integrating=q.integrate(r*v*v);dm=np.zeros_like(r,dtype=LD);history=[]
    for iteration in range(8):
        b=bref-2*dm/r;logw=np.log1p(-2*dm/(r*bref))/2;w=np.exp(logw);ga=s['gamma'];gm=ga-1
        thermal=np.divide(np.expm1(gm*logw),gm,out=logw.copy(),where=abs(gm)>1e-14)
        eg=w*(s['eg']+s['pg']*thermal);pg=s['pg']*np.exp(ga*logw)
        source=4*LD(np.pi)*r*r*A**4*LD(G)/LD(C)**4*(eg+s['er'])+r*r*v*v/2
        integ,faces=q.integrate(np.exp(integrating)*source);mass=np.exp(-integrating)*integ
        newdm=mass-mref;error=float(np.max(abs((newdm-dm)/(r*bref))));history.append(error);dm=newdm
        if error<1e-17:break
    assert history[-1]<1e-17,history
    b=1-2*mass/r;mass_faces=np.exp(-face_integrating)*faces
    nr=mass/(r*r*b)+4*LD(np.pi)*r*A**4*LD(G)/LD(C)**4*(pg+s['pr'])/b+r*v*v/2
    nuint,nuf=q.integrate(nr)
    # At infinity the fixed old scalar Cauchy gradient continues. Its mass
    # correction obeys d(delta_m)/dr=-r*Phi^2*delta_m, exactly at Pi_phi=0.
    qo=Quadrature(np.linspace(0,1,17),order);zz=qo.r;rrg=data.rend/zz
    mg,ng,bg,cg,vg=data.bg.metric(rrg.ravel()/data.R);phig=(vg/data.R).reshape(zz.shape).astype(LD)
    pot=LD(data.rend)**2*phig**2/zz**3;potint,potf=qo.integrate(pot)
    old_outer=data.bg.metric(np.array([data.rend/data.R]));delta_outer=mass_faces[-1]-LD(old_outer[0][0])*LD(data.R)
    delta=delta_outer*np.exp(-(potf[-1]-potint));oldmass=mg.reshape(zz.shape).astype(LD)*LD(data.R)
    oldb=1-2*oldmass/rrg;newb=oldb-2*delta/rrg
    _,lapse_integral=qo.integrate(delta/(LD(data.rend)*oldb*newb))
    nu_outer=np.log(LD(old_outer[1][0]))-lapse_integral[-1]
    nu=nu_outer-(nuf[-1]-nuint);N=np.exp(nu)
    kappa=np.zeros_like(r,dtype=LD);known=r>=data.re[0]
    F=data.flux_interp(r[known]).astype(LD)
    kappa[known]=4*LD(np.pi)*r[known]*A[known]**4*LD(G)*F/(LD(C)**5*np.sqrt(b[known]))
    mt=-LD(C)*N*r*b*kappa
    inv=w*np.sqrt(bref/b)-1
    out={key:q.sample(val,data.re) for key,val in dict(mass=mass,delta_mass=dm,nu=nu,N=N,b=b,phi=phi,Phi=v,
        nu_prime=nr,gas_energy=eg,gas_pressure=pg,density=s['rho']*w,density_log_ratio=logw,Krr=kappa,mass_time_derivative=mt).items()}
    # Scalar acceleration is free evolution, not an Einstein constraint.
    mr=source-r*v*v*mass;lam=(mr/r-mass/r**2)/b
    trace=eg-3*pg;acc=LD(C)**2*N*N*(b*(vp+(2/r+nr-lam)*v)-4*LD(np.pi)*A**4*LD(G)/LD(C)**4*(-4*phi)*trace)
    out['scalar_acceleration']=q.sample(acc,data.re)
    out.update(radius_E=data.re,panel_radius=r,panel_mass=mass,panel_lapse=N,
        panel_scalar_acceleration=acc,panel_Krr=kappa,edges=q.edges,mass_faces=mass_faces,nu_faces=nu_outer-(nuf[-1]-nuf),
        ADM_mass=mass_faces[-1]+(LD(data.bg.M)*LD(data.R)-LD(old_outer[0][0])*LD(data.R))+delta_outer*np.expm1(-potf[-1]),
        picard_history=history,maximum_inventory_defect=float(np.max(abs(inv))),order=order)
    return out


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,previous.old.optical.timeout);signal.alarm(120)
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    data=Data();rows=[]
    for order in [12,20]:
        row=solve(data,order);np.savez_compressed(OUT/f'initial-{order}.npz',**row);rows.append(row)
    a,b=rows;error={key:float(np.max(abs(a[key]-b[key]))/max(np.max(abs(b[key])),LD('1e-100'))) for key in ['mass','nu_prime','density','N']}
    old=data.bg.sample(data.re/data.R);oldm=old['m']*data.R;oldnu=np.log(old['N']);oldb=1-2*oldm/data.re
    np.savez_compressed(OUT/'corrected-initial-state.npz',**{k:b[k] for k in ['radius_E','mass','delta_mass','nu','N','b','phi','Phi','nu_prime','gas_energy','gas_pressure','density','density_log_ratio','Krr','mass_time_derivative','scalar_acceleration']},
        reference_density=data.rho,reference_gas_energy=data.eg,reference_gas_pressure=data.p,reference_gamma=data.gamma,
        initial_photon_energy=data.er,initial_photon_radial_pressure=data.pr,initial_photon_flux=data.flux,scalar_normal_time_derivative=np.zeros_like(data.re),
        ADM_mass=b['ADM_mass'],free_scalar_Cauchy_data=True,global_native_EOS_certificate=False,full_goal_complete=False)
    result=dict(classification='Counterexample candidate',passed=max(error.values())<1e-11 and b['maximum_inventory_defect']<1e-14,
        controls=error,inventory_relative=b['maximum_inventory_defect'],picard_histories=[x['picard_history'] for x in rows],
        ADM_mass_cm=float(b['ADM_mass']),old_ADM_mass_cm=float(data.bg.M*data.R),
        maximum_lapse_log_change=float(max(abs(b['nu']-oldnu))),maximum_radial_metric_log_change=float(max(abs(np.log(oldb/b['b'])/2))),
        maximum_density_log_change=float(max(abs(b['density_log_ratio']))),maximum_Krr_per_cm=float(max(abs(b['Krr']))),
        maximum_mass_rate_cm_s=float(max(abs(b['mass_time_derivative']))),maximum_scalar_acceleration_per_s2=float(max(abs(b['scalar_acceleration']))),
        scalar_Cauchy_stationarity=False,initial_constraint_projection=True,full_native_continuum_EOS=False,deep_core_momentum_source_known=False,
        full_GR_evolution=False,final_charge_solved=False,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
