"""Counterexample candidate: one material flux between evolving interior and atmosphere.

Interior mechanics are conservative cell momenta with acoustic perturbation
fluxes on the existing native volumes. The actual shared face uses the existing
nonlinear atmospheric HLL solver and native cold EOS. No new spatial resolution
or full nonlinear interior constitutive/metric accuracy is asserted.
"""
from pathlib import Path
from types import MethodType
import inspect
import json
import os
import signal
import sys
import textwrap
import time
import numpy as np
import sympy as sp
import def_native_interior_feedback as previous

old=previous.old;C=previous.C;write=previous.write;sha=previous.sha
OUT=previous.OUT.parent/'def-native-material-join'
replace=old.interface.old.replace


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b024c099c',
        claim='Replace the copied material ghost and subsequent debit by a simultaneous shared mass, momentum, Killing energy and neutral-H flux between evolving interior and atmosphere. Determine whether the positive conditional direct charge survives.',
        model='Retain16 native interior volumes,152 frequencies,8 angles,512 actual atmospheric cells and3.434ms. Evolve interior cell momentum conservatively with acoustic perturbation fluxes; retain the analytic initial gas/gravity imbalance in a background-balanced quadrature. At the-400m interface use actual interior fractional density, log-temperature, relative neutral fraction and velocity mapped onto the saved native face background, then the existing nonlinear HLL flux. Apply this identical flux with opposite signs in each SSP stage.',
        limits='First-order density/inventory and acoustic interior mechanics; face thermochemical perturbations are constant over the last34km center-to-face gap. This is a declared spatial reconstruction, not spatial convergence or a native inventory continuum certificate. Fixed metric, finite H/Thomson microphysics and transient initial radiation remain.',
        reuse='All Phase116 native thermodynamic/spectral banks and photon solver, same atmospheric native ColdEOS and primitive inverse. Old trajectories preserved as a model comparison; no unchanged trajectory repeated.',
        budget=dict(check_seconds=25,pilot_seconds=35,production_seconds=450,readout_seconds=60,CPU_threads=1,memory_GB=3,native_calls=100,paths=[64,128]),
        forecast='Phase116 actual64/128 paths took134/176s. Joint SSP adds16-cell material stages, with unmeasured late primitive/source cost. Measure two coupled steps per path, continue their saved prefixes only if1.7x projected remainder plus10s fits450s.',
        gates=dict(energy=1e-8,baryon=1e-10,source=1e-9,native=.002,time_direct=.02,time_trace=.02,history_direct=.02,cadence_direct=.02,quadrature_direct=.02,density=.001,speed=.00001),
        stop='Stop on positivity, support, conservation, measured budget or time-readout failure; no automatic extra paths, finer grid, longer horizon or weakened gate.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),Path(old.__file__),previous.OUT/'bank.npz',previous.motion.OUT/'native.npz',previous.OUT/'coupled-128.npz']}))


# Retain the real atmospheric reconstruction, HLL flux and stage limiter.
reconstruction=textwrap.dedent(inspect.getsource(old.interface.Flow.reconstruct))
reconstruction=replace(reconstruction,'ghost_delta=delta[:,0].copy()',
    'ghost_delta=self.join_state[:3]-self.background_left[:,0]')
hydro=replace(old.hydro_source,'    rate=-C*np.diff(flux)/m.vol',
    '    flux[:,0]=self.join_flux(R[:,0],y[0])\n    rate=-C*np.diff(flux)/m.vol')
namespace=dict(old.ns,thermo=old.interface.thermo,OUT=OUT)
exec(compile(reconstruction+'\n'+hydro,__file__,'exec'),namespace)


class Coupled(previous.Coupled):
    def __init__(self):
        super().__init__();b=self.bulk;d=b.d;f=self.flow;m=self.m
        self.Pi=np.zeros(b.n);self.join_mass=0.;self.join_momentum=0.;self.join_energy=0.;self.join_neutral=0.
        self.join0=np.r_[f.background_left[:,0],f.eos.y0]
        f.join_state=self.join0.copy();f.join_flux=self.join_flux
        f.reconstruct=MethodType(namespace['reconstruct'],f);f.hydro=MethodType(namespace['rhs'],f)
        f.incoming_y=lambda y:f.join_state[3]
        self.edge=d['edges'];self.area_gas=4*np.pi*self.edge**2*self.f0['af']
        self.face_rho=self.f0['rho_face'];self.face_K=np.interp(self.edge,d['r'],self.mech.K[0])
        self.face_p=np.interp(self.edge,d['r'],self.f0['p0'])
        f.eos.y=np.array([f.eos.y0]);pp=f.eos(np.array([self.join0[0]]),np.array([self.join0[2]]))[0][0]
        self.face_p[-1]=pp*f.eos.rho0*C*C
        self.base_momentum_flux=self.area_gas*self.face_p
        self.mflux=np.zeros(4)

    def velocity(self):return self.Pi/(self.cx*self.mass*C)

    def join_flux(self,right,donor_y):
        f=self.flow;m=self.m;L=f.join_state[:,None];R=right[:,None];a=m.af[:1]
        UL,FL,tl=f.conserved(*L,a);UR,FR,tr=f.conserved(*R,a)
        sl=min(0.,((L[1]-tl[-1])/(1-L[1]*tl[-1])).item(),((R[1]-tr[-1])/(1-R[1]*tr[-1])).item())
        sr=max(0.,((L[1]+tl[-1])/(1+L[1]*tl[-1])).item(),((R[1]+tr[-1])/(1+R[1]*tr[-1])).item())
        flux=((sr*FL-sl*FR+sl*sr*(UR-UL))/(sr-sl))[:,0]*a[0]*m.area[0]
        flux[3]=flux[0]*(float(L[3,0]) if flux[0]>=0 else donor_y)
        scale=4*np.pi*m.RJ**2*f.eos.rho0*C
        self.mflux=flux*scale*np.array([1.,C,C*C,f.eos.nH])
        return flux

    def material_state(self,u,eta):
        d=self.bulk.d
        return [self.h.copy(),self.Pi.copy(),self.energy(u),self.mass*d['thermo'][:,4]*d['y0']*(1+eta)]

    def recover_material(self,state,t,theta):
        self.h,self.Pi=state[0].copy(),state[1].copy();self.mass=self.mass0-self.h[1:]+self.h[:-1]
        assert self.mass.min()>0;self.set_material(t);b=self.bulk;d=b.d;a=d['a']
        eta=state[3]/(self.mass*d['thermo'][:,4]*d['y0'])-1
        u=(state[2]-(a-self.m.a0)*self.cx*C*C*self.mass)/a/self.mass-self.kinetic()/self.mass
        theta=b.eos.recover(u,eta,theta)
        self.max_frame=max(self.max_frame,float(max(abs(self.velocity()))));assert self.max_frame<1e-5
        f=self.flow;f.join_state=self.join0.copy();f.join_state[0]*=1+b.eos.x[-1]
        f.join_state[1]=self.velocity()[-1];f.join_state[2]+=theta[-1];f.join_state[3]*=1+eta[-1]
        assert f.join_state[0]>0 and 0<f.join_state[3]<1
        return u,theta,eta

    def material_rhs(self,u,theta,eta):
        b=self.bulk;d=b.d;v=self.velocity()*C;p=b.eos.gas(theta,eta)[0];dp=p-self.f0['p0']
        vl=np.r_[0.,v];vr=np.r_[v,0.];pl=np.r_[dp[0],dp];pr=np.r_[dp,dp[-1]]
        sound=np.sqrt(self.face_K/(self.cx*self.face_rho));Z=self.cx*self.face_rho*sound
        self.deep_dt=float(.35*np.min(np.diff(self.edge)*self.f0['B']/d['a']/(np.maximum(sound[:-1],sound[1:])+abs(v))))
        vf=(vl+vr)/2+(pl-pr)/(2*Z);ps=(pl+pr)/2+Z*(vl-vr)/2
        vf[0]=0.;ps[0]=dp[0]
        mass=self.area_gas*self.face_rho*vf
        momentum=self.area_gas*(self.face_p+ps+self.cx*self.face_rho*vf*vf)
        donor_u=np.where(mass>=0,np.r_[u[0],u],np.r_[u,u[-1]])
        nh=d['thermo'][:,4];y=d['y0']*(1+eta)
        donor_y=np.where(mass>=0,np.r_[nh[0]*y[0],nh*y],np.r_[nh*y,nh[-1]*y[-1]])
        energy=mass*((self.f0['af']-self.m.a0)*self.cx*C*C+self.f0['af']*(donor_u+(self.face_p+ps)/self.face_rho))
        neutral=mass*donor_y
        mass[-1],momentum[-1],energy[-1],neutral[-1]=self.mflux
        de=(self.mass-self.mass0)/b.volume*self.cx*C*C+(self.mass*u-self.mass0*b.u0)/b.volume
        geometry=b.volume/self.f0['B']*(-de*self.f0['ap']+2*d['a']*dp/d['r'])
        # Background quadrature retains the saved nonzero gas/gravity force;
        # it does not remove the radiation-driven physical initial imbalance.
        force=-np.diff(momentum-self.base_momentum_flux)+b.volume*self.f0['initial_support']+geometry
        return [mass,force,-np.diff(energy),-np.diff(neutral)]

    def local_joint(self,U,I,u,theta,eta,t,duration):
        elapsed=0.;ledger=np.zeros(6);discard=np.zeros(4);steps=0;state=self.material_state(u,eta);m=self.m
        while elapsed<duration:
            u,theta,eta=self.recover_material(state,t+elapsed,theta)
            k,p,l,dt=self.local_rhs(U,I.sum(0),t+elapsed);first=self.material_rhs(u,theta,eta);port=self.mflux.copy()
            dt=min(dt,self.deep_dt,duration-elapsed);assert dt>0
            for _ in range(12):
                trial=U+dt*k;phot=I[1]+dt*p;assert (I[0]+phot).min()>=0
                trialstate=[a+dt*b for a,b in zip(state,first)]
                uu,tt,ee=self.recover_material(trialstate,t+elapsed+dt,theta)
                k2,p2,l2,limit=self.local_rhs(trial,I[0]+phot,t+elapsed+dt);second=self.material_rhs(uu,tt,ee);port2=self.mflux.copy()
                limit=min(limit,self.deep_dt)
                if dt<=limit*(1+1e-12):break
                dt=min(dt/2,limit)
            else:raise AssertionError('Joint SSP stage time cap')
            U=(U+trial+dt*k2)/2;I[1]=(I[1]+phot+dt*p2)/2
            state=[(a+b+dt*c)/2 for a,b,c in zip(state,trialstate,second)]
            u,theta,eta=self.recover_material(state,t+elapsed+dt,tt)
            assert I.sum(0).min()>=0
            tiny=U[0]<self.flow.eos.floor;discard+=np.sum(U[:,tiny]*m.vol[tiny],axis=1);U[:,tiny]=0
            ledger+=dt*(l+l2)/2
            impulse=dt*(port+port2)/2
            for key,value in zip(['join_mass','join_momentum','join_energy','join_neutral'],impulse):setattr(self,key,getattr(self,key)+float(value))
            elapsed+=dt;steps+=1;assert steps<20000
        self.mechanical_steps+=steps;self.j=first[0]
        return U,I,u,theta,eta,ledger,discard,steps


# Reuse the nonlinear photon solve and exact opposite moving-source work.
radiation=textwrap.dedent(inspect.getsource(previous.Coupled.radiate))
radiation=replace(radiation,'next_j=start_j.copy();next_j[1:-1]+=h*self.mech.pref*(self.mech.face_average@force)',
    'next_j=start_j+h*b.volume*force')
radiation=radiation.replace('self.j','self.Pi')
runner=textwrap.dedent(inspect.getsource(previous.Coupled.run))
runner=runner.replace("'deep_escape']","'deep_escape','join_mass','join_momentum','join_energy','join_neutral']")
runner=replace(runner,"self.h=z['h'];self.j=z['j'];", "self.Pi=z['Pi'];self.h=z['h'];self.j=z['j'];")
runner=replace(runner,"['U','I','bulk_I','u','theta','eta','h','j','mass','t']", "['U','I','bulk_I','u','theta','eta','h','j','Pi','mass','t']")
runner=replace(runner,'j=self.j.copy(),mass=', 'j=self.j.copy(),Pi=self.Pi.copy(),mass=')
runner=replace(runner,'j=self.j,mass=', 'j=self.j,Pi=self.Pi,mass=')
runner=replace(runner,'U,I,ll,dd,ss=self.local(U,I,k*h,h/2);ledger+=ll;discard+=dd;ssp+=ss\n            u,theta,eta=self.mechanics(u,theta,eta,k*h,h/2,ll)',
    'U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h,h/2);ledger+=ll;discard+=dd;ssp+=ss')
runner=replace(runner,'U,I,ll,dd,ss=self.local(U,I,k*h+h/2,h/2);ledger+=ll;discard+=dd;ssp+=ss\n            u,theta,eta=self.mechanics(u,theta,eta,k*h+h/2,h/2,ll)',
    'U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h+h/2,h/2);ledger+=ll;discard+=dd;ssp+=ss')
run_ns=dict(vars(previous),OUT=OUT)
exec(compile(radiation+'\n'+runner,__file__,'exec'),run_ns)
Coupled.radiate=run_ns['radiate'];Coupled.run=run_ns['run']


def checks():
    assert not (OUT/'checks.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(25)
    m=Coupled();f=m.flow;b=m.bulk;n=b.n;zero=np.zeros(n);U=f.initial.copy()
    state=m.material_state(b.u0,zero);u,theta,eta=m.recover_material(state,0.,zero)
    V=f.primitive(U);L,R=f.reconstruct(V,0.);base=m.join_flux(R[:,0],V[3,0]);baseport=m.mflux.copy()
    assert abs(base[0])<1e-20 and abs(base[2])<1e-20
    rates=m.material_rhs(u,theta,eta);imbalance=np.max(abs(rates[1]-b.volume*m.f0['initial_support']))/max(np.max(abs(rates[1])),1.)
    assert imbalance<1e-10
    f.join_state=m.join0.copy();f.join_state[2]+=.001;m.join_flux(R[:,0],V[3,0]);left=m.mflux.copy()
    f.join_state=m.join0.copy();right=R[:,0].copy();right[2]+=.001;m.join_flux(right,V[3,0]);rightport=m.mflux.copy()
    assert left[1]>baseport[1] and rightport[1]>baseport[1]
    # HLL equal-density states initially carry zero mass at zero velocity;
    # a pressure jump first acts on momentum. Test advective signs separately.
    f.join_state=m.join0.copy();f.join_state[1]=1e-8;m.join_flux(R[:,0],V[3,0]);outward=float(m.mflux[0])
    f.join_state=m.join0.copy();right=R[:,0].copy();right[1]=-1e-8;m.join_flux(right,V[3,0]);inward=float(m.mflux[0])
    assert outward>0 and inward<0
    # Independent native state at the actual face mapping, not the deep center.
    native=old.optical.ex.Native(cap=100);errs=[]
    for dt,ratio,x in [(0.,1.,0.),(.017,.65,1e-7),(-.01,1.2,-1e-7)]:
        rho=m.join0[0]*np.exp(x);lt=m.join0[2]+dt;y=m.join0[3]*ratio
        s=native.state(np.log(rho),lt,y);f.eos.y=np.array([y]);p,u,*_=f.eos(np.array([rho]),np.array([lt]))
        errs.append(max(abs(float(p[0])*f.eos.rho0*C*C/s['raw'][1]-1),abs(float(u[0])*C*C/s['raw'][2]-1)))
    assert max(errs)<.002
    F,h=sp.symbols('F h');assert sp.expand(-h*F+h*F)==0
    d=sp.symbols('d0:4');assert sp.expand(sum(d[i]-d[i+1] for i in range(3))-(d[0]-d[-1]))==0
    row=dict(classification='Counterexample candidate',passed=True,zero_state_physical_force_relative=float(imbalance),native_face_relative=max(errs),native_calls=native.ion.calls,
        left_heating_pressure_force_dyn=float(left[1]-baseport[1]),right_heating_pressure_force_dyn=float(rightport[1]-baseport[1]),
        left_outward_g_s=outward,right_inward_g_s=inward,shared_four_flux_opposite_sign=True,seconds=time.monotonic()-start)
    write(OUT/'checks.json',row);write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='One identical numerical face flux cancels in summed adjacent conservation equations, including both SSP RK2 stages. No continuum or EOS theorem.'))
    for name,src in [('hydro',hydro),('reconstruction',reconstruction),('radiation',radiation),('run',runner)]:(OUT/f'expanded-{name}.py').write_text(src)
    print(json.dumps(row),flush=True);signal.alarm(0)


def pilot():
    assert json.loads((OUT/'checks.json').read_text())['passed'];assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(35);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'pilot-{steps}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)!=2 or not all(r['passed'] for r in rows) else sum(r['seconds']/2*(r['steps']-2) for r in rows)
    upper=None if forecast is None else 1.7*forecast+10
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,upper_seconds=upper,eligible=bool(upper is not None and upper<450),seconds=time.monotonic()-start))
    signal.alarm(0)


def production():
    plan=json.loads((OUT/'plan.json').read_text());assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'production.json').exists()
    for p,h in plan['bindings'].items():
        if p!=str(Path(__file__)):assert sha(p)==h,p
    assert sha(__file__)==json.loads((OUT/'implementation-receipt.json').read_text())['source_sha256']
    started=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(450);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'coupled-{steps}',restart=f'pilot-{steps}');rows.append(row)
        if not row['passed']:break
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=bool(len(rows)==2 and all(r['passed'] for r in rows)),paths=rows,seconds=time.monotonic()-started,
        actual_shared_material_flux=True,full_GR_feedback=False,spatial_continuum_certified=False,final_charge_solved=False))
    signal.alarm(0)


def receipt():
    write(OUT/'implementation-receipt.json',dict(classification='Counterexample candidate',source_sha256=sha(__file__),
        preserved='first-producer.py',failure='Preflight incorrectly expected an immediate HLL mass flux from an equal-density zero-velocity thermal pressure jump. It correctly produces a momentum flux first; no trajectory had started.',
        correction='Check the pressure-force response and independent incoming/outgoing velocity signs. Explicitly retain the deep acoustic CFL and speed gate. No physical trajectory, resolution or acceptance gate was changed.',
        remaining_check_seconds=22,remaining_pilot_seconds=35,production_seconds=450))


if __name__=='__main__':globals()[sys.argv[1]]()
