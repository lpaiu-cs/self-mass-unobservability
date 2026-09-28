"""Counterexample candidate: evolving native interior mechanics/radiation.

Reuse the nonlinear H/heat/photon solve and atmosphere. Interior displacement,
mass and compression evolve online; their density/inventory and moving-frame
corrections are first order on the retained native16-cell background.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
import def_native_interior_motion as motion

prior=motion.prior;old=prior.prior;cold=prior.cold;chem=motion.chemistry
C=motion.C;write=motion.write;sha=motion.sha
OUT=motion.OUT.parent/'def-native-interior-feedback'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='39e0bae7e',
        claim='Evolve the actual interior pressure, heat, hydrogen and photons together with baryon/momentum motion, including material compression/advection and opposite photon momentum/work. Decide whether the Phase115 charge cancellation survives this feedback.',
        reuse='Existing16 native volumes,152 frequencies,8 angles,448 atmosphere path and3.434ms. Native EOS/rho derivatives from Phase115; new missing spectral derivatives only. No fine896 path or duration increase.',
        equations='Existing nonlinear unsplit radiation/H/heat stage and moving atmospheric SSP stages. Interior mechanical half stages use actual current pressure and conservative face mass/enthalpy/species fluxes. Radiation impulse and moving-frame work use the opposite moments of the same collision. Native density and advected initial-inventory gradients are retained to first order.',
        boundary='Same shared photon face. Actual atmosphere inner mass/species flux supplies the opposite deep boundary flux. The mechanical pressure boundary and full scalar/metric evolution are not yet free matched.',
        native=dict(missing_spectral_derivatives='16 cells, initial/final thermo-chemical points, rho offsets-1e-4/0/+1e-4 and frequency offsets+-1e-4',call_cap=900,seconds=45),
        pilot=dict(steps=2,paths=[64,128],seconds=30),
        production=dict(paths=[64,128],seconds=360,CPU_threads=1,memory_GB=3,forecast='Measure actual full coupled pilot seconds per step; require2x projected remainder within360s. Late cold primitive recovery and outer iterations remain unmeasured.'),
        gates=dict(native=.002,source_exchange=1e-9,energy=1e-8,baryon=1e-10,species=1e-9,time_trace=.02,time_charge=.02,linear_density=.001,linear_velocity=.00001),
        limits='First-order interior perturbative density/inventory/moving-frame closure on the original16 cells, finite H and Thomson microphysics, fixed GR metric. Numerical convergence of this model is not final physical charge or a continuum EOS/GR certificate.',
        stop='Stop on support, positivity, conservation, pilot forecast or time budget failure. Preserve failed paths. No automatic extra paths, grid, horizon or weaker gate.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(motion.__file__),Path(old.__file__),Path(cold.__file__),
            motion.OUT/'native.npz',motion.OUT/'forcing-896.npz',old.OUT/'cells-448-steps-128.npz',old.OUT/'cells-448-steps-64.npz']}))


def bank():
    assert not (OUT/'bank.npz').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(45)
    b=old.Coupled(448,8).bulk;d=b.d;z=np.load(motion.OUT/'forcing-896.npz');native=chem.old.Native(cap=900)
    h=1e-4;density=[];frequency=[];err=[]
    for it in [0,-1]:
        dr=[];df=[]
        for j in range(b.n):
            chem.setup(native,d,j);lt=np.log(d['T'][j])+z['theta'][it,j];y=d['y0'][j]*(1+z['eta'][it,j]);E=d['Einf']/d['a'][j]
            values=[]
            for x in [-h,0,h]:
                s=native.state(x,lt,y);chi,em=chem.prior.coefficients(native,s,E);values.append(np.array([chi+em,em]))
                if x==0:
                    freq=[]
                    for f in [-h,h]:
                        cc,ee=chem.prior.coefficients(native,s,E*np.exp(f));freq.append(np.array([cc+ee,ee]))
                    center=s
            base=values[1];mask=base>1e-280
            dr.append(np.divide(values[2]-values[0],2*h*base,out=np.zeros_like(base),where=mask))
            df.append(np.divide(freq[1]-freq[0],2*h*base,out=np.zeros_like(base),where=mask))
            aa,ee,*_=b.eos.radiation(z['theta'][it],z['eta'][it]);weight=d['num']*E
            for true,estimate in [(base[0],aa[j]),(base[1],ee[j])]:err.append(float(np.sum(abs(estimate-true)*weight)/max(float(true@weight),1e-280)))
        density.append(dr);frequency.append(df)
    np.savez_compressed(OUT/'bank.npz',density=density,frequency=frequency,native_calls=native.ion.calls)
    row=dict(classification='Counterexample candidate',passed=max(err)<.002,coefficient_anchor_relative=max(err),native_calls=native.ion.calls,seconds=time.monotonic()-start)
    write(OUT/'bank.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert row['passed']


class Table:
    def __init__(self,base):
        self.base=base;self.d=base.d;self.ratios=base.ratios;self.temperature_bounds=base.temperature_bounds
        self.n=base.n;self.saved=np.load(motion.OUT/'native.npz');self.spectral=np.load(OUT/'bank.npz');d=self.d
        bg=np.load(motion.OUT/'forcing-896.npz');r=d['r'];self.r=r;self.x=np.zeros(self.n);self.xi=np.zeros(self.n);self.f=0.
        rho_grad=bg['rho_prime']/d['rho'];T_grad=np.gradient(np.log(d['T']),r);y_grad=np.gradient(d['y0'],r)/d['y0']
        self.initial=base.gas(np.zeros(self.n),np.zeros(self.n));p,u,ut,uy,pt,py,ne,net,ney=self.initial
        raw=self.saved['raw'];self.nr=(raw[:,:,1,13]-raw[:,:,2,13])/(2e-4*1.66053906660e-24)
        self.inventory=np.array([bg['pgprime']-self.saved['pr'][0]*rho_grad-pt*T_grad-py*y_grad,
            bg['uprime']-self.saved['ur'][0]*rho_grad-ut*T_grad-uy*y_grad,
            np.gradient(ne,r)-self.nr[0]*rho_grad-net*T_grad-ney*y_grad])
        ab,em,at,et,ay,ey=base.radiation(np.zeros(self.n),np.zeros(self.n));self.rinventory=[]
        for k,(v,vt,vy) in enumerate([(ab,at,ay),(em,et,ey)]):
            dv=np.gradient(v,r,axis=0)-v*self.spectral['density'][0,:,k]*rho_grad[:,None]-vt*T_grad[:,None]-vy*y_grad[:,None]
            dv+=v*self.spectral['frequency'][0,:,k]*(bg['ap']/bg['a'])[:,None]
            self.rinventory.append(np.divide(dv,v,out=np.zeros_like(v),where=v>1e-280))
        self.rinventory=np.array(self.rinventory)

    def coeff(self,key):return self.saved[key][0]*(1-self.f)+self.saved[key][1]*self.f

    def gas(self,theta,eta):
        p,u,ut,uy,pt,py,ne,net,ney=self.base.gas(theta,eta)
        p=p+self.coeff('pr')*self.x-self.inventory[0]*self.xi
        u=u+self.coeff('ur')*self.x-self.inventory[1]*self.xi
        ne=ne+(self.nr[0]*(1-self.f)+self.nr[1]*self.f)*self.x-self.inventory[2]*self.xi
        assert np.all(p>0) and np.all(ne>0)
        return p,u,ut,uy,pt,py,ne,net,ney

    def radiation(self,theta,eta):
        v=list(self.base.radiation(theta,eta));dr=self.spectral['density'][0]*(1-self.f)+self.spectral['density'][1]*self.f
        for k in [0,1]:
            factor=1+dr[:,k]*self.x[:,None]-self.rinventory[k]*self.xi[:,None]
            assert factor.min()>0,'Linear native opacity domain'
            for index in [k,k+2,k+4]:v[index]=v[index]*factor
        return v

    def recover(self,u,eta,theta):
        theta=theta.copy()
        for _ in range(8):
            _,value,ut,*_=self.gas(theta,eta);dx=(u-value)/ut;theta-=dx
            if max(abs(dx))<2e-13:return theta
        raise AssertionError('Interior conservative temperature recovery')


class Bulk(old.deep.Model):
    def __init__(self,angles=8):
        super().__init__(angles);self.eos=Table(self.eos)
        self.extra=np.zeros_like(self.initial);self.extra_u=np.zeros(self.n);self.extra_y=np.zeros(self.n)

    def evaluate(self,theta,eta,xref,uref,yref,h,jacobian=True):
        return super().evaluate(theta,eta,xref+h*self.extra/self.scale,uref+h*self.extra_u,yref+h*self.extra_y,h,jacobian)


class Coupled(old.Coupled):
    def __init__(self,feedback=True):
        super().__init__(448,8);original=self.bulk;bulk=Bulk(8)
        # Reuse the exact geometry shortened to the shared atmosphere face.
        for key in ['d','W','volume','gas_weight','photon_energy_weight','area','A','stream']:setattr(bulk,key,getattr(original,key))
        self.bulk=b=bulk;self.flow.eos=cold.ColdEOS();self.spectrum=cold.ColdSpectrum();prior.tight_primitive(self.flow)
        self.mech=motion.Mechanics(448);self.h=np.zeros(b.n+1);self.j=np.zeros(b.n+1);self.mass0=b.d['rho']*b.volume
        self.mass=self.mass0.copy();self.feedback=feedback;self.mechanical_energy=0.;self.radiation_work=0.;self.boundary_species=0.;self.max_frame=0.;self.max_density=0.
        self.f0=np.load(motion.OUT/'forcing-896.npz');self.cx=float(b.d['cx']);self.maximum_inner_iterations=0;self.mechanical_steps=0

    def set_material(self,t):
        b=self.bulk;tab=b.eos;d=b.d;tab.x=(self.mass-self.mass0)/self.mass0;tab.xi=self.mech.xi@self.h;tab.f=t/old.END
        assert max(abs(tab.x))<.001,'Interior perturbative density gate'
        rho=self.mass/b.volume;b.factor=d['a']*C*rho*d['thermo'][:,4]
        b.en_scaled=d['num'][None,:]*d['Einf']*b.scale/(d['a']**4*rho)[:,None]
        b.num_scaled=d['num'][None,:]*b.scale/(d['a']**3*rho*d['thermo'][:,4]*d['y0'])[:,None]
        b.gas_weight=self.mass*d['a'];self.max_density=max(self.max_density,float(max(abs(tab.x))))

    def velocity(self):
        d=self.f0;face=self.j/(4*np.pi*d['edges']**2*d['af']*d['rho_face'])
        return (face[:-1]+face[1:])/2/C

    def kinetic(self):return .5*self.cx*self.mass*(self.velocity()*C)**2

    def energy(self,u):
        a=self.bulk.d['a'];a0=self.m.a0
        return a*(self.mass*u+self.kinetic())+(a-a0)*self.cx*C*C*self.mass

    def hydro_force(self,theta,eta):
        b=self.bulk;d=self.f0;p,u,*_=b.eos.gas(theta,eta);dp=p-d['p0']
        de=(self.mass/b.volume-d['rho'])*self.cx*C*C+(self.mass/b.volume*u-d['rho']*d['u0'])
        f=-d['af'][1:-1]/d['Bf'][1:-1]*(self.mech.grad@dp)-self.mech.gravity*(self.mech.face_average@(de+dp))+self.mech.force0
        return self.mech.pref*f

    def mechanics(self,u,theta,eta,t,dt,atmo_ledger):
        if not self.feedback:return u,theta,eta
        b=self.bulk;d=b.d;scale=self.gas_scale/C**2;dm=float(atmo_ledger[0]*scale)
        jy=float(atmo_ledger[2]*scale*self.flow.eos.nH)
        p=b.eos.gas(theta,eta)[0];j=self.j.copy();j[1:-1]+=dt/2*self.hydro_force(theta,eta);j[-1]=dm/dt
        rho_face=self.f0['rho_face'];udonor=np.where(j>=0,np.r_[u[0],u],np.r_[u,u[-1]])
        y=d['y0']*(1+eta);nh=d['thermo'][:,4];ydonor=np.where(j>=0,np.r_[nh[0]*y[0],nh*y],np.r_[nh*y,nh[-1]*y[-1]])
        species=j*ydonor;species[-1]=jy/dt
        # Evolve the same Killing total-energy flux as the atmosphere. This
        # determines compression heat after removing actual rest and kinetic
        # energy; a second independent p*dV debit would double count work.
        pface=np.interp(d['edges'],d['r'],p);af=self.f0['af'];a=d['a'];a0=self.m.a0
        energy_flux=j*((af-a0)*self.cx*C*C+af*(udonor+pface/rho_face))
        energy_flux[-1]=float(atmo_ledger[1]*self.gas_scale)/dt
        energy=self.energy(u)-dt*np.diff(energy_flux)
        number=self.mass*nh*y-dt*np.diff(species)
        before=float(self.energy(u).sum());self.h+=dt*j;self.mass=self.mass0-self.h[1:]+self.h[:-1]
        self.set_material(t+dt);eta=number/(self.mass*nh*d['y0'])-1;self.j=j.copy()
        for _ in range(3):
            u=(energy-(a-a0)*self.cx*C*C*self.mass)/a/self.mass-self.kinetic()/self.mass
            theta=b.eos.recover(u,eta,theta);new=j.copy();new[1:-1]+=dt/2*self.hydro_force(theta,eta);self.j=new
        u=(energy-(a-a0)*self.cx*C*C*self.mass)/a/self.mass-self.kinetic()/self.mass;theta=b.eos.recover(u,eta,theta)
        self.mechanical_energy+=float(self.energy(u).sum())-before;self.boundary_species+=jy;self.mechanical_steps+=1
        return u,theta,eta

    def collision(self,I,theta,eta,beta):
        b=self.bulk;d=b.d;a=d['a'];rho=self.mass/b.volume;ne=b.eos.gas(theta,eta)[6]
        ab,em,*_=b.eos.radiation(theta,eta);factor=b.factor[:,None,None]
        static=factor*(em[:,None,:]-(ab-em)[:,None,:]*I)
        sc=b.scfactor[:,None,None]*ne[:,None,None]*(np.einsum('qk,ikf->iqf',b.S,I)-I)
        if not self.feedback:return static+sc,np.zeros_like(I),np.zeros((b.n,3)),np.zeros_like(I)
        de=b.eos.spectral['frequency'][0]*(1-b.eos.f)+b.eos.spectral['frequency'][1]*b.eos.f
        boost=-beta[:,None,None]*b.mu[None,:,None]
        bound_delta=factor*boost*((em*(1+de[:,1]))[:,None,:]*(1+I)-(ab*(1+de[:,0]))[:,None,:]*I)
        kap=ne*6.6524587321e-25/rho
        moving,escape=self.scattering(I,rho,beta,kap,a)
        extra=bound_delta+moving-sc
        # Frequency-box losses are a separate radiation port, as in atmosphere.
        return static+sc+extra,extra,escape,bound_delta

    def radiate(self,xb,I,u,theta,eta,h):
        if not self.feedback:return super().radiate(xb,I,u,theta,eta,h)
        b=self.bulk;beta=self.velocity();start_j=self.j.copy();initial_K=self.kinetic();trial_K=initial_K.copy();best=None;guess=xb*b.scale;gt=theta.copy();gy=eta.copy()
        for iteration in range(4):
            total,extra,escape,bound_extra=self.collision(guess,gt,gy,beta)
            en=b.d['num']*b.d['Einf'];a=b.d['a'];rho=self.mass/b.volume
            photon_force=np.einsum('iqf,q,f->i',total,b.w*b.mu,en)/a**4+escape[:,2]/a
            work=-(trial_K-initial_K)/(b.volume*h)
            b.extra=extra
            b.extra_u=-np.einsum('iqf,q,f->i',extra,b.w,en)/(a**4*rho)-escape[:,1]/(a*rho)+work/rho
            b.extra_y=np.einsum('iqf,q,f->i',bound_extra,b.w,b.d['num'])/(a**3*rho*b.d['thermo'][:,4]*b.d['y0'])
            result=super().radiate(xb,I,u,theta,eta,h);bb,aa,uu,tt,yy,port,it,err=result
            current,_,esc,_=self.collision(bb*b.scale,tt,yy,beta)
            force=-np.einsum('iqf,q,f->i',current,b.w*b.mu,en)/(a**4*C)-esc[:,2]/(a*C)
            next_j=start_j.copy();next_j[1:-1]+=h*self.mech.pref*(self.mech.face_average@force)
            self.j=next_j;new_K=self.kinetic();self.j=(start_j+next_j)/2;new_beta=self.velocity();change=float(max(abs(new_beta-beta)));self.j=start_j
            best=result,next_j,work,escape
            if iteration>0 and change<1e-13 and max(abs(gt-tt))<1e-11:break
            beta=new_beta;guess=bb*b.scale;gt=tt;gy=yy;trial_K=new_K
        else:raise AssertionError(('Interior moving-source iteration',change))
        self.maximum_inner_iterations=max(self.maximum_inner_iterations,iteration+1);self.j=best[1];self.max_frame=max(self.max_frame,float(max(abs(beta))));assert self.max_frame<1e-5
        self.radiation_work+=float(h*np.sum(b.volume*a*best[2]));self.deep_escape=getattr(self,'deep_escape',0.)+float(h*(b.volume@best[3][:,1]))
        return best[0]

    def run(self,steps,label,stop_after=None,restart=None):
        assert not (OUT/(label+'.npz')).exists();started=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=old.END/steps
        U=f.initial.copy();I=np.stack([self.initial_I,np.zeros_like(self.initial_I)]);xb=b.initial/b.scale;u=b.u0.copy();theta=np.zeros(b.n);eta=np.zeros(b.n)
        ledger=np.zeros(6);discard=np.zeros(4);boundary=np.zeros(4);balance=0.;snapshots=[];history=[];begin=0;failure=None;ssp=0;residual=0.
        scalar_names=['mechanical_energy','radiation_work','boundary_species','max_frame','max_density','maximum_inner_iterations','mechanical_steps','deep_escape']
        if restart is not None:
            z=np.load(OUT/(restart+'.npz'));U=z['U'];I=z['I'];xb=z['bulk_I']/b.scale;u=z['u'];theta=z['theta'];eta=z['eta'];begin=int(z['completed_steps'])
            self.h=z['h'];self.j=z['j'];self.mass=self.mass0-self.h[1:]+self.h[:-1];ledger=z['ledger'];discard=z['discard'];boundary=z['boundary']
            for key in scalar_names:setattr(self,key,float(z['scalar_'+key]))
            self.set_material(begin*h);balance=float(z['balance']);ssp=int(z['ssp']);residual=float(z['residual'])
            history=list(z['history']);snapshots=[{key:z['snapshot_'+key][i] for key in ['U','I','bulk_I','u','theta','eta','h','j','mass','t']} for i in range(len(z['snapshot_t']))]
        p0=b.eos.base.gas(np.zeros(b.n),np.zeros(b.n))[0];count=steps if stop_after is None else stop_after
        def record(t):
            p=b.eos.gas(theta,eta)[0];trace=self.mass*(self.cx*C*C+u)-self.mass0*(self.cx*C*C+b.u0)-3*(p-p0)*b.volume
            nonrest=self.mass*(u-b.u0)+(self.mass-self.mass0)*b.u0-3*(p-p0)*b.volume
            history.append([t,float(nonrest@b.d['a']),float(max(abs(self.velocity()))),float(max(abs(b.eos.x))),float(np.sum(self.mass-self.mass0))])
            snapshots.append(dict(U=U.copy(),I=I.copy(),bulk_I=xb*b.scale,u=u.copy(),theta=theta.copy(),eta=eta.copy(),h=self.h.copy(),j=self.j.copy(),mass=self.mass.copy(),t=t))
        if begin==0:record(0.)
        try:
            for k in range(begin,count):
                U,I,ll,dd,ss=self.local(U,I,k*h,h/2);ledger+=ll;discard+=dd;ssp+=ss
                u,theta,eta=self.mechanics(u,theta,eta,k*h,h/2,ll)
                xb,I,u,theta,eta,port,it,err=self.radiate(xb,I,u,theta,eta,h);boundary+=port;residual=max(residual,err)
                U,I,ll,dd,ss=self.local(U,I,k*h+h/2,h/2);ledger+=ll;discard+=dd;ssp+=ss
                u,theta,eta=self.mechanics(u,theta,eta,k*h+h/2,h/2,ll)
                photon=float(np.sum((xb*b.scale-b.initial)*b.photon_energy_weight)+np.sum((I.sum(0)-self.initial_I)*self.energy_weight))
                a=b.d['a'];delta=self.mass-self.mass0
                deep=float(np.sum(a*(self.mass*(u-b.u0)+delta*b.u0+self.kinetic())+(a-m.a0)*self.cx*C*C*delta))
                atmosphere=float((np.sum((U[2]-f.initial[2])*m.vol)+discard[2])*self.gas_scale)
                expected=boundary[0]-ledger[5]-getattr(self,'deep_escape',0.)
                balance=max(balance,abs(photon+deep+atmosphere-expected))
                if (k+1)%max(1,steps//16)==0 or k+1==count:record((k+1)*h)
                if (k+1)%16==0:print(label,'STEP',k+1,'SECONDS',time.monotonic()-started,flush=True)
        except Exception as exc:failure=repr(exc)
        completed=k+1 if failure is None else k
        response=max(abs(boundary[0]),max(abs(np.asarray(history)[:,1])),1.)
        baryon=float(abs(np.sum(self.mass-self.mass0,dtype=np.longdouble)+(np.sum((U[0]-f.initial[0])*m.vol,dtype=np.longdouble)+discard[0])*self.gas_scale/C**2))
        initial_atmo=float(np.sum(f.initial[0]*m.vol)*self.gas_scale/C**2)
        # The small port gets a separate test, beyond the total initial mass.
        portmass=abs(float(ledger[0]*self.gas_scale/C**2));joint=baryon/max(portmass,1.)
        np.savez_compressed(OUT/(label+'.npz'),U=U,I=I,bulk_I=xb*b.scale,u=u,theta=theta,eta=eta,h=self.h,j=self.j,mass=self.mass,ledger=ledger,discard=discard,boundary=boundary,
            completed_steps=completed,history=history,balance=balance,ssp=ssp,residual=residual,
            **{'scalar_'+key:getattr(self,key,0.) for key in scalar_names},
            **{'snapshot_'+key:np.array([r[key] for r in snapshots]) for key in snapshots[0]})
        row=dict(classification='Counterexample candidate',passed=bool(failure is None and balance/response<1e-8 and baryon/initial_atmo<1e-10 and residual<1e-9),failure=failure,
            steps=steps,completed_steps=completed,seconds=time.monotonic()-started,energy_relative=balance/response,joint_baryon_relative_to_initial_atmosphere=baryon/initial_atmo,
            joint_baryon_relative_to_actual_port=joint,maximum_deep_residual=residual,maximum_density=self.max_density,maximum_frame_velocity=self.max_frame,
            maximum_source_iterations=self.maximum_inner_iterations,local_SSP_steps=ssp,deep_radiation_work_erg=self.radiation_work,deep_spectral_escape_erg=getattr(self,'deep_escape',0.),
            native_density_and_inventory_order=1,fixed_metric=True,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/(label+'.json'),row);print(json.dumps(row),flush=True);return row


def pilot():
    assert json.loads((OUT/'bank.json').read_text())['passed'];assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.optical.timeout);signal.alarm(30);rows=[]
    for steps in [64,128]:
        row=Coupled().run(steps,f'pilot-{steps}',2);rows.append(row)
        if not row['passed']:break
    forecast=None if len(rows)<2 or not all(r['passed'] for r in rows) else sum(r['seconds']/2*(r['steps']-2) for r in rows)
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,forecast_upper_seconds=None if forecast is None else forecast*2,
        eligible=bool(forecast is not None and forecast*2<360),seconds=time.monotonic()-start));signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
