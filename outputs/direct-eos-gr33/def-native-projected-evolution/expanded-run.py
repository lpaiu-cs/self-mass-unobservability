def run(self,steps,label,stop_after=None,restart=None):
    assert not (OUT/(label+'.npz')).exists();started=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=old.END/steps
    U=f.initial.copy();I=np.stack([self.initial_I,np.zeros_like(self.initial_I)]);xb=b.initial/b.scale;u=b.u0.copy();theta=np.zeros(b.n);eta=np.zeros(b.n)
    ledger=np.zeros(6);discard=np.zeros(4);boundary=np.zeros(4);balance=0.;snapshots=[];history=[];begin=0;failure=None;ssp=0;residual=0.
    scalar_names=['mechanical_energy','radiation_work','boundary_species','max_frame','max_density','maximum_inner_iterations','mechanical_steps','deep_escape','join_mass','join_momentum','join_energy','join_neutral']
    if restart is not None:
        z=np.load(OUT/(restart+'.npz'));U=z['U'];I=z['I'];xb=z['bulk_I']/b.scale;u=z['u'];theta=z['theta'];eta=z['eta'];begin=int(z['completed_steps'])
        self.Pi=z['Pi'];self.h=z['h'];self.j=z['j'];self.mass=self.mass0-self.h[1:]+self.h[:-1];ledger=z['ledger'];discard=z['discard'];boundary=z['boundary']
        for key in scalar_names:setattr(self,key,float(z['scalar_'+key]))
        self.set_material(begin*h);balance=float(z['balance']);ssp=int(z['ssp']);residual=float(z['residual'])
        history=list(z['history']);snapshots=[{key:z['snapshot_'+key][i] for key in ['U','I','bulk_I','u','theta','eta','h','j','Pi','mass','t']} for i in range(len(z['snapshot_t']))]
    p0=self.f0['p0'];count=steps if stop_after is None else stop_after
    def record(t):
        p=b.eos.gas(theta,eta)[0];trace=self.mass*(self.cx*C*C+u)-self.mass0*(self.cx*C*C+b.u0)-3*(p-p0)*b.volume
        nonrest=self.mass*(u-b.u0)+(self.mass-self.mass0)*b.u0-3*(p-p0)*b.volume
        history.append([t,float(nonrest@b.d['a']),float(max(abs(self.velocity()))),float(max(abs(b.eos.x))),float(np.sum(self.mass-self.mass0))])
        snapshots.append(dict(U=U.copy(),I=I.copy(),bulk_I=xb*b.scale,u=u.copy(),theta=theta.copy(),eta=eta.copy(),h=self.h.copy(),j=self.j.copy(),Pi=self.Pi.copy(),mass=self.mass.copy(),t=t))
    if begin==0:record(0.)
    def save(path,completed,compressed=False):
        temporary=path.with_suffix('.tmp')
        with temporary.open('wb') as fp:
            (np.savez_compressed if compressed else np.savez)(fp,U=U,I=I,bulk_I=xb*b.scale,u=u,theta=theta,eta=eta,h=self.h,j=self.j,Pi=self.Pi,mass=self.mass,ledger=ledger,discard=discard,boundary=boundary,
                completed_steps=completed,history=history,balance=balance,ssp=ssp,residual=residual,resumable=failure is None,
                **{'scalar_'+key:getattr(self,key,0.) for key in scalar_names},
                **{'snapshot_'+key:np.array([r[key] for r in snapshots]) for key in snapshots[0]})
            fp.flush();os.fsync(fp.fileno())
        os.replace(temporary,path)
    k=begin-1
    try:
        for k in range(begin,count):
            U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h,h/2);ledger+=ll;discard+=dd;ssp+=ss
            xb,I,u,theta,eta,port,it,err=self.radiate(xb,I,u,theta,eta,h);boundary+=port;residual=max(residual,err)
            U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h+h/2,h/2);ledger+=ll;discard+=dd;ssp+=ss
            photon=float(np.sum((xb*b.scale-b.initial)*b.photon_energy_weight)+np.sum((I.sum(0)-self.initial_I)*self.energy_weight))
            a=b.d['a'];delta=-np.diff(self.h)
            deep=float(np.sum(a*(self.mass*(u-b.u0)+delta*b.u0+self.kinetic())+(a-m.a0)*self.cx*C*C*delta))
            atmosphere=float((np.sum((U[2]-f.initial[2])*m.vol)+discard[2])*self.gas_scale)
            expected=boundary[0]-ledger[5]-getattr(self,'deep_escape',0.)
            balance=max(balance,abs(photon+deep+atmosphere-expected))
            if (k+1)%max(1,steps//16)==0 or k+1==count:record((k+1)*h)
            if (k+1)%16==0:
                save(OUT/(label+'-checkpoint.npz'),k+1)
                print(label,'STEP',k+1,'SECONDS',time.monotonic()-started,flush=True)
    except Exception as exc:failure=repr(exc)
    completed=k+1 if failure is None else k
    response=max(abs(boundary[0]),max(abs(np.asarray(history)[:,1])),1.)
    baryon=float(abs(np.sum(-np.diff(self.h),dtype=np.longdouble)+(np.sum((U[0]-f.initial[0])*m.vol,dtype=np.longdouble)+discard[0])*self.gas_scale/C**2))
    initial_atmo=float(np.sum(f.initial[0]*m.vol)*self.gas_scale/C**2)
    # The small port gets a separate test, beyond the total initial mass.
    portmass=abs(float(ledger[0]*self.gas_scale/C**2));joint=baryon/max(portmass,1.)
    save(OUT/(label+'.npz'),completed,True)
    row=dict(classification='Counterexample candidate',passed=bool(failure is None and balance/response<1e-8 and baryon/initial_atmo<1e-10 and residual<1e-9),failure=failure,
        steps=steps,completed_steps=completed,seconds=time.monotonic()-started,energy_relative=balance/response,joint_baryon_relative_to_initial_atmosphere=baryon/initial_atmo,
        joint_baryon_relative_to_actual_port=joint,maximum_deep_residual=residual,maximum_density=self.max_density,maximum_frame_velocity=self.max_frame,
        maximum_source_iterations=self.maximum_inner_iterations,local_SSP_steps=ssp,deep_radiation_work_erg=self.radiation_work,deep_spectral_escape_erg=getattr(self,'deep_escape',0.),
        native_density_and_inventory_order=1,fixed_metric=True,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/(label+'.json'),row);print(json.dumps(row),flush=True);return row
