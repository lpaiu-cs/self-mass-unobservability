def run(self,label):
    assert not (OUT/f'{label}.json').exists();m=self.base;start=time.monotonic();U=self.initial.copy();t=0.;end=float(m.hist_t[-1]);ledger=np.zeros(5);discard=np.zeros(4);steps=0;next_dump=0.;history=[];snapshots=[]
    while t<end:
        k,l,dt=self.rhs(U,t);dt=min(dt,end-t);trial=U+dt*k
        assert min(trial[0])>-1e-13,'Density positivity'
        k2,l2,_=self.rhs(trial,t+dt);nxt=(U+trial+dt*k2)/2
        tiny=nxt[0]<self.eos.floor;discard+=np.sum(nxt[:,tiny]*m.vol[tiny],axis=1);nxt[:,tiny]=0
        U=nxt;ledger+=dt*(l+l2)/2;t+=dt;steps+=1
        assert steps<20000,'Evolution step cap'
        if t>=next_dump or t>=end:
            outside=np.sum(U[0,m.x>=0]*m.vol[m.x>=0]);energy=np.sum((U[2]-self.initial[2])*m.vol)
            history.append([t,outside,energy,*ledger,*discard]);snapshots.append(U.copy());next_dump+=end/32
            np.savez_compressed(OUT/f'{label}-progress.npz',U=U,initial=self.initial,time_seconds=t,ledger=ledger,conserved_discard=discard,
                history=history,snapshots=snapshots,steps=steps,next_dump=next_dump,seed=self.seed)
    rho,v,sigma,y=self.primitive(U);p,u,gamma,T,kap=self.eos(rho,sigma)
    ir,iv,isent=self.base.primitive(self.base.initial);iT=self.base.eos(ir,isent)[3];self.eos.y=self.eos.y0;ip,iu,*_=self.eos(ir,np.log(iT));self.eos.y=y
    trace=-self.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*u-3*p-(self.base.initial[0]*iu-3*ip)
    scale=4*np.pi*m.RJ**2*self.eos.rho0
    energy_residual=np.sum((U[2]-self.initial[2])*m.vol)+discard[2]-ledger[1]
    response=max(np.sum(m.vol*(rho*v*v+p)),abs(history[-1][2]),1e-100)
    baryon=abs(np.sum((U[0]-self.initial[0])*m.vol)+discard[0]-ledger[0])/np.sum(self.initial[0]*m.vol)
    energy_error=abs(energy_residual)/response
    entropy=(self.eos.evaluate(rho,sigma)[-1]-float(self.eos.d['s0']))/self.eos.sunit
    np.savez_compressed(OUT/f'{label}.npz',U=U,initial=self.initial,x_cm=m.x,volume=m.vol,rho=rho*self.eos.rho0,velocity_cm_s=v*C,logT=sigma,sigma=entropy,T=T,pressure=p*self.eos.rho0*C*C,
        history=history,snapshots=snapshots,conserved_discard=discard,initial_entropy=self.base.primitive(self.base.initial)[2])
    species_residual=np.sum((U[3]-self.initial[3])*m.vol)+discard[3]-ledger[3]-ledger[4]
    species_error=abs(species_residual)/max(np.sum(self.initial[3]*m.vol),1e-100)
    passed=species_error<1e-9 and self.max_absorption<.001 and baryon<1e-10 and energy_error<1e-8 and self.max_recovery<1e-8
    row=dict(classification='Counterexample candidate',passed=bool(passed),cells=self.n,steps=steps,seconds=time.monotonic()-start,
        baryon_ledger_relative=float(baryon),energy_ledger_over_response=float(energy_error),energy_ledger_residual_erg=float(energy_residual*scale*C*C),
        maximum_primitive_energy_relative=self.max_recovery,scalar_root_fallbacks=self.scalar_roots,minimum_logT=self.minimum_sigma,maximum_logT=self.maximum_sigma,
        final_minimum_sigma=float(min(entropy[rho>=self.eos.floor])),final_maximum_sigma=float(max(entropy[rho>=self.eos.floor])),
        gas_outside_original_radius_g=float(history[-1][1]*scale),integrated_trace_energy_erg=float(np.sum(trace*m.vol)*scale*C*C),
        dilute_baryon_g=float(discard[0]*scale),dilute_Killing_nonrest_energy_erg=float(discard[2]*scale*C*C),
        inner_baryon_into_layer_g=float(ledger[0]*scale),inner_Killing_nonrest_energy_into_layer_erg=float((ledger[1]-ledger[2])*scale*C*C),
        photon_energy_into_gas_erg=float(ledger[2]*scale*C*C),maximum_scattering_optical_depth=self.max_optical,
        neutral_fraction_minimum=self.min_y,maximum_speed_over_c=self.max_speed,
        species_ledger_relative=float(species_error),maximum_absorption_energy_optical_depth=self.max_absorption,
        finite_H_reactions=True,photon_Killing_energy_paired=True,
        static_frame_bath_approximation=True,full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/f'{label}.json',row);print(label,json.dumps(row),flush=True);assert passed
    return row
