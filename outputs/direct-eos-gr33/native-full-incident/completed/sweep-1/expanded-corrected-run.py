def run(self,steps,label,limit=None,restart=None):
    assert not (OUT/f'{label}.npz').exists();start=time.monotonic();h=self.t[-1]/steps;base_h=h;gamma=1/4;flags,front_cells,front_times=split_flags(self,steps)
    lu_by_parts={parts:[splu(sparse.eye(self.n*self.q,format='csc')-v*base_h/parts*self.A) for v in [5/12,1/4]] for parts in [1,2]}
    x=np.zeros_like(self.I[0]);g=np.zeros((self.n,2));ledger=np.zeros(2);escape=np.zeros(3);impulse=np.zeros(self.n);transfer=np.zeros((self.n,2));ports=np.zeros((2,2));port_history=[];photon_history=[];gas_history=[];transfer_history=[]
    error=0.;species_error=0.;records=[];times=[];begin=0;moment=0.;stage_weights=[];actual_edges=[0.]
    def record(t):
        physical_x=x+self.lift(t)[0]
        k=int(np.argmin(abs(self.t-t))); pmap=self.point(k)['pressure_map']
        assert abs(self.t[k]-t)<1e-18 or t==count*base_h,'Only canonical pressure snapshots'
        # Noncanonical pilot endpoint gets an interpolated pressure map.
        if abs(self.t[k]-t)>1e-18:
            j=max(0,min(np.searchsorted(self.t,t)-1,15));f=(t-self.t[j])/(self.t[j+1]-self.t[j]);pmap=(1-f)*self.point(j)['pressure_map']+f*self.point(j+1)['pressure_map']
        pressure=(np.einsum('nj,nj->n',pmap,g)+self.local(t)['pressure_source'])*self.volume
        port_history.append(ports.copy()*AMPLITUDE);photon_history.append(physical_x*self.scale*AMPLITUDE);gas_history.append(g.copy()*AMPLITUDE);transfer_history.append(transfer.copy()*AMPLITUDE)
        times.append(t);records.append(np.stack([np.sum(physical_x*self.Eweight,axis=(1,2)),g[:,0]*self.eu+np.array([np.interp(t,self.t,v) for v in self.energy_offset.T]),g[:,1]*self.nu,impulse,np.sum(abs(physical_x)*self.Eweight,axis=(1,2)),np.sum(physical_x*self.Eweight*self.model.bulk.mu2[None,:,None],axis=(1,2)),pressure])*AMPLITUDE)
    count=steps if limit is None else limit
    if restart is None:record(0.)
    else:
        z=np.load(OUT/f'{restart}.npz');previous=json.loads((OUT/f'{restart}.json').read_text());begin=previous['completed_steps']
        x=z['delta_packet_scaled_occupation']/(self.scale*AMPLITUDE);g=z['delta_material']/AMPLITUDE
        ledger=z['ledger']/AMPLITUDE;escape=z['escape']/AMPLITUDE;impulse=z['moments'][-1,3]/AMPLITUDE
        ports=z['radial_ports'][-1]/AMPLITUDE;port_history=list(z['radial_ports']);photon_history=list(z['photon_history_scaled_occupation']);gas_history=list(z['material_history']);transfer_history=list(z['collision_transfer']);transfer=z['collision_transfer'][-1]/AMPLITUDE;self.max_residual=previous['linear_relative'];self.max_iterations=previous['max_Krylov_iterations'];self.angular_times=list(z['accepted_angular_times']);self.angular=list(z['accepted_angular_luminosity']);stage_weights=list(z['accepted_angular_quadrature_weights']);actual_edges=list(z['actual_step_edges']);times=z['t'].tolist();records=list(z['moments']);error=previous['energy_balance_relative'];species_error=previous['species_balance_relative']
    if restart is not None:
        H,_,e,_,_=self.lift(begin*base_h)
        x-=H;ledger-=self.moments(H*self.scale);escape-=e/AMPLITUDE
    def stream(xx):return (self.A@xx.reshape(self.n*self.q,self.nf)).reshape(xx.shape)
    for k in range(begin,count):
        parts=2 if flags[k] else 1;h=base_h/parts;lus=lu_by_parts[parts]
        for sub in range(parts):
            t=k*base_h+sub*h
            pair,mechanical=stages(self,t,h,x,g,lus)
            (y,gy,p1,g1,es1,l1,e1),(z,gz,p2,g2,es2,l2,e2)=pair
            transfer+=h*((1-gamma)*g1+gamma*g2-mechanical)*np.stack([self.eu,self.nu],axis=-1)
            ledger+=h*np.array([-np.sum(mechanical[:,1]*self.nu),np.sum(mechanical[:,0]*self.eu)])
            ports+=h*((1-gamma)*self.boundary_ports(t+h/3,y)+gamma*self.boundary_ports(t+h,z))
            esc=h*((1-gamma)*es1+gamma*es2);escape+=esc.sum(1)+h*((1-gamma)*l1[:3]+gamma*l2[:3])
            ledger+=h*((1-gamma)*(self.port(y*self.scale)+l1[3:])+gamma*(self.port(z*self.scale)+l2[3:]))-esc[:2].sum(1)
            impulse-=h*np.sum(((1-gamma)*p1+gamma*p2)*self.Eweight*self.mu[None,:,None],axis=(1,2))+esc[2]
            x,g=z,gz;moment=max(moment,e1,e2)
            H=self.lift(t+h)[0];physical_x=x+H;physical_ledger=ledger+self.moments(H*self.scale)
            photon=self.moments(physical_x*self.scale);total=np.array([photon[0]-np.sum(g[:,1]*self.nu),photon[1]+np.sum(g[:,0]*self.eu)])
            norm=np.array([max(np.sum(abs(physical_x)*self.Nweight),np.sum(abs(g[:,1])*self.nu),1.),max(np.sum(abs(physical_x)*self.Eweight),np.sum(abs(g[:,0])*self.eu),1.)])
            defect=abs(total-physical_ledger)/np.maximum(norm,abs(physical_ledger));species_error=max(species_error,float(defect[0]));error=max(error,float(defect[1]))
            stage_weights.extend(h*RK_B);actual_edges.append(t+h)
        if (k+1)%max(1,steps//16)==0 or k+1==count:
            record((k+1)*base_h)
            path=OUT/f'{label}-checkpoint.npz';tmp=path.with_suffix('.tmp')
            with tmp.open('wb') as handle:np.savez_compressed(handle,x=x,g=g,ledger=ledger,escape=escape,impulse=impulse,completed_steps=k+1,t=times,moments=records,ports=ports,photon_history=photon_history,gas_history=gas_history,port_history=port_history,transfer_history=transfer_history,transfer=transfer,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular,accepted_angular_quadrature_weights=stage_weights,actual_step_edges=actual_edges)
            tmp.replace(path);write(OUT/f'{label}-progress.json',dict(completed_steps=k+1,planned_steps=steps,seconds=time.monotonic()-start))
    H,_,e,_,_=self.lift(count*base_h)
    x+=H;ledger+=self.moments(H*self.scale);escape+=e/AMPLITUDE
    solve_seconds=time.monotonic()-start
    np.savez_compressed(OUT/f'{label}.npz',t=times,moments=records,delta_packet_scaled_occupation=x*self.scale*AMPLITUDE,delta_material=g*AMPLITUDE,material_energy_units=self.eu,material_neutral_units=self.nu,ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r,radial_ports=port_history,photon_history_scaled_occupation=photon_history,material_history=gas_history,collision_transfer=transfer_history,energy_offset_reference=AMPLITUDE*self.energy_offset,energy_offset_t=self.t,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular,accepted_angular_quadrature_weights=stage_weights,actual_step_edges=actual_edges,split_macro_steps=flags,front_cells=front_cells,front_times=front_times,time_integrator='RadauIIA2',RK_A=RK_A,RK_b=RK_B,RK_c=RK_C)
    row=dict(classification='Counterexample candidate',reference=self.reference,steps=steps,completed_steps=count,new_steps=count-begin,seconds=time.monotonic()-start,operator_point_seconds=self.point_seconds,operator_points=self.point_count,stepping_seconds=solve_seconds-self.point_seconds,energy_balance_relative=error,species_balance_relative=species_error,linear_relative=self.max_residual,max_Krylov_iterations=self.max_iterations,frequency_moment_relative=moment,endpoint_photon_reference_energy_erg=float(np.sum(x*self.Eweight)*AMPLITUDE),endpoint_material_reference_energy_erg=float(np.sum(g[:,0]*self.eu)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(g[:,0])*self.eu)*AMPLITUDE),additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
    actual=g[:,0]*self.eu+np.array([np.interp(count*base_h,self.t,v) for v in self.energy_offset.T])
    row.update(endpoint_material_nonrest_energy_erg=row['endpoint_material_reference_energy_erg'],endpoint_material_reference_energy_erg=float(np.sum(actual)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(actual))*AMPLITUDE),velocity_jet_relative=self.velocity_jet_error,primitive_mapping_relative=self.mapping_error,conserved_material_input_applied=True)
    row.update(maximum_extended_stage_residual=max(v['extended_residual'] for v in LINEAR),maximum_refinement_steps=max(v['corrections'] for v in LINEAR),initial_GMRES_stagnations=sum(v['initial_info']!=0 for v in LINEAR))
    row.update(actual_completed_steps=len(stage_weights)//2,actual_new_steps=len(stage_weights)//2-sum(1+int(v) for v in flags[:begin]),split_completed_macro_steps=int(sum(flags[:count])),base_steps=steps)
    row['passed']=self.velocity_jet_error<1e-4 and self.mapping_error<1e-10 and max(error,species_error)<1e-8 and self.max_residual<1e-12
    write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);return row
