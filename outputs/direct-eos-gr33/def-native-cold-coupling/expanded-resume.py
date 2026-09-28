def run(self,steps,label,stop_after=None):
    assert not (OUT/(label+'.npz')).exists();start=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=prior.END/steps;count=steps
    saved=np.load(prior.OUT/'cells-896-steps-128.npz');assert saved['snapshot_t'][-1]==104*h
    U=saved['snapshot_U'][-1].copy();I=saved['snapshot_I'][-1].copy();xb=saved['snapshot_bulk_I'][-1]/b.scale
    theta=saved['snapshot_theta'][-1].copy();eta=saved['snapshot_eta'][-1].copy();u=b.eos.gas(theta,eta)[1]
    refU=U.copy();refI=I.copy();refb=xb.copy();refu=u.copy()
    ledger=np.zeros(6);discard=np.zeros(4);boundary=np.zeros(4);balance=0.;gas_balance=0.;local_balance=0.;iterations=0;residual=0.;substeps=0;snapshots=[];failure=None
    keys=['t','bulk_trace','atmosphere_trace','outside_mass','maximum_speed','surface_luminosity_per_frequency']
    history=[{k:saved[k][j] for k in keys} for j in range(104)]
    np.savez_compressed(OUT/'restart104.npz',U=refU,I=refI,bulk_I=refb*b.scale,u=refu,theta=theta,eta=eta,t=104*h)
    p0=b.eos.gas(theta,eta)[0];f.eos.y=f.eos.y0;iu=f.eos(U[0],f.initial_temperature)[1];ip=f.eos(U[0],f.initial_temperature)[0]
    def record(t):
        rho,v,lt,y=f.primitive(U);p,uu,*_=f.eos(rho,lt)
        trace=-f.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*uu-3*p-(f.initial[0]*iu-3*ip)
        bulk_trace=b.d['rho']*(u-b.u0)-3*(b.eos.gas(theta,eta)[0]-p0)
        history.append(dict(t=t,bulk_trace=float(bulk_trace@(b.volume*b.d['a'])),atmosphere_trace=float(trace@m.vol*self.gas_scale),
            outside_mass=float(np.sum(U[0,m.x>=0]*m.vol[m.x>=0])*self.gas_scale/(C*C)),maximum_speed=float(max(abs(v))),
            surface_luminosity_per_frequency=4*np.pi*C*self.area[-1]*(self.w*self.mu@I.sum(0)[-1])*self.number*self.E))
    record(104*h)
    try:
        for j in range(104,count):
            U,I,ll,dd,ss=self.local(U,I,j*h,h/2);ledger+=ll;discard+=dd;substeps+=ss
            xb,I,u,theta,eta,port,it,err=self.radiate(xb,I,u,theta,eta,h);boundary+=port;iterations=max(iterations,it);residual=max(residual,err)
            U,I,ll,dd,ss=self.local(U,I,j*h+h/2,h/2);ledger+=ll;discard+=dd;substeps+=ss
            record((j+1)*h)
            photon=float(np.sum(((xb-refb)*b.scale)*b.photon_energy_weight)+np.sum((I.sum(0)-refI.sum(0))*self.energy_weight))
            material=float(b.gas_weight@(u-refu)+(np.sum((U[2]-refU[2])*m.vol)+discard[2])*self.gas_scale)
            expected=boundary[0]+ledger[1]*self.gas_scale-ledger[5];balance=max(balance,abs(photon+material-expected))
            gas_res=(np.sum((U[2]-refU[2])*m.vol)+discard[2]-ledger[1]-ledger[3])*self.gas_scale;gas_balance=max(gas_balance,abs(gas_res))
            atmosphere_res=float(np.sum((I[1]-refI[1])*self.energy_weight))+(np.sum((U[2]-refU[2])*m.vol)+discard[2]-ledger[1])*self.gas_scale-boundary[1]+ledger[5]
            local_balance=max(local_balance,abs(atmosphere_res))
            if j==108:
                elapsed=time.monotonic()-start;forecast=elapsed+(127-j)*(elapsed/5)*1.5+15
                write(OUT/'remaining-measured-budget.json',dict(first_five_seconds=elapsed,forecast_total_seconds=forecast,cap_seconds=190,eligible=forecast<190))
                assert forecast<190,'Remaining measured trajectory budget'
            if (j+1)%4==0 or j==108:print(json.dumps(dict(step=j+1,Tmin_K=float(np.exp(f.primitive(U)[2][U[0]>=f.eos.floor]).min()),seconds=time.monotonic()-start)),flush=True)
            if (j+1)%max(1,steps//16)==0 or j+1==count:snapshots.append(dict(U=U.copy(),I=I.copy(),bulk_I=xb*b.scale,theta=theta.copy(),eta=eta.copy(),t=(j+1)*h))
    except Exception as exc:failure=repr(exc)
    arrays={k:np.array([row[k] for row in history]) for k in history[0]};response=max(abs(arrays['bulk_trace']).max(),abs(arrays['atmosphere_trace']).max(),abs(boundary[0]),1.)
    gas_response=max(abs(arrays['atmosphere_trace']).max(),abs(ledger[3]*self.gas_scale),1.)
    mass=abs(np.sum((U[0]-refU[0])*m.vol)+discard[0]-ledger[0])/np.sum(f.initial[0]*m.vol)
    species=abs(np.sum((U[3]-refU[3])*m.vol)+discard[3]-ledger[2]-ledger[4])/np.sum(f.initial[3]*m.vol)
    np.savez_compressed(OUT/(label+'.npz'),**arrays,U=U,I=I,bulk_I=xb*b.scale,theta=theta,eta=eta,initial_U=f.initial,initial_I=self.initial_I,
        ledger=ledger,discard=discard,boundary_energy=boundary,**{f'snapshot_{k}':np.array([row[k] for row in snapshots]) for k in ['U','I','bulk_I','theta','eta','t']})
    result=dict(classification='Counterexample candidate',passed=bool(failure is None and balance/response<1e-8 and gas_balance/gas_response<1e-8 and local_balance/gas_response<1e-8 and mass<1e-10 and species<1e-9 and residual<1e-9),
        failure=failure,atmosphere_cells=f.atmosphere_cells,actual_fluid_cells=f.n,angles=self.q,steps=steps,completed_steps=len(history)-1,local_SSP_steps=substeps,
        seconds=time.monotonic()-start,total_energy_relative=balance/response,gas_energy_relative=gas_balance/gas_response,atmosphere_photon_energy_relative=local_balance/gas_response,baryon_relative=float(mass),species_relative=float(species),
        maximum_deep_residual=residual,maximum_deep_Newton=iterations,maximum_velocity_over_c=self.maximum_frame_speed,scatter_number_relative=self.scatter_number_error,
        spectral_escape_energy_erg=float(ledger[5]),atmosphere_trace_erg=float(arrays['atmosphere_trace'][-1]),bulk_trace_erg=float(arrays['bulk_trace'][-1]),
        atmosphere_generated_inward_energy_into_deep_erg=float(-boundary[2]),shared_interface_net_outward_photon_energy_erg=float(boundary[3]),
        photon_energy_deposited_in_moving_gas_erg=float(ledger[3]*self.gas_scale),
        restart_step=104,ledger_scope="Segment104-128 only; earlier conservation remains separately certified",coupled_photons_and_moving_atmosphere=failure is None,space_time_comparison_passed=False,full_GR_scalar_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
    write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result
