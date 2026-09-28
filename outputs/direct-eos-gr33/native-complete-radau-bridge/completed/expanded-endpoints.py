def source():
    initialize();pressure=FunctionType(precision.pressure_function.__code__,dict(precision.pressure_function.__globals__,OUT=OUT))()
    model=base.gr.Response();rows=[];outputs=[]
    for n in [64,128]:
        m=run.owner.Model(n);p=dict(np.load(base.saved(n)));count,mom,ports=recovered(n)
        times=p['actual_step_edges'][:count+1];assert abs(times[-1]-read(OUT/'plan.json')['original_return_horizon_seconds'])<1e-18
        d=base.retained.template()
        for key,v in list(d.items()):
            if v.ndim and len(v)==17:d[key]=np.zeros((len(times),*v.shape[1:]),dtype=v.dtype)
        d['t']=times.copy();material=m.material;coeff=model.coeff(material.rE)
        Eg0,Pg0,Kg0,Er0,Pr0=[coeff[key]*C**4/base.gr.base.G*material.V for key in ['Eg','Pg','Kg','Er','Pr']]
        contrasts=background_contrasts(material,Eg0,Pg0)
        b=model.model.bulk;weights=4*np.pi*np.r_[b.W,model.model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
        gas=[];photons=[];baryons=[];probes=[];mapping=[];references=[];ledgers=[];discard=np.zeros((m.n,4),LD)
        cumulative=np.concatenate([np.zeros((1,2,2),LD),np.cumsum(np.sum((ports*p['joint_stage_weights'][:2*count,None,None]).reshape(count,2,2,2),axis=1),axis=0)],axis=0)
        for i,t in enumerate(times):
            q=np.zeros((4,m.n),LD) if i==0 else p['joint_stage_conserved_scaled'][2*i-1].copy()
            g=restored_gas(m,q);off=~material.active(t)
            discard[off]+=g[off]*m.units[off];g[off]=0
            z=np.array([g[:,2]*m.bu,g[:,3]*m.su,g[:,0]*m.eu,g[:,1]*m.nu],LD)
            field=m.geometry(float(t))[0];vol=(3*field[0]+field[2])*AMP
            k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k])
            pp=[];qs=[]
            for j,v in [(k,1-w),(k+1,w)]:
                bank=dict(np.load(base.feedback.old.OUT/f'bank-{m.reference}/point-{j}.npz'))
                pp.append(v*np.array(pressure(material,j,z,field,bank)));qs.append(v*material.point(j)['Q'])
            p0,delta,error=sum(pp);Q=sum(qs);backgroundE=(Q[2]+material.rest*Q[0])/material.a
            nonrest=z[2]*AMP/material.a+((1-w)*contrasts[k]+w*contrasts[k+1])*vol
            pg,pr=delta*AMP+(Kg0[None]-p0)*vol
            B=z[0]*AMP;gas.append([nonrest,pg,pr]);baryons.append(B);probes.append(error*AMP)
            local=m.local(float(t));reference=(m.pressure(float(t),g)+local['pressure_source'])*m.volume*AMP
            references.append(reference);mapping.append(delta[0]*AMP-reference)
            I=(1-w)*m.I[k]+w*m.I[k+1];Ebg=np.einsum('nqf,nqf->n',I,weights)/material.a;Pbg=np.einsum('nqf,nqf,q->n',I,weights,b.mu2)/material.a
            em,pm=np.zeros((2,m.n),LD) if i==0 else mom[2*i-1,:2]
            phi=field[0]*AMP;lam=field[2]*AMP
            photons.append([em/material.a-Ebg*vol+4*Er0*phi+(Er0+Pr0)*lam,pm/material.a-Pbg*vol+4*Pr0*phi+(3*Pr0-model.ratio4*Er0)*lam])
            if i:
                expected=np.sum(p['joint_stage_weights'][:2*i,None,None]*(p['joint_native_rates_scaled'][:2*i]+p['joint_collision_rates_scaled'][:2*i]),axis=0,dtype=LD)
                actual=g*m.units+discard;ledgers.append((np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))).astype(float).tolist())
        gas=np.array(gas);photons=np.array(photons);B=np.array(baryons);rest=B*LD(m.model.cx)*LD(C)**2
        d.update(baryon_g=B,gas_nonrest_energy_erg=gas[:,0],nonrest_trace_erg=gas[:,0]-gas[:,2]-2*gas[:,1],nonrest_stress_erg=gas[:,0]-gas[:,2],pressure_volume_erg=gas[:,1],photon_energy_erg=photons[:,0],photon_radial_pressure_erg=photons[:,1],metric_stress_erg=rest+gas[:,0]-gas[:,2]+photons[:,0]-photons[:,1],inner_cumulative_energy_erg=cumulative[:,0,1],outer_cumulative_energy_erg=cumulative[:,1,1])
        norm=max(np.max(np.sum(abs(np.array(references)),axis=-1)),LD('1e-290'))
        row=dict(clock=n,accepted_recovery_steps=count,pressure_probe=float(np.max(np.sum(abs(np.array(probes)),axis=-1))/norm),pressure_mapping=float(np.max(np.sum(abs(np.array(mapping)),axis=-1))/norm),local_material_ledger=ledgers)
        write(OUT/f'source-{n}-check.json',row);np.savez_compressed(OUT/'gr'/f'source-{n}.npz',**d)
        assert row['pressure_probe']<.002 and row['pressure_mapping']<1e-12 and np.max(ledgers)<1e-8,row
        rows.append(row);outputs.append(d);del m,p;gc.collect()
    keys=['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg','metric_stress_erg']
    errors={k:aligned(*outputs,k) for k in keys}
    result=dict(classification='Counterexample candidate',passed=max(errors.values())<.02,source_time=errors,rows=rows,strict_prefix_only=True,old_longer_recovery_failure_preserved=True,final_charge_conclusion='unadjudicated')
    write(OUT/'sources.json',result);print(json.dumps(result),flush=True)
