def source():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    pressure=FunctionType(prior.precision.pressure_function.__code__,dict(prior.precision.pressure_function.__globals__,OUT=OUT))()
    model=base.gr.Response();checks=[]
    for n in [64,128]:
        m=run.owner.Model(n);p=dict(np.load(base.saved(n)));count=len(p['actual_step_edges'])-1
        photon=dict(np.load(INPUT/f'recovered-{n}.npz'))
        d=dict(np.load(OUT/'gr'/f'endpoint-{n}.npz'));assert len(d['t'])==count+1
        material=m.material;cf=model.coeff(material.rE)
        Eg0,Pg0,Kg0,Er0,Pr0=[cf[key]*C**4/base.gr.base.G*material.V for key in ['Eg','Pg','Kg','Er','Pr']]
        contrasts=background_contrasts(material,Eg0,Pg0)
        b=model.model.bulk;weights=4*np.pi*np.r_[b.W,model.model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
        banks={k:dict(np.load(base.feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz')) for k in range(17)}
        probes=[];mapping=[];references=[];values={k:[] for k in KEYS};dense_errors=[];poly_errors=[];endpoint_errors=[]
        def readout(t,g,mom,ports,field_override=None,record=True):
            z=np.array([g[:,2]*m.bu,g[:,3]*m.su,g[:,0]*m.eu,g[:,1]*m.nu],LD)
            field=m.geometry(float(t))[0] if field_override is None else field_override;vol=(3*field[0]+field[2])*AMP
            k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k])
            pp=[];qs=[]
            for j,v in [(k,1-w),(k+1,w)]:
                pp.append(v*np.array(pressure(material,j,z,field,banks[j])));qs.append(v*material.point(j)['Q'])
            p0,delta,error=sum(pp);Q=sum(qs);background=(Q[2]+material.rest*Q[0])/material.a
            nonrest=z[2]*AMP/material.a+((1-w)*contrasts[k]+w*contrasts[k+1])*vol;pg,pr=delta*AMP+(Kg0[None]-p0)*vol;B=z[0]*AMP
            local=m.local(float(t));reference=(m.pressure(float(t),g)+local['pressure_source'])*m.volume*AMP
            if record:probes.append(error*AMP);references.append(reference);mapping.append(delta[0]*AMP-reference)
            I=(1-w)*m.I[k]+w*m.I[k+1];Ebg=np.einsum('nqf,nqf->n',I,weights)/material.a;Pbg=np.einsum('nqf,nqf,q->n',I,weights,b.mu2)/material.a
            phi,lam=field[0]*AMP,field[2]*AMP;em,pm=mom[:2]
            pe=em/material.a-Ebg*vol+4*Er0*phi+(Er0+Pr0)*lam
            pp=pm/material.a-Pbg*vol+4*Pr0*phi+(3*Pr0-model.ratio4*Er0)*lam
            return dict(zip(KEYS,[B,nonrest,nonrest-pr-2*pg,nonrest-pr,pg,pe,pp,B*LD(m.model.cx)*LD(C)**2+nonrest-pr+pe-pp,ports[0,1],ports[1,1]]))
        gas=lambda q:restored_gas(m,q)
        cumulative=np.zeros((2,2),LD)
        unit=np.zeros((5,m.n),LD);unit[0]=m.driver.zc['alpha'];unit[2]=m.driver.centers*m.driver.zc['Phi']
        zero=np.zeros((m.n,4),LD);zero_m=np.zeros((3,m.n),LD);zero_p=np.zeros((2,2),LD);T=d['t'][-1]
        gt=np.unique(np.r_[m.t[m.t<T],T]);gg=[readout(v,zero,zero_m,zero_p,unit,False) for v in gt]
        geometry={k:np.array([[a[k] for a in gg[:-1]],[(b[k]-a[k])/(v-u) for a,b,u,v in zip(gg[:-1],gg[1:],gt[:-1],gt[1:])]]) for k in KEYS}
        for step in range(count):
            times=p['joint_stage_times'][2*step:2*step+2];h=LD(4)*p['joint_stage_weights'][2*step+1];t=LD(p['actual_step_edges'][step])
            initial=np.zeros((m.n,4),LD) if step==0 else gas(p['joint_stage_conserved_scaled'][2*step-1]);initial[~material.active(t)]=0
            pair=np.array([gas(q) for q in p['joint_stage_conserved_scaled'][2*step:2*step+2]])
            g0=initial.copy();m0=np.zeros((3,m.n),LD) if step==0 else photon['photon_moments'][2*step-1]
            mp=photon['photon_moments'][2*step:2*step+2];port=photon['radial_ports'][2*step:2*step+2]
            def at(theta,zero_field=False):
                u=LD(theta);q=np.array([LD('1.5')*u-LD('.75')*u*u,LD('.75')*u*u-LD('.5')*u])
                return readout(t+h*u,quadratic(g0,pair,u),quadratic(m0,mp,u),cumulative+h*np.einsum('j,jab->ab',q,port),np.zeros((5,m.n),LD) if zero_field else None,not zero_field)
            rates=(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled'])[2*step:2*step+2]/m.units
            for j,u in enumerate([LD(1)/3,LD(1)]):
                q=np.array([LD('1.5')*u-LD('.75')*u*u,LD('.75')*u*u-LD('.5')*u])
                defect=pair[j]-g0-h*np.einsum('j,jnk->nk',q,rates)
                dense_errors.append((np.sum(abs(defect)*m.units,axis=0)/np.maximum(np.sum(abs(pair[j])*m.units,axis=0),LD('1e-290'))).astype(float).tolist())
            rows=[at(LD(j)/3,True) for j in range(4)];co={k:cubic([r[k] for r in rows]) for k in KEYS}
            for u in [LD(1)/6,LD(1)/2,LD(5)/6]:
                actual=at(u)
                phi=m.driver.wave(float(t+h*u),m.driver.xc)[0]/(m.driver.centers*AMP)
                prediction={k:evaluate(co[k],u)+geometry_at(geometry,gt,np.array([t+h*u]),k)[0]*phi if k not in KEYS[-2:] else evaluate(co[k],u) for k in KEYS}
                poly_errors.append({k:float(np.sum(abs(prediction[k]-actual[k]))/max(np.sum(abs(actual[k])),max(np.sum(abs(r[k])) for r in rows),LD('1e-290'))) for k in KEYS})
            for k in KEYS:values[k].append([r[k] for r in rows])
            cumulative+=np.sum(p['joint_stage_weights'][2*step:2*step+2,None,None]*port,axis=0,dtype=LD)
            closing=pair[-1].copy();closing[~material.active(p['actual_step_edges'][step+1])]=0
            end=readout(p['actual_step_edges'][step+1],closing,mp[-1],cumulative)
            endpoint_errors.append({k:float(np.sum(abs(end[k]-d[k][step+1]))/max(np.max(np.sum(abs(d[k]),axis=-1)) if d[k].ndim>1 else np.max(abs(d[k])),LD('1e-290'))) for k in KEYS})
        norm=max(np.max(np.sum(abs(np.array(references)),axis=-1)),LD('1e-290'))
        row=dict(clock=n,dense_stage_max=float(np.max(dense_errors)),polynomial_max=max(max(r.values()) for r in poly_errors),endpoint_max=max(max(r.values()) for r in endpoint_errors),
            pressure_probe=float(np.max(np.sum(abs(np.array(probes)),axis=-1))/norm),pressure_mapping=float(np.max(np.sum(abs(np.array(mapping)),axis=-1))/norm),dense_stage=dense_errors,polynomial=poly_errors,endpoint=endpoint_errors)
        row['passed']=row['dense_stage_max']<1e-12 and row['polynomial_max']<1e-12 and row['endpoint_max']<1e-12 and row['pressure_probe']<.002 and row['pressure_mapping']<1e-12
        d.update({'state_coeff_'+k:cubic(np.array(v).swapaxes(0,1)) for k,v in values.items()});d.update({'geometry_coeff_'+k:v for k,v in geometry.items()})
        d.update(geometry_times=gt,drive_x=m.driver.xc,drive_duration=m.driver.D,drive_centers=m.driver.centers,drive_radius=m.driver.r0,drive_amplitude=incident.ETA,drive_times=m.driver.times,drive_born=np.array([np.interp(m.driver.xc,m.driver.tx,v) for v in m.driver.born['U']]))
        knots,polys=coefficients(d);row['source_polynomial_intervals']=len(knots)-1
        del knots,polys;np.savez_compressed(OUT/'gr'/f'source-{n}.npz',**d);write(OUT/f'source-{n}-check.json',row)
        assert row['passed'],row;checks.append(row);del m,p;gc.collect()
    old=read(OUT/'endpoint-sources.json')
    write(OUT/'sources.json',dict(classification='Counterexample candidate',representation_controls_passed=True,source_time_passed=old['passed'],unchanged_source_time=old['source_time'],rows=checks,original_failure_preserved=True,GR_return_admitted=False,final_charge_conclusion='unadjudicated'))
