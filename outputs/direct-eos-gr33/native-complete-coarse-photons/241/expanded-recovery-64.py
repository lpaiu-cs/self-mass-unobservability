def run(n):
    assert read(seed_folder(n)/f'recovered-{n}.json')['passed']
    initialize=FunctionType(base.prior.initialize.__code__,dict(base.prior.initialize.__globals__,OUT=OUT));initialize()
    z=dict(np.load(saved(n)));m=interval_model(z,n,int(np.load(seed_path(n))['step']))
    target=float(z['actual_step_edges'][117]);edges=z['actual_step_edges'];stop=int(np.argmin(abs(edges-target)))
    assert stop>0 and abs(edges[stop]-target)<1e-18
    x=np.asarray(z['photon_history_scaled_occupation'][0]/(m.scale*AMP),LD);assert not np.any(x)
    zero=np.zeros((m.n,4),LD);zeroJ=sparse.csr_matrix((4*m.n,4*m.n));timings=[]
    begin,x,logs,moments,collisions,packets,ports,seed_error=seed(m,z,n)
    def stream(v):return (m.A@v.reshape(m.n*m.q,m.nf)).reshape(v.shape)
    def relative(v,reference):return float(np.sum(abs(v-reference),dtype=LD)/max(np.sum(abs(reference),dtype=LD),LD('1e-290')))
    for step in range(begin,stop):
        if step!=begin and step in interval_starts(z):
            units=m.units.copy();scale=np.array(m.scale,copy=True);del m;gc.collect();m=interval_model(z,n,step)
            assert np.array_equal(m.units,units) and np.array_equal(m.scale,scale)
        started=time.monotonic();t=edges[step];h=4*z['joint_stage_weights'][2*step+1];times=z['joint_stage_times'][2*step:2*step+2]
        assert np.array_equal(h*B,z['joint_stage_weights'][2*step:2*step+2])
        checkpoint(n,x=x,step=step,t=t,moments=moments,collisions=collisions,packets=packets,ports=ports,logs=np.array(json.dumps(logs)))
        ids=np.arange(2*step,2*step+2);assert np.max(abs(z['joint_stage_times'][ids]-times))<1e-18
        gas=[]
        for q in z['joint_stage_conserved_scaled'][ids]:gas.append(restored_gas(m,q))
        cs=[m.local(now) for now in times];ss=[m.source(now) for now in times]
        force=np.array([m.collision(c,np.zeros_like(x),g,True)[0]+s[0]/(m.scale*AMP) for c,g,s in zip(cs,gas,ss)])
        rhs=(x[None]+h*np.einsum('ij,jnqf->inqf',A,force)).ravel()
        shape=(2,*x.shape)
        def mat(value):
            v=value.reshape(shape);rates=np.array([stream(vv)+m.collision(c,vv,zero)[0] for c,vv in zip(cs,v)])
            return (v-h*np.einsum('ij,jnqf->inqf',A,rates)).ravel()
        inv=[];lus=[]
        for j,c in enumerate(cs):
            # Reuse the existing photon rank-two angular inverse, with gas
            # feedback disabled only in this photon-only preconditioner.
            bare=dict(c,**{k:np.zeros_like(c[k]) for k in ['B','Bb','Be']})
            inv.append(m.inverse_pair(bare,h*A[j,j],zeroJ))
            lus.append(splu(sparse.eye(m.n*m.q,format='csc')-h*A[j,j]*m.A))
        def pre(value):
            rows=[]
            for j,v in enumerate(value.reshape(shape)):
                streamed=lus[j].solve(np.asarray(v,float).reshape(m.n*m.q,m.nf)).reshape(x.shape)
                rows.append(m.unpack(inv[j](streamed,np.zeros((m.n,4))))[0])
            return np.asarray(rows,float).ravel()
        op=LinearOperator((len(rhs),)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
        sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]
        if step==begin:
            proposal=np.load(OLD/'rejected-original-64.npz');assert np.array_equal(x,proposal['x_initial'])
            sol=proposal['photon_stage_solution'].copy().ravel()
        if n==128 and step==begin and not RESUME:
            proposal=dict(np.load(OLD/'rejected-original-128.npz'));assert int(proposal['step'])==step and np.array_equal(proposal['x_initial'],x)
            sol=proposal['photon_stage_solution'].copy().ravel()
        for attempt in range(12):
            residual=rhs-op.matvec(sol);err=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),LD('1e-290')))
            r=residual.reshape(shape);s=sol.reshape(shape);b=rhs.reshape(shape)
            physical=[float(np.sum(abs(r)*w)/max(np.sum(abs(b)*w),np.sum(abs(s)*w),1.)) for w in [m.Nweight,m.Eweight]]
            if err<1e-14 and max(physical)<1e-13:break
            if attempt>=11:
                np.savez_compressed(OUT/f'rejected-linear-{n}.npz',step=step,time=t,step_size=h,x_initial=x,solution=sol,rhs=rhs,residual=residual)
                raise AssertionError(('Conditional photon residual',err,physical))
            history=[];start=time.monotonic()
            delta,info=gmres(op,np.asarray(residual,float),M=P,rtol=1e-12,atol=0.,restart=40,maxiter=10,callback=history.append,callback_type='pr_norm')
            sol+=np.asarray(delta,LD);calls.append(dict(info=int(info),iterations=len(history),seconds=time.monotonic()-start))
            write(OUT/f'linear-{n}.json',dict(step=step,calls=calls))
        pairs=sol.reshape(shape);errors=[]
        for j,(xx,g,c) in enumerate(zip(pairs,gas,cs)):
            _,q,*_=m.collision(c,xx,g,True);actual=q*m.units;stored=z['joint_collision_rates_scaled'][ids[j]]
            errors.append([relative(actual[:,k],stored[:,k]) for k in range(4)])
            collisions.append(actual)
            moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*AMP)
            ports.append(m.boundary_ports(float(times[j]),xx)*AMP);packets.append(m.angular[-1])
        packet_error=relative(np.array(packets[-2:]),z['accepted_angular_luminosity'][ids])
        row=dict(step=step,linear_relative=err,physical_relative=physical,collision_relative=errors,packet_relative=packet_error,calls=calls)
        row['original_joint_equation']=original_equation(m,z,step,t,h,x,gas,pairs,times,cs,ss)
        row['strict_archival_rate_identity_passed']=bool(np.max(errors)<1e-12)
        row['snapshot']=snapshot(n,m,z,step,pairs[-1],collisions,ports,packets)
        logs.append(row);write(OUT/f'progress-{n}.json',dict(completed=step+1,total=stop,last=row))
        if packet_error>=1e-12:
            np.savez_compressed(OUT/f'rejected-{n}.npz',x_initial=x,photon_stage_solution=pairs,gas=gas,time=t,step_size=h,stage_times=times,collision_rates=collisions[-2:],photon_moments=moments[-2:],ports=ports[-2:],angular=packets[-2:])
            raise AssertionError(row)
        x=pairs[-1].copy();timings.append(time.monotonic()-started)
        checkpoint(n,x=x,step=step+1,t=edges[step+1],moments=moments,collisions=collisions,packets=packets,ports=ports,logs=np.array(json.dumps(logs)))
    expected=np.load(ANCHOR117)['x']*m.scale*AMP;actual=x*m.scale*AMP
    endpoint=relative(actual,expected)
    result=dict(classification='Counterexample candidate',passed=endpoint<1e-12,clock=n,recovered_steps=stop,
        horizon_seconds=target,endpoint_relative=endpoint,rows=logs,step_seconds=timings,
        same_material_history_unchanged=True,new_physical_steps=0,reused_steps=begin,new_conditional_steps=stop-begin,seed_endpoint_relative=seed_error,strict_archival_identity_admitted=False,original_equation_audited_steps=list(range(begin,stop)),final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/f'recovered-{n}.npz',times=z['joint_stage_times'][:2*stop],weights=z['joint_stage_weights'][:2*stop],
        photon_moments=np.array(moments),collision_rates=np.array(collisions),angular=np.array(packets),radial_ports=np.array(ports),endpoint_occupation=actual)
    write(OUT/f'recovered-{n}.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);assert result['passed'],result
