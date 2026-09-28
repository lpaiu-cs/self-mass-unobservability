def run(self,steps,name,limit=None,restart=None):
    start=time.monotonic();z=np.zeros((4,self.n));ledger=np.zeros(4,dtype=LD);discard=ledger.copy();norm=ledger.copy();t=0.;substeps=0
    history=[z.copy()];times=[t];ledgers=[ledger.copy()];discards=[discard.copy()];norms=[norm.copy()];first=0;directional=[]
    if restart:
        data=np.load(OUT/f'{restart}.npz');z=data['delta_scaled'].copy();ledger=data['ledger_scaled'].astype(LD);discard=data['discard_scaled'].astype(LD);norm=data['norm_scaled'].astype(LD)
        t=float(data['time']);first=int(data['completed']);substeps=int(data['substeps']);history=list(data['history_scaled']);times=list(data['t']);ledgers=list(data['ledgers_scaled']);discards=list(data['discards_scaled']);norms=list(data['norms_scaled'])
    for k in range(first,min(steps,limit or steps)):
        end=self.t[-1]*(k+1)/steps
        while t<end-1e-18:
            r,l,cfl=self.rhs(t,z);dt=min(end-t,cfl*64/steps)
            for attempt in range(10):
                trial=z+dt*r;r2,l2,cfl2=self.end_rhs(t+dt,trial)
                middle=.75*z+.25*(trial+dt*r2);r3,l3,cfl3=self.rhs(t+dt/2,middle)
                if dt<=min(cfl2,cfl3)*64/steps*(1+1e-10):break
                dt=min(dt/2,cfl2*64/steps,cfl3*64/steps)
            else:raise AssertionError('Material SSP time cap')
            z=(z+2*(middle+dt*r3))/3;ledger+=LD(dt)*(l.astype(LD)/6+l2.astype(LD)/6+LD(2)/3*l3.astype(LD));norm+=LD(dt)*(np.sum(abs(r),axis=1,dtype=LD)/6+np.sum(abs(r2),axis=1,dtype=LD)/6+LD(2)/3*np.sum(abs(r3),axis=1,dtype=LD))
            t+=dt;active=self.active(t);discard+=np.sum(z[:,~active],axis=1,dtype=LD);z[:,~active]=0.;substeps+=1
            assert np.isfinite(z).all() and substeps<10000
        # Save all macro endpoints; canonical17 are a subset.
        t=end;history.append(z.copy());times.append(t);ledgers.append(ledger.copy());discards.append(discard.copy());norms.append(norm.copy())
        if k in [1,steps//2-1,steps-1]:
            r1=self.rhs(t,z,1.)[0];r2=self.rhs(t,z,.5)[0]
            error=np.sum(abs(r1-r2),axis=1)/np.maximum(np.sum(abs(r2),axis=1),1.)
            directional.append(dict(time=t,relative=error.tolist()))
        balance=abs(np.sum(z,axis=1,dtype=LD)+discard-ledger)/np.maximum(norm,1.)
        assert max(balance)<1e-8,('Material conservation',balance.tolist())
        write(OUT/f'{name}-progress.json',dict(completed=k+1,steps=steps,time=t,substeps=substeps,seconds=time.monotonic()-start))
    completed=len(times)-1
    np.savez_compressed(OUT/f'{name}.npz',delta_scaled=z,ledger_scaled=ledger,discard_scaled=discard,norm_scaled=norm,time=t,completed=completed,substeps=substeps,t=times,history_scaled=history,ledgers_scaled=ledgers,discards_scaled=discards,norms_scaled=norms,radius_E=self.rE,amplitude=AMP)
    j,f,_,_=self.fields(t);p=self.point(j);q=self.point(j+1);Q=(1-f)*p['Q']+f*q['Q'];active=self.active(t)
    units=np.maximum(abs(Q),1.);units[1]=np.maximum(Q[0]*C*C,1.)
    relative=float(np.max(abs(z[:,active])*AMP/units[:,active]));directional_error=max((max(x['relative']) for x in directional),default=0.)
    result=dict(classification='Counterexample candidate',passed=bool(max(balance)<1e-8 and relative<1e-6 and directional_error<.002),steps=steps,reference=self.reference,completed=completed,substeps=substeps,
        balance_relative=np.asarray(balance,float).tolist(),directional=directional,directional_relative=directional_error,forward_probe_indicator=self.probe_error,
        maximum_true_relative_state=relative,endpoint_sum=np.sum(z,axis=1).tolist(),endpoint_L1=np.sum(abs(z),axis=1).tolist(),endpoint_sum_physical=(AMP*np.sum(z,axis=1)).tolist(),endpoint_L1_physical=(AMP*np.sum(abs(z),axis=1)).tolist(),
        discard_physical=np.asarray(AMP*discard,float).tolist(),raw_owner_calls=self.raw_calls,seconds=time.monotonic()-start,full_GR_feedback=False,final_charge_solved=False)
    write(OUT/f'{name}.json',result);return result
