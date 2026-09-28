def hydro_sources(pilot,selected=None,label=None):
    start=time.monotonic();m=Material(driven=False);model=m.model;f=model.flow;by={(r['kind'],r['it'],r['cell']):r for r in read(EOS/'production-samples.json')}
    import verify_native_conserved_history_charge as cached
    high=cached.CachedNative(16000);high.prefix=dict(high.prefix);low=current.support.fast_native(16000);low.prefix=dict(low.prefix)
    rows=[];durations=[];states={};new=[];maxp=0.;maxu=0.;maxrho=0.
    for r in by.values():
        if r['kind']=='atmosphere':states[(float(np.log(r['rho'])),r['lt'],r['y'])]=np.array(r['raw'])
    for p in FACE_INPUTS:
        for r in read(p):states[tuple(r['coordinates'])]=np.array(r['raw'])
    label=label or ('pilot' if pilot else 'production')
    native_calls0=high.ion.calls+low.ion.calls
    acoustic=dict(face_calls=0,center_calls=0,max_speed=0.,max_gamma_change=0.,max_deep_K_change=0.)
    baselineK=model.face_K.copy()
    # ponytail: exact-coordinate cache is local to this saved history; no EOS interpolation is inferred from it.
    class NativeEOS:
        def __init__(self,parent):self.parent=parent
        def __getattr__(self,key):return getattr(self.parent,key)
        @property
        def y(self):return self.parent.y
        @y.setter
        def y(self,value):self.parent.y=value
        def evaluate(self,rho,lt):
            nonlocal maxp,maxu,maxrho
            values=np.array(self.parent.evaluate(rho,lt));y=np.broadcast_to(self.y,rho.shape)
            for j in np.flatnonzero(rho>=self.floor):
                key=(float(np.log(rho[j])),float(lt[j]),float(y[j]))
                if key not in states:
                    n=low if rho[j]<self.parent.original.floor else high;calls=n.ion.calls;s=n.state(*key);states[key]=s['raw'].copy()
                    new.append(dict(coordinates=key,raw=s['raw'].tolist(),owner='low' if n is low else 'high',population_error=float(s['population_error']),native_calls=n.ion.calls-calls))
                    assert s['population_error']<1e-12
                raw=states[key];p=raw[1]/(self.rho0*C*C);u=raw[2]/C**2
                maxp=max(maxp,abs(p/values[0,j]-1));maxu=max(maxu,abs(u-values[1,j])/max(abs(u),1e-300));maxrho=max(maxrho,abs(raw[0]/(rho[j]*self.rho0)-1))
                assert maxp<.002 and maxu<.002 and maxrho<1e-10,(key,maxp,maxu,maxrho)
                values[0,j]=p;values[1,j]=u
            assert high.ion.calls+low.ion.calls-native_calls0<=16000
            for j in np.flatnonzero(rho>=self.floor):
                key=(float(np.log(rho[j])),float(lt[j]),float(y[j]))
                if len(rho)==m.n-m.nb:
                    acoustic['center_calls']+=1
                    continue
                gamma=lookup[key];acoustic['face_calls']+=1
                acoustic['max_gamma_change']=max(acoustic['max_gamma_change'],abs(gamma/values[2,j]-1))
                values[2,j]=gamma
                cs=np.sqrt(gamma*values[0,j]/(rho[j]*(self.cx+values[1,j])+values[0,j]))
                acoustic['max_speed']=max(acoustic['max_speed'],float(cs))
            return tuple(values)
        def __call__(self,rho,lt):return self.evaluate(rho,lt)[:5]
    indices=selected if selected is not None else ([0,8,16] if pilot else [k for k in range(17) if not (HYDRO/f'point-{k}.json').exists()])
    record_conserved(m)
    setup=time.monotonic()-start
    if not pilot and selected is None:
        p=read(HYDRO/'pilot.json');assert p['eligible']
        upper=2*max(p['point_seconds'])*len(indices)+setup+20
        write(HYDRO/'dispatch.json',dict(classification='Counterexample candidate',upper_seconds=upper,remaining=indices,eligible=upper<CAPS['hydro_production']))
        assert upper<CAPS['hydro_production']
    try:
        for k in indices:
            tick=time.monotonic();f.acoustic_records=[];p=m.point(k);before=m.raw(k,np.zeros_like(p['Q']),np.zeros((5,m.n)),0.)
            theta,eta=before[3]['theta'],before[3]['eta'];prim=before[3]['primitive'].copy();pg=np.array([by['deep',k,j]['pressure'] for j in range(m.nb)])
            nt=theta+np.array([by['deep',k,j]['lt']-by['deep',k,j]['old_lt'] for j in range(m.nb)])
            for j in np.flatnonzero(p['active'][m.nb:]):
                r=by['atmosphere',k,int(j)];prim[:,j]=[r['rho'],r['v'],r['lt'],r['y']]
            oldprimitive=f.primitive;oldeos=f.eos;gas=model.bulk.eos.gas
            model.face_K=np.interp(model.edge,model.bulk.d['r'],np.array([deep[j]['K'][k] for j in range(m.nb)]))
            acoustic['max_deep_K_change']=max(acoustic['max_deep_K_change'],float(np.max(abs(model.face_K/baselineK-1))))
            def native_gas(t,e):
                v=list(gas(t,e));v[0]=pg.copy();return v
            model.bulk.eos.gas=native_gas;f.eos=NativeEOS(oldeos);f.primitive=lambda U:tuple(prim.copy())
            f.join_state[2]+=nt[-1]-theta[-1]
            try:
                af,ag,dt,_=m.atmosphere(m.d['snapshot_U'][m.ids[k]],m.t[k]);df,dg,ddt=m.deep(m.d['snapshot_u'][m.ids[k]],nt,eta)
                factor=4*np.pi*model.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH]);af=af*factor[:,None];ag=ag*m.V[m.nb:]*f.eos.rho0*C*C
                shared=float(abs(df[1,-1]+C*model.base_momentum_flux[-1]-af[1,0])/max(abs(af[1,0]),1.));assert shared<1e-12
                af[1,0]=df[1,-1];ag[0]+=C*model.base_momentum_flux[-1]
                allshared=float(np.max(abs(df[:,-1]-af[:,0])/np.maximum(abs(af[:,0]),1.)));assert allshared<1e-12
                F=np.c_[df[:,:-1],af];G=np.r_[dg,ag]
                stable_flux,stable_rounding=difference(m,m.d['snapshot_u'][m.ids[k]],nt,eta,baselineK)
            finally:f.primitive=oldprimitive;f.eos=oldeos;model.bulk.eos.gas=gas;model.face_K=baselineK.copy()
            previous=np.load(OLDHYDRO/f'point-{k}.npz')
            dF=stable_flux
            dG=np.zeros_like(G,dtype=np.longdouble)
            dt=min(dt,.35*np.min(np.diff(model.m.rf)/(C*model.m.a/model.m.B*(np.max(abs(prim[1]))+acoustic['max_speed']+1e-100))))
            rate=-np.diff(dF,axis=1);rate[1]+=dG;ledger=dF[:,0]-dF[:,-1];ledger[1]+=dG.sum()
            balance=float(np.max(abs(rate.sum(1)-ledger)/np.maximum(np.sum(abs(rate),axis=1),1.)));assert balance<1e-8
            epsilon=np.finfo(F.dtype).eps;assert epsilon==np.finfo(np.longdouble).eps
            np.savez_compressed(HYDRO/f'point-{k}.npz',t=m.t[k],flux=dF,gravity=dG,rate=np.asarray(rate,float),ledger=np.asarray(ledger,float),dt=min(dt,ddt),arithmetic_epsilon=epsilon)
            duration=time.monotonic()-tick;durations.append(duration)
            rounding=stable_rounding
            resolution=float(np.max(rounding/np.maximum(np.sum(abs(rate),axis=1),1.)))
            row=dict(classification='Counterexample candidate',k=k,seconds=duration,shared=allshared,balance=balance,owner=p['owner_error'],flux_resolution=resolution,arithmetic_epsilon=float(epsilon),rate_L1=np.asarray(np.sum(abs(rate),axis=1),float).tolist());assert row['owner']<1e-8 and resolution<.002
            write(HYDRO/f'point-{k}.json',row);rows.append(row)
            write(HYDRO/f'{label}-native-states.json',new)
        if selected is not None:result=dict(complete_subset=True)
        elif pilot:
            upper=2*max(durations)*14+setup+20;result=dict(eligible=upper<CAPS['hydro_production'],upper_seconds=upper,point_seconds=durations)
        else:result=dict(complete=all((HYDRO/f'point-{k}.json').exists() for k in range(17)))
        result.update(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start,setup_seconds=setup,point_seconds=durations,
            new_native_states=len(new),native_calls=high.ion.calls+low.ion.calls-native_calls0,maximum_pressure_relative=maxp,maximum_energy_relative=maxu,maximum_density_relative=maxrho,native_acoustic_derivatives_certified=False)
        result['acoustic']=acoustic
        write(HYDRO/f'{label}.json',result);return result
    finally:
        write(HYDRO/f'{label}-native-states.json',new)
        write(HYDRO/f'{label}-calls.json',dict(native_calls=high.ion.calls+low.ion.calls-native_calls0,states=len(new),seconds=time.monotonic()-start))
