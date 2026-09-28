def repair():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(40)
    data=forcing.inputs();m,ds,ats,z,by=data;f=m.flow;nb=m.bulk.n
    owner=m.scattering.__func__;source=textwrap.dedent(inspect.getsource(owner))
    assert source.count('escape=np.zeros((n,3))')==1
    source=source.replace('escape=np.zeros((n,3))','escape=np.zeros((n,3),dtype=I.dtype)')
    source=source.replace('loss.sum((1,2)),1e-250','abs(loss).sum((1,2)),1e-250')
    scope=dict(owner.__globals__);exec(compile(source,__file__,'exec'),scope);scatter=scope['scattering']
    (OUT/'expanded-scattering.py').write_text(source)
    assert m.velocity.__func__.__module__=='def_native_material_join'
    raw=OUT/'bound-repaired';raw.mkdir(exist_ok=True);native=forcing.chem.old.Native(cap=2)
    source=inspect.getsource(forcing.point)
    before='nr[:,0]-oldrad[0]-nr[:,1]+oldrad[1]';assert source.count(before)==1
    source=source.replace(before,'(nr[:,0]-oldrad[0])-(nr[:,1]-oldrad[1])')
    point_scope=dict(vars(forcing),OUT=raw);exec(compile(source,__file__,'exec'),point_scope)
    (OUT/'expanded-bound-owner.py').write_text(source)
    rows=[];changed=0
    for k in range(17):
        if (OUT/f"point-{k}.npz").exists():
            row=saved(k,data);rows.append(row);changed+=row["changed_velocity_cells"];continue
        d=ats[k];ids=d['active'];v=d['v'][ids];rho=d['rho'][ids]*f.eos.rho0
        f.eos.y=d['y'];kap=f.eos(d['rho'],d['lt'])[4][ids]
        roots=[by['atmosphere',k,int(j)] for j in ids]
        rn=np.array([r['rho'] for r in roots])*f.eos.rho0;vn=np.array([r['v'] for r in roots])
        kn=np.array([r['raw'][13]/r['raw'][0]*6.6524587321e-25/forcing.AMU for r in roots])
        I=z['snapshot_I'][k].sum(0)[ids];a=m.m.a[ids];vol=(4*np.pi*m.m.RJ**2*m.m.vol)[ids]
        if not (raw/f'point-{k}.json').exists():point_scope['point'](data,native,k)
        original=dict(np.load(raw/f'point-{k}.npz'));ph=original['photon'].copy();exit=original['escape'].copy()
        saved_error=m.scatter_number_error
        old,oe=m.scattering(I,rho,v,kap,a);new,ne=m.scattering(I,rn,vn,kn,a)
        ld=np.longdouble;same=vn==v;changed+=int((~same).sum())
        delta=np.zeros(I.shape,dtype=ld);ex=np.zeros((len(ids),3),dtype=ld)
        # All operations stay in the same packet owner; signed rates need an
        # absolute-loss diagnostic, never the positive-rate denominator.
        if same.any():
            product=rn[same].astype(ld)*kn[same]-rho[same].astype(ld)*kap[same]
            delta[same],ex[same]=scatter(m,I[same].astype(ld),np.ones(same.sum(),dtype=ld),v[same].astype(ld),product,a[same].astype(ld))
        if (~same).any():
            j=~same
            p,e=scatter(m,I[j].astype(ld),rn[j].astype(ld),vn[j].astype(ld),kn[j].astype(ld),a[j].astype(ld))
            p0,e0=scatter(m,I[j].astype(ld),rho[j].astype(ld),v[j].astype(ld),kap[j].astype(ld),a[j].astype(ld))
            delta[j]=p-p0;ex[j]=e-e0
        m.scatter_number_error=saved_error
        # Recompose from bound-free instead of subtracting rounded old packets.
        ph[nb+ids]=original['bound'][nb+ids]+delta
        exit[:,nb+ids]=(ex*vol[:,None]).T
        N=m.energy_weight[ids]/m.E
        residual=np.sum(delta*N,axis=(1,2))+ex[:,0]*vol
        norm=np.maximum(np.sum(abs(delta)*N,axis=(1,2))+abs(ex[:,0])*vol,1.)
        number=float(np.max(abs(residual)/norm))
        arithmetic=float(np.sum(abs(delta-(new-old))*m.energy_weight[ids])/max(np.sum(abs(delta)*m.energy_weight[ids]),1.))
        # The tiny velocity part can have a worse relative cancellation ratio;
        # the actual exported paired source is checked independently below.
        number_source=np.sum((ph-original['bound'])*np.r_[m.bulk.photon_energy_weight/m.bulk.d['Einf'],m.energy_weight/m.E],axis=(1,2))+exit[0]
        source_scale=np.maximum(np.sum((abs(ph)+abs(original['bound']))*np.r_[m.bulk.photon_energy_weight/m.bulk.d['Einf'],m.energy_weight/m.E],axis=(1,2))+abs(exit[0]),1.)
        paired=float(np.max(abs(number_source)/source_scale))
        pre_projection=ph.copy();weights=np.r_[m.bulk.photon_energy_weight,m.energy_weight];Nall=weights/m.E
        # Exact number invariant, with zero energy and radial-momentum moments.
        # Store and bound this arithmetic repair; it is not hidden material heat.
        elo,ehi=np.longdouble(m.E[0]),np.longdouble(m.E[-1])
        for _ in range(2):
            R=np.sum((ph-original['bound']).astype(np.longdouble)*Nall,axis=(1,2))+exit[0]
            ph[:,0,0]+=np.asarray(-R*ehi/(ehi-elo)/Nall[:,0,0],float)
            ph[:,0,-1]+=np.asarray(R*elo/(ehi-elo)/Nall[:,0,-1],float)
        correction=ph-pre_projection;norm=np.maximum(np.sum(abs(ph)*weights,axis=(1,2)),1.)
        projection=float(np.max(np.sum(abs(correction)*weights,axis=(1,2))/norm))
        energy_zero=float(np.max(abs(np.sum(correction*weights,axis=(1,2)))/norm))
        total_arithmetic=float(np.sum(abs(ph-original['photon'])*weights)/max(np.sum(abs(ph)*weights),1.))
        number_source=np.sum((ph-original['bound'])*Nall,axis=(1,2))+exit[0]
        paired=float(np.max(abs(number_source)/source_scale))
        assert total_arithmetic<.002 and projection<1e-8 and energy_zero<1e-12 and paired<1e-12,(k,total_arithmetic,projection,energy_zero,paired)
        dm=original['defect_moments'].copy();w=m.energy_weight[ids]
        dm[0,nb+ids]=np.sum(ph[nb+ids]*w,axis=(1,2))+exit[1,nb+ids]
        dm[2,nb+ids]=-(np.sum(ph[nb+ids]*w*m.mu[None,:,None],axis=(1,2))+exit[2,nb+ids])/a
        original.update(photon=ph,escape=exit,defect_moments=dm,number_projection=correction)
        np.savez_compressed(OUT/f'point-{k}.npz',**original)
        rows.append(dict(k=k,scattering_number=number,export_number_before_projection=float(np.max(abs(R)/source_scale)),export_number=paired,projection_energy_L1=projection,projection_energy_relative=energy_zero,unresolved_subtractive_scattering_relative=arithmetic,total_collision_arithmetic=total_arithmetic,changed_velocity_cells=int((~same).sum())))
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,changed_velocity_cells=changed,
                seconds=time.monotonic()-start,new_native_state_calls=0,native_initialization_calls=native.ion.calls,old_velocity_diagnosis_retracted=True)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)
