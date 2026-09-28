def run():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(40)

    rows=json.loads((prior.OUT/'production-samples.json').read_text());by={(r['it'],r['cell']):r for r in rows}
    m=Coupled();b=m.bulk;e=b.eos;native=chem.old.Native(cap=2);initial_calls=native.ion.calls
    points,z,ids=prior.inputs(m,128);coefficients=[];moments=[];checks=[]
    weight=b.d['num']*b.d['Einf'];energy_weight=b.volume[:,None,None]*b.w[None,:,None]*weight[None,None,:]/b.d['a'][:,None,None]**3
    number_weight=b.volume[:,None,None]*b.w[None,:,None]*b.d['num'][None,None,:]/b.d['a'][:,None,None]**3
    force_weight=energy_weight*b.mu[None,:,None]/(C*b.d['a'][:,None,None])
    for it,k in enumerate(ids):
        m.Pi=z['snapshot_Pi'][k].copy();m.h=z['snapshot_h'][k].copy();m.j=z['snapshot_j'][k].copy()
        m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(points[it]['t'])
        eta=z['snapshot_eta'][k];theta=e.recover(z['snapshot_u'][k],eta,z['snapshot_theta'][k]);I=z['snapshot_bulk_I'][k];beta=m.velocity()
        oldgas=e.gas(theta,eta);oldrad=e.radiation(theta,eta);base=e.base.radiation(theta+e.theta0,eta)
        nr=[];ne=[]
        for j in range(b.n):
            r=by[it,j];s=dict(r,fraction=np.array(r['fraction']))
            chi,emit=chem.prior.coefficients(native,s,b.d['Einf']/b.d['a'][j]);nr.append([chi+emit,emit]);ne.append(r['raw'][13]/AMU)
        nr=np.array(nr);ne=np.array(ne)-e.inventory[2]*e.xi
        df=e.spectral['frequency'][0]*(1-e.f)+e.spectral['frequency'][1]*e.f
        for channel in [0,1]:
            # Preserve the original additive advected-inventory term. Other
            # density and frequency effects are evaluated natively above.
            nr[:,channel]-=base[channel]*e.rinventory[channel]*e.xi[:,None]*(1+df[:,channel]*e.frequency_shift[:,None])
        assert nr.min()>=0 and ne.min()>0
        positive=[]
        for channel,field in [(0,I),(1,np.ones_like(I)),(1,I)]:
            weighted=field*energy_weight*b.factor[:,None,None]
            true=nr[:,channel,None,:]*weighted;old=oldrad[channel][:,None,:]*weighted
            positive.append(float(np.sum(abs(true-old))/max(np.sum(true),1e-300)))
        before=m.collision(I,theta,eta,beta)
        gas,rad=e.gas,e.radiation
        def native_gas(t,y):
            v=list(gas(t,y));v[6]=ne;return v
        def native_rad(t,y):
            v=list(rad(t,y));v[0]=nr[:,0];v[1]=nr[:,1];return v
        e.gas=native_gas;e.radiation=native_rad
        try:after=m.collision(I,theta,eta,beta)
        finally:e.gas=gas;e.radiation=rad
        # Independent stationary coefficient residual, plus the moving terms
        # returned by the same owner; unlike a positive-rate comparison this
        # retains absorption/emission cancellation in the net source.
        raw=b.factor[:,None,None]*((nr[:,1]-oldrad[1])[:,None,:]-(nr[:,0]-oldrad[0]-nr[:,1]+oldrad[1])[:,None,:]*I)
        ds=b.scfactor[:,None,None]*(ne-oldgas[6])[:,None,None]*(np.einsum('qk,ikf->iqf',b.S,I)-I)
        subtractive=after[0]-before[0]
        def delta_gas(t,y):
            v=list(gas(t,y));v[6]=ne-oldgas[6];return v
        def delta_rad(t,y):
            v=list(rad(t,y));v[0]=nr[:,0]-oldrad[0];v[1]=nr[:,1]-oldrad[1];return v
        e.gas=delta_gas;e.radiation=delta_rad;scatter_error=m.scatter_number_error
        try:delta=m.collision(I,theta,eta,beta)
        finally:e.gas=gas;e.radiation=rad;m.scatter_number_error=scatter_error
        defect=delta[0];estimate=raw+delta[1]+ds
        owner=float(np.sum(abs(defect-estimate)*energy_weight)/max(np.sum(abs(defect)*energy_weight),1.))
        def integrated(v):
            full,extra,escape,bound_extra=v
            # Scattering conserves photon count, including declared spectral
            # exits. Only bound-free creation changes neutral H.
            ab,em=(oldrad[:2] if v is before else ((nr[:,0]-oldrad[0],nr[:,1]-oldrad[1]) if v is delta else (nr[:,0],nr[:,1])))
            bound=b.factor[:,None,None]*(em[:,None,:]-(ab-em)[:,None,:]*I)+bound_extra
            energy=np.sum(full*energy_weight,axis=(1,2))+b.volume*escape[:,1]
            neutral=np.sum(bound*number_weight,axis=(1,2))
            force=-np.sum(full*force_weight,axis=(1,2))-b.volume*escape[:,2]/(b.d['a']*C)
            return np.array([energy,neutral,force])
        original=integrated(before);change=integrated(delta);actual=original+change;moments.append([original,actual,change])
        coefficients.append(np.stack([oldrad[0],oldrad[1],nr[:,0],nr[:,1]]))
        checks.append(dict(time=points[it]['t'],positive_energy_weighted_errors=positive,
            electron_relative=float(np.max(abs(ne/oldgas[6]-1))),owner_relative=owner))
    assert native.ion.calls==initial_calls
    times=np.array([r['t'] for r in points]);moments=np.array(moments)
    integrated=np.trapezoid(np.sum(abs(moments),axis=-1),times,axis=0)
    net_relative=integrated[2]/np.maximum(integrated[0],1.)
    maxcoef=max(max(r['positive_energy_weighted_errors']) for r in checks);owner=max(r['owner_relative'] for r in checks)
    np.savez_compressed(OUT/'native-collision.npz',t=times,coefficients=coefficients,moments=moments,
        moments_order=np.array(['Killing_photon_energy_per_s','neutral_creation_per_s','gas_force_dyn']),
        integrated_absolute=integrated)
    result=dict(classification='Counterexample candidate',passed=bool(maxcoef<.002 and max(net_relative)<.02 and owner<1e-10),
        maximum_positive_channel_relative=maxcoef,net_exchange_integrated_relative=net_relative.tolist(),owner_relative=owner,
        checks=checks,native_initialization_calls=initial_calls,new_native_state_calls=0,seconds=time.monotonic()-start,
        native_coefficients_applied_to_collision=True,constitutive_repair_required=bool(maxcoef>=.002 or max(net_relative)>=.02),
        frozen_trajectory_only=True,finite_coupled_response_bounded=False,coupled_EOS_evolution=False,
        full_source_error_enclosed=False,final_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)
