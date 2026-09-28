def sources():
    assert json.loads((OUT/'production.json').read_text())['passed'];assert not (OUT/'sources.json').exists();start=time.monotonic()
    signal.signal(signal.SIGALRM,old.base.flow.old.optical.timeout);signal.alarm(75);gr=wave.Response();histories=[];allstress=[];rows=[]
    b=gr.model.bulk;model=gr.model
    weights=4*np.pi*np.r_[b.W,model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
    for steps,ref in [[64,128],[128,128]]:
        m=Material(ref,steps);d=np.load(material_path(steps,ref));ids=[int(np.argmin(abs(d['t']-t))) for t in m.t];z=d['history_scaled'][ids];histories.append(z)
        c=gr.coeff(m.rE);V=m.V;factor=C**4/wave.base.G
        Eg0,Pg0,Kg0,Er0,Pr0=[c[key]*factor*V for key in ['Eg','Pg','Kg','Er','Pr']]
        phi=m.metric['delta_u'];lam=m.metric['delta_lambda'];s=3*phi+lam
        p=np.load(photons.path(steps,ref));pid=[int(np.argmin(abs(p['t']-t))) for t in m.t];mom=p['moments'][pid]
        background=np.load(old.base.flow.OUT/f'coupled-{ref}.npz');bid=m.ids
        I=np.concatenate([background['snapshot_bulk_I'][bid],background['snapshot_I'][bid].sum(1)],axis=1)
        Ebg=np.einsum('tnqf,nqf->tn',I,weights)/m.a;Pbg=np.einsum('tnqf,nqf,q->tn',I,weights,b.mu2)/m.a
        mu4=(b.edges_mu[1:]**5-b.edges_mu[:-1]**5)/(5*np.diff(b.edges_mu));I0=np.concatenate([b.initial,model.initial_I]);E0=np.sum(I0*weights,axis=(1,2));ratio4=np.sum(I0*weights*mu4[None,:,None],axis=(1,2))/E0;R40=ratio4*Er0
        photonE=mom[:,0]/m.a-Ebg*s+4*Er0*phi+(Er0+Pr0)*lam
        photonP=mom[:,5]/m.a-Pbg*s+4*Pr0*phi+(3*Pr0-R40)*lam
        source=[];pressure_errors=[];balance=float(np.max(abs(np.sum(d['history_scaled'],axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        for k,t in enumerate(m.t):
            point=m.point(k);field=m.fields(t)[2];p0,delta,err=stress.pressure(m,k,z[k],field)
            q=point['Q'];total=(z[k,2]+m.rest*z[k,0])/m.a*AMP;backgroundE=(q[2]+m.rest*q[0])/m.a
            eF=total+(Eg0+Pg0-backgroundE)*s[k];pF=delta*AMP+(Kg0[None]-p0)*s[k]
            source.append(np.array([eF,pF[1],eF-pF[1]-2*pF[0],pF[0]]));pressure_errors.append(err*AMP)
        source=np.array(source);combined=np.concatenate([source,np.stack([photonE,photonP],axis=1)],axis=1);assert np.isfinite(combined).all()
        allstress.append(combined);pressure_error=float(np.max(np.sum(abs(np.asarray(pressure_errors)),axis=2))/max(np.max(np.sum(abs(source[:,[3,1]]),axis=2)),1.))
        # This is a measured waveform mismatch, not a contraction certificate.
        joint=mom[:,[1,2]];actual=z[:,[2,3]]*AMP;residual=np.max(np.sum(abs(joint-actual),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1.)
        base=template();rest=z[:,0].astype(LD)*LD(AMP)*LD(m.model.cx)*LD(C)**2
        base.update(baryon_g=z[:,0]*AMP,gas_nonrest_energy_erg=source[:,0].astype(LD)-rest,nonrest_trace_erg=source[:,2].astype(LD)-rest,
            nonrest_stress_erg=source[:,0].astype(LD)-source[:,1]-rest,pressure_volume_erg=source[:,3],photon_energy_erg=photonE,photon_radial_pressure_erg=photonP,metric_stress_erg=source[:,0]+photonE-source[:,1]-photonP,
            inner_cumulative_energy_erg=p['radial_ports'][pid,0,1],outer_cumulative_energy_erg=p['radial_ports'][pid,1,1])
        label=f'{steps}-reference-{ref}';np.savez_compressed(GR/f'source-{label}.npz',**base)
        np.savez_compressed(OUT/f'stress-{label}.npz',t=m.t,radius_E=m.rE,material=source,photon_energy=photonE,photon_radial_pressure=photonP)
        rows.append(dict(steps=steps,reference=ref,conservation=balance,pressure_probe=pressure_error,legacy_material_sweep_comparison=None,energy_H_waveform_residual=residual.tolist()))
    def compare(a,b):return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    comparisons=dict(time=compare(histories[0],histories[1]),stress_time=compare(allstress[0],allstress[1]))
    result=dict(classification='Counterexample candidate',passed=max(v for row in comparisons.values() for v in row)<.02 and max(max(r['conservation']/1e-8,r['pressure_probe']/.002) for r in rows)<1,
        comparisons=comparisons,paths=rows,seconds=time.monotonic()-start,returned_photon_transfer_applied_to_material=True,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(OUT/'sources.json',result);print(json.dumps(result));signal.alarm(0)
