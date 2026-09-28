def block(sweep):
    plan=read(OUT/'block-plan.json');assert plan['sweep']==sweep
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    run.initialize(sweep);rows=[];gamma=1-1/np.sqrt(2);started=time.monotonic()
    for n in [64,128]:
        setup=time.monotonic();m=run.coupled.Response(n);photon,material=run.paths(sweep)
        p=np.load(photon/f'steps-{n}-reference-128.npz');d=np.load(material/f'steps-{n}-reference-128.npz')
        ids=[int(np.argmin(abs(d['t']-t))) for t in m.t];motion=d['history_scaled'][ids].copy()
        offset=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*LD(run.C)**2*motion[:,0];motion[:,2]-=offset
        new=dict(motion=np.asarray(motion,float),energy_offset=np.asarray(offset,float),
            mechanical=np.asarray(motion[:,[2,3]]-p['collision_transfer'].astype(LD).transpose(0,2,1)/LD(run.AMP),float),
            xi=np.array([m.material.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:m.nb])] for z in motion]))
        old={k:getattr(m,k) for k in new}
        def relative(a,b):return float(np.max(np.sum(abs(a-b),axis=-1))/max(np.max(np.sum(abs(b),axis=-1)),1e-290))
        inputs={k:relative(old[k],new[k]) for k in ['xi','energy_offset']}
        for j,name in enumerate(['B','S']):inputs[name]=relative(old['motion'][:,j],new['motion'][:,j])
        for j,name in enumerate(['Etilde','H']):
            inputs['M_'+name]=relative(old['mechanical'][:,j],new['mechanical'][:,j])
            inputs['dM_'+name]=relative(np.diff(old['mechanical'][:,j],axis=0),np.diff(new['mechanical'][:,j],axis=0))
        def select(v):
            for k,w in v.items():setattr(m,k,w)
        def pair(t,mechanical=None):
            select(old);a=m.local(t);select(new);b=m.local(t)
            for k in ['loss','sc','B','Bb','Be','pressure_map']:assert np.array_equal(a[k],b[k]),k
            for k in ['S','esc']:
                diff=a[k]-b[k]
                assert not (np.any(diff.data) if hasattr(diff,'nnz') else np.any(diff)),k
            if mechanical is not None:a['mechanical'],b['mechanical']=mechanical
            s=m.source(t)[0]/(m.scale*run.AMP)
            x=a['q']+s;y=b['q']+s
            ga=m.gas(a['q'],a['qb'],a['qe'])+a['mechanical']
            gb=m.gas(b['q'],b['qb'],b['qe'])+b['mechanical']
            difference=[np.sum(abs(x-y)*m.Eweight),np.sum(abs(x-y)*m.Nweight),
                np.sum(abs(ga[:,0]-gb[:,0])*m.eu),np.sum(abs(ga[:,1]-gb[:,1])*m.nu)]
            norm=[np.sum(abs(y)*m.Eweight),np.sum(abs(y)*m.Nweight),np.sum(abs(gb[:,0])*m.eu),np.sum(abs(gb[:,1])*m.nu)]
            return np.asarray(difference,LD),np.asarray(norm,LD),(a['mechanical'].copy(),b['mechanical'].copy())
        times,_,port=packets(photon/f'steps-{n}-reference-128.npz');pair(times[0]);mark=time.monotonic();pair(times[0]);warm=time.monotonic()-mark
        point=m.point_seconds/m.point_count;forecast=2*(2*17*point+stage_count*warm+2*(mark-setup-m.point_seconds)+10)
        assert forecast<300,('block audit forecast',forecast)
        maximum=np.zeros(4,LD);scale=np.zeros(4,LD)
        for first,last in times.reshape(-1,2):
            a,b,mech=pair(first);maximum=np.maximum(maximum,a);scale=np.maximum(scale,b)
            a,b,_=pair(last,mech);maximum=np.maximum(maximum,a);scale=np.maximum(scale,b)
        source=(maximum/np.maximum(scale,LD('1e-290'))).astype(float).tolist()
        paired=read(OUT/'result.json')['rows'][[64,128].index(n)]['photon_material_energy_H_residual']
        rows.append(dict(steps=n,lagged_inputs=inputs,stage_forcing_defect=source,paired_E_H=paired,
            stage_count=len(times),angular_port_relative=port,stage_operator_identical=True,forecast_seconds=forecast))
        print(json.dumps(rows[-1]),flush=True)
    maximum=max(max(list(r['lagged_inputs'].values())+r['stage_forcing_defect']+r['paired_E_H']) for r in rows)
    write(OUT/'block-result.json',dict(classification='Counterexample candidate',passed=maximum<.002,rows=rows,
        maximum_block_defect=maximum,actual_Radau_stage_clock_used=True,no_extra_sweep_authorized=True,
        no_new_evolution_steps=True,scope=plan['scope'],seconds=time.monotonic()-started))
    print(json.dumps(read(OUT/'block-result.json')),flush=True)

    assert maximum<.002,('Actual finite reciprocal block',maximum)
