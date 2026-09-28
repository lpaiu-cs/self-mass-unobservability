def analyze():
    assert not (OUT/'result.json').exists();signal.alarm(60);begin=time.monotonic()
    catalog=json.loads((a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    st=np.load(OUT/'station.npz');T=float(st['T']);rho=float(st['rho'])
    saved=np.load(a.OUT/'eos-state.npz');optical=np.load(a.OUT/'optical-levels.npz')
    inv=saved['number_fractions'][0];nav=6.02214076e23
    # Compile-time parameter (not a dynamic-library symbol); bind its source.
    constants=a.p.SOURCE/'mod_free_eos_constants.f90'
    assert 'avogadro = 6.02214076e23_fp_kind' in constants.read_text()
    conversion=rho*nav*float(st['mass_scale'])
    coeff=[];ratio_errors=[];level_map={};rows=[];coefficient_rows=[];parent_fractions=[]
    for stage,row in enumerate(catalog['rows']):
        element=row['element_index']-1;q=row['charge'];i=row['ion_index']-1
        lower=-st['zero'][element] if q==0 else st['ce'][i-1]+st['dv'][i-1]-st['binding'][i-1]*st['tc2'][0]
        upper=st['ce'][i]+st['dv'][i]-st['binding'][i]*st['tc2'][0]
        ln_ratio=upper-lower
        ratio_errors.append(float(abs(np.log(inv[element,q+1]/inv[element,q])-ln_ratio)))
        # This uses native free-energy stationarity coefficients, independently
        # of the saved ion population ratio. Electron activity is in dv.
        parent_match=np.flatnonzero(np.all(saved['ids'][0]==[element+1,q+1],axis=1))
        assert len(parent_match)<=1
        parent_L=0. if not len(parent_match) else saved['value'][0,parent_match[0],3]/saved['scale'][0,parent_match[0]]
        parent_fraction=np.exp(-parent_L);parent_fractions.append(parent_fraction)
        # Native photoionization channels end in the next ground state.
        # Make its LTE fraction explicit instead of treating the whole stage
        # as a ground state. The factor cancels only under this declared LTE.
        K=np.exp(-np.longdouble(ln_ratio))*optical[f'fraction_{stage}']/parent_fraction
        coeff.append(K)
        for local,index in enumerate(np.flatnonzero(row['included'])):
            level=row['global_levels'][index];level_map[level]=(stage,local)
            coefficient_rows.append([level,stage,local,float(K[local]),row['binding_cm_inverse'][index]*erg,row['anchor_shift_erg'],float(parent_fraction)])
        rows.append(dict(Z=row['Z'],charge=q,log_parent_over_stage=float(ln_ratio),parent_ground_fraction=float(parent_fraction),chemical_log_ratio_error=ratio_errors[-1]))
    raw=np.loadtxt(OUT/'fort.92');lines=[];excluded=[];bb_errors=[]
    echarge=1.602176634e-19*.1*c
    me=ryd*(h/echarge**2)*(h/echarge)**2*c/(2*np.pi**2)
    sigma_integral=np.pi*echarge**2/(me*c)
    for index,line in enumerate(raw):
        lo,hi,mode=int(line[3]),int(line[4]),int(line[5]);f=float(line[7])
        reason=None
        if abs(mode)!=1:reason='disabled_or_unsupported_mode'
        elif lo not in level_map or hi not in level_map:reason='outside_common_retained_levels'
        elif f<=0:reason='nonpositive_oscillator_strength'
        if reason:excluded.append(dict(row=index,reason=reason));continue
        stage,l=level_map[lo];stage2,u=level_map[hi];assert stage==stage2
        weight=optical[f'weight_{stage}'];lw=optical[f'logw_{stage}'];pop=optical[f'populations_{stage}']*conversion
        energy=(optical[f'excitation_cm_inverse_{stage}'][u]-optical[f'excitation_cm_inverse_{stage}'][l])*erg
        if energy<=0:excluded.append(dict(row=index,reason='nonpositive_photon_energy'));continue
        nu=energy/h;ratio=np.exp(lw[u]-lw[l]);up=min(np.longdouble(1),ratio);down=min(np.longdouble(1),1/ratio)
        amplitude=sigma_integral*f
        abs0=pop[l]*amplitude*up
        stim=pop[u]*amplitude*weight[l]/weight[u]*down
        spontaneous=stim*(2*h*nu**3/c**2)
        planck=2*h*nu**3/c**2/np.expm1(energy/(k*T))
        net=abs0-stim
        err=abs(spontaneous-net*planck)/max(abs(spontaneous),abs(net*planck),1e-300)
        bb_errors.append(float(err));assert net>0
        lines.append([stage,l,u,float(nu),f,float(abs0),float(stim),float(spontaneous),float(net),float(up),float(down)])
    cross=np.load(OUT/'cross-sections.npz');meta=cross['metadata'];sigma=cross['sigma'];forward=np.zeros(len(sigma));reverse=np.zeros(len(sigma))
    continuum=np.loadtxt(OUT/'fort.98');channels={}
    for record in continuum:
        level=int(record[1])
        if level not in level_map:continue
        assert level not in channels,'Multiple photoionization channels need an explicit sum'
        assert record[2]==record[3],'An excited parent channel needs its own threshold'
        channels[level]=record
    assert set(channels)==set(level_map),'Never invent a hydrogenic channel for an unlisted level'
    for stage,row in enumerate(catalog['rows']):
        mask=meta[:,0]==stage;local=meta[mask,1].astype(int)
        pop=optical[f'populations_{stage}']*conversion
        parent=inv[row['element_index']-1,row['charge']+1]*conversion*parent_fractions[stage]
        forward[mask],reverse[mask],_,_,_=continuum_rates(pop[local],parent,coeff[stage][local],sigma[mask],meta[mask,2],T,h,c,k)
    u=meta[:,2];stim=reverse*np.exp(-u);net=forward-stim
    nu=u*k*T/h;spontaneous=stim*(2*h*nu**3/c**2)
    planck=2*h*nu**3/c**2/np.expm1(u)
    positive=forward>0
    bf_error=float(np.max(abs(spontaneous[positive]-net[positive]*planck[positive])/np.maximum(abs(spontaneous[positive]),abs(net[positive]*planck[positive]))))
    assert np.all(sigma[meta[:,3]<0]==0) and np.all(net>=0)
    assert all(np.isfinite(x).all() for x in [forward,reverse,stim,net,spontaneous])
    np.savez_compressed(OUT/'rates.npz',lines=np.array(lines),continuum_coefficients=np.array(coefficient_rows),bf_forward=forward,bf_reverse_prefactor=reverse,bf_stimulated=stim,bf_spontaneous=spontaneous,bf_net=net)
    write(OUT/'coverage.json',dict(classification='Counterexample candidate',atomic_stages=rows,parsed_lines=len(raw),retained_lines=len(lines),excluded=excluded,
        threshold_zero_samples=int(np.sum(meta[:,3]<0)),continuum_positive_samples=int(positive.sum()),sampled_levels=len(level_map),
        unresolved=['Missing H/He II provider amplitudes in this metal/He I extraction','Missing ion stages and source model transitions','Dissolved bound states and oscillator-strength redistribution','Line broadening, finite-profile redistribution and resonant continuum integration','Microscopic nonideal recombination activity correction','Opacity derivatives and off-equilibrium population kinetics']))
    passed=bool(max(ratio_errors)<1e-9 and max(bb_errors)<1e-9 and bf_error<1e-9)
    result=dict(classification='Counterexample candidate',passed=passed,decision='COMMON_FINITE_ATOMIC_POINTWISE_RATES_IMPLEMENTED_PHYSICAL_COMPLETENESS_OPEN',
        same_state_bitwise=True,native_EOS_calls=1,atomic_stages=len(rows),common_levels=len(level_map),bound_bound_lines=len(lines),bound_free_samples=len(sigma),
        chemical_log_ratio_max=max(ratio_errors),bound_bound_balance_relative=max(bb_errors),bound_free_balance_relative=bf_error,
        line_net_strength_sum_cm_inverse_Hz=float(np.array(lines)[:,8].sum()),seconds=time.monotonic()-begin,avogadro=nav,
        new_stellar_steps=0,coupled_runs=0,physical_opacity_certified=False,full_dynamic_charge_solved=False,
        bindings={str(p):digest(p) for p in [a.OUT/'eos-state.npz',a.OUT/'optical-levels.npz',a.OUT/'catalog.json',OUT/'station.npz',OUT/'cross-sections.npz',OUT/'fort.92',OUT/'fort.93',OUT/'fort.98',CACHE/'atomic-rates',a.LIB,a.CACHE/'gas.so',constants]})
    write(OUT/'result.json',result)
    n,K,s,z,b=sp.symbols('n K s z b',positive=True)
    A=n*K*s;R=n*K*s*z;J=R*b;B=b*z/(1-z)
    assert sp.simplify(J-(A-R)*B)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='If the native thermodynamic coefficient K independently implies n_bound=n_parent_ground K, absorption A=n_bound sigma, stimulated R=n_parent_ground K sigma exp(-u), and spontaneous J=R 2h nu^3/c^2 obey J=(A-R) Bnu. Actual sigma is retained in both directions.',
        limitation='Conditional fixed-state equilibrium closure. It does not derive microscopic nonideal recombination or off-equilibrium activities.'))
    signal.alarm(0);print(json.dumps(result,indent=2),flush=True);assert passed
