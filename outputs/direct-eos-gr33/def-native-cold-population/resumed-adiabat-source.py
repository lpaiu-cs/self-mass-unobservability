def resumed():
    from scipy.interpolate import PchipInterpolator
    from scipy.integrate import cumulative_simpson
    begin=time.monotonic();signal.alarm(60)
    ion=Ions(cap=800);fan=ion.fan;lr=np.log(fan.rho);lt=np.log(fan.T)
    base=ion.snapshot(lr,lt,np.zeros(318));target=base['number_fractions'];s0=base['eos'][3]
    d=np.load(old.OUT/'accepted-adiabat-prefix.npz');density=np.linspace(0,-6,49)
    start=len(d['T']);assert start==34 and np.array_equal(d['log_density_ratio'],density[:start])
    values=list(d['raw']);affinities=list(d['fields']);counts=list(d['number_fractions']);solvedT=list(d['T'])
    errors=[[(row[3]-s0)*T/(1.5*row[1]/row[0]),float(np.max(abs(q-target))/target.sum())] for row,T,q in zip(values,solvedT,counts)]
    fields=affinities[-1].copy();lastT=np.log(solvedT[-1]);lastx=density[start-1]
    raw=(SOURCE/'mod_ionization_data.f90').read_text()
    body=raw[raw.index('monatomic_ip(nions) ='):];body=body[:body.index(']')]
    body='\n'.join(line.split('!')[0] for line in body.splitlines())
    potentials=np.array([float(v.replace('d','e')) for v in re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[edED][-+]?[0-9]+)?)_fp_kind',body)])
    assert len(potentials)==316,('Native ionization constants',len(potentials))
    c2=float(json.loads((native.OUT.parent/'def-photon-shared-atomic/catalog.json').read_text())['constants'][2])
    energy=np.concatenate([np.cumsum(potentials[ion.starts[e]:ion.starts[e+1]]) for e in range(24)])*c2
    charges=np.concatenate([np.arange(1,z+1) for z in ion.Z])
    active=np.concatenate([target[e,1:z+1]>target[e].sum()*1e-18 for e,z in enumerate(ion.Z)])
    diss=float(re.search(r'h2diss = ([0-9.]+)_fp_kind',raw).group(1))*c2
    molecule_ion=diss+potentials[0]*c2-21375.95*c2
    try:
        for i,x in enumerate(density):
            if i<start:continue
            theta=lastT+2*(x-lastx)/3
            fields[:316]+=np.where(active,energy*(np.exp(-theta)-np.exp(-lastT))+charges*((x-lastx)-1.5*(theta-lastT)),0.)
            fields[316]+=-diss*(np.exp(-theta)-np.exp(-lastT))-(x-lastx)+1.5*(theta-lastT)
            fields[317]+=molecule_ion*(np.exp(-theta)-np.exp(-lastT))+(x-lastx)-1.5*(theta-lastT)
            for iteration in range(8):
                assert np.exp(theta)>=100,'Registered cold domain boundary'
                a,n,error=ion.constrain(lr+x,theta,target,fields,target_molecules=base['molecular_H_fractions'],tolerance=1e-12);fields=ion.fields.copy()
                row=a['eos'];T=np.exp(theta);scale=1.5*row[1]/row[0]
                residual=(row[3]-s0)*T/scale
                if abs(residual)<2e-10:break
                assert abs(residual)<.15,'Entropy root left bounded local step'
                theta-=residual
            else:raise AssertionError('Fixed-ion entropy solve')
            values.append(row);affinities.append(fields.copy());counts.append(a['number_fractions']);errors.append([residual,error])
            solvedT.append(T)
            lastT=theta;lastx=x
            write(OUT/'adiabat-progress.json',dict(nodes=i+1,EOS_calls=ion.calls,seconds=time.monotonic()-begin,T=T))
        rows=np.array(values)
        # Retain the solved temperatures independently of any EOS output index.
        T=np.array(solvedT)
        curves=[]
        for stride in [2,1]:
            d=-density[::stride];a=rows[::stride]
            gamma=-PchipInterpolator(d,np.log(a[:,1])).derivative()(d)
            enthalpy=fan.cx*fan.c**2+a[:,2]+a[:,1]/a[:,0]
            cs=np.sqrt(gamma*a[:,1]/a[:,0]/enthalpy)
            rapidity=cumulative_simpson(cs,x=d,initial=0);v=np.tanh(rapidity);xi=(v-cs)/(1-v*cs)
            assert np.all(gamma>1) and np.all(np.diff(xi)>0)
            curves.append(dict(gamma=gamma,cs=cs,rapidity=rapidity,velocity=v,xi=xi))
        contrast=float(max(abs(curves[0]['rapidity']-curves[1]['rapidity'][::2]))/max(curves[1]['rapidity']))
        eq=dict(np.load(native.OUT/'fine.npz'));comparisons=[]
        for x in [-1.,-2.,-4.,-6.]:
            i=int(np.argmin(abs(density-x)));j=int(np.argmin(abs(eq['log_density_ratio']-x)))
            comparisons.append(dict(log_density_ratio=x,fixed_T=float(T[i]),LTE_T=float(eq['T'][j]),
                pressure_ratio=float(rows[i,1]/eq['raw'][j,1]),fixed_gamma=float(curves[1]['gamma'][i]),
                LTE_gamma=float(eq['raw'][j,4]),fixed_H_ion=float(rows[i,14]),LTE_H_ion=float(eq['raw'][j,14])))
        np.savez_compressed(OUT/'adiabat.npz',log_density_ratio=density,T=T,raw=rows,fields=affinities,
            number_fractions=counts,errors=errors,**curves[1])
        result=dict(classification='Counterexample candidate',passed=bool(contrast<.002 and np.max(abs(errors),axis=0)[0]<1e-8),
            nodes=len(rows),EOS_calls=ion.calls,seconds=time.monotonic()-begin,
            coarse_fine_rapidity_relative=contrast,maximum_entropy_scaled=float(np.max(abs(errors),axis=0)[0]),
            maximum_population_error=float(np.max(abs(errors),axis=0)[1]),comparisons=comparisons,
            physical_frozen_ion_limit_certified=False,full_fluid_evolution=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/'adiabat.json',result);assert result['passed'];print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'adiabat-failure.json',dict(error=repr(exc),EOS_calls=ion.calls,seconds=time.monotonic()-begin,nodes=len(values)))
        raise
    finally:
        ion.save('adiabat-native-states.npz')
        write(OUT/'density-fallbacks.json',dict(rows=ion.density_fallbacks,source_sha256=sha(__file__)))
    signal.alarm(0)
