def bank():
    assert not (OUT/'bank-result.json').exists(); start=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout); signal.alarm(90)
    errors=[]; files=[]
    for reference in [128]:
        m=Response(reference); folder=OUT/f'bank-{reference}'; folder.mkdir()
        Eunit=[]; Nunit=[]
        for k in range(17):
            m.restore(k); c=m.coefficients(); jets={}
            for param in ['dr','dt','dy']:
                plus=m.coefficients(**{param:1e-5}); minus=m.coefficients(**{param:-1e-5})
                jets[param]={key:(plus[key]-minus[key])/2e-5 for key in ['emit','loss','sc','energy','neutral','u','p']}
                if reference==128 and k in [0,8,16]:
                    plus=m.coefficients(**{param:5e-6}); minus=m.coefficients(**{param:-5e-6})
                    for key in ['emit','loss','sc','energy','neutral','u','p']:
                        fine=(plus[key]-minus[key])/1e-5; coarse=jets[param][key]
                        if key in ['emit','loss']:
                            w=m.weights*m.E*(1+m.I[k])
                            err=np.sum(abs(fine-coarse)*w)/max(np.sum(abs(fine)*w),1.)
                        else: err=np.max(abs(fine-coarse))/max(np.max(abs(fine)),1.)
                        errors.append(dict(reference=reference,k=k,parameter=param,quantity=key,relative=float(err)))
            active=c['rho']>0; assert np.all(jets['dt']['energy'][active]>0)
            gammaT=np.divide(c['p']/np.maximum(c['rho'],1e-300)-jets['dr']['u'],jets['dt']['u'],out=np.zeros(m.n),where=active)
            record=dict(c,adiabatic_logT=gammaT)
            record.update({param+'_'+key:value for param,v in jets.items() for key,value in v.items()})
            p=folder/f'point-{k}.npz'; np.savez_compressed(p,**record); files.append(str(p))
            Eunit.append(jets['dt']['energy']); Nunit.append(c['neutral'])
        np.savez_compressed(folder/'units.npz',energy=np.maximum(np.max(Eunit,axis=0),1.),neutral=np.maximum(np.max(Nunit,axis=0),1.))
    result=dict(classification='Counterexample candidate',passed=max(e['relative'] for e in errors)<1e-4,derivative_checks=errors,seconds=time.monotonic()-start,
        new_native_bank_calls=0,coefficient_model='existing corrected native interpolation, centered relative jets',full_uniform_derivative_bound=False)
    write(OUT/'bank-result.json',result); signal.alarm(0); print(json.dumps({k:v for k,v in result.items() if k!='derivative_checks'}),flush=True); assert result['passed']
