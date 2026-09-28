def pilot():
    pass;folder=OUT/'sweep-1/photons';rows=[];then=time.monotonic()
    for n in [64,128]:
        m=prior.c.Response(n);row=m.run(n,f'pilot-{n}',n//16)
        with np.load(folder/f'pilot-{n}.npz') as p:
            tt=p['accepted_angular_times'];w=p['accepted_angular_quadrature_weights']
            h=m.t[-1]/n;expected=(h*(np.arange(n//16)[:,None]+RK_C)).ravel()
            assert tt.shape==w.shape==expected.shape and np.max(abs(tt-expected))<1e-18
            flux=(p['accepted_angular_luminosity']+p['accepted_angular_quadrature_correction'])@(np.arange(1,8,2)/32)
            error=float(abs(w@flux-p['radial_ports'][-1,1,1])/max(w@abs(flux),1e-290))
            assert error<1e-12,error
        row['actual_angular_quadrature_relative']=error;row['angular_quadrature_includes_separate_known_correction']=True;row['time_integrator']='RadauIIA2'
        write(folder/f'pilot-{n}.json',row);rows.append(row);assert row['passed'],row
        del m;gc.collect()
    data=[np.load(folder/f'pilot-{n}.npz')['moments'][:,[0,1,2,3,5,6]] for n in [64,128]]
    errors=np.max(np.sum(abs(data[0]-data[1]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(data[1]),axis=2),axis=0),1e-290)
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02),rows=rows,time_comparison=errors.astype(float).tolist(),
        seconds=time.monotonic()-then,full_horizon_completed=False,physical_final_charge_solved=False,full_goal_complete=False,
        remaining='A passing prefix is not a full-path result. Adapt all angular/time consumers and measure a justified production budget before continuing these exact prefixes.')
    write(folder/'pilot.json',result);print(json.dumps(result),flush=True);assert result['passed'],result
