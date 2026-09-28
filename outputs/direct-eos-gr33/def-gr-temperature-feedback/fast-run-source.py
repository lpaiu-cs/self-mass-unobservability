def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for path,h in plan['bindings'].items():assert go.task.digest(Path(path))==h,path
    pilotrow=json.loads((OUT/'fast-pilot.json').read_text());assert pilotrow['forecast_seconds']<600
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    started=time.monotonic();patch.install();p=Problem();m=p.model;nr=len(m.original.native)
    half=np.zeros((go.COUNT//2+1,2*nr),np.clongdouble);L=go.laguerre(12)
    ell=np.zeros(go.COUNT,np.longdouble);ell[:2048]=L[-1];kernel=go.fft(ell)[:go.COUNT//2+1]/go.COUNT
    base_end=np.zeros(m.size,np.longdouble);corr_end=base_end.copy();dEend=np.zeros(len(m.heat.edges),np.longdouble)
    pilot=np.load(OUT/'pilot.npz');saved={int(k):i for i,k in enumerate(pilot['ids'])}
    for k in range(1,go.COUNT//2+1):
        z,factor=go.contour(k,12)
        if k in saved:
            i=saved[k];q,dq,E,dE=[pilot[name][i] for name in ['base_q','correction_q','base_E','correction_E']]
        else:q,dq,E,dE=p.transform_pair(z)
        half[k,:nr]=factor*z*p.speed*(m.nativeV[0]@dq);half[k,nr:]=factor*(m.nativeV[1]@dq)
        mult=(1 if k==go.COUNT//2 else 2)*kernel[k]*factor
        base_end+=(mult*q).real;corr_end+=(mult*dq).real;dEend+=(mult*dE).real
        if k%512==0:print('FEEDBACK',k,round(time.monotonic()-started,2),flush=True)
    half[-1]=half[-1].real;native=np.empty((3,65,2*nr),np.longdouble);coarse=np.empty((65,2*nr),np.longdouble)
    for lo in range(0,2*nr,128):
        hi=min(lo+128,2*nr);a=go.coefficients(half[:,lo:hi]);b=go.coefficients(half[::2,lo:hi])
        for j,n in enumerate(go.DEGREES):native[j,:,lo:hi]=L[:,:n]@a[:n]
        coarse[:,lo:hi]=L@b[:2048]
    initial_raw=float(abs(native[:,0]).max());native[:,0]=0;coarse[0]=0
    baseline=np.load(patch.OUT/'p4-2048.npz');w=baseline['weights'];masks=baseline['masks']
    def norm(values):
        v=values[...,:nr];f=values[...,nr:]
        return np.stack([np.sqrt(np.sum(w*v*v,axis=-1)),np.sqrt(np.sum(w*f*f,axis=-1)),
            *[np.sqrt(np.sum(w[mask]*v[...,mask]**2,axis=-1)/w[mask].sum()) for mask in masks]],axis=-1)
    base_native=np.c_[baseline['native_velocity'],baseline['native_scalar']]
    scale=norm(base_native).max(0); effect=norm(native[-1]).max(0)/scale
    time_error=norm(native[-1]-native[-2]).max(0)/scale;contour_error=norm(native[-1]-coarse).max(0)/scale
    endpoint_relative=float(abs(base_end-baseline['q']).max()/max(abs(baseline['q']).max(),1e-100))
    balance=float(abs(np.sum(-np.diff(dEend),dtype=np.longdouble))/max(abs(dEend).max(),1e-100))
    np.savez_compressed(OUT/'response.npz',times=go.TIMES,radius=m.original.native,weights=w,masks=masks,
        correction=native,coarse_correction=coarse,base_q=base_end,correction_q=corr_end,correction_heat_energy=dEend,
        total_native_velocity=baseline['native_velocity']+native[-1,:,:nr],total_native_scalar=baseline['native_scalar']+native[-1,:,nr:])
    passed=bool(max(time_error)<1e-6 and max(contour_error)<1e-6 and endpoint_relative<1e-10 and balance<2e-13 and p.loop_residual<1e-10)
    row=dict(classification='Counterexample candidate',actual_temperature_feedback_evolved=True,numerical_target_passed=passed,
        temperature_only_feedback_below_target=bool(passed and max(effect)<1e-4),
        comparisons={key:dict(correction_relative=float(effect[j]),time_absolute_to_baseline=float(time_error[j]),contour_absolute_to_baseline=float(contour_error[j])) for j,key in enumerate(go.task.FIELDS)},
        baseline_endpoint_replay=endpoint_relative,heat_balance=balance,linear_residual=p.error,feedback_residual=p.loop_residual,
        maximum_corrections=max(p.iterations,pilotrow['corrections']),sampled_second_to_first=max(p.second_to_first,pilotrow['second_to_first']),
        sampled_heat_feedback_relative=max(p.feedback_relative,pilotrow['sampled_heat_feedback_relative']),raw_initial_correction_max=initial_raw,
        seconds=time.monotonic()-started,total_compute_seconds=time.monotonic()-started+pilotrow['seconds'],memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        full_temperature_metric_coefficient_feedback=False,full_GR_photon_feedback_evolved=False,full_nonlinear_evolution=False,full_dynamic_charge_solved=False,
        scope='Temperature-only feedback on the fixed transport geometry/rates/coefficients. This new loop is solved at p4 with absolute correction tolerance relative to the accepted baseline, not a new spatial/physical certification or rigorous continuum error bound.')
    write(OUT/'result.json',row);signal.alarm(0);print('FEEDBACK RESULT',json.dumps(row),flush=True)
