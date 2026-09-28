"""Counterexample candidate: independent saved-flux, port and source audit."""
import json,resource,signal,time
import numpy as np
import def_retained_native_return as run


def diagnose():
    start=time.monotonic()
    def timeout(*_):raise TimeoutError('30s saved-flux diagnosis cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(30)
    try:
        run.initialize();m=run.Material(driven=False);rows=[]
        for k in [0,8,14,16]:
            d=np.load(run.HYDRO/f'point-{k}.npz');p=m.point(k);df=d['flux'].astype(np.longdouble);dg=d['gravity'].astype(np.longdouble)
            rate=-np.diff(df,axis=1);rate[1]+=dg;norm=np.maximum(np.sum(abs(rate),axis=1),1.)
            f=p['flux'].astype(np.longdouble);g=p['gravity'].astype(np.longdouble)
            floor=16*np.finfo(float).eps*np.sum(abs(f+df)+abs(f),axis=1);floor[1]+=16*np.finfo(float).eps*np.sum(abs(g+dg)+abs(g))
            row=dict(k=k,owner=p['owner_error'],rate_L1=np.asarray(norm,float).tolist(),floor=np.asarray(floor,float).tolist(),ratio=np.asarray(floor/norm,float).tolist())
            rows.append(row)
        result=dict(classification='Counterexample candidate',rows=rows,new_native_calls=0,physical_time_steps=0,seconds=time.monotonic()-start)
        run.write(run.OUT/'hydro-arithmetic-diagnosis.json',result);print(json.dumps(result),flush=True)
    finally:
        signal.alarm(0);run.write(run.OUT/'diagnosis-receipt.json',dict(seconds=time.monotonic()-start,cap_seconds=30))


def main(seconds):
    start=time.monotonic()
    def timeout(*_):raise TimeoutError('Remaining readout audit budget')
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,seconds)
    run.initialize();m=run.Material(driven=False);hydro=[];rho0=m.model.flow.eos.rho0
    for k in range(17):
        d=np.load(run.HYDRO/f'point-{k}.npz');p=m.point(k);df=d['flux'].astype(np.longdouble);dg=d['gravity'].astype(np.longdouble)
        rate=-np.diff(df,axis=1);rate[1]+=dg;ledger=df[:,0]-df[:,-1];ledger[1]+=dg.sum()
        norm=np.maximum(np.sum(abs(rate),axis=1),1.)
        error=float(max(np.max(np.sum(abs(rate-d['rate']),axis=1)/norm),np.max(abs(ledger-d['ledger'])/norm),np.max(abs(rate.sum(1)-ledger)/norm)))
        f=p['flux'].astype(np.longdouble);g=p['gravity'].astype(np.longdouble)
        eps=float(d['arithmetic_epsilon']);assert eps==np.finfo(np.longdouble).eps
        floor=16*eps*np.sum(abs(f+df)+abs(f),axis=1);floor[1]+=16*eps*np.sum(abs(g+dg)+abs(g))
        resolution=float(np.max(floor/norm));assert error<1e-8 and resolution<.002
        hydro.append(dict(k=k,conservative_export_relative=error,arithmetic_indicator=resolution))
    faces={};density=0.;duplicates=0
    for p in sorted(run.HYDRO.glob('*-native-states.json')):
        for row in run.read(p):
            key=tuple(row['coordinates']);raw=np.asarray(row['raw'])
            density=max(density,abs(raw[0]/(rho0*np.exp(key[0]))-1))
            assert row['population_error']<1e-12 and np.isfinite(raw).all() and raw[1]>0
            if key in faces:
                assert np.max(abs(raw-faces[key])/np.maximum(abs(raw),1.))<1e-10;duplicates+=1
            faces[key]=raw
    assert density<1e-10 and len(faces)<=16000
    angular=[];sources=[];gamma=1-1/np.sqrt(2)
    for n in [64,128]:
        p=np.load(run.PHOTON/f'steps-{n}-reference-128.npz');t=p['accepted_angular_times'];L=p['accepted_angular_luminosity']
        h=float(p['t'][-1])/n;expected=np.ravel(h*(np.arange(n)[:,None]+np.array([gamma,1.])))
        assert len(t)==len(L)==2*n and np.max(abs(t-expected))<1e-18
        pairs=(L@((2*np.arange(4)+1)/32)).reshape(n,2)
        cumulative=np.r_[0.,np.cumsum(h*(pairs@np.array([1-gamma,gamma])))];ids=np.rint(p['t']/h).astype(int)
        error=float(np.max(abs(cumulative[ids]-p['radial_ports'][:,1,1]))/max(h*np.sum(abs(pairs)),1.));assert error<1e-12
        angular.append(dict(steps=n,stages=len(t),relative=error))
        d=np.load(run.GR/f'source-{n}-reference-128.npz');s=np.load(run.MATERIAL/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(np.longdouble)*np.longdouble(d['cx'])*np.longdouble(run.C)**2
        E=rest+d['gas_nonrest_energy_erg']
        errors=[E-s[:,0],rest+d['nonrest_trace_erg']-s[:,2],rest+d['nonrest_stress_erg']-(s[:,0]-s[:,1]),
            d['pressure_volume_erg']-s[:,3],d['metric_stress_erg']-(E+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        error=float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(s)),1.));assert error<1e-12
        sources.append(dict(steps=n,relative=error))
    assert run.read(run.OUT/'result.json')['passed']
    import def_native_global_scalar_closure as exterior
    d=np.load(run.OUT/'infinity/charge-128-a8-r8.npz');s=np.load(run.GR/'applied-source.npz')
    eps=exterior.G/run.C**4*d['arrived_energy_erg']/float(s['M_cm']);alpha=-float(s['K_cm'])/float(s['M_cm'])
    rebuilt=(d['compact']+d['normalized_exterior']+alpha*eps)/(1-eps)
    infinity_error=float(np.max(abs(rebuilt-d['normalized']))/max(np.max(abs(rebuilt)),1e-300));assert infinity_error<1e-12
    for j,source in enumerate(m.extended_sources):(run.OUT/f'extended-owner-{j}.py').write_text(source)
    result=dict(classification='Counterexample candidate',passed=True,hydro=hydro,unique_native_face_states=len(faces),
        repeated_exact_coordinates=duplicates,native_density_relative=density,accepted_angular_ports=angular,source_identities=sources,
        infinity_normalization_relative=infinity_error,infinity_controls_passed=run.read(run.OUT/'infinity/result.json')['passed'],
        seconds=time.monotonic()-start,scope='Saved conservative flux and angular-stage/source identities;16eps is a diagnostic, not a uniform EOS or nonlinear remainder bound.',
        coupled_fixed_point_verified=False,final_charge_solved=False)
    run.write(run.OUT/'verification.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    r=run.read(run.OUT/'readout-receipt.json');p=run.read(run.OUT/'infinity-receipt.json')
    remaining=r['cap_seconds']-(p['total_charged_wall_seconds']-r['total_charged_wall_seconds']+r['seconds']);assert remaining>0
    start=time.monotonic();cpu=time.process_time()
    try:main(remaining)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);elapsed=time.monotonic()-start
        run.write(run.OUT/'verification-receipt.json',dict(seconds=elapsed,CPU_seconds=time.process_time()-cpu,cap_seconds=remaining,
            readout_group_seconds=r['cap_seconds']-remaining+elapsed,
            total_charged_wall_seconds=p['total_charged_wall_seconds']+elapsed,
            total_charged_CPU_seconds=p['total_charged_CPU_seconds']+time.process_time()-cpu))
