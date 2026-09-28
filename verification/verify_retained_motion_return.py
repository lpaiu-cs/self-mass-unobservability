"""Independent saved-state checks for the one finite material/photon return."""
import json, resource, signal, time
import numpy as np
import def_retained_motion_return as run


def audit():
    run.initialize();m=run.Response();rows=[];LD=np.longdouble
    direct=run.read(run.OUT/'direct-motion/result.json');assert direct['passed']
    assert not run.read(run.OUT/'infinity/result.json')['passed']
    for n in [64,128]:
        p=np.load(run.TOTAL/f'steps-{n}-reference-128.npz')
        a=np.load(run.BEFORE/f'photons/steps-{n}-reference-128.npz')
        d=np.load(run.PHOTON/f'steps-{n}-reference-128.npz')
        old=np.load(run.BEFORE/f'material/steps-{n}-reference-128.npz')
        new=np.load(run.MATERIAL/f'steps-{n}-reference-128.npz')
        ids=np.rint(np.arange(17)*n/16).astype(int)
        assert np.array_equal(old['t'],new['t']) and np.max(abs(old['t'][ids]-p['t']))<1e-18
        state=old['history_scaled'][ids]*run.AMP
        exact_g=np.moveaxis(state[:,[2,3]],1,2)/np.stack([m.eu,m.nu],axis=-1)+d['material_history']
        assert np.array_equal(p['material_history'],exact_g)
        assert np.array_equal(p['collision_transfer'],a['collision_transfer']+d['collision_transfer'])
        assert np.array_equal(p['moments'][:,[1,2]],state[:,[2,3]]+d['moments'][:,[1,2]])
        assert np.array_equal(p['photon_history_scaled_occupation'],a['photon_history_scaled_occupation']+d['photon_history_scaled_occupation'])
        ports=[]
        for z in [p,d]:
            gamma=1-1/np.sqrt(2);h=float(z['t'][-1])/n
            angular=z['accepted_angular_luminosity'].astype(LD);weights=np.arange(1,8,2,dtype=LD)/32
            cumulative=np.r_[LD(0),np.cumsum(h*(angular@weights).reshape(n,2)@np.array([1-gamma,gamma],LD))]
            error=float(np.max(abs(cumulative[ids]-z['radial_ports'][:,1,1]))/max(h*np.sum(abs(angular)@weights),1.))
            assert error<1e-12;ports.append(error)
        source=np.load(run.GR/f'source-{n}-reference-128.npz');stress=np.load(run.MATERIAL/f'stress-{n}-reference-128.npz')['material']
        rest=source['baryon_g'].astype(LD)*LD(source['cx'])*LD(run.C)**2
        energy=rest+source['gas_nonrest_energy_erg']
        errors=[energy-stress[:,0],rest+source['nonrest_trace_erg']-stress[:,2],
                rest+source['nonrest_stress_erg']-(stress[:,0]-stress[:,1]),source['pressure_volume_erg']-stress[:,3]]
        identity=float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(stress)),1.));assert identity<1e-12
        actual=new['history_scaled'][ids]*run.AMP
        mismatch=(np.max(np.sum(abs(actual[:,[2,3]]-p['moments'][:,[1,2]]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(actual[:,[2,3]]),axis=2),axis=0),1.)).tolist()
        change=(np.max(np.sum(abs(actual-state),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1.)).tolist()
        inc=np.load(run.OUT/f'direct-motion/source-{n}.npz')
        dq=(new['history_scaled'][ids].astype(LD)-old['history_scaled'][ids].astype(LD))*LD(run.AMP)
        rest=inc['baryon_g']*LD(source['cx'])*LD(run.C)**2
        # Response does not own a material model; the saved fixed energy reference does.
        reference=run.prior.Material(driven=False)
        energy=(dq[:,2]+LD(reference.rest)*dq[:,0])/reference.a
        error=float(np.max(abs(energy-rest-inc['gas_nonrest_energy_erg']))/max(np.max(abs(energy)),1.));assert error<1e-12
        assert np.array_equal(inc['baryon_g'],dq[:,0])
        assert np.array_equal(inc['photon_energy_erg'],d['moments'][:,0]/reference.a)
        rows.append(dict(steps=n,actual_material_baseline_exact=True,collision_only_ledger_exact=True,
            accepted_port_relative=ports,source_identity_relative=identity,energy_H_waveform_residual=mismatch,
            material_update_relative=change,direct_increment_energy_identity=error))
    import def_native_global_scalar_closure as exterior
    z=np.load(run.OUT/'direct-motion/charge-128-a8-r8.npz');s=np.load(run.GR/'applied-source.npz')
    eps=exterior.G/run.C**4*z['arrived_energy_erg']/float(s['M_cm']);alpha=-float(s['K_cm'])/float(s['M_cm'])
    rebuilt=(z['compact']+z['normalized_exterior']+alpha*eps)/(1-eps)
    error=float(np.max(abs(rebuilt-z['normalized']))/max(np.max(abs(rebuilt)),1e-300));assert error<1e-12
    run.write(run.OUT/'verification.json',dict(classification='Counterexample candidate',passed=True,paths=rows,
        normalized_identity_relative=error,actual_reciprocal_sweep_completed=True,
        direct_motion_controls=direct['motion_return_controls'],first_readout_failure_preserved=True,
        pressure_increment_arithmetic_diagnostic=direct['pressure'],
        no_native_defect_or_mechanical_double_count=True,coupled_fixed_point_verified=False,
        uniform_EOS_derivative_bound=False,final_charge_solved=False,full_goal_complete=False))


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    def timeout(*_):raise TimeoutError('60s shared combination/audit budget')
    signal.signal(signal.SIGALRM,timeout)
    previous=run.read(run.OUT/'combine-receipt.json')['seconds'];assert previous<60
    signal.setitimer(signal.ITIMER_REAL,60-previous);start=time.monotonic();cpu=time.process_time();error=None
    try:audit()
    except Exception as exc:error=repr(exc);raise
    finally:
        signal.alarm(0);run.write(run.OUT/'audit-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=run.sha(__file__)))
