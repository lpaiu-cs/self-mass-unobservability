"""Physical launch-energy input on the existing actual Radau clocks.

No source-work diagnostic is added and no physical trajectory is integrated.
This is the input to the next physical exterior/metric connection, not its
completed response or a final charge certificate.
"""
from pathlib import Path
import fcntl,json,resource,time
import numpy as np

lock=Path('.native-reader-initialization.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
import reconcile_native_mass_energy as old
old.inf.prior.initialize()
fcntl.flock(lock,fcntl.LOCK_UN)

OUT=Path('native-physical-port255-work');ACTUAL=Path('native-returned-krylov254-work')
PRIMARY=Path('native-common-arithmetic239-work');FIELD=Path('native-retarded-extension248-work')
read,write,sha,LD=old.read,old.write,old.sha,old.LD
assert not OUT.exists();OUT.mkdir()
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2);old.inf.incident.native.deadline(300)
start=time.monotonic();error=None
try:
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Construct the physical instantaneous outer launch-energy increment from the same primary and applied-return geometry on the exact existing Radau clocks.',
        method='Lphys=Lref+(u+nu)_face*L0. L0 is the SAME128background outgoing angular spectrum interpolated on the original17canonical knots. No second volume/speed factor. High uses the actual155input; low uses the actual applied128return metric, also for coarse. Independent4/8geometry controls remain separate.',
        decision='Feed these current inputs into exterior energy/work and the same-history GR boundary; do not relabel the current254reference-energy result as complete physical infinity charge.',
        budget_seconds=300,virtual_GiB=8,CPU_affinity=2,
        forecast='166same boundary owner measured4.279s with746MB RSS on384stages. Current700stage evaluations and two geometry orders estimated10..30s, later cost unmeasured; five minutes allowed.',
        physical_steps=0,scientific_gates_changed=False,final_charge_conclusion='unadjudicated'))
    background=old.inf.prior.EV/'coupled-128.npz';clock=old.inf.incident.METRIC/'corrected/metric-128-g8.npz'
    d=np.load(background);t=np.load(clock)['t'];ids=np.array([np.argmin(abs(d['snapshot_t']-v)) for v in t])
    assert np.max(abs(d['snapshot_t'][ids]-t))<1e-18
    outer=d['snapshot_I'][ids].sum(1)[:,-1];assert outer.shape==(17,8,152)
    files=[Path(__file__),Path(old.__file__),Path(old.inf.incident.__file__),background,clock]
    values={};rows=[]
    for order in [4,8]:
        driver=old.inf.incident.Driver(order);model=driver.model;b=model.bulk
        L=LD(2*np.pi*old.C)*model.area[-1]*(outer[:,b.mu>0]@(b.d['num']*b.d['Einf']))
        assert L.shape==(17,4)
        metric_path=ACTUAL/f'metric/metric-128-g{order}.npz';field_path=FIELD/f'gr/fields-128-g{order}.npz'
        metric=np.load(metric_path);field=np.load(field_path)
        assert np.array_equal(metric['t'],field['t'])
        assert abs(float(field['radius_E'][-1])/driver.r0-1)<1e-12
        low_ell=metric['delta_nu_faces'][:,-1].astype(LD)+LD(driver.z0['alpha'][0])*field['U'][:,-1].astype(LD)/LD(driver.r0)
        files += [metric_path,field_path,old.inf.incident.FIELDS/f'born-g{order}.npz']
        for n in [64,128]:
            path=PRIMARY/f'sweep-1/photons/complete-{n}.npz';p=np.load(path);stage=p['joint_stage_times'];weights=p['joint_stage_weights']
            assert np.array_equal(p['accepted_angular_times'][:len(stage)],stage)
            j=np.clip(np.searchsorted(t,stage,side='left')-1,0,len(t)-2);w=((stage-t[j])/(t[j+1]-t[j])).astype(LD)
            lum=(1-w[:,None])*L[j]+w[:,None]*L[j+1]
            high=[]
            for now in stage:
                u0=driver.wave(now,np.array([0.]))[0][0];uo=driver.wave(now,driver.xout)[0]
                ext=np.sum(driver.eq.h*((2*driver.zout['Phi']*uo/(driver.rout*driver.zout['b'])).reshape(driver.eq.r.shape)@driver.eq.w),dtype=LD)
                high.append((driver.z0['alpha'][0]/driver.r0+driver.z0['Phi'][0])*u0-ext)
            k=np.array([np.argmin(abs(metric['t']-v)) for v in stage]);assert np.max(abs(metric['t'][k]-stage))<1e-18
            for label,ell in [('high',np.asarray(high,LD)),('low',low_ell[k])]:
                conversion=lum*ell[:,None];energy=np.sum(weights[:,None]*conversion,axis=0,dtype=LD)@(np.arange(1,8,2,dtype=LD)/32)
                np.savez_compressed(OUT/f'{label}-{n}-g{order}.npz',stage_t=stage,stage_weights=weights,background_angular_luminosity=lum,face_log_lapse=ell,physical_energy_conversion=conversion)
                values[label,n,order]=energy;rows.append(dict(component=label,clock=n,geometry_order=order,physical_launch_energy_conversion_erg=float(energy)))
            files.append(path)
    controls={label:dict(time=float(abs(values[label,64,8]-values[label,128,8])/abs(values[label,128,8])),
        quadrature=float(abs(values[label,128,4]-values[label,128,8])/abs(values[label,128,8]))) for label in ['high','low']}
    result=dict(classification='Counterexample candidate',rows=rows,controls=controls,
        passed=all(v['time']<.02 and v['quadrature']<.002 for v in controls.values()),
        actual_applied_low_metric='128-g8 for both physical clocks',
        propagation_work_complete=False,corrected_metric_applied_to_actual_return=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)})
    write(OUT/'input.json',result);print(json.dumps({k:v for k,v in result.items() if k!='bindings'}),flush=True);assert result['passed'],controls
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'input-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
