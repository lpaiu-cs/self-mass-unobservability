"""Continue the accepted corrected solution, checking each canonical interval.

Counterexample candidate. Reuse the actual6/12-step prefixes. Stop promptly
on the original time/physical gates; never run an unchecked full fine path.
"""
from pathlib import Path
import gc,json,os,resource,shutil,sys,time
import numpy as np
import resolve_native_pressure_front as prior

OUT=Path('native-pressure-return173-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
PREPARE_CAP=30;RUN_CAP=4800
KEYS=[0,1,2,3,5,6]


def status(value):
    path=OUT/'status.json';tmp=path.with_suffix('.tmp')
    write(tmp,value);tmp.replace(path)


def prepare():
    assert not OUT.exists();accepted=read(OLD/'result.json');assert accepted['passed'];OUT.mkdir()
    files=[Path(__file__),Path(prior.__file__),Path(prior.prior.__file__),Path(prior.radau.__file__),
        Path(prior.original.__file__),OLD/'result.json',OLD/'grid-check.json',OLD/'plan.json',OLD/'grid-sum-repair.json']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    for n in [64,128]:
        for suffix in ['.npz','.json']:
            p=OLD/f'sweep-1/photons/pilot-{n}{suffix}';q=OUT/p.relative_to(OLD)
            shutil.copyfile(p,q);assert sha(p)==sha(q);files.append(p)
    rows=accepted['rows'];grid=read(OLD/'grid-check.json')['rows'];pilot=read(OLD/'sweep-1/photons/pilot.json')
    step=max(r['stepping_seconds']/r['actual_new_steps'] for r in rows)
    setup=max(0.,pilot['seconds']-sum(r['seconds'] for r in rows))/2
    overhead=setup+max(r['seconds']-r['stepping_seconds'] for r in rows)
    remaining=sum(g['full_planned_actual_steps']-r['actual_completed_steps'] for g,r in zip(grid,rows))
    nominal=remaining*step+30*overhead;upper=1.75*nominal;assert upper<RUN_CAP
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d2361a1ad6adf24732250736b9650faadcc2010c',
        claim='Finish the SAME corrected direct-source photon/thermal/H response over the declared3.434431ms interval, retaining all physical ledgers for the subsequent free-material/GR return.',
        decision='Complete both paths only while every successive shared canonical interval passes the original six2percent time gates and physical checks. This is an actual continuation toward final charge, not a separate diagnostic trajectory. No final charge can be accepted from photons alone.',
        reuse='Copy the accepted1726/12-step states byte-for-byte; restart rather than recompute them. Same source, fixed local bisection rule, spatial grid, physical input, background, stage solver and amplitudes.',
        actual_steps=dict(base=[64,128],full=[g['full_planned_actual_steps'] for g in grid],completed=[r['actual_completed_steps'] for r in rows],remaining=remaining),
        budget=dict(run_seconds=RUN_CAP,prepare_seconds=PREPARE_CAP,CPU_threads=1,virtual_GiB=3,pairlets=15,max_sweeps=1),
        forecast=dict(measured_seconds_per_actual_step=step,measured_constructor_seconds=setup,measured_per_path_overhead_seconds=overhead,
            nominal_remaining_seconds=nominal,margin_factor=1.75,upper_remaining_seconds=upper,
            limits='Later Krylov/assembly/IO rates are extrapolations, not guarantees. Rebuild each clock model serially to remain within3GiB; update the cost admission after every actual pairlet.'),
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12),
        stop='Stop at the first failed physical/time gate, cost admission or80minute hard cap. Preserve last accepted checkpoints and the failed interval. No new grid, extra bisection, automatic restart, extra sweep or source change.',
        consumers='All packet checks use saved actual stage times, step edges and quadrature weights. Legacy uniform SDIRK/Radau packet consumers are not admitted for this path. Free material and retarded GR must later use this same accepted solution and boundary history.',
        limits='One simultaneous photon/thermal/H increment with zero previous free-motion input and extra metric. Free-material reciprocity, additional GR, same-solution energy/mass closure and final charge remain required; full EOS/derivative/spatial/boundary/nonlinear/static/observational limits stay open.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    status(dict(state='prepared',classification='Counterexample candidate',completed_shared_intervals=1,total_shared_intervals=16,
        final_charge_conclusion='unadjudicated',forecast=read(OUT/'plan.json')['forecast']))


def verify_packet(path,prefix):
    with np.load(path) as p,np.load(prefix) as old:
        edges=p['actual_step_edges'];dt=np.diff(edges);tt=p['accepted_angular_times'];w=p['accepted_angular_quadrature_weights']
        assert np.all(dt>0) and np.max(abs(tt-(edges[:-1,None]+dt[:,None]*prior.radau.RK_C).ravel()))<1e-18
        assert np.max(abs(w-(dt[:,None]*prior.radau.RK_B).ravel()))<1e-18
        packets=w[:,None]*p['accepted_angular_luminosity'];aw=np.arange(1,8,2)/32
        cumulative=np.r_[0.,np.cumsum(packets@aw,dtype=np.longdouble)]
        ids=np.array([np.argmin(abs(edges-t)) for t in p['t']]);assert np.max(abs(edges[ids]-p['t']))<1e-18
        err=float(np.max(abs(cumulative[2*ids]-p['radial_ports'][:,1,1]))/max(np.sum(abs(packets)@aw),1e-290));assert err<1e-12,err
        for key in ['t','moments','material_history','photon_history_scaled_occupation','collision_transfer','radial_ports',
                    'accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','actual_step_edges']:
            assert np.array_equal(p[key][:len(old[key])],old[key]),('Accepted prefix changed',key)
        return err,p['moments'][:,KEYS].copy(),len(edges)-1


def run():
    plan=read(OUT/'plan.json')
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    lock=OUT/'run.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
    prior.OUT=OUT;prior.initialize();folder=OUT/'sweep-1/photons';start=time.monotonic()
    identity=dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],started_unix=time.time(),source_sha256=sha(__file__))
    labels={n:f'pilot-{n}' for n in [64,128]};completed={64:6,128:12};comparisons=[];timings=[]
    maximum_step=plan['forecast']['measured_seconds_per_actual_step'];maximum_overhead=plan['forecast']['measured_per_path_overhead_seconds']
    try:
        for k in range(2,17):
            pair=[];data=[]
            for n in [64,128]:
                status(dict(identity,state='running',classification='Counterexample candidate',active_base_clock=n,
                    target_shared_interval=k,completed_shared_intervals=k-1,total_shared_intervals=16,
                    actual_completed_steps=completed,elapsed_seconds=time.monotonic()-start,cap_seconds=RUN_CAP,
                    final_charge_conclusion='unadjudicated'))
                mark=time.monotonic();m=prior.original.c.Response(n)
                label=f'steps-{n}-reference-128' if k==16 else f'interval-{k:02d}-{n}'
                row=m.run(n,label,k*(n//16),restart=labels[n]);wall=time.monotonic()-mark
                assert row['passed'],row
                err,moments,count=verify_packet(folder/f'{label}.npz',folder/f'pilot-{n}.npz')
                row.update(actual_angular_quadrature_relative=err,accepted_prefix_preserved=True,continuation_wall_seconds=wall)
                write(folder/f'{label}.json',row);pair.append(row);data.append(moments);completed[n]=count;labels[n]=label
                maximum_step=max(maximum_step,row['stepping_seconds']/row['actual_new_steps'])
                maximum_overhead=max(maximum_overhead,wall-row['stepping_seconds']);timings.append(wall)
                del m;gc.collect()
            errors=np.max(np.sum(abs(data[0]-data[1]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(data[1]),axis=2),axis=0),1e-290)
            row=dict(classification='Counterexample candidate',shared_interval=k,passed=bool(max(errors)<.02),time_comparison=errors.astype(float).tolist(),rows=pair,
                elapsed_seconds=time.monotonic()-start,actual_completed_steps=completed.copy())
            write(OUT/f'comparison-{k:02d}.json',row);comparisons.append(row);print(json.dumps({key:value for key,value in row.items() if key!='rows'}),flush=True)
            assert row['passed'],('Time gate',k,row['time_comparison'])
            remaining=sum(v-completed[n] for n,v in zip([64,128],plan['actual_steps']['full']))
            upper=1.75*(remaining*maximum_step+2*(16-k)*maximum_overhead)
            write(OUT/f'cost-{k:02d}.json',dict(elapsed_seconds=time.monotonic()-start,remaining_actual_steps=remaining,
                upper_remaining_seconds=upper,cap_seconds=RUN_CAP,measured_max_seconds_per_step=maximum_step,measured_max_overhead=maximum_overhead))
            assert time.monotonic()-start+upper<RUN_CAP,('Cost admission',k,time.monotonic()-start,upper)
        result=dict(classification='Counterexample candidate',passed=True,full_horizon_photon_thermal_H_completed=True,
            actual_steps=completed,rows=comparisons[-1]['rows'],time_comparison=comparisons[-1]['time_comparison'],
            accepted_prefixes_reused=True,all_canonical_time_checks_passed=True,actual_packet_weights_verified=True,
            free_material_response_completed=False,physical_final_charge_solved=False,full_goal_complete=False,
            final_charge_conclusion='unadjudicated until the SAME free-material/GR solution and energy-boundary closure are completed',seconds=time.monotonic()-start)
        write(folder/'result.json',result);write(OUT/'result.json',result);status(dict(identity,state='completed',result=result))
    except BaseException as exc:
        failure=dict(classification='Counterexample candidate',passed=False,error=repr(exc),last_labels=labels,
            actual_completed_steps=completed,last_comparison=comparisons[-1] if comparisons else None,elapsed_seconds=time.monotonic()-start,
            full_horizon_photon_thermal_H_completed=False,physical_final_charge_solved=False,full_goal_complete=False,final_charge_conclusion='unadjudicated')
        write(OUT/'failure.json',failure);status(dict(identity,state='failed',result=failure));raise


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','run'];receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));prior.original.inf.incident.native.deadline(PREPARE_CAP if action=='prepare' else RUN_CAP)
    start=time.monotonic();cpu=time.process_time();error=None
    try:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
