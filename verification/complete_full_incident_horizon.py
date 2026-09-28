"""Counterexample candidate: bounded full-horizon continuation of the same solve.

Reuse184's exact restart and173's actual packet audit. The model, equations,
local fronts and gates are unchanged. No new source or GR sweep is dispatched.
"""
from pathlib import Path
import fcntl,gc,json,os,resource,shutil,sys,time
import numpy as np
import continue_full_material_front as prior
import continue_native_pressure_response as packet_owner

OUT=Path('native-full-horizon185-work');OLD=prior.OUT
owner=prior.owner;LD=prior.LD;AMP=prior.AMP
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=30,check=15,run=16200)
HISTORY_KEYS=list(prior.HISTORY.values())+['t','moments','radial_ports',
    'photon_history_scaled_occupation','material_history','collision_transfer',
    'accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights',
    'actual_step_edges','energy_offset_reference','energy_offset_t']


def path(label):return OUT/f'sweep-1/photons/{label}.npz'


def status(value):
    target=OUT/'status.json';temp=target.with_suffix('.tmp');write(temp,value);temp.replace(target)


def prepare():
    assert not OUT.exists();OUT.mkdir();reuse={};files=[]
    assert read(OLD/'time-result.json')['passed'] and read(OLD/'prefix-audit.json')['passed']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/pilot-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    for src in inputs:
        dest=OUT/src.relative_to(OLD)
        # Immutable input names never serve as output names in this module.
        os.link(src,dest);reuse[str(dest.relative_to(OUT))]=dict(path=str(src),sha256=sha(src));files.append(src)
    steps={};rates={};overheads=[]
    for n,name in [(64,'coarse'),(128,'fine')]:
        row=read(OLD/f'pilot-{n}-result.json');receipt=read(OLD/f'{name}-receipt.json')
        with np.load(path(f'pilot-{n}')) as p:
            full=int(sum(1+p['split_macro_steps']));done=len(p['actual_step_edges'])-1
            assert row['passed'] and row['completed_steps']==n//16 and 'restart_x' in p
            steps[n]=dict(full=full,completed=done,remaining=full-done)
        rates[n]=row['stepping_seconds']/row['actual_new_steps']
        overheads.append(receipt['seconds']-row['stepping_seconds'])
    overhead=max(overheads)+3.;nominal=sum(steps[n]['remaining']*rates[n] for n in steps)+30*overhead
    upper=1.75*nominal+120;assert upper<CAPS['run']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    files += [Path(__file__),OLD/'time-result.json',OLD/'prefix-audit.json',OLD/'audit-result.json',
        OLD/'restart-check.json',OLD/'dispatch-admission.json',owner.drive.FIELDS/'born-g8.npz']
    write(OUT/'reuse.json',reuse)
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='d4c5fd443',
        previous_turn='Progress: exact joint-state/history restart was repaired and applied to actual continuation. Both paths reachedT/16and passed original time, stage, constitutive and same-solution ledger gates.',
        claim='Finish the SAME corrected-EOS full incident photon/B/S/Etilde/H response over the declared3.434431ms, so its actual stress and energy/boundary history can feed its own GR and charge readout.',
        decision='Proceed only while all six photon/thermal and four material2percent time channels and same-solution ledgers pass at every canonical endpoint. Failure stops the pair, rather than spending the full interval first.',
        reuse='Hard-link184accepted states and original constructor inputs without changing their bytes. Preserve original0.1073ms extra output and the canonical history. Restart atT/16; never recompute the accepted prefix. No old directional response or charge is added.',
        equations='Frozen181full primary+firstBorn incident field,183input-defined one-level front split,184four-variable exact restart, same EOS/background/atmosphere, amplitude, native branch rule, Radau equation and max3Newton. No additional geometry or self-generated GR yet.',
        actual_steps=steps,
        forecast=dict(measured_seconds_per_step=rates,measured_per_path_overhead_seconds=overhead,
            nominal_remaining_seconds=nominal,assumed_range_seconds=[nominal,CAPS['run']],margin_factor=1.75,
            upper_remaining_seconds=upper,
            caveat='Later Newton, Krylov, assembly and I/O costs are unmeasured extrapolations. Recompute using the maximum observed continued cost for each clock after every completed pair. Starting-from-zero183costs are retained separately; current exact restart needs fewer iterations.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,new_sweeps=1,new_pairs=15,
        gates=dict(time=.02,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,max_Newton_solves=3),
        restart='Restore actual internal variables and all four gas-stage/floor/conserved/angular histories. Clear only the current segment LINEAR diagnostic list after its previous values have been retained in the restart row; do not count old solves twice.',
        stop='First physical, time, constitutive, prefix, ledger, cost-admission or4.5hour wall cap failure. Preserve the rejected interval and last accepted pair. No automatic restart, clock/refinement/split/method change, extension beyondT, second sweep or GR run.',
        limits='Retained directional response to declared external input. Full native EOS derivative certification, nonlinear Einstein evolution, boundary/spatial controls, physical static comparison and observational matching remain open. Full-horizon time success alone is not final scientific closure.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    status(dict(state='prepared',completed_shared_intervals=1,total_shared_intervals=16,final_charge_conclusion='unadjudicated'))


def verify(file,prefix):
    port,moments,count=packet_owner.verify_packet(file,prefix)
    with np.load(file) as p,np.load(prefix) as old:
        for key in HISTORY_KEYS:
            assert prior.prior.fingerprint(p[key][:len(old[key])])==prior.prior.fingerprint(old[key]),key
        assert np.array_equal(p['joint_stage_times'],p['accepted_angular_times'])
        assert np.array_equal(p['joint_stage_weights'],p['accepted_angular_quadrature_weights'])
        q=p['conserved_material_history'][-1]/AMP
        final=np.column_stack([q[2]-p['energy_offset_reference'][-1]/AMP,q[3],q[0],q[1]])
        expected=np.sum(p['joint_stage_weights'][:,None,None].astype(LD)*(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled']),axis=0,dtype=LD)
        actual=final+p['material_floor_discard_scaled']
        error=np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))
        assert max(error)<1e-8,('Same-solution local material ledger',error.tolist())
        current=json.loads(str(p['joint_checks_json']));previous=json.loads(str(old['joint_checks_json']))
        for key in ['newton','stages']:assert current[key][:len(previous[key])]==previous[key]
        stage=max(a[-1]['relative'] for a in current['newton']);physical=max(max(a[-1]['moments']) for a in current['newton'])
        assert stage<1e-12 and physical<1e-13
        return dict(angular_port_relative=port,local_material_balance=error.astype(float).tolist(),
            maximum_true_stage=stage,maximum_true_physical_stage=physical,
            maximum_newton_iterations=max(map(len,current['newton'])),actual_completed_steps=count,
            accepted_prefix_preserved=True),moments,p['conserved_material_history'].copy(),p['t'].copy()


def check():
    rows=[]
    for n in [64,128]:
        report,_,_,_=verify(path(f'pilot-{n}'),path(f'pilot-{n}'))
        rows.append(dict(clock=n,**report))
    from fractions import Fraction as F
    for power in range(3):assert F(3,4)*F(1,3)**power+F(1,4)==F(1,power+1)
    write(OUT/'admission-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        new_physical_steps=0,symbolic=dict(classification='Proven',passed=True,scope='Exact Radau moments0..2 only.'),
        final_charge_conclusion='unadjudicated'))


def run():
    assert read(OUT/'admission-check.json')['passed']
    # Native advisory lock prevents concurrent producers and releases on exit.
    lock=(OUT/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (OUT/'failure.json').exists() and not (OUT/'result.json').exists()
    plan=read(OUT/'plan.json');prior.OUT=OUT;prior.initialize()
    identity=dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
        started_unix=time.time(),source_sha256=sha(__file__),plan_sha256=sha(OUT/'plan.json'))
    labels={n:f'pilot-{n}' for n in [64,128]};completed={n:plan['actual_steps'][str(n)]['completed'] for n in labels}
    rates={int(n):v for n,v in plan['forecast']['measured_seconds_per_step'].items()}
    overhead=plan['forecast']['measured_per_path_overhead_seconds'];start=time.monotonic();comparisons=[]
    try:
        for k in range(2,17):
            rows=[];ph=[];gas=[];clocks=[]
            for n in [64,128]:
                status(dict(identity,state='running',active_base_clock=n,target_shared_interval=k,
                    completed_shared_intervals=k-1,total_shared_intervals=16,actual_completed_steps=completed,
                    elapsed_seconds=time.monotonic()-start,cap_seconds=CAPS['run'],final_charge_conclusion='unadjudicated'))
                mark=time.monotonic();m=owner.Model(n);label=f'interval-{k:02d}-{n}'
                m.run.__globals__['LINEAR'].clear()
                row=m.run(n,label,k*(n//16),restart=labels[n])
                assert row['passed'],row
                audit,photons,material,tt=verify(path(label),path(labels[n]))
                wall=time.monotonic()-mark
                row.update(audit,continuation_wall_seconds=wall,full_incident_input_applied=True,
                    self_GR_return_closed=False,full_horizon_completed=k==16,final_charge_conclusion='unadjudicated')
                write(path(label).with_suffix('.json'),row)
                rows.append(row);ph.append(photons);gas.append(material);clocks.append(tt)
                completed[n]=audit['actual_completed_steps'];labels[n]=label
                rates[n]=max(rates[n],row['stepping_seconds']/row['actual_new_steps'])
                overhead=max(overhead,wall-row['stepping_seconds'])
                del m;gc.collect()
            assert np.array_equal(*clocks)
            photon=owner.joint.previous.run.c.relative(*ph);material=owner.joint.previous.run.c.relative(*gas)
            result=dict(classification='Counterexample candidate',passed=max(photon+material)<.02,shared_interval=k,
                photon_time=photon,material_time=material,same_horizon_seconds=float(clocks[0][-1]),rows=rows,
                elapsed_seconds=time.monotonic()-start,actual_completed_steps=completed.copy(),
                final_charge_conclusion='unadjudicated',full_goal_complete=False)
            write(OUT/f'comparison-{k:02d}.json',result);comparisons.append(result)
            print(json.dumps({a:b for a,b in result.items() if a!='rows'}),flush=True)
            assert result['passed'],('Original time gate',k,photon,material)
            remaining={n:plan['actual_steps'][str(n)]['full']-completed[n] for n in labels}
            upper=1.75*(sum(remaining[n]*rates[n] for n in labels)+2*(16-k)*overhead)+120
            elapsed=time.monotonic()-start
            write(OUT/f'cost-{k:02d}.json',dict(classification='Conjectural',elapsed_seconds=elapsed,
                remaining_actual_steps=remaining,maximum_seconds_per_step=rates,maximum_overhead_seconds=overhead,
                assumed_upper_remaining_seconds=upper,eligible=elapsed+upper<CAPS['run']))
            assert elapsed+upper<CAPS['run'],('Cost admission',k,elapsed,upper)
        result.update(full_horizon_completed=True,all_canonical_time_checks_passed=True,
            same_solution_photon_material_history_preserved=True,self_GR_return_closed=False,
            physical_final_charge_solved=False,full_goal_complete=False)
        write(OUT/'result.json',result);status(dict(identity,state='completed',result=result))
    except BaseException as exc:
        failure=dict(classification='Counterexample candidate',passed=False,error=repr(exc),
            last_labels=labels,actual_completed_steps=completed,last_comparison=comparisons[-1] if comparisons else None,
            elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False)
        write(OUT/'failure.json',failure);status(dict(identity,state='failed',result=failure));raise
    finally:lock.close()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            for p,v in read(OUT/'reuse.json').items():assert sha(OUT/p)==v['sha256'],p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
