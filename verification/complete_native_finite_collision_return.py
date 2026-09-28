"""Finish the original finite collision return within the reallocated total cap."""
from pathlib import Path
import gc,hashlib,json,math,signal,sys,time

BASE=Path('outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/conserved-history-charge/collision-forcing/response/material-charge/finite/zero-exact/compensated/motion-feedback')
SOURCE=BASE/'finite-collision'/'extended';APPLIED=SOURCE/'applied';OUT=APPLIED/'completed-return'
CAPS=dict(source=500,photon_pilot=90,photon_production_and_sum=662,material_pilot=60,material_production=340,sources=40,charge=90,audit=49)
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def prepare():
    assert not OUT.exists();OUT.mkdir()
    previous=json.loads(Path('outputs/direct-eos-gr33/native-finite-collision-remainder-manifest.json').read_text())
    pilot=json.loads((APPLIED/'cached-input-pilot.json').read_text());failure=json.loads((APPLIED/'cached_sample-failure.json').read_text())
    assert not pilot['eligible'] and not (SOURCE/'production.json').exists()
    saved=previous['available_knots'];missing=previous['missing_knots'];point=pilot['point']['seconds']
    forecast=2*point*len(missing)+(failure['seconds']-point)+60+10
    old_total=2421.523589098;charged=previous['charged_readout_action_seconds']+previous['seconds']
    assert forecast<CAPS['source'] and charged+sum(CAPS.values())<old_total
    files=[Path(__file__),Path('verification/propagate_native_finite_collision_remainder.py'),Path('verification/repair_native_finite_collision_remainder.py'),
           Path('verification/verify_native_finite_collision_return.py'),APPLIED/'cached-input-pilot.json',APPLIED/'cached_sample-failure.json',BASE/'result.json',BASE/'return/gr/result.json']
    files+=[SOURCE/f'point-{k}.{suffix}' for k in saved for suffix in ['json','npz']]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1f94ddc48',previous_goal_turn='progress: sampled finite remainder and resource failure changed the execution plan',
        claim='Complete the remaining actual finite collision sources, propagate the added source through simultaneous photons/E/H, apply the new TOTAL photon transfers to finite material, then read the actual combined GR charge.',
        decision='Determine the charge change produced by the previously unpropagated finite collision remainder. Preserve the interpolation, native-derivative, continuous residual and nonlinear/exterior limitations.',
        reassessment='The previous source attempt is closed and its failed forecasts remain unchanged. Reallocate unused downstream time within the same aggregate allowance; do not repeat completed states or another independent diagnostic pilot.',
        reused_knots=saved,remaining_knots=missing,source_forecast_seconds=forecast,
        forecast_assumptions='Twice the latest cached point cost for each remaining point, the measured last initialization/writing overhead,60s for imports and10s additional archive margin. Host IO remains variable. One process avoids repeated imports between all stages.',
        budgets=CAPS,earlier_charged_action_seconds=charged,total_previous_allowance_seconds=old_total,new_aggregate_allowance_seconds=charged+sum(CAPS.values()),
        allocations='Source500s; photon production including superposition662s instead of900; material production340s instead of450. Keep the90/60s pilots. Remaining source/audit89s accounts for the previous saved audit. GR90s unchanged. No added total time.',
        execution='Start the source deadline before heavy imports. Clamp every legacy alarm to the current stage deadline, including legacy disarm calls; only the outer owner can finish the stage. Reuse one Python process and immutable-byte cache. CPU1, virtual memory3GiB.',
        photon_dispatch='Require the existing equal-horizon pilot and original numerical gates. Add50s to its2x production forecast for total-photon superposition; require the sum below662s.',
        material_dispatch='Use existing late-call-count forecast and two actual prefix steps; require its2x forecast below340s. Reuse the successful prefixes.',
        stop='One pipeline attempt. Stop on any original scientific gate, stage deadline, or forecast. No smaller amplitude, extra horizon/path/grid, weakened tolerance or automatic extra waveform iteration.',
        unchanged_model='Same531 cells,8 angles,152 frequencies,17 original source knots,3.434431ms,64/128 clocks on background128. One deferred collision correction, not fully nonlinear radiation or a uniform EOS error bound.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in files}))
    print(json.dumps(dict(prepared=True,remaining=missing,source_forecast_seconds=forecast,remaining_stage_caps=sum(CAPS.values()))),flush=True)


def reallocate():
    assert not (OUT/'reallocation-plan.json').exists()
    failure=json.loads((OUT/'failure.json').read_text());assert failure['stage']=='photon_pilot'
    p=json.loads((OUT/'photon-pilot-662.json').read_text());assert not p['eligible']
    assert all(r['passed'] for r in p['rows']) and max(p['equal_horizon'])<.02
    caps={k:v for k,v in CAPS.items() if k not in ['source','photon_pilot']};caps['photon_production_and_sum']=900
    plan=json.loads((OUT/'plan.json').read_text())
    assert failure['seconds']+sum(caps.values())<sum(CAPS.values())
    required=p['upper_seconds']+50+60;assert required<900
    p.update(eligible=True,previous_662s_eligible=False,production_and_sum_cap_seconds=900,
             production_and_sum_forecast=required,restart_import_allowance_seconds=60)
    write(APPLIED/'pilot.json',p)
    paths=[Path(__file__),OUT/'plan.json',OUT/'second-producer.py',OUT/'failure.json',OUT/'photon-pilot-662.json',
           SOURCE/'production.json',APPLIED/'pilot.json']+[APPLIED/f'pilot-{n}.npz' for n in [64,128]]
    plan.update(checkpoint='1f94ddc48',budgets=caps,earlier_pipeline_seconds=failure['seconds'],
        reassessment='All17 sources completed in224.461726827s, leaving275.538273173s of their500s allowance. The two photon prefixes passed every numerical gate, but their2x forecast739.5652834s plus50s sum allowance exceeded662s. Preserve that failure and both prefixes. Transfer238s of the unused source allowance to restore a900s production+sum cap, including60s restart import margin. No original state or prefix replay; same one physical correction.',
        remaining_forecast_seconds=required,remaining_caps_seconds=sum(caps.values()),
        aggregate_action_cap_seconds=plan['earlier_charged_action_seconds']+failure['seconds']+sum(caps.values()),
        stop='Resume the accepted photon prefixes once, then the previously registered material/GR stages. No extra correction, grid, horizon, tolerance change or repeat of completed stages.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths})
    assert plan['aggregate_action_cap_seconds']<plan['total_previous_allowance_seconds']
    write(OUT/'reallocation-plan.json',plan);print(json.dumps(dict(prepared=True,forecast=required,cap=900,remaining_caps=sum(caps.values()))),flush=True)


def execute(resume=False):
    import resource
    result_path=OUT/('resume-run.json' if resume else 'run.json');failure_path=OUT/('resume-failure.json' if resume else 'failure.json')
    assert not result_path.exists() and not failure_path.exists()
    plan=json.loads((OUT/('reallocation-plan.json' if resume else 'plan.json')).read_text());caps=plan['budgets']
    for path,h in plan['bindings'].items():assert sha(path)==h,path
    first=next(iter(caps));started=time.monotonic();real_alarm=signal.alarm;deadline=started+caps[first];rows=[]
    def timeout(*_):raise TimeoutError('Phase144 registered stage deadline')
    def bounded_alarm(seconds):
        left=deadline-time.monotonic()
        if left<=0:raise TimeoutError('Phase144 registered stage deadline')
        return real_alarm(max(1,math.ceil(min(left,seconds) if seconds else left)))
    signal.signal(signal.SIGALRM,timeout);signal.alarm=bounded_alarm
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));real_alarm(caps[first])
    stage=first
    try:
        import propagate_native_finite_collision_remainder as run
        import_seconds=time.monotonic()-started
        for stage,cap in caps.items():
            begin=started if stage==first else time.monotonic();deadline=begin+cap;bounded_alarm(cap)
            print(json.dumps(dict(stage=stage,event='begin',cap_seconds=cap)),flush=True)
            if stage=='source':
                s=run.repair.State()
                for k in plan['remaining_knots']:
                    assert not (SOURCE/f'point-{k}.json').exists();assert s.point(k)['passed']
                points=[json.loads((SOURCE/f'point-{k}.json').read_text()) for k in range(17)];assert all(p['passed'] for p in points)
                write(SOURCE/'production.json',dict(classification='Counterexample candidate',passed=True,rows=points,
                    saved_points_reused=plan['reused_knots'],new_states=len(plan['remaining_knots']),seconds=time.monotonic()-begin,
                    import_seconds=import_seconds,maximum_remainder_over_linear=max(p['rows'][0]['remainder_over_linear'] for p in points),uniform_remainder_bound=False))
                del s
            elif stage=='photon_pilot':
                run.photon_pilot();p=json.loads((APPLIED/'pilot.json').read_text())
                p.update(production_and_sum_forecast=p['upper_seconds']+50,production_and_sum_cap_seconds=662)
                p['eligible']=p['eligible'] and p['production_and_sum_forecast']<662;write(APPLIED/'pilot.json',p);assert p['eligible']
            elif stage=='photon_production_and_sum':
                run.photon_production();p=json.loads((APPLIED/'result.json').read_text())
                p.update(actual_finite_motion_in_photons=False,increment_only=True,prescribed_finite_motion_in_saved_base=True)
                write(APPLIED/'result.json',p);run.combine()
            elif stage=='material_pilot':
                run.material_pilot();p=json.loads((run.RETURN/'pilot.json').read_text());p.update(eligible=p['upper_seconds']<340,production_cap_seconds=340)
                write(run.RETURN/'pilot.json',p);assert p['eligible']
            elif stage=='audit':
                import verify_native_finite_collision_return as audit
                original=run.write
                def separate(path,value):
                    p=Path(path)
                    if p.resolve().parent==APPLIED.resolve() and p.name in ['symbolic.json','audit-plan.json','audit.json']:p=OUT/p.name
                    return original(p,value)
                run.write=separate
                try:audit.main()
                finally:run.write=original
            else:getattr(run,stage)()
            gc.collect();row=dict(stage=stage,seconds=time.monotonic()-begin,cap_seconds=cap,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
            assert row['seconds']<=cap;rows.append(row);write(OUT/f'{stage}-receipt.json',row);print(json.dumps(row),flush=True)
        write(result_path,dict(classification='Counterexample candidate',passed=True,stages=rows,seconds=time.monotonic()-started,all_stage_imports_charged=True,
            finite_remainder_propagated_to_GR=True,uniform_nonlinear_remainder_bound=False,final_charge_solved=False,full_goal_complete=False))
    except Exception as exc:
        write(failure_path,dict(stage=stage,error=repr(exc),completed_stages=rows,seconds=time.monotonic()-started));raise
    finally:real_alarm(0);signal.alarm=real_alarm


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    elif action=='reallocate':reallocate()
    else:
        assert action in ['run','resume'];execute(action=='resume')
