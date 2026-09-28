"""Resume accepted material prefixes with an unchanged physical-step budget."""
from pathlib import Path
import json,resource,sys,time
import resolve_native_material_flux_tangent as prior

run=prior.run;OUT=prior.OUT;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=10,production=450,residual=30)
ORIGINAL_READ=run.read;ORIGINAL_WRITE=run.write


def prepare():
    pilot=read(prior.paths(1)[1]/'pilot.json')
    assert pilot['passed'] and not pilot['eligible'];forecasts=[]
    for row in pilot['rows']:
        old=read(Path('native-mixed-return162-work/sweep-1/material')/f"steps-{row['steps']}-reference-128.json")
        # The old directional RHS calls raw twice per point. The analytic
        # flux tangent calls it once; raw-call count is no longer a work unit.
        cost=max(old['seconds']/old['substeps'],row['seconds']/row['substeps'])
        setup=max(old['worker_wall_seconds']-old['seconds'],row['worker_wall_seconds']-row['seconds'])
        remain=old['substeps']-row['substeps'];assert remain>0
        forecasts.append(dict(clock=row['steps'],old_full_substeps=old['substeps'],accepted_substeps=row['substeps'],
                              seconds_per_substep=cost,remaining_seconds=remain*cost+setup+5))
    upper=2*sum(r['remaining_seconds'] for r in forecasts)
    assert upper<450,upper
    files=[Path(__file__),Path(prior.__file__),Path(prior.face.__file__),Path(run.__file__),
           OUT/'flux-tangent-plan.json',OUT/'flux-tangent-check.json',OUT/'flux-tangent-symbolic.json',
           OUT/'flux-tangent-pilot-receipt.json',prior.paths(1)[1]/'pilot.json']
    for n in [64,128]:files.extend([prior.paths(1)[1]/f'pilot-{n}.npz',prior.paths(1)[0]/f'steps-{n}-reference-128.npz'])
    write(OUT/'material-continuation-plan.json',dict(classification='Counterexample candidate',eligible=True,
        preserved_budget_rejection_seconds=pilot['upper_remaining_seconds'],upper_remaining_seconds=upper,forecasts=forecasts,
        correction='Use actual SSP substeps as the common work unit across the old two-raw and new analytic-flux operators. Retain the larger old/new measured cost, old full-period late CFL count, setup and2x margin. This changes neither450s production cap nor640s total cap.',
        assumptions='The current later physical CFL and per-step cost can differ. Hard wall/memory caps and unchanged numerical gates remain active; stop rather than enlarge them.',
        claim='Apply the accepted same173photon history to the full free-material trajectory, preserving its new accepted4/8prefix arrays.',
        arithmetic='Nominal constitutive partial step is1e-5 with half/double controls. Legacy field nominal_arithmetic_probe16described the replaced global conservative probe and is corrected at the metadata write boundary.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in files}))
    print(json.dumps(dict(eligible=True,forecast_seconds=upper/2,upper_seconds=upper,cap=450)),flush=True)


def current_read(p):
    value=ORIGINAL_READ(p)
    if Path(p)==prior.paths(1)[1]/'pilot.json':
        assert value['passed'] and not value['eligible']
        value=dict(value,eligible=read(OUT/'material-continuation-plan.json')['eligible'])
    return value


def metadata(p,value):
    if 'nominal_arithmetic_probe' in value:
        value=dict(value);value['legacy_global_probe_label']=value.pop('nominal_arithmetic_probe')
        value['constitutive_log_step']=1e-5;value['gross_flux_difference_used']=False
        value['forward_probe_indicator']=None
    ORIGINAL_WRITE(p,value)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'material-continuation-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+read(OUT/'probe-localization.json')['seconds']
        assert spent+CAPS[action]<=640,(spent,CAPS[action])
        if action=='prepare':prepare()
        else:
            for p,h in read(OUT/'material-continuation-plan.json')['bindings'].items():assert sha(p)==h,p
            for p,h in read(OUT/'flux-tangent-plan.json')['bindings'].items():assert sha(p)==h,p
            run.initialize=prior.initialize;run.paths=prior.paths;run.read=current_read;run.write=metadata
            if action=='production':run.material(False)
            else:run.residual()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
