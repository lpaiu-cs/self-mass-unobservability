"""Finish the unfinished preregistered direct ray, then run the registered phase-263 propagation.

The phase-263 execution plan, its failed pipeline status, the one-hour reference
receipt and the saved packet 0 stay unchanged. Gates, tolerances, grids, period
and path count are the registered ones. Production still starts only after the
independent reference passes.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
import numpy as np
out=Path('native-primitive-photon263-work');old=out/'direct-reference';reference=out/'direct-reference-264'
status=out/'pipeline-status-264.json';plan_path=out/'execution-plan-264.json'
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
REFERENCE_SECONDS=21600;PATH_SECONDS=43200;ACTIONS=['fine','angular','temporal','geometry'];CPUS=[4,6,10,14]
previous=read(out/'execution-plan.json');failed=read(out/'pipeline-status.json')
assert failed['state']=='failed' and failed['completed']==[] and 'TimeoutError' in failed['error']
assert 'TimeoutError' in read(old/'receipt.json')['error'] and not (old/'result.json').exists()
for f,h in previous['bindings'].items():assert sha(f)==h,f
pilot=read(out/'pilot.json');assert read(out/'pilot-receipt.json')['error'] is None
assert read(out/'check.json')['passed'] and read(out/'symbolic.json')['passed'] and read(out/'method-comparison.json')['passed']
assert pilot['forecast_upper_seconds_per_full_path']<previous['maximum_path_seconds']<PATH_SECONDS
assert not plan_path.exists() and not status.exists() and not reference.exists() and not (out/'result.json').exists()
assert not any((out/a).exists() or (out/f'{a}-receipt.json').exists() for a in ACTIONS+['collect'])
keys=['delta_radius_cm','delta_direction','delta_log_H','delta_arrival_seconds','integrated_log_H_work','stress_weak_moments_erg']
for folder,group in read(out/'method-comparison.json')['relative'].items():
    for cell,expected in group.items():
        x,y=np.load(Path(folder)/f'pilot-{cell}.npz'),np.load(out/f'pilot-{cell}.npz')
        for key in keys:
            error=float(np.max(abs(x[key]-y[key]))/max(np.max(abs(y[key])),1e-290))
            assert error==expected[key] and error<.002
files=[Path(__file__),Path('.phase264-direct-reference.py'),Path('.phase264-production.py'),Path('.phase263-production.py'),
    out/'execution-plan.json',out/'pipeline-status.json',out/'controller-start.json',
    old/'plan.json',old/'receipt.json',old/'progress.json',old/'packet-0.json',old/'packet-0.npz']
bindings=dict(previous['bindings']);bindings.update({str(p):sha(p) for p in files})
boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
write(plan_path,dict(classification='Conjectural',claim=previous['claim'],decision=previous['decision'],
    continuation_of='Phase263 execution-plan.json. Its failed pipeline, one-hour reference receipt and saved packet 0 are preserved unchanged.',
    failure='The registered 3600 s direct-reference budget stopped packet 31 after packet 0 completed in 1820.997 s. No production path had started.',
    reference_budget_seconds=REFERENCE_SECONDS,maximum_path_seconds=PATH_SECONDS,maximum_concurrent_paths=4,CPU_affinities=CPUS,
    forecast_upper_seconds_per_path=pilot['forecast_upper_seconds_per_full_path'],
    reason_for_budget='The one-hour reference cap was not loosened with the eight-hour path cap. Allow 6 hours for the single unfinished ray (packet 0 cost 1821 s; packet 31 exceeded 1771 s). Production cannot resume partial snapshots, so a stop forces a full rerun: raise the cap to 12 hours, above the unchanged 6.09 hour conservative forecast. No grid, horizon, path count, tolerance or gate change.',
    numerical_gates=previous['numerical_gates'],source_time_control=previous['source_time_control'],
    stop_rule='If the reference fails its gate or cap, do not start production. If a path fails, stop all paths and re-plan; do not extend automatically.',
    bindings=bindings,physical_final_charge_solved=False))
start=time.monotonic();workers=[];done=[]
write(out/'controller-start-264.json',dict(pid=os.getpid(),start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=boot,
    started_unix=time.time(),execution_plan_sha256=sha(plan_path)))
def launch(script,args,cpu,label):
    for f,h in bindings.items():assert sha(f)==h,f
    with (out/f'{label}.stdout.log').open('w') as stdout,(out/f'{label}.stderr.log').open('w') as stderr:
        child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,script,*args],stdout=stdout,stderr=stderr)
    workers.append(child);return child
try:
    child=launch('.phase264-direct-reference.py',[],CPUS[0],'direct-reference-264')
    write(status,dict(state='running',action='independent_reference',pid=child.pid,started_unix=time.time()))
    code=child.wait();done.append(dict(action='independent_reference',returncode=code))
    receipt=read(reference/'receipt.json');assert code==0 and receipt['error'] is None,receipt
    assert read(reference/'result.json')['passed']
    group=[(a,launch('.phase264-production.py',[a],c,a)) for a,c in zip(ACTIONS,CPUS)]
    write(status,dict(state='running',action='propagation',workers={a:w.pid for a,w in group},completed=done,started_unix=time.time()))
    pending=dict(group)
    while pending:
        for action,child in list(pending.items()):
            code=child.poll()
            if code is not None:
                done.append(dict(action=action,returncode=code));del pending[action];assert code==0,(action,code)
        if pending:time.sleep(5)
    child=launch('.phase264-production.py',['collect'],CPUS[1],'collect');code=child.wait()
    done.append(dict(action='collect',returncode=code));assert code==0
    write(status,dict(state='completed',completed=done,seconds=time.monotonic()-start,physical_final_charge_solved=False))
except Exception as exc:
    for child in workers:
        if child.poll() is None:child.terminate()
    for child in workers:child.wait()
    write(status,dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start,physical_final_charge_solved=False))
    raise
