"""Measured continuation with exact-transform and independent ray gates."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
import numpy as np
out=Path('native-primitive-photon263-work');reference=out/'direct-reference'
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
pilot=read(out/'pilot.json');assert read(out/'pilot-receipt.json')['error'] is None
assert read(out/'check.json')['passed'] and read(out/'symbolic.json')['passed']
assert pilot['forecast_upper_seconds_per_full_path']<28800
assert not (out/'execution-plan.json').exists()
files=[Path(__file__),Path('.phase263-production.py'),Path('.phase262-production.py'),Path('.phase263-direct-reference.py'),
    out/'plan.json',out/'check.json',out/'pilot.json',out/'pilot-receipt.json',out/'method-comparison.json',reference/'plan.json']
keys=['delta_radius_cm','delta_direction','delta_log_H','delta_arrival_seconds','integrated_log_H_work','stress_weak_moments_erg']
for folder,group in read(out/'method-comparison.json')['relative'].items():
    for cell,expected in group.items():
        a=Path(folder)/f'pilot-{cell}.npz';b=out/f'pilot-{cell}.npz';files.extend([a,b])
        x,y=np.load(a),np.load(b)
        for key in keys:
            error=float(np.max(abs(x[key]-y[key]))/max(np.max(abs(y[key])),1e-290))
            assert error==expected[key] and error<.002
bindings=dict(read(out/'plan.json')['bindings']);bindings.update({str(p):sha(p) for p in dict.fromkeys(files)})
boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
pid=14897;proc=Path(f'/proc/{pid}')
if proc.exists():
    assert proc.joinpath('cmdline').read_bytes().split(bytes([0]))[1]==b'.phase263-direct-reference.py'
    ticks=proc.joinpath('stat').read_text().split()[21]
else:
    assert (reference/'receipt.json').exists();ticks=None
write(out/'execution-plan.json',dict(classification='Conjectural',
    claim='Complete actual returned-metric exterior photon stress through the full original period using the verified exact primitive characteristic equations.',
    decision='Only an independent original time-jet trajectory pass and the unchanged0.2percent whole-source controls admit the new exterior stress. No physical charge claim follows before full boundary feedback.',
    original_four_hour_cost_verdict=False,forecast_upper_seconds_per_path=pilot['forecast_upper_seconds_per_full_path'],
    maximum_path_seconds=28800,maximum_concurrent_paths=4,CPU_affinities=[4,6,10,14],
    reason_for_budget='Measured primitive worst cohort80.19s gives a conservative2x full-path bound6.09hours. Preserve the old4hour rejection; allow8hours to finish once with margin under the user instruction to avoid repeated tight-budget stops. No physical grid, horizon, path count or accuracy gate increase.',
    numerical_gates=dict(rtol=2e-8,atol=2e-11,actual_direct_comparison=.002,angular_invariant=1e-10,source_and_boundary_quadrature=.002),
    source_time_control='Not a separate64/128exterior geometry trajectory comparison; this remains required before final physical acceptance.',
    reference_handle=dict(pid=pid,start_ticks=ticks,boot_id=boot),bindings=bindings,
    physical_final_charge_solved=False))
start=time.monotonic();workers=[];done=[]
write(out/'controller-start.json',dict(pid=os.getpid(),start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=boot,started_unix=time.time(),execution_plan_sha256=sha(out/'execution-plan.json')))
try:
    write(out/'pipeline-status.json',dict(state='running',action='independent_reference',pid=pid))
    while ticks is not None:
        try:
            if proc.joinpath('stat').read_text().split()[21]!=ticks:break
        except FileNotFoundError:break
        time.sleep(5)
    receipt=read(reference/'receipt.json');assert receipt['error'] is None,receipt
    assert read(reference/'result.json')['passed']
    done.append(dict(action='independent_reference',passed=True))
    def launch(action,cpu):
        for f,h in bindings.items():assert sha(f)==h,f
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,'.phase263-production.py',action],stdout=stdout,stderr=stderr)
        workers.append(child);return child
    group=[(a,launch(a,c)) for a,c in zip(['fine','angular','temporal','geometry'],[4,6,10,14])]
    write(out/'pipeline-status.json',dict(state='running',action='propagation',workers={a:p.pid for a,p in group},completed=done,started_unix=time.time()))
    pending=dict(group)
    while pending:
        for action,child in list(pending.items()):
            code=child.poll()
            if code is not None:
                done.append(dict(action=action,returncode=code));del pending[action];assert code==0,(action,code)
        if pending:time.sleep(5)
    child=launch('collect',6);code=child.wait();done.append(dict(action='collect',returncode=code));assert code==0
    write(out/'pipeline-status.json',dict(state='completed',completed=done,seconds=time.monotonic()-start,physical_final_charge_solved=False))
except Exception as exc:
    for child in workers:
        if child.poll() is None:child.terminate()
    for child in workers:child.wait()
    write(out/'pipeline-status.json',dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start,physical_final_charge_solved=False))
    raise
