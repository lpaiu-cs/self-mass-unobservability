"""Continue the identified live pilot, then run only its eligible four paths."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-returned-moments262-work')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
pilot=494; ticks='39789'; boot='83309146-1468-49d9-928a-350c68b74e90'
assert Path('/proc/sys/kernel/random/boot_id').read_text().strip()==boot
assert Path(f'/proc/{pilot}/stat').read_text().split()[21]==ticks
assert Path(f'/proc/{pilot}/cmdline').read_bytes().split(bytes([0]))[1:3]==[b'verification/propagate_returned_moments.py',b'pilot']
assert not (out/'execution-plan.json').exists()
files=[Path(__file__),Path('.phase262-production.py'),Path('.phase262-check-repair.py'),
    out/'plan.json',out/'check.json',out/'arithmetic-check.json',out/'check-repair-receipt.json']
bindings={str(p):sha(p) for p in files};bindings.update(read(out/'plan.json')['bindings'])
write(out/'execution-plan.json',dict(classification='Conjectural',
    claim='Complete actual returned-metric exterior photon stress on the original period, after the existing three-cell pilot passes identity gates and a four-hour per-path forecast.',
    pilot=dict(pid=pilot,start_ticks=ticks,boot_id=boot),
    paths=['fine','angular','temporal','geometry'],CPU_affinities=[4,6,10,14],
    maximum_concurrent_paths=4,maximum_path_seconds=14400,maximum_pilot_seconds=1800,
    expected_wall_time='Pending measured pilot. No production if forecast_upper_seconds_per_full_path >=14400.',
    gates=dict(energy_identity=.002,angular_invariant=1e-10,quadrature=.002),
    exclusions='Separate returned source-clock convergence, reciprocal scalar and background operator, complete boundary feedback and final physical charge remain open.',
    bindings=bindings))
workers=[];done=[];start=time.monotonic()
write(out/'controller-start.json',dict(pid=os.getpid(),start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=boot,started_unix=time.time(),execution_plan_sha256=sha(out/'execution-plan.json')))
try:
    write(out/'pipeline-status.json',dict(state='running',action='pilot',pid=pilot))
    while Path(f'/proc/{pilot}').exists():
        if Path(f'/proc/{pilot}/stat').read_text().split()[21]!=ticks:break
        time.sleep(5)
    receipt=read(out/'pilot-receipt.json');assert receipt['error'] is None,receipt
    done.append(dict(action='pilot',error=None))
    if not read(out/'pilot.json')['eligible']:
        write(out/'pipeline-status.json',dict(state='cost_review_required',completed=done,pilot=read(out/'pilot.json')))
        sys.exit(0)
    def launch(action,cpu):
        for f,h in bindings.items():assert sha(f)==h,f
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,'.phase262-production.py',action],stdout=stdout,stderr=stderr)
        workers.append(child);return child
    group=[(a,launch(a,c)) for a,c in zip(['fine','angular','temporal','geometry'],[4,6,10,14])]
    write(out/'pipeline-status.json',dict(state='running',action='propagation',workers={a:p.pid for a,p in group},completed=done))
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
