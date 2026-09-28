"""Pilot and, only within the registered cost and gates, apply existing controls."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-dynamic-vacuum261-work');read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
producer=Path('verification/propagate_exterior_vacuum.py');collector=Path('.phase261-vacuum-collect.py')
assert read(out/'check.json')['passed']
assert read(out/'check-receipt.json')['source_sha256']==sha(producer)
assert not (out/'controller-start.json').exists()
bound={str(p):sha(p) for p in [producer,collector,Path(__file__),out/'plan.json',out/'check.json']}
start=time.monotonic();workers=[];done=[]
write(out/'controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),bindings=bound,
    scope='Actual incident-metric exterior stress propagation, not full physical charge',CPU_affinities=[4,6,10,14]))
def launch(action,cpu):
    for p,h in bound.items():assert sha(p)==h,p
    cmd=['taskset','-c',str(cpu),sys.executable,str(collector)] if action=='collect' else ['taskset','-c',str(cpu),sys.executable,str(producer),action]
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr)
    workers.append(child);return child
try:
    child=launch('pilot',4)
    write(out/'pipeline-status.json',dict(state='running',action='pilot',child_pid=child.pid))
    code=child.wait();done.append(dict(action='pilot',returncode=code));assert code==0
    assert read(out/'pilot.json')['eligible']
    group=[(a,launch(a,c)) for a,c in zip(['fine','angular','temporal','geometry'],[4,6,10,14])]
    write(out/'pipeline-status.json',dict(state='running',action='propagation',workers={a:p.pid for a,p in group},completed=done))
    pending=dict(group)
    while pending:
        for action,child in list(pending.items()):
            code=child.poll()
            if code is not None:
                done.append(dict(action=action,returncode=code));del pending[action];assert code==0,(action,code)
        if pending:time.sleep(5)
    child=launch('collect',6);write(out/'pipeline-status.json',dict(state='running',action='collect',child_pid=child.pid,completed=done))
    code=child.wait();done.append(dict(action='collect',returncode=code));assert code==0
    write(out/'pipeline-status.json',dict(state='completed',completed=done,seconds=time.monotonic()-start,physical_final_charge_solved=False))
except BaseException as exc:
    for child in workers:
        if child.poll() is None:child.terminate()
    for child in workers:
        if child.poll() is None:child.wait()
    write(out/'pipeline-status.json',dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start,physical_final_charge_solved=False))
    raise
