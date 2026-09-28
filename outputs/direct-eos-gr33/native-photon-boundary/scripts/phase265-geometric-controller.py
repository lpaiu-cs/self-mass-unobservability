"""Run the 575-time geometric photon lapse: prepare, pilot, parallel work, collect."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

module=Path('verification/propagate_geometric_boundary_clock.py')
out=Path('native-geometric-clock265-work')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='/home/lpaiu/work/nutimo_pilot/request13_deps:verification')
status=Path('native-geometric-clock265-status.json')
def write(p,v):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
source=sha(module);start=time.monotonic();done=[];workers=[]
def launch(action,cpu):
    assert sha(module)==source
    log=Path(f'.phase265-geometric-{action}')
    with log.with_suffix('.stdout.log').open('w') as o,log.with_suffix('.stderr.log').open('w') as e:
        child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(module),action],stdout=o,stderr=e,env=env)
    workers.append(child);return child
try:
    for action in ['prepare','pilot']:
        write(status,dict(state='running',action=action,completed=done,seconds=time.monotonic()-start))
        code=launch(action,1).wait();done.append(dict(action=action,returncode=code));assert code==0,(action,code)
    plan=read(out/'execution-plan.json')
    group=[(f'work-{k}',launch(f'work-{k}',cpu)) for k,cpu in zip(range(plan['workers']),plan['CPU_affinities'])]
    write(status,dict(state='running',action='work',workers={a:c.pid for a,c in group},completed=done,seconds=time.monotonic()-start))
    pending=dict(group)
    while pending:
        for action,child in list(pending.items()):
            code=child.poll()
            if code is not None:done.append(dict(action=action,returncode=code));del pending[action];assert code==0,(action,code)
        if pending:time.sleep(5)
    code=launch('collect',1).wait();done.append(dict(action='collect',returncode=code));assert code==0
    write(status,dict(state='completed',completed=done,seconds=time.monotonic()-start,result=read(out/'result.json')))
except Exception as exc:
    for child in workers:
        if child.poll() is None:child.terminate()
    for child in workers:child.wait()
    write(status,dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start))
    raise
