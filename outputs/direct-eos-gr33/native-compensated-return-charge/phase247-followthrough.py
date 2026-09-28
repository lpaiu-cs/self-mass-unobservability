"""Parallel independent short returned-charge fields, then the original audit."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-compensated-charge247-work/check');producer=Path('verification/read_compensated_return_charge.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
assert read(out/'prepare-receipt.json')['error'] is None;assert not (out/'controller-start.json').exists()
start=time.monotonic();bound=sha(producer);completed=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,
    controller_sha256=sha(__file__),parallel_fields={'field1288':0,'field648':8,'field1284':12},sequential_actions=['collect','audit','compare']))
def launch(action,cpu):
    assert sha(producer)==bound
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(producer),'check_'+action],stdout=stdout,stderr=stderr)
try:
    workers={a:launch(a,c) for a,c in [('field1288',0),('field648',8),('field1284',12)]}
    write('controller-status.json',dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed))
    for action,child in workers.items():completed.append(dict(action=action,returncode=child.wait()))
    assert not any(v['returncode'] for v in completed),completed
    for action in ['collect','audit','compare']:
        child=launch(action,0);write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
        code=child.wait();completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
    write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
