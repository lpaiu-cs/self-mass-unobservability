"""Parallel independent GR fields, then gated actual coupled return."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-complete-return236-work');producer=Path('verification/return_complete_history_gr.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert json.loads((out/'prepare-receipt.json').read_text())['error'] is None
assert not (out/'controller-start.json').exists();start=time.monotonic();bound=sha(producer);completed=[]
assert set([4,6,8]).issubset(os.sched_getaffinity(0))
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,
    controller_sha256=sha(__file__),parallel_fields={'field1288':4,'field648':6,'field1284':8},
    sequential_actions=['collect','metric','coarse','fine','audit']))


def launch(action,cpu):
    assert sha(producer)==bound
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(producer),action],stdout=stdout,stderr=stderr)


workers={action:launch(action,cpu) for action,cpu in [('field1288',4),('field648',6),('field1284',8)]}
write('controller-status.json',dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed))
for action,child in workers.items():completed.append(dict(action=action,returncode=child.wait()))
if any(v['returncode'] for v in completed):
    write('controller-status.json',dict(state='failed',action='parallel_fields',completed=completed,elapsed_seconds=time.monotonic()-start));sys.exit(1)
for action in ['collect','metric','coarse','fine','audit']:
    child=launch(action,4)
    write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
    code=child.wait();completed.append(dict(action=action,returncode=code))
    if code:
        write('controller-status.json',dict(state='failed',action=action,completed=completed,elapsed_seconds=time.monotonic()-start));sys.exit(code)
write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
