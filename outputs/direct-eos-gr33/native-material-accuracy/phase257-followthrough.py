"""Complete the original two returned paths with bounded inner refinement."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-material-accuracy257-work');producer=Path('verification/finish_returned_material_accuracy.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
assert read(out/'prepare-receipt.json')['error'] is None
assert not (out/'controller-start.json').exists()
os.sched_setaffinity(0,{3});bound=sha(producer);started=time.monotonic();done=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    producer_sha256=bound,controller_sha256=sha(__file__),CPU_affinity=3,virtual_GiB=16,
    actions=['check','fine','audit'],budgets=dict(check=600,fine=21600,audit=600),scientific_gates_changed=False))
try:
    for action in ['check','fine','audit']:
        assert sha(producer)==bound
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen([sys.executable,'.phase251-reader-launch.py',str(producer),action],stdout=stdout,stderr=stderr)
            write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=done,elapsed_seconds=time.monotonic()-started))
            code=child.wait();done.append(dict(action=action,returncode=code));assert code==0,(action,code)
    write('controller-status.json',dict(state='completed',completed=done,elapsed_seconds=time.monotonic()-started,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('controller-status.json',dict(state='failed',completed=done,error=repr(exc),elapsed_seconds=time.monotonic()-started));raise
