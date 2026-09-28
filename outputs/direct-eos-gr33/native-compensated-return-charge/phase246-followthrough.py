"""Use the same dense adapter after the actual long returned source completes."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-dense-returned246-work');dependency=Path('native-returned-source245-work')
producer=Path('verification/read_dense_returned_source.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
os.sched_setaffinity(0,{8});bound=sha(producer)
assert read(out/'check/result.json')['representation_controls_passed']
assert read(out/'check/source-receipt.json')['source_sha256']==bound
assert not (out/'controller-start.json').exists();start=time.monotonic();completed=[];entry=read(dependency/'controller-start.json')
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,
    controller_sha256=sha(__file__),dependency=entry,CPU_affinity=8,actions=['prepare','geometry','source'],max_dependency_wait_seconds=32400))
try:
    while True:
        assert sha(producer)==bound;state=read(dependency/'controller-status.json');assert state['state']!='failed',state
        if state['state']=='completed':break
        p=Path('/proc')/str(entry['pid'])
        assert p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
        assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert time.monotonic()-start<32400,'Dependency wait cap'
        write('controller-status.json',dict(state='waiting',action='same245long_return_source',dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(60)
    for action in ['prepare','geometry','source']:
        assert sha(producer)==bound
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
            write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start));code=child.wait()
        completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
    write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
