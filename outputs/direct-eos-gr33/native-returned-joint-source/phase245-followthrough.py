"""Consume the actual completed GR return, never a provisional snapshot."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-returned-source245-work');dependency=Path('native-complete-return236-work')
producer=Path('verification/read_returned_joint_source.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())


def write(name,value):
    p=out/name;temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');os.replace(temp,p)


os.sched_setaffinity(0,{6});bound=sha(producer)
assert read(out/'check/result.json')['representation_controls_passed']
assert read(out/'check/source-receipt.json')['source_sha256']==bound
assert not (out/'controller-start.json').exists();start=time.monotonic();completed=[]
entry=read(dependency/'controller-start.json')
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    producer_sha256=bound,controller_sha256=sha(__file__),dependency=entry,CPU_affinity=6,
    actions=['prepare','source'],max_dependency_wait_seconds=32400,
    final_charge_conclusion='unadjudicated',full_goal_complete=False))
try:
    while True:
        assert sha(producer)==bound;state=read(dependency/'controller-status.json')
        assert state['state']!='failed',state
        if state['state']=='completed':break
        p=Path('/proc')/str(entry['pid'])
        assert p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
        assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert time.monotonic()-start<32400,'Dependency wait cap'
        write('controller-status.json',dict(state='waiting',action='actual236paired_return',
            dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start,completed=completed))
        time.sleep(60)
    for action in ['prepare','source']:
        assert sha(producer)==bound;out.joinpath('full').mkdir(exist_ok=True) if action=='source' else None
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
            write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,
                elapsed_seconds=time.monotonic()-start,completed=completed));code=child.wait()
        completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
    write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('controller-status.json',dict(state='failed',error=repr(exc),completed=completed,elapsed_seconds=time.monotonic()-start,
        final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
