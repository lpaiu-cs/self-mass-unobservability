"""Apply the full-period metric and continue the actual coupled return."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-full-return249-work');producer=Path('verification/align_return_restart_clocks.py')
dependencies=[Path('native-complete-return236-work'),Path('native-retarded-extension248-work')]
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
os.sched_setaffinity(0,{3});bound=sha(producer)
assert read(out/'restart-regression.json')['passed'] and read(out/'check-receipt.json')['source_sha256']==sha('verification/complete_returned_period.py')
assert read(out/'clock-alignment.json')['passed'] and read(out/'clock-check-receipt.json')['source_sha256']==bound
assert not (out/'controller-start.json').exists();entries={str(p):read(p/'controller-start.json') for p in dependencies}
start=time.monotonic();completed=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    controller_sha256=sha(__file__),producer_sha256=bound,dependencies=entries,max_dependency_wait_seconds=43200,
    CPU_affinity=3,virtual_GiB=16,actions=['prepare','metric','coarse','fine','audit'],
    final_charge_conclusion='unadjudicated',full_goal_complete=False))
try:
    while True:
        assert sha(producer)==bound;states={}
        for p in dependencies:
            state=read(p/'controller-status.json');states[str(p)]=state['state'];assert state['state']!='failed',state
            if state['state']=='completed':continue
            entry=entries[str(p)];proc=Path('/proc')/str(entry['pid'])
            assert proc.exists() and proc.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
            assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        if all(s=='completed' for s in states.values()):break
        assert time.monotonic()-start<43200,'Dependency wait cap'
        write('controller-status.json',dict(state='waiting',action='actual236return_and248full_field',dependencies=states,elapsed_seconds=time.monotonic()-start));time.sleep(60)
    for action in ['prepare','metric','coarse','fine','audit']:
        assert sha(producer)==bound
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            child=subprocess.Popen([sys.executable,'.phase251-reader-launch.py',str(producer),action],stdout=stdout,stderr=stderr)
            write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start));code=child.wait()
        completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
        os.link(out/f'{action}-receipt.json',out/f'full/{action}-receipt.json')
    final=dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write('full/controller-status.json',final);write('controller-status.json',final)
except BaseException as exc:
    write('controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
