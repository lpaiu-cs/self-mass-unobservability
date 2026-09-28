"""Complete corrected fine stages, then gate the original paired time audit."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
os.sched_setaffinity(0,{10})
out=Path('native-common-arithmetic239-work');producer=Path('verification/complete_common_arithmetic_fine.py')
dependency=Path('native-true-momentum238-work')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert json.loads((out/'restart-check.json').read_text())['passed']
assert not (out/'controller-start.json').exists();start=time.monotonic();completed=[];bound=sha(producer)
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,
    controller_sha256=sha(__file__),actions=['fine','audit'],CPU_affinity=10))
for action in ['fine','audit']:
    assert sha(producer)==bound
    if action=='audit':
        wait_start=time.monotonic();entry=json.loads((dependency/'controller-start.json').read_text())
        while not (dependency/'prefix-receipt.json').exists():
            state=json.loads((dependency/'controller-status.json').read_text());p=Path('/proc')/str(entry['pid'])
            live=p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
            live=live and entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
            if state['state']!='running' or not live or time.monotonic()-wait_start>3600:
                write('controller-status.json',dict(state='failed',action='prefix_dependency',completed=completed,
                    dependency_state=state,live=live,elapsed_seconds=time.monotonic()-start));sys.exit(1)
            write('controller-status.json',dict(state='running',action='wait_for_live_prefix',completed=completed,
                dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(30)
    write('controller-status.json',dict(state='running',action=action,completed=completed,elapsed_seconds=time.monotonic()-start))
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        code=subprocess.call([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
    completed.append(dict(action=action,returncode=code))
    if code:
        write('controller-status.json',dict(state='failed',action=action,completed=completed,elapsed_seconds=time.monotonic()-start));sys.exit(code)
write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,
    final_charge_conclusion='unadjudicated',full_goal_complete=False))
