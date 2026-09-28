"""Run the accepted same-history GR return through both original clocks."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-independent-fine232-work');producer=Path('verification/continue_independent_fine_native.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert json.loads((out/'restart-check.json').read_text())['passed']
assert not (out/'controller-start.json').exists();start=time.monotonic();completed=[];bound=sha(producer)
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,controller_sha256=sha(__file__),actions=['fine']))
for action in ['fine']:
    assert sha(producer)==bound
    write('controller-status.json',dict(state='running',action=action,completed=completed,elapsed_seconds=time.monotonic()-start))
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        code=subprocess.call([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
    completed.append(dict(action=action,returncode=code))
    if code:
        write('controller-status.json',dict(state='failed',action=action,completed=completed,elapsed_seconds=time.monotonic()-start));sys.exit(code)
write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
