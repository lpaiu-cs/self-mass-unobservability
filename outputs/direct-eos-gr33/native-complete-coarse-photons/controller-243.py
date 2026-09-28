"""Recover only the missing coarse photons, then join actual saved captures."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
os.sched_setaffinity(0,{2})
out=Path('native-preimage-photon243-work');producer=Path('verification/finish_exact_preimage_photons.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert json.loads((out/'prepare-receipt.json').read_text())['error'] is None
assert not (out/'controller-start.json').exists();start=time.monotonic();completed=[];bound=sha(producer)
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),producer_sha256=bound,
    controller_sha256=sha(__file__),actions=['check','recover','assemble'],CPU_affinity=2))
for action in ['check','recover','assemble']:
    assert sha(producer)==bound
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        child=subprocess.Popen([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
        write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed))
        code=child.wait()
    completed.append(dict(action=action,returncode=code))
    if code:
        write('controller-status.json',dict(state='failed',action=action,completed=completed,elapsed_seconds=time.monotonic()-start));sys.exit(code)
write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,
    final_charge_conclusion='unadjudicated',full_goal_complete=False))
