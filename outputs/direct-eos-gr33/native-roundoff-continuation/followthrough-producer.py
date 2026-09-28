"""Run the two authorized remaining paths once, with their saved gates."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-roundoff-polish197-work')
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert not (out/'controller-start.json').exists()
assert json.loads((out/'restart-check.json').read_text())['passed']
start=time.time();receipts=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=start,
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),actions=['coarse','fine','audit']))
with (out/'controller.stdout.log').open('w') as stdout,(out/'controller.stderr.log').open('w') as stderr:
    for action in ['coarse','fine','audit']:
        write('controller-status.json',dict(state='running',action=action,completed=receipts,elapsed_seconds=time.time()-start))
        child=subprocess.run([sys.executable,'verification/polish_full_native_roundoff.py',action],stdout=stdout,stderr=stderr)
        receipts.append(dict(action=action,returncode=child.returncode))
        if child.returncode:
            write('controller-status.json',dict(state='failed',action=action,completed=receipts,elapsed_seconds=time.time()-start));sys.exit(child.returncode)
write('controller-status.json',dict(state='completed',completed=receipts,elapsed_seconds=time.time()-start,
    final_charge_conclusion='unadjudicated',full_goal_complete=False))


