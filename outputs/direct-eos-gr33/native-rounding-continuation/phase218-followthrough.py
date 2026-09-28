"""Finish the original pair from saved116; read its own GR only on acceptance."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-rounding-continuation218-work')
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert not (out/'controller-start.json').exists()
assert json.loads((out/'restart-check.json').read_text())['passed']
producer=Path('verification/continue_rounding_aware_native.py');source_sha=sha(producer)
reader=Path('verification/read_completed_joint_gr.py');reader_sha=sha(reader)
start=time.monotonic();receipts=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),source_sha256=sha(__file__),producer_sha256=source_sha,reader_sha256=reader_sha,actions=['coarse','fine','audit']))
with (out/'controller.stdout.log').open('w') as stdout,(out/'controller.stderr.log').open('w') as stderr:
    for action in ['coarse','fine','audit']:
        assert sha(producer)==source_sha
        write('controller-status.json',dict(state='running',action=action,completed=receipts,elapsed_seconds=time.monotonic()-start))
        child=subprocess.run([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
        receipts.append(dict(action=action,returncode=child.returncode))
        if child.returncode:
            write('controller-status.json',dict(state='failed',action=action,completed=receipts,elapsed_seconds=time.monotonic()-start));sys.exit(child.returncode)
    write('controller-status.json',dict(state='completed',completed=receipts,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    for action in ['prepare','source','fields']:
        assert sha(reader)==reader_sha
        write('gr-controller-status.json',dict(state='running',action=action))
        child=subprocess.run([sys.executable,str(reader),action]+([str(out)] if action=='prepare' else []),stdout=stdout,stderr=stderr)
        if child.returncode:
            write('gr-controller-status.json',dict(state='failed',action=action,returncode=child.returncode));sys.exit(child.returncode)
    write('gr-controller-status.json',dict(state='completed',final_charge_conclusion='unadjudicated',full_goal_complete=False))
