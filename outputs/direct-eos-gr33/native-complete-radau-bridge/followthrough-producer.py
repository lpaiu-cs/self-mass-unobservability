"""Consume the completed accepted photon history after its live producer ends."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-complete-radau224-control');up=Path('native-captured-photon223-work')
check=Path('native-complete-radau224-check-work');producer=Path('verification/read_complete_radau_history.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def write(name,value):
    p=out/name;tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(p)

assert read(check/'result.json')['numerical_controls_passed']
assert read(check/'sources.json')['representation_controls_passed']
assert read(check/'driver-polynomial-audit.json')['passed']
assert not out.exists();out.mkdir();start=time.monotonic();bound=sha(producer);owner=read(up/'controller-start.json')
boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip();assert boot==owner['boot_id']
write('start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=boot,started_unix=time.time(),
    producer_sha256=bound,controller_sha256=sha(__file__),upstream=owner,
    actions=['wait_for_accepted223','prepare','endpoints','source','fields'],wait_cap_seconds=10800,
    scope='Same accepted common15/16history only. Source/field time failures remain failures. No physical return or final charge admission. No worker is restarted automatically.'))
state='waiting';completed=[]
try:
    while True:
        try:status=read(up/'controller-status.json')
        except json.JSONDecodeError:status={}
        if status.get('state')=='completed':break
        if status.get('state')=='failed':raise RuntimeError('Original photon history recovery failed; no GR dispatch')
        assert time.monotonic()-start<10800,'Upstream wait cap'
        stat=Path(f"/proc/{owner['pid']}/stat");assert stat.exists(),'Upstream controller handle missing'
        assert stat.read_text().split()[21]==owner['process_start_ticks'],'Upstream controller handle changed'
        assert Path('/proc/sys/kernel/random/boot_id').read_text().strip()==owner['boot_id'],'Upstream boot changed'
        write('status.json',dict(state='waiting',upstream_pid=owner['pid'],elapsed_seconds=time.monotonic()-start))
        time.sleep(60)
    assert read(up/'result.json')['same_solution_accepted_history_recovered']
    for action in ['prepare','endpoints','source','fields']:
        assert sha(producer)==bound,'Consumer source changed'
        write('status.json',dict(state='running',action=action,completed=completed,elapsed_seconds=time.monotonic()-start))
        with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
            code=subprocess.call([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
        completed.append(dict(action=action,returncode=code))
        assert code==0,(action,code)
    write('status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('status.json',dict(state='stopped',error=repr(exc),completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
