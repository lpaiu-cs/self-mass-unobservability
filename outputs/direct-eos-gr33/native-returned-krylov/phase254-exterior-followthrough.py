"""Read the completed full-period charge with its own exterior and mass ports."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
root=Path('native-krylov-exterior254-work');dep=Path('native-krylov-charge254-work/full')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=root/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
assert read(root/'common/audit.json')['passed']
assert not (root/'full-controller-start.json').exists()
files=[Path('verification/read_krylov_return_exterior.py'),Path('verification/read_returned_exterior.py'),Path('verification/verify_returned_exterior.py'),
       Path('.phase251-reader-launch.py'),root/'common/result.json',root/'common/audit.json',dep/'controller-start.json']
bound={str(p):sha(p) for p in files};entry=read(dep/'controller-start.json');start=time.monotonic();done=[]
write('full-controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),controller_sha256=sha(__file__),
    bindings=bound,dependency=entry,wait_cap_seconds=43200,action_cap_seconds=1800,CPU_affinity=2,virtual_GiB=16,
    measured_common_seconds=read(root/'common-receipt.json')['seconds'],new_physical_steps=0))
try:
    while True:
        for p,h in bound.items():assert sha(p)==h,p
        state=read(dep/'controller-status.json');assert state['state']!='failed',state
        if state['state']=='completed':break
        proc=Path('/proc')/str(entry['pid'])
        assert proc.exists() and proc.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
        assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert time.monotonic()-start<43200,'Dependency wait cap'
        write('full-controller-status.json',dict(state='waiting',dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(60)
    for module in ['read_krylov_return_exterior.py']:
        for p,h in bound.items():assert sha(p)==h,p
        with (root/(module+'.stdout.log')).open('w') as out,(root/(module+'.stderr.log')).open('w') as err:
            child=subprocess.Popen(['taskset','-c','2',sys.executable,'.phase251-reader-launch.py','verification/'+module,'full'],stdout=out,stderr=err)
            write('full-controller-status.json',dict(state='running',action=module,child_pid=child.pid,completed=done,elapsed_seconds=time.monotonic()-start))
            rc=child.wait();done.append(dict(action=module,returncode=rc));assert rc==0,(module,rc)
    write('full-controller-status.json',dict(state='completed',completed=done,elapsed_seconds=time.monotonic()-start,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    write('full-controller-status.json',dict(state='failed',completed=done,error=repr(exc),elapsed_seconds=time.monotonic()-start));raise
