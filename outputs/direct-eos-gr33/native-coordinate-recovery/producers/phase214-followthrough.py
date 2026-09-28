"""Complete two independent saved-history recoveries within their frozen caps."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib,json,os,subprocess,sys,time

out=Path('native-coordinate-photon214-work');producer=Path('verification/resume_coordinate_exact_photons.py')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
assert not (out/'controller-start.json').exists()
assert json.loads((out/'restart-check.json').read_text())['passed']
assert int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))>18*1024**2
source_sha=sha(producer);start=time.monotonic();jobs={};handles=[];completed=[]
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),source_sha256=sha(__file__),producer_sha256=source_sha,actions=['coarse','fine','audit'],parallel_paths=2,stop_policy='First failed path cancels the other owned reconstruction; both preserve their last accepted atomic checkpoint.'))
try:
    for action in ['coarse','fine']:
        assert sha(producer)==source_sha
        stdout=(out/f'{action}.stdout.log').open('w');stderr=(out/f'{action}.stderr.log').open('w');handles.extend([stdout,stderr])
        jobs[action]=subprocess.Popen([sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
    write('controller-status.json',dict(state='running',workers={a:p.pid for a,p in jobs.items()},completed=[]))
    failed=False
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(p.wait):action for action,p in jobs.items()}
        for future in as_completed(pending):
            action=pending[future];code=future.result();completed.append(dict(action=action,returncode=code))
            if code and not failed:
                failed=True
                for other,p in jobs.items():
                    if other!=action and p.poll() is None:p.terminate()
            write('controller-status.json',dict(state='failed' if failed else 'running',workers={a:p.pid for a,p in jobs.items()},completed=completed,elapsed_seconds=time.monotonic()-start))
    if failed:sys.exit(1)
    assert sha(producer)==source_sha
    with (out/'audit.stdout.log').open('w') as stdout,(out/'audit.stderr.log').open('w') as stderr:
        code=subprocess.call([sys.executable,str(producer),'audit'],stdout=stdout,stderr=stderr)
    completed.append(dict(action='audit',returncode=code))
    write('controller-status.json',dict(state='completed' if code==0 else 'failed',completed=completed,elapsed_seconds=time.monotonic()-start,full_declared_period_recovered=False,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    sys.exit(code)
finally:
    for h in handles:h.close()
