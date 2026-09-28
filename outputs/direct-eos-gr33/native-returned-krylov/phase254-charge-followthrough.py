"""Continue the same actual full return all the way through its charge reader."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
import numpy as np
root=Path('native-krylov-charge254-work');check=sys.argv[1]=='check';prefix='check_' if check else ''
out=root/('check-retry' if check else 'full');out.mkdir(parents=True,exist_ok=True)
producer=Path('verification/read_krylov_return_charge.py');dependency=Path('native-returned-krylov254-work')
launcher=Path('.phase251-reader-launch.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
assert not (out/'controller-start.json').exists()
bound={str(p):sha(p) for p in [producer,launcher,Path('.phase247-polynomial-audit.py')]}
if not check:
    assert read(root/'regression.json')['passed'] and read(root/'check-retry/controller-status.json')['state']=='completed'
    assert read(root/'regression-receipt.json')['source_sha256']==sha("verification/read_full_return_charge.py")
    for p,h in read(root/'regression.json')['bindings'].items():assert sha(p)==h,p
    bound[str(root/'regression.json')]=sha(root/'regression.json')
    entry=read(dependency/'controller-start.json');bound[str(dependency/'controller-start.json')]=sha(dependency/'controller-start.json')
start=time.monotonic();completed=[];workers={}
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    controller_sha256=sha(__file__),bindings=bound,dependency=None if check else entry,
    max_dependency_wait_seconds=43200,field_cap_seconds=10800,virtual_GiB_per_process=16,
    parallel_field_CPUs=[4,10,14],physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False))
def launch(action,cpu=6):
    for p,h in bound.items():assert sha(p)==h,p
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(launcher),str(producer),prefix+action],stdout=stdout,stderr=stderr)
def sequential(action):
    child=launch(action);workers[action]=child
    write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
    code=child.wait();completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
try:
    if not check:
        while True:
            state=read(dependency/'controller-status.json');assert state['state']!='failed',state
            if state['state']=='completed':break
            p=Path('/proc')/str(entry['pid'])
            assert p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
            assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
            assert time.monotonic()-start<43200,'Dependency wait cap'
            write('controller-status.json',dict(state='waiting',action='same254full_actual_return',dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(60)
    for action in ['endpoint_prepare','endpoint_source','dense_prepare','dense_geometry','dense_source','charge_prepare','charge_polynomial']:sequential(action)
    assert read(out/'charge/polynomial-audit.json')['passed']
    if not check:
        reference=Path('native-compensated-charge247-work');dense=Path('native-dense-returned246-work')
        suffix='full' if all((reference/f'full/{a}-receipt.json').exists() for a in ['field1288','field648','field1284']) else 'check'
        reference/=suffix;dense/=suffix
        def cuts(folder,n):
            with np.load(folder/f'gr/source-{n}.npz') as d:
                knots=np.unique(np.r_[d['t'],d['geometry_times'],d['metric_times']]);return int(np.count_nonzero((knots>=0)&(knots<=d['t'][-1]))-1)
        output_ratio=read(out/'charge/plan.json')['output_times']/read(reference/'plan.json')['output_times']
        cut_ratio=max(cuts(out/'dense',n)/cuts(dense,n) for n in [64,128])
        receipts=[read(reference/f'{a}-receipt.json') for a in ['field1288','field648','field1284']]
        assert all(r['error'] is None for r in receipts)
        nominal=max(r['seconds'] for r in receipts)*output_ratio*cut_ratio;upper=2*nominal+300
        rss=1.5*max(r['peak_RSS_bytes'] for r in receipts)*max(output_ratio,cut_ratio)
        forecast=dict(classification='Conjectural',reference=str(reference),receipts=receipts,output_ratio=output_ratio,source_cut_ratio=cut_ratio,
            nominal_longest_worker_seconds=nominal,planning_upper_seconds=upper,estimated_peak_RSS_bytes_per_worker=rss,
            assumption='Measured identical operator scaled by output samples times source cuts; later cost is unmeasured, not a guaranteed ETA.',
            passed=upper<10800 and rss<16*1024**3)
        write('dispatch-forecast.json',forecast);assert forecast['passed'],forecast
    workers={a:launch(a,c) for a,c in [('charge_field1288',4),('charge_field648',10),('charge_field1284',14)]}
    write('controller-status.json',dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed,elapsed_seconds=time.monotonic()-start))
    while any(v.poll() is None for v in workers.values()):
        failed={k:v.returncode for k,v in workers.items() if v.poll() not in [None,0]};assert not failed,failed;time.sleep(5)
    completed.extend(dict(action=k,returncode=v.returncode) for k,v in workers.items())
    assert not any(v.returncode for v in workers.values()),completed
    for action in ['charge_collect','charge_audit','charge_compare']+(['regression'] if check else []):sequential(action)
    write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,
        conditional_compact_sign_survives_one_return=read(out/'charge/result.json')['conditional_compact_sign_survives_one_return'],
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    for child in workers.values():
        if child.poll() is None:child.terminate()
    for child in workers.values():
        if child.poll() is None:child.wait()
    write('controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
