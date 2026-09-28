"""Complete full-period primary GR by extending the verified causal prefix."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
out=Path('native-retarded-extension248-work');dependency=Path('native-full-captured244-work')
producer=Path('verification/extend_unchanged_source_prefix.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=out/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
assert read(out/'regression.json')['passed'] and read(out/'prefix-check.json')['passed'];assert read(out/'prefix-check-receipt.json')['source_sha256']==sha(producer)
assert not (out/'controller-start.json').exists();entry=read(dependency/'controller-start.json')
start=time.monotonic();bound=sha(producer);completed=[];workers={}
write('controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    controller_sha256=sha(__file__),producer_sha256=bound,dependency=entry,max_dependency_wait_seconds=43200,
    parallel_fields={'field1288':1,'field648':5,'field1284':9},field_cap_seconds=7200,virtual_GiB_per_worker=16,
    final_charge_conclusion='unadjudicated',full_goal_complete=False))
def launch(action,cpu=1):
    assert sha(producer)==bound
    with (out/f'{action}.stdout.log').open('w') as stdout,(out/f'{action}.stderr.log').open('w') as stderr:
        return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(producer),action],stdout=stdout,stderr=stderr)
def sequential(action):
    child=launch(action);workers[action]=child
    write('controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
    code=child.wait();completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
try:
    while True:
        state=read(dependency/'controller-status.json');assert state['state']!='failed',state
        if state['state']=='completed':break
        p=Path('/proc')/str(entry['pid'])
        assert p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
        assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert time.monotonic()-start<43200,'Dependency wait cap'
        write('controller-status.json',dict(state='waiting',action='same244full_source',dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(60)
    sequential('prepare')
    import numpy as np
    def cuts(folder,n):
        with np.load(folder/f'gr/source-{n}.npz') as d:
            C=2.99792458e10;knots=np.unique(np.r_[d['t'],d['geometry_times'],d['drive_times'],-d['drive_x']/C,d['drive_duration']-d['drive_x']/C])
            return int(np.count_nonzero((knots>=0)&(knots<=d['t'][-1]))-1)
    plan=read(out/'plan.json');reg=read(out/'regression.json')
    cut_ratio=max(cuts(dependency,n)/cuts(Path('native-complete-radau224-work'),n) for n in [64,128])
    query_ratio=(plan['new_output_times']+2)/5
    nominal=max(r['seconds'] for r in reg['rows'])*query_ratio*cut_ratio
    upper=2*nominal+300;memory=read(out/'check-receipt.json')['peak_RSS_bytes']*1.5*cut_ratio
    forecast=dict(classification='Conjectural',measured_five_query_seconds=[r['seconds'] for r in reg['rows']],
        new_output_times=plan['new_output_times'],source_cut_ratio=cut_ratio,query_ratio=query_ratio,
        nominal_longest_worker_seconds=nominal,planning_upper_seconds=upper,estimated_peak_RSS_bytes_per_worker=memory,
        assumption='Conservative scaling includes fixed setup in every query; full tail speed is unmeasured.',passed=upper<7200 and memory<16*1024**3)
    write('dispatch-forecast.json',forecast);assert forecast['passed'],forecast
    workers={a:launch(a,c) for a,c in [('field1288',1),('field648',5),('field1284',9)]}
    write('controller-status.json',dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed,elapsed_seconds=time.monotonic()-start))
    while any(v.poll() is None for v in workers.values()):
        failed={k:v.returncode for k,v in workers.items() if v.poll() not in [None,0]};assert not failed,failed;time.sleep(5)
    completed.extend(dict(action=k,returncode=v.returncode) for k,v in workers.items());assert not any(v.returncode for v in workers.values()),completed
    for action in ['collect','audit']:sequential(action)
    write('controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    for child in workers.values():
        if child.poll() is None:child.terminate()
    for child in workers.values():
        if child.poll() is None:child.wait()
    write('controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
