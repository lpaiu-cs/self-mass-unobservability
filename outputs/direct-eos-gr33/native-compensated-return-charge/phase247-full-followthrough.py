"""Continue the accepted charge readout after its own long returned solution."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

root=Path('native-compensated-charge247-work');out=root/'full';short=root/'check'
dependency=Path('native-dense-returned246-work');producer=Path('verification/read_compensated_return_charge.py')
polynomial=Path('.phase247-polynomial-audit.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    p=root/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)

assert not (root/'full-controller-start.json').exists()
assert read(short/'controller-status.json')['state']=='completed'
assert read(short/'result.json')['charge_comparison_admitted']
assert read(short/'polynomial-audit.json')['passed']
bound={str(p):sha(p) for p in [producer,polynomial,dependency/'controller-start.json',short/'result.json',short/'polynomial-audit.json']}
for p,h in read(short/'polynomial-audit.json')['bindings'].items():assert sha(p)==h,p
entry=read(dependency/'controller-start.json');start=time.monotonic();completed=[];workers={}
write('full-controller-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),
    controller_sha256=sha(__file__),bindings=bound,dependency=entry,parallel_fields={'field1288':0,'field648':8,'field1284':12},
    max_dependency_wait_seconds=43200,field_cap_seconds=7200,virtual_GiB_per_worker=16,
    decision='Continue the same long returned solution once actual source and independent polynomial checks pass. No new material steps or physical refinement. Forecast must fit the existing generous field caps; otherwise preserve evidence and reassess.',
    final_charge_conclusion='unadjudicated',full_goal_complete=False))

def launch(action,cpu=0):
    for p,h in bound.items():assert sha(p)==h,p
    command=[sys.executable,str(polynomial),'full'] if action=='polynomial' else [sys.executable,str(producer),action]
    with (root/f'full-{action}.stdout.log').open('w') as stdout,(root/f'full-{action}.stderr.log').open('w') as stderr:
        return subprocess.Popen(['taskset','-c',str(cpu),*command],stdout=stdout,stderr=stderr)

def sequential(action):
    child=launch(action);workers[action]=child
    write('full-controller-status.json',dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
    code=child.wait();completed.append(dict(action=action,returncode=code));assert code==0,(action,code)

try:
    while True:
        state=read(dependency/'controller-status.json');assert state['state']!='failed',state
        if state['state']=='completed':break
        p=Path('/proc')/str(entry['pid'])
        assert p.exists() and p.joinpath('stat').read_text().split()[21]==entry['process_start_ticks']
        assert entry['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert time.monotonic()-start<43200,'Dependency wait cap'
        write('full-controller-status.json',dict(state='waiting',action='same246long_dense_source',dependency_pid=entry['pid'],elapsed_seconds=time.monotonic()-start));time.sleep(60)
    sequential('prepare');sequential('polynomial')
    assert read(out/'polynomial-audit.json')['passed']
    import numpy as np
    def cuts(folder,n):
        with np.load(folder/f'gr/source-{n}.npz') as d:
            knots=np.unique(np.r_[d['t'],d['geometry_times'],d['metric_times']])
            return int(np.count_nonzero((knots>=0)&(knots<=d['t'][-1]))-1)
    plan=read(out/'plan.json');short_plan=read(short/'plan.json')
    output_ratio=plan['output_times']/short_plan['output_times']
    cut_ratio=max(cuts(dependency/'full',n)/cuts(dependency/'check',n) for n in [64,128])
    receipts=[read(short/f'{action}-receipt.json') for action in ['field1288','field648','field1284']]
    nominal=max(r['seconds'] for r in receipts)*output_ratio*cut_ratio
    upper=2*nominal+300
    rss_upper=1.5*max(r['peak_RSS_bytes'] for r in receipts)*max(output_ratio,cut_ratio)
    forecast=dict(classification='Conjectural',measured_short_receipts=receipts,output_ratio=output_ratio,source_cut_ratio=cut_ratio,
        nominal_longest_worker_seconds=nominal,planning_upper_seconds=upper,estimated_peak_RSS_bytes_per_worker=rss_upper,
        assumption='Work scales with output samples times source intervals; memory with the larger count. Unmeasured full-path scaling, not a guaranteed ETA.',
        passed=upper<7200 and rss_upper<16*1024**3)
    write('full-dispatch-forecast.json',forecast);assert forecast['passed'],forecast
    workers={action:launch(action,cpu) for action,cpu in [('field1288',0),('field648',8),('field1284',12)]}
    write('full-controller-status.json',dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed,elapsed_seconds=time.monotonic()-start))
    while any(v.poll() is None for v in workers.values()):
        failed={k:v.returncode for k,v in workers.items() if v.poll() not in [None,0]}
        assert not failed,failed
        time.sleep(5)
    completed.extend(dict(action=k,returncode=v.returncode) for k,v in workers.items())
    assert not any(v.returncode for v in workers.values()),completed
    for action in ['collect','audit','compare']:sequential(action)
    write('full-controller-status.json',dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,
        conditional_compact_sign_survives_one_return=read(out/'result.json')['conditional_compact_sign_survives_one_return'],final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    for child in workers.values():
        if child.poll() is None:child.terminate()
    for child in workers.values():
        if child.poll() is None:child.wait()
    write('full-controller-status.json',dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
