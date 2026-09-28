"""Finish the corrected original pair and its same-solution charge/mass audit."""
from pathlib import Path
import hashlib,json,os,signal,subprocess,sys,time

out=Path('native-short-return259-work')
charge=Path('native-short-return-charge259-work/full')
exterior=Path('native-short-return-exterior259-work')
producer=Path('verification/continue_returned_short_krylov.py')
launcher=Path('.phase251-reader-launch.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,value):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,p)

assert read(out/'prepare-receipt.json')['error'] is None
assert read(out/'check-receipt.json')['error'] is None
assert not (out/'controller-start.json').exists()
bound={str(p):sha(p) for p in [producer,launcher,out/'charge-reader.py',out/'exterior-reader.py',Path(__file__),Path('.phase247-polynomial-audit.py')]}
start=time.monotonic();done=[];workers=[]
entry=dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),bindings=bound,
    CPU_affinity=4,physical_paths=[119,231],same_solution_readout_automatic=True,
    final_charge_conclusion='unadjudicated',full_goal_complete=False)
write(out/'controller-start.json',entry)
def launch(module,action,folder,cpu):
    for p,h in bound.items():assert sha(p)==h,p
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/f'{action}.stdout.log').open('w') as stdout,(folder/f'{action}.stderr.log').open('w') as stderr:
        child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(launcher),str(module),action],stdout=stdout,stderr=stderr)
    workers.append(child);return child
def run(module,action,folder,cpu=4):
    child=launch(module,action,folder,cpu)
    row=dict(state='running',action=action,child_pid=child.pid,completed=done,elapsed_seconds=time.monotonic()-start)
    write(out/'pipeline-status.json',row);write(folder/'controller-status.json',row)
    if action=='coarse':
        old=Path('native-left-metric258-work');identity=read(old/'controller-start.json')
        old_status=read(old/'pipeline-status.json');assert old_status['state']=='running' and old_status['action']=='coarse'
        old_pid=old_status['child_pid'];stat=Path(f'/proc/{old_pid}/stat').read_text().split()
        old_ticks=stat[21];assert old_pid==161206 and old_ticks=='990419'
        assert identity['boot_id']==entry['boot_id']=='abb8c3bc-36fe-4e31-929d-1d66b6717368'
        superseded=False
        while child.poll() is None:
            if not superseded:
                rows=[read(p) for p in (out/'sweep-1/photons').glob('interval-*-64.json')]
                latest=max((r['actual_completed_steps'] for r in rows if r['passed']),default=0)
                old_count=read(old/'capture-64.json')['actual_steps']
                if latest>=old_count and latest>read(out/'plan.json')['reused_coarse_actual_steps']:
                    current=read(old/'pipeline-status.json')
                    assert current['state']=='running' and current['action']=='coarse' and current['child_pid']==old_pid
                    assert Path(f'/proc/{old_pid}/stat').read_text().split()[21]==old_ticks
                    command=Path(f'/proc/{old_pid}/cmdline').read_bytes().replace(bytes([0]),b' ').decode()
                    assert 'verification/repair_returned_metric_endpoint.py coarse' in command
                    record=dict(classification='Counterexample candidate',reason='Measured identical-system15.68x speedup; the new accepted canonical prefix has caught the old actual progress. Intentional cost supersession, not a physical gate failure.',old_pid=old_pid,old_start_ticks=old_ticks,boot_id=entry['boot_id'],old_actual_steps=old_count,new_canonical_steps=latest,old_source_and_plan_unchanged=True,scientific_gates_changed=False,time_unix=time.time())
                    write(out/'supersession-intent.json',record)
                    os.kill(old_pid,signal.SIGTERM);superseded=True
                    write(out/'supersession.json',record)
            time.sleep(5)
    code=child.wait();done.append(dict(module=str(module),action=action,returncode=code));assert code==0,(module,action,code)
try:
    for action in ['coarse','fine','audit']:run(producer,action,out)
    write(out/'controller-status.json',dict(state='completed',completed=list(done),elapsed_seconds=time.monotonic()-start))
    write(charge/'controller-start.json',entry)
    reader=out/'charge-reader.py'
    for action in ['endpoint_prepare','endpoint_source','dense_prepare','dense_geometry','dense_source','charge_prepare','charge_polynomial']:
        run(reader,action,charge,6)
    assert read(charge/'charge/polynomial-audit.json')['passed']
    group=[(action,launch(reader,action,charge,cpu)) for action,cpu in [('charge_field1288',4),('charge_field648',10),('charge_field1284',14)]]
    write(out/'pipeline-status.json',dict(state='running',action='parallel_charge_fields',workers={a:p.pid for a,p in group},elapsed_seconds=time.monotonic()-start))
    for action,child in group:
        code=child.wait();done.append(dict(action=action,returncode=code));assert code==0,(action,code)
    for action in ['charge_collect','charge_audit','charge_compare']:run(reader,action,charge,6)
    write(charge/'controller-status.json',dict(state='completed',elapsed_seconds=time.monotonic()-start))
    run(out/'exterior-reader.py','full',exterior,2)
    result=read(exterior/'full/audit.json');assert result['passed']
    write(exterior/'controller-status.json',dict(state='completed',elapsed_seconds=time.monotonic()-start))
    write(out/'pipeline-status.json',dict(state='completed',completed=done,elapsed_seconds=time.monotonic()-start,
        compact_sign_survives=read(charge/'charge/result.json')['conditional_compact_sign_survives_one_return'],
        frozen_exterior_sign_survives=result['conditional_compact_sign_survives_frozen_exterior'],
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
    for child in workers:
        if child.poll() is None:child.terminate()
    for child in workers:
        if child.poll() is None:child.wait()
    failure=dict(state='failed',completed=done,error=repr(exc),elapsed_seconds=time.monotonic()-start,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'pipeline-status.json',failure)
    if not (out/'result.json').exists():write(out/'controller-status.json',failure)
    raise
