"""Evolve the one-return photon-boundary pair with the user-approved H material gate, then its same-solution readers."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-photon-boundary-hgate265-work')
charge=Path('native-photon-boundary-hgate-charge265-work/full')
exterior=Path('native-photon-boundary-hgate-exterior265-work')
producer=Path('verification/apply_photon_geometric_boundary_hgate.py')
launcher=Path('.phase251-reader-launch.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,value):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,p)

assert read(out/'prepare-receipt.json')['error'] is None
assert read(out/'check-receipt.json')['error'] is None and read(out/'boundary-check.json')['passed']
assert not (out/'controller-start.json').exists()
bound={str(p):sha(p) for p in [producer,launcher,out/'charge-reader.py',out/'exterior-reader.py',Path(__file__),Path('.phase247-polynomial-audit.py')]}
start=time.monotonic();done=[];workers=[]
entry=dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),bindings=bound,
    CPU_affinities=dict(coarse=4,fine=6,readers=6,fields=[4,12,14],exterior=2),physical_paths=[119,231],
    concurrent_physical_paths=True,same_solution_readout_automatic=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
write(out/'controller-start.json',entry)
def launch(module,action,folder,cpu):
    for p,h in bound.items():assert sha(p)==h,p
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/f'{action}.stdout.log').open('w') as stdout,(folder/f'{action}.stderr.log').open('w') as stderr:
        child=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(launcher),str(module),action],stdout=stdout,stderr=stderr)
    workers.append(child);return child
def run(module,action,folder,cpu):
    child=launch(module,action,folder,cpu)
    row=dict(state='running',action=action,child_pid=child.pid,completed=done,elapsed_seconds=time.monotonic()-start)
    write(out/'pipeline-status.json',row);write(folder/'controller-status.json',row)
    code=child.wait();done.append(dict(module=str(module),action=action,returncode=code));assert code==0,(module,action,code)
def group(module,actions,folder):
    children=[(action,launch(module,action,folder,cpu)) for action,cpu in actions]
    write(out/'pipeline-status.json',dict(state='running',action='parallel',workers={a:c.pid for a,c in children},completed=done,elapsed_seconds=time.monotonic()-start))
    pending=dict(children)
    while pending:
        for action,child in list(pending.items()):
            code=child.poll()
            if code is not None:
                done.append(dict(module=str(module),action=action,returncode=code));del pending[action];assert code==0,(module,action,code)
        if pending:time.sleep(5)
try:
    group(producer,[('coarse',4),('fine',6)],out)
    run(producer,'audit',out,6)
    write(out/'controller-status.json',dict(state='completed',completed=list(done),elapsed_seconds=time.monotonic()-start))
    write(charge/'controller-start.json',entry)
    reader=out/'charge-reader.py'
    for action in ['endpoint_prepare','endpoint_source','dense_prepare','dense_geometry','dense_source','charge_prepare','charge_polynomial']:
        run(reader,action,charge,6)
    assert read(charge/'charge/polynomial-audit.json')['passed']
    group(reader,[('charge_field1288',4),('charge_field648',12),('charge_field1284',14)],charge)
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
