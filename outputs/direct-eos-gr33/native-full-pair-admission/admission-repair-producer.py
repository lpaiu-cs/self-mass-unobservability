"""Repair a missing audit input; rerun the unchanged verifier, never physics."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time
out=Path('native-full-admission250-work');target=Path('native-common-arithmetic239-work');coarse=Path('native-true-momentum238-work')
producer=Path('verification/complete_common_arithmetic_fine.py')
read=lambda p:json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()
def write(p,v):
    p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,indent=2)+'\n');os.replace(q,p)
def stopped(folder):
    state=read(folder/'controller-status.json');entry=read(folder/'controller-start.json')
    assert state['state']=='failed' and not state.get('completed'),state
    proc=Path('/proc')/str(entry['pid'])
    assert not proc.exists() or proc.joinpath('stat').read_text().split()[21]!=entry['process_start_ticks'],'Original controller still live'
    return state,entry
assert not (out/'plan.json').exists();out.mkdir(exist_ok=True);start=time.monotonic()
failed=read(target/'controller-status.json');assert failed['state']=='failed' and failed['action']=='audit'
assert read(target/'fine-receipt.json')['error'] is None
assert "interval-15-64.npz" in (target/'audit.stderr.log').read_text()
assert read(target/'audit-receipt.json')['source_sha256']==sha(producer)
files=[target/f'sweep-1/photons/complete-{n}.npz' for n in [64,128]]
physical={str(p):sha(p) for p in files}
for name in ['controller-status.json','controller-start.json','audit-receipt.json','audit.stdout.log','audit.stderr.log','result.json','prefix-result.json','prefix-admission.json']:
    dst=out/('failed-'+name)
    if dst.exists():assert sha(dst)==sha(target/name),name
    else:shutil.copyfile(target/name,dst)
bindings={};plan=read(coarse/'plan.json')
for ext in ['.npz','.json']:
    src=coarse/f'sweep-1/photons/interval-15-64{ext}';dst=target/src.relative_to(coarse)
    upstream=[p for p in plan['bindings'] if p.endswith('/interval-15-64'+ext)]
    assert upstream and all(sha(src)==plan['bindings'][p] for p in upstream)
    assert not dst.exists();os.link(src,dst);bindings[str(dst)]=sha(src)
# These are audit-owned names, not accepted physics. Preserve their failed
# versions above before recreating the unchanged verifier's output contract.
for name in ['audit-receipt.json','prefix-result.json','prefix-admission.json']:
    (target/name).rename(out/('original-'+name))
write(out/'plan.json',dict(classification='Conjectural',claim='Run the original full paired audit using its omitted but already saved coarse15/16prefix.',
    new_physical_steps=0,new_photon_solves=0,source_unchanged=sha(producer),physical_bindings=physical,
    linked_missing_audit_inputs=bindings,max_seconds=600,stop='Any original verification gate fails; never rerun physical paths or relax gates.'))
with (out/'audit.stdout.log').open('w') as stdout,(out/'audit.stderr.log').open('w') as stderr:
    code=subprocess.call([sys.executable,str(producer),'audit'],stdout=stdout,stderr=stderr)
for p,h in physical.items():assert sha(p)==h,('Physical history changed',p)
result=dict(classification='Counterexample candidate',returncode=code,seconds=time.monotonic()-start,
    original_producer_unchanged=sha(producer)==read(out/'plan.json')['source_unchanged'],
    completed_physical_histories_unchanged=True,new_physical_steps=0,new_photon_solves=0,
    failed_original_audit_preserved=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
write(out/'result.json',result);assert code==0,read(target/'audit-receipt.json')
accepted=read(target/'result.json');assert accepted['passed'] and accepted['same_corrected_native_arithmetic'] and accepted['common_saved_prefix_material_ledger_revalidated']
result.update(passed=True,time_relative=accepted['time_relative'],full_horizon_seconds=accepted['full_horizon_seconds']);write(out/'result.json',result)
write(target/'controller-status.json',dict(state='completed',completed=[dict(action='fine',returncode=0),dict(action='audit',returncode=0)],
    elapsed_seconds=time.time()-read(target/'controller-start.json')['started_unix'],
    repaired_after_terminal_audit_input_failure=True,repair_result=str(out/'result.json'),repair_producer_sha256=sha(__file__),
    original_failed_controller=str(out/'failed-controller-status.json'),final_charge_conclusion='unadjudicated',full_goal_complete=False))
rearmed=[]
for folder in [Path('native-full-captured244-work'),Path('native-retarded-extension248-work'),Path('native-full-return249-work')]:
    state,entry=stopped(folder);assert not (folder/'plan.json').exists() and not (folder/'full/plan.json').exists()
    archive=out/folder.name;archive.mkdir()
    for name in ['controller-start.json','controller-status.json']:(folder/name).rename(archive/name)
    rearmed.append(dict(folder=str(folder),old_entry=entry,failed_state=state))
write(out/'rearmed-waiters.json',dict(classification='Counterexample candidate',only_terminal_waiters_with_no_dispatched_actions=True,rows=rearmed))
print(json.dumps(result))
