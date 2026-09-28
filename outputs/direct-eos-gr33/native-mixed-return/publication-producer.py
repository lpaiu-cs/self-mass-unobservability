from pathlib import Path
import hashlib,json,py_compile,shutil,subprocess,sys

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-mixed-return162-work'
out=root/'outputs/direct-eos-gr33/native-mixed-return'
manifest=out.parent/'native-mixed-return-manifest.json'
master=root/'paper/revision-manifest.json'
modules=['return_native_mixed_gr','resolve_native_mixed_return','repair_native_mixed_material','audit_native_mixed_return']
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())

def package(sweep):
    assert not manifest.exists()
    final=read(work/'result.json');audit=read(work/'audit.json');block=read(work/f'sweep-{sweep}/block-result.json')
    assert final['passed'] and audit['passed'] and block['passed']
    assert final['actual_mixed_field_applied_to_photons_and_free_material'] and not final['full_goal_complete']
    compact=read(work/f'sweep-{sweep}/gr/result.json')
    photon=read(work/f'sweep-{sweep}/photons/result.json');material=read(work/f'sweep-{sweep}/material/production.json')
    sources=read(work/f'sweep-{sweep}/material/sources.json')
    assert all(r['passed'] for r in [compact,photon,material,sources])
    data=read(master);before=sha(master);prefixes=read(root/'.phase162-doc-prefixes.json')
    previous=read(out.parent/'native-mixed-gr-manifest.json');preserved={}
    for name,h in previous['sha256'].items():
        if not name.startswith('docs/'):
            assert sha(root/name)==h,name;preserved[name]=h
    for name,v in prefixes.items():
        assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],name
    for module in modules:
        p=root/f'verification/{module}.py';assert sha(p)==sha(runtime/f'verification/{module}.py')
        py_compile.compile(str(p),doraise=True)
    versions={sha(p):p.name for p in work.glob('*.py')}
    versions.update({sha(root/f'verification/{m}.py'):f'verification/{m}.py' for m in modules})
    receipts={p.relative_to(work).as_posix():read(p) for p in work.rglob('*-receipt.json')}
    for name,v in receipts.items():assert v['source_sha256'] in versions,(name,v['source_sha256'])
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),
        CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),
        CPU_threads=1,virtual_GiB=3,registered_total_seconds=4200,
        scope='Recorded action bodies including failures and read-only budget reassessment. Excludes imports, preflight, inspection, writing, publication and Git. Not total task walltime.')
    assert cost['recorded_action_seconds']<4200
    dest=out/'completed';dest.mkdir(parents=True);copies={};omitted={};size=0
    for p in work.rglob('*'):
        if not p.is_file():continue
        assert not p.is_symlink();rel=p.relative_to(work)
        if p.name.endswith('-checkpoint.npz'):
            assert (p.parent/p.name.replace('-checkpoint.npz','.npz')).exists()
            omitted[rel.as_posix()]=sha(p);continue
        q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(p)==sha(q);copies[rel.as_posix()]=sha(q);size+=q.stat().st_size
    shutil.copyfile(__file__,out/'publication-producer.py')
    write(out/'publication.json',dict(copies=copies,bytes=size,omitted_terminal_checkpoints=omitted,
        omission_scope='Only duplicate terminal restart checkpoints omitted from Git; hashes retained and runtime files left intact.',
        receipt_source_versions={k:versions[v['source_sha256']] for k,v in receipts.items()}))
    write(out/'preservation.json',dict(previous_phase161_nondoc=preserved,document_prefixes=prefixes,
        before_task_checkpoint='b27561096b6dda642fcb2439f3a003bd509e0375',previous_master_sha256=before))
    result=dict(classification='Counterexample candidate',passed=True,actual_mixed_field_returned=True,
        finite_block_coupling_passed=True,corrected_mixed_charge=audit['corrected_mixed_endpoint'],
        transport_return_charge=audit['physical_transport_return_endpoint'],
        transport_over_corrected_mixed=audit['transport_over_corrected_mixed'],
        selected_endpoint_without_tiny_return=audit['selected_endpoint_without_tiny_return'],
        compact_raw=compact,infinity=final,audit=audit,block=block,photon=photon,material=material,sources=sources,cost=cost,
        inherited_compact_metadata=audit['inherited_readout_metadata'],
        exterior_mixed_source_closed=False,physical_mixed_ADM_closed=False,physical_final_charge_solved=False,
        full_goal_complete=False,contribution='Loophole progress: the computed compact mixed field is actually returned through photons and free matter, its finite coupled inputs pass, and actual signed emission is read at fixed-operator infinity. No universal or observational conclusion.')
    write(out/'final-result.json',result)
    files=[root/f'verification/{m}.py' for m in modules]+[root/'notes/REQUEST162_MIXED_GR_TRANSPORT_KO.md']
    files += [root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(classification='Counterexample candidate',passed=True,actual_mixed_field_returned=True,
        finite_block_coupling_passed=True,actual_signed_emission_read_at_fixed_operator_infinity=True,
        corrected_mixed_charge=audit['corrected_mixed_endpoint'],transport_return_charge=audit['physical_transport_return_endpoint'],
        transport_over_corrected_mixed=audit['transport_over_corrected_mixed'],
        exterior_mixed_source_closed=False,physical_mixed_ADM_closed=False,physical_final_charge_solved=False,full_goal_complete=False,
        accepted_scope=final.get('scope',audit['scope']),next_decisive_lever='Close the omitted exterior mixed GR source and physical mixed mass/normalization, then apply them to the same final charge. Preserve physical-error and static/observational obligations.',
        cost=cost,previous_master_commit='b27561096b6dda642fcb2439f3a003bd509e0375',previous_master_sha256=before,
        preserved_phase161_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,
        raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_mixed_return']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],name
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-mixed-return-162-paths').write_bytes(b'\0'.join(x.encode() for x in paths)+b'\0')
    if mode in ['staged','head']:
        ref=':' if mode=='staged' else 'HEAD:'
        if mode=='staged':
            actual=subprocess.check_output(['git','diff','--cached','--name-only','-z'],cwd=root).decode().split('\0')
            assert set(filter(None,actual))==set(paths)
        for name in paths:assert hashlib.sha256(subprocess.check_output(['git','show',ref+name],cwd=root)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),prefixes=6,full_goal_complete=False,cost=m['cost'])))

if __name__=='__main__':
    if sys.argv[1]=='package':package(int(sys.argv[2]))
    check(sys.argv[1])
