from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-collision-source169-work';out=root/'outputs/direct-eos-gr33/native-collision-source'
manifest=out.parent/'native-collision-source-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert not final['passed'] and not final['full_goal_complete']
    final['actual_coupled_Radau_prefix_completed']=True;final['time_accuracy_accepted']=False
    final['error_localization']=read(work/'error-profile.json');final['pressure_decomposition']=read(work/'pressure-profile.json')
    final['failure_reason']='Original2percent material pressure time gate failed; inherited export_repair field does not describe this failure.'
    final['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    prefixes=read(root/'.phase169-doc-prefixes.json');prior=read(out.parent/'native-radau-transfer-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/integrate_native_collision_source.py';assert sha(module)==sha(runtime/'verification'/module.name)
    hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==5 and all(v['source_sha256'] in hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'pilot-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),read_only_profile_seconds=final['pressure_decomposition']['seconds'],CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),original_total_seconds=450,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including failed JSON export; imports, derivation, writing, publication and Git excluded. One4/8-step pair, no full-horizon or free-material production.')
    assert cost['recorded_action_seconds']<450
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.rglob('*'):
        if not p.is_file():continue
        q=dest/p.relative_to(work);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p);copies[p.relative_to(work).as_posix()]=sha(q)
    shutil.copyfile(module,dest/'final-producer.py');copies['final-producer.py']=sha(module)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bound_to_preserved_versions=True,scientific_time_gate_passed=False))
    before=sha(master);data=read(master)
    tails={
      'model-definition':'분류: Proven. 기존 아핀 C(t)와 정확한 H의 시간 모멘트로 알려진 충돌 원천을 누적 적분하고 Radau 단계에 맞추는 항등식을 확인했다. 광자와 물질에 같은 원천을 사용한다. 전체 시간 오차 보증은 아니다.',
      'observable-targets':'분류: Counterexample candidate. 물질 에너지·수소·충격량의 시간 차이는 원2% 기준에 들어왔으나 압력 차이는6.2624%로 실패했다. 수정된 동일 결합 해의 최종 전하 결론은 판정 불가다.',
      'adiabatic-limit':'분류: Conjectural. 원천의 정확한 적분은 정적 비교를 벗어나는 신호의 증거가 아니다. 최종 전하는 지배 오차를 고친 동일 결합 해에서 판정해야 한다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 실제4/8단계 결합 해에 C(t)H 적분을 적용했다. 보존과 실제 각도 출구는 통과했으나 물질 압력 시간 기준은 실패해 전체 기간과 GR 반환을 실행하지 않았다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 압력은 에너지·수소의 큰 기여가 상쇄된 잔차이고 새 압력 차이94.42%가 셀15·16에 집중된다. 충돌 항만 적분한 수정은 전체 시간 병목을 해소하지 못했다. 원 실패를 보존하고 수송·출구까지 일관된 알려진 원천 적분을 검토한다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. 사용자 기준은 지배 오차를 고친 동일 결합 해에서 최종 전하 결론이 유지되는가다. actual_coupled_prefix_completed=true이나 time_accuracy_accepted,full_horizon_completed,physical_final_charge_solved,full_goal_complete는false다. 원천·보존·예비 통과는 중간 검증이다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계169 — 충돌 원천 적분의 실제 예비 대조\n\n'+line+' [단계169 보고서](../notes/REQUEST169_COLLISION_SOURCE_PREFIX_KO.md).\n').encode())
    note=root/'notes/REQUEST169_COLLISION_SOURCE_PREFIX_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:5개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. import·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase168_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note,root/'AGENTS.md']+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='39341e99f70723f4e66671546853d4676201db80',previous_master_sha256=before,
        preserved_phase168_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_collision_known_source_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase168_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-collision-source-169-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=False,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
