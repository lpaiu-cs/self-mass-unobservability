from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-direct-radau171-work';out=root/'outputs/direct-eos-gr33/native-direct-radau'
manifest=out.parent/'native-direct-radau-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert not final['passed'] and not final['full_goal_complete']
    final['actual_coupled_Radau_prefix_completed']=True;final['time_accuracy_accepted']=False
    final['error_localization']=read(work/'error-profile.json')
    final['failure_reason']='Original2percent material pressure time gate failed; other five time gates and physical stage conservation passed.'
    final['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    prefixes=read(root/'.phase171-doc-prefixes.json');prior=read(out.parent/'native-known-forcing-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/apply_native_direct_radau.py';assert sha(module)==sha(runtime/'verification'/module.name)
    hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==5 and all(v['source_sha256'] in hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'pilot-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v.get('CPU_seconds',0.) for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v.get('peak_RSS_bytes',0) for v in receipts.values()),original_total_seconds=300,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including failed time comparison; CPU/RSS excludes the unmeasured read-only profile; imports, derivation, writing, publication and Git excluded. One4/8-step pair, no full-horizon or free-material production.')
    assert cost['recorded_action_seconds']<300
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
      'model-definition':'분류: Counterexample candidate. H=0인 기존 경로에서 실제 단계 원천을 그대로 넣은 두 단계 Radau 결합 해를 계산했다. 원 면 적색편이 원천과 동일하고 0·부호·배율 및 원 물리 장부를 확인했다. 별도 유효 원천률을 사용하지 않는다.',
      'observable-targets':'분류: Counterexample candidate. 실제 직접 원천 Radau의 다섯 시간 성분은 통과했으나 물질 압력11.4298%가 실패했다. 수정된 동일 결합 해에서 최종 전하 결론은 판정 불가다.',
      'adiabatic-limit':'분류: Imported from prior work. 단계170의 스칼라 stiff 극한은 원시함수 제거의 유한 단계 한계를 보였다. 원천 직접 적용은 그 특정 한계를 피하지만 실제 시간 수렴이나 정적 관측 식별성을 대신하지 않는다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 실제4/8단계와 각도 출구·보존 장부를 완료했다. 압력 차이97.93%가 셀15에 집중된다. 같은 셀의 압력 모드와 구동 도달 구간의 오차를 해결하기 전 전체 진화·GR 반환을 수락하지 않는다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 원천 직접 적용만으로 압력 시간 병목을 해소하지 못했다. 서로 다른 미세 해의 작은 차이나 성분별 좋은 결과를 골라 수락하지 않는다. 압력 잔차 생성의 실제 결합 동역학을 먼저 특정하고 원 기준을 유지한다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. 지배 오차를 고친 동일 결합 해의 최종 전하 유지 여부는 아직 unadjudicated다. time_accuracy_accepted,full_horizon_completed,physical_final_charge_solved,full_goal_complete는false다. 보존·다섯 성분 통과를 최종 성과로 세지 않는다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계171 — 실제 단계 원천 직접 적용\n\n'+line+' [단계171 보고서](../notes/REQUEST171_DIRECT_SOURCE_RADAU_KO.md).\n').encode())
    note=root/'notes/REQUEST171_DIRECT_SOURCE_RADAU_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:5개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. import·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase170_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='d40b5c21fa57d904f0ada4194862f1842fc8f1fa',previous_master_sha256=before,
        preserved_phase170_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_direct_source_Radau_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase170_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-direct-radau-171-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=False,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
