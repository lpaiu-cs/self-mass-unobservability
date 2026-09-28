from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-radau-transfer168-work';out=root/'outputs/direct-eos-gr33/native-radau-transfer'
manifest=out.parent/'native-radau-transfer-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert not final['passed'] and not final['full_goal_complete']
    final['actual_coupled_Radau_prefix_completed']=True;final['time_accuracy_accepted']=False
    final['error_localization']=read(work/'error-profile.json')
    prefixes=read(root/'.phase168-doc-prefixes.json');prior=read(out.parent/'native-conservative-redshift-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/solve_native_radau_transfer.py';assert sha(module)==sha(runtime/'verification'/module.name)
    hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==4 and all(v['source_sha256'] in hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'pilot-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),original_total_seconds=300,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including failed JSON export; imports, derivation, writing, publication and Git excluded. One4/8-step pair, no full-horizon or free-material production.')
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
      'model-definition':'분류: Proven. 두 collocation 절점1/3·1의 Radau 행렬[[5/12,-1/12],[3/4,1/4]]과 가중치[3/4,1/4]의 고전적3차 조건을 확인했다. 실제 stiff 광자·물질 방정식의 균일 시간 오차 보증은 아니다.',
      'observable-targets':'분류: Counterexample candidate. 같은 수정 원천의 실제 두 단계 동시 Radau 계산에서 물질 압력 대조는0.885%로 통과했지만 물질 에너지·수소·충격량은3.8–4.3%로 실패했다. 새 GR 전하를 수락하지 않는다.',
      'adiabatic-limit':'분류: Conjectural. 시간 적분기 변경은 정적 재정의로 흡수되지 않는 신호의 증거가 아니다. 같은 입력의 실제 응답 정확도와 원 질량·정적 비교 경계를 유지한다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 같은64/128 시계의4/8단계 Radau 예비 구간을 완료했다. 보존·단계 잔차·실제 가중 각도 출구는 통과했지만 원2% 시간 기준은 실패했다. 저장 packet에 새 가중치를 명시했으며 이전 SDIRK 판독기로 GR 반환하지 않았다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 고전적3차 Radau만으로 물질 시간 정확도를 해결하지 못했다. 세 실패 관측량 차이의94.3–94.9%가 셀15에 집중된다. 구동 도달·원천 적분을 먼저 특정해야 하며 전체 시계나 방법 차수를 자동 확대하지 않는다. 완료 뒤 JSON 변환 실패도 보존하고 단계 재실행 없이 복구했다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. actual_coupled_Radau_prefix_completed=true다. time_accuracy_accepted,full_horizon_completed,physical_final_charge_solved,full_goal_complete는false다. 실제 면 적색편이 수정의 최종 GR·질량 폐쇄 효과는 아직 판정하지 못했다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계168 — 동시 Radau 적분의 실제 예비 대조\n\n'+line+' [단계168 보고서](../notes/REQUEST168_RADAU_TRANSFER_PREFIX_KO.md).\n').encode())
    note=root/'notes/REQUEST168_RADAU_TRANSFER_PREFIX_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:4개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. import·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase167_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='a06c4fa043faaeff42258118b3547b911bf86fc2',previous_master_sha256=before,
        preserved_phase167_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_Radau_coupled_time_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase167_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-radau-transfer-168-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=False,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
