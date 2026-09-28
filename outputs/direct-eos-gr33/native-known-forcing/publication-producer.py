from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-known-forcing170-work';out=root/'outputs/direct-eos-gr33/native-known-forcing'
manifest=out.parent/'native-known-forcing-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert not final['passed'] and not final['full_goal_complete']
    final['actual_coupled_Radau_prefix_completed']=True;final['time_accuracy_accepted']=False
    final['stiff_limit']=read(work/'stiff-limit.json')
    final['failure_reason']='All six original2percent time gates failed; physical stage equations and conservation passed.'
    final['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    prefixes=read(root/'.phase170-doc-prefixes.json');prior=read(out.parent/'native-collision-source-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/integrate_native_known_forcing.py';assert sha(module)==sha(runtime/'verification'/module.name)
    hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==5 and all(v['source_sha256'] in hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'pilot-receipt.json','failed-symbolic-prepare-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),original_total_seconds=390,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including failed symbolic preparation; imports, derivation, writing, publication and Git excluded. One4/8-step pair, no full-horizon or free-material production.')
    assert cost['recorded_action_seconds']<390
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
      'model-definition':'분류: Proven. 스칼라 stiff 모형에서 알려진 LH의 정확한 누적 적분률을 단계에 넣으면 실제 해가0으로 가는 극한에도 H-H_eff가 남는다. 연속 변수 변환·보존 항등식과 유한 단계의 stiff 정확성을 구분한다.',
      'observable-targets':'분류: Counterexample candidate. 수송·충돌·출구를 함께 적분한 실제 대조는 여섯 시간 기준을 모두 실패했다. 수정된 동일 결합 해에서 최종 전하 결론은 아직 판정 불가다.',
      'adiabatic-limit':'분류: Proven. 스칼라 xprime=lambda*x+t^3에서 실제 원천을 직접 넣은 Radau는 lambda*x_j=-S_j 극한을 보존하지만 정확한 H 제거의 유한 단계는 일반적으로 그러지 않는다. 이 수치적 극한을 물리적 정적 no-go와 혼동하지 않는다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 실제 방출 표본과 알려진 원천의 적분 보정량을 따로 보존했다. 원 각도 출구·보존은 통과했으나 광자 에너지·압력 차이는8.71/10.75%로 커졌다. 전체 진화·GR 반환은 시작하지 않았다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 원천 적분 정확도만 높이는 경로는 채택하지 않는다. 같은 실제 방정식의 단계 평형을 지키는 원천 직접 적용을 검토하며, 시계·기간·차수를 자동 확대하지 않는다. 원 실패와 첫 기호 준비 실패를 보존했다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. 현재 최종 전하 결론은 unadjudicated다. 지배 오차를 고친 동일 결합 해의 최종 전하 유지 여부만 성과 기준으로 삼으며 원천·보존·예비 통과로 대신하지 않는다. full_goal_complete=false를 유지한다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계170 — 수송 원천 적분 실패와 강한 감쇠 극한\n\n'+line+' [단계170 보고서](../notes/REQUEST170_KNOWN_FORCING_STIFF_LIMIT_KO.md).\n').encode())
    note=root/'notes/REQUEST170_KNOWN_FORCING_STIFF_LIMIT_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:5개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. import·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase169_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note,root/'AGENTS.md']+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='873a09ad13d0dd9df748d8ab6ea871e119f21068',previous_master_sha256=before,
        preserved_phase169_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_known_forcing_stiff_limit_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase169_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-known-forcing-170-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=False,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
