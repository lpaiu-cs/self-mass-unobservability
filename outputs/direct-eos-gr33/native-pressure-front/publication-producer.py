from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-pressure-front172-work';out=root/'outputs/direct-eos-gr33/native-pressure-front'
manifest=out.parent/'native-pressure-front-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert final['passed'] and not final['full_goal_complete']
    final['actual_coupled_Radau_prefix_completed']=True;final['time_accuracy_accepted']=True;final['acceptance_scope']='0.214651945ms prefix only; no full coupled solution or charge'
    final['pressure_mode_inspection']=read(work/'inspection.json')
    final['pressure_prefix_bottleneck_resolved']=True
    final['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    prefixes=read(root/'.phase172-doc-prefixes.json');prior=read(out.parent/'native-direct-radau-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/resolve_native_pressure_front.py';assert sha(module)==sha(runtime/'verification'/module.name)
    hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==6 and all(v['source_sha256'] in hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'check-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v.get('CPU_seconds',0.) for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v.get('peak_RSS_bytes',0) for v in receipts.values()),original_total_seconds=420,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including preserved grid-sum check failure; imports, derivation, writing, publication and Git excluded. One6/12-step pair, no full-horizon or free-material production.')
    assert cost['recorded_action_seconds']<420
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.rglob('*'):
        if not p.is_file():continue
        q=dest/p.relative_to(work);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p);copies[p.relative_to(work).as_posix()]=sha(q)
    shutil.copyfile(module,dest/'final-producer.py');copies['final-producer.py']=sha(module)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bound_to_preserved_versions=True,scientific_time_gate_passed=True))
    before=sha(master);data=read(master)
    tails={
      'model-definition':'분류: Counterexample candidate. 같은 직접 원천·두 단계 Radau에서 느린 수송 핵부의 입력 도달 구간만 한 번 이분했다. 실제 하위 단계의 충돌·물질·일·출구와 명시적 각도 가중치를 함께 유지했다.',
      'observable-targets':'분류: Counterexample candidate. 같은 실제 결합 해의 압력 시간 차이가11.4298%에서1.6284%로 줄어 여섯 성분 모두 원2% 기준을 통과했다. 수락 범위는 예비 구간이며 최종 전하 결론은 미판정이다.',
      'adiabatic-limit':'분류: Proven. 영 초기값의 선두 모형 xprime=t^3,Pprime=x는 P=t^5/20을 준다. 이 국소 성장식이나 예비 수렴을 물리적 정적 no-go 또는 전체 오차 상계로 확대하지 않는다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 원0.214651945ms 구간의 실제6/12단계 대조와 보존·출구가 통과했다. 기본64/128을 유지하되 실제 단계를 숨기지 않으며, 같은 규칙의 전체117/227단계는 아직 실행하지 않았다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 압력 사상 변화나 극단적인 국소 물질 감쇠만으로 실패를 설명하지 못했다. 실제 구동 도달 이후의 시간 해상으로 원 압력 예비 기준을 통과했으며, 원 실패·격자 합산 실패를 보존했다. 전체 해의 수락 여부는 별도다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. prefix_time_accuracy_accepted=true이고 full_horizon_completed,physical_final_charge_solved,full_goal_complete는false다. 다음 판정은 이 수정의 동일 결합 해를 끝까지 연결해 최종 전하 결론을 다시 읽는 것이다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계172 — 같은 결합 해의 압력 도달 구간 해상\n\n'+line+' [단계172 보고서](../notes/REQUEST172_PRESSURE_ONSET_RESOLUTION_KO.md).\n').encode())
    note=root/'notes/REQUEST172_PRESSURE_ONSET_RESOLUTION_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:6개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. import·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase171_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='4f469fdfe71e6eeb2fee53dd02cd716664284963',previous_master_sha256=before,
        preserved_phase171_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_pressure_onset_resolved_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase171_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-pressure-front-172-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=True,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
