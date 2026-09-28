from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-conservative-redshift167-work';out=root/'outputs/direct-eos-gr33/native-conservative-redshift'
manifest=out.parent/'native-conservative-redshift-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json')
    assert not final['passed'] and final['actual_source_applied_to_coupled_prefix'] and not final['full_goal_complete']
    prefixes=read(root/'.phase167-doc-prefixes.json');prior=read(out.parent/'native-mass-energy-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/apply_native_conservative_redshift.py';assert sha(module)==sha(runtime/'verification'/module.name)
    producer_hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==11 and all(v['source_sha256'] in producer_hashes for v in receipts.values())
    assert {k for k,v in receipts.items() if v['error']}=={'photon_pilot-1-receipt.json','photon_retry-1-receipt.json','photon_exact-1-receipt.json','finish-1-receipt.json'}
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),original_total_seconds=3600,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies including failed actions; imports, read-only exploration, derivation, writing, publication and Git excluded. Three4/8-step coupled prefixes only; no full-horizon or free-material production.')
    assert cost['recorded_action_seconds']<3600
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.rglob('*'):
        if not p.is_file():continue
        q=dest/p.relative_to(work);q.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,q);assert sha(q)==sha(p);copies[p.relative_to(work).as_posix()]=sha(q)
    shutil.copyfile(module,dest/'final-producer.py');copies['final-producer.py']=sha(module)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bound_to_preserved_versions=True,scientific_time_gate_passed=False))
    before=sha(master);data=read(master)
    tails={
      'model-definition':'분류: Proven. 같은 선형 광자·물질 방정식에서 x=y+H는 yprime=A*y+S+A*H-Hprime를 준다. 실제 적색편이 원천의 적분을 제거해도 연속 방정식은 동일하며, 유한 단계 결과의 동일성은 주장하지 않는다.',
      'observable-targets':'분류: Counterexample candidate. 면 적색편이 수정 원천을 동시 광자·열·수소 응답에 실제 적용했으나 원2% 시간 게이트가 실패했다. 새 GR 전하는 아직 없으며 기존 선택 전하 부호를 전체 수정 결과로 확대하지 않는다.',
      'adiabatic-limit':'분류: Conjectural. 원천 원시함수의 해석적 계산은 정적 비교 모형의 흡수 가능성이나 궤도 구동을 판정하지 않는다. 실제 동시 결합 전달의 시간 정확도를 먼저 해결해야 한다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 같은0.214652ms의4/8단계 예비 계산에서 실제 원천 적분을 사용한 광자 에너지·압력은2% 시간 대조를 통과했지만 물질 에너지·수소·충격량·압력은4.4–10.2%로 실패했다. 보존 장부 통과와 시간 정확도를 구분한다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 원천 직접 적용·선두항 제거·실제 원시함수 제거의 세 예비 결합 계산 모두 원2% 시간 대조에서 실패했다. 각 원 파일과 생산자를 보존하고 전 기간 실행은 하지 않았다. 다음 병목은 광자→물질 동시 전달의 시간 적분이다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. actual_source_applied_to_coupled_prefix와 simultaneous_photon_thermal_H_prefix_completed는true다. time_accuracy_accepted,free_material_production_completed,full_horizon_coupled_response_completed,actual_mass_flux_balance_verified,exact_ADM_conservation_verified,physical_final_charge_solved,full_goal_complete는false이며 GR_charge_correction은null이다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계167 — 실제 적색편이 수정 결합 예비 계산\n\n'+line+' [단계167 보고서](../notes/REQUEST167_CONSERVATIVE_REDSHIFT_COUPLED_PREFIX_KO.md).\n').encode())
    note=root/'notes/REQUEST167_CONSERVATIVE_REDSHIFT_COUPLED_PREFIX_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 집계:11개 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. 원3600초 한도 안이다. import·읽기 전용 탐사·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase166_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(final,previous_master_commit='027adc70c4575def0a7a326f21e1c7623f2c66cf',previous_master_sha256=before,
        preserved_phase166_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_conservative_redshift_coupled_prefix']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,h in read(out/'preservation.json')['previous_phase166_nondoc'].items():assert sha(root/p)==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-conservative-redshift-167-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,time_gate_passed=False,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
