from pathlib import Path
import hashlib,json,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-mass-energy166-work';out=root/'outputs/direct-eos-gr33/native-mass-energy'
manifest=out.parent/'native-mass-energy-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)

def package():
    assert not manifest.exists();final=read(work/'result.json');assert final['passed'] and not final['full_goal_complete']
    prefixes=read(root/'.phase166-doc-prefixes.json');prior=read(out.parent/'native-generated-scalar-manifest.json')
    preserved={k:h for k,h in prior['sha256'].items() if not k.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    module=root/'verification/reconcile_native_mass_energy.py';assert sha(module)==sha(runtime/'verification'/module.name)
    producer_hashes={sha(p) for p in work.glob('*producer.py')}|{sha(module)}
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==6 and all(v['source_sha256'] in producer_hashes for v in receipts.values())
    assert [k for k,v in receipts.items() if v['error']]==['ledger-receipt.json']
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),original_total_seconds=270,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies only; imports, bounded read-only exploratory commands, derivation, writing, publication and Git excluded. No new fluid steps, roots or rays.')
    assert cost['recorded_action_seconds']<270
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.iterdir():
        assert p.is_file();q=dest/p.name;shutil.copyfile(p,q);assert sha(q)==sha(p);copies[p.name]=sha(q)
    shutil.copyfile(module,dest/'final-producer.py');copies['final-producer.py']=sha(module)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bound_to_preserved_versions=True))
    before=sha(master);data=read(master)
    tails={
      'model-definition':'분류: Proven. 기준 광자 수의 공간 적색편이 일은 같은 upwind 면 유량과 lapse 섭동 차이의 곱으로 구성할 수 있다. 중심의 연속 미분을 따로 사용하면 이산 에너지 곱셈 법칙은 자동 성립하지 않는다.',
      'observable-targets':'분류: Counterexample candidate. 실제 입력에서 누락된 공간 일0.3853807705erg는 기존 paired 질량 잔여 절댓값의1.253845배다. 앞선 선택 전하 부호 구간은 이 수송 수정의 결합 응답을 포함하지 않는다. 실제 전하를 다시 읽기 전 갱신하지 않는다.',
      'adiabatic-limit':'분류: Conjectural. 유한 셀 적색편이 일의 수정은 정적 비교 모형이나 궤도 구동을 확정하지 않는다. 먼저 같은 수송의 에너지 장부와 GR 질량을 연결해야 한다.',
      'nonadiabatic-regime':'분류: Counterexample candidate. 기존64/128단계 입력만 재사용하여 중심 미분 일과 면 유량 일을 대조했다. 이 차이와 각도 평균 보정을 실제 광자·물질 결합 증분으로 반환하는 일이 다음 단계다.',
      'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 질량 잔여를 외부 상계만으로 설명하는 접근에 이어 공간 수송의 곱셈 법칙 불일치를 특정했다. 원1e-12 ledger 등가 검사 실패, 출구 net-cell 판독, 단독 공간 counterterm 변경은 각각 보존하며 수락 보정으로 대체하지 않는다.',
      'dynamic-charge-completion':'분류: Counterexample candidate. discrete_redshift_work_mismatch_identified=true다. correction_applied_to_coupled_transport,actual_mass_flux_balance_verified,exact_ADM_conservation_verified,physical_final_charge_solved,full_goal_complete는false다. 다음은 면 commutator 원천을 실제 결합 진화에 적용하고 GR 전하·동일 경계 flux를 다시 읽는 것이다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';k=p.relative_to(root).as_posix();v=prefixes[k]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
        with p.open('ab') as f:f.write(('\n\n## 단계166 — 질량 장부와 이산 적색편이 일\n\n'+line+' [단계166 보고서](../notes/REQUEST166_MASS_ENERGY_RECONCILIATION_KO.md).\n').encode())
    note=root/'notes/REQUEST166_MASS_ENERGY_RECONCILIATION_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 기록: 여섯 action 본문{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트. 원270초 한도 안이다. import·읽기 전용 탐사·유도·발행 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase165_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary={k:final[k] for k in ['classification','passed','root','source_work_correction_erg','source_work_over_old_residual','rows','correction_applied_to_coupled_transport','actual_mass_flux_balance_verified','exact_ADM_conservation_verified','physical_final_charge_solved','full_goal_complete','next']}
    summary.update(cost=cost,previous_master_commit='b9ae9ab14',previous_master_sha256=before,
        preserved_phase165_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_mass_energy_reconciliation']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for p,h in m['sha256'].items():assert sha(root/p)==h and data['sha256'][p]==h,p
    for p,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-mass-energy-166-paths').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        for p in paths:assert hashlib.sha256(git('show',(':' if mode=='staged' else 'HEAD:')+p)).hexdigest()==sha(root/p),p
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
