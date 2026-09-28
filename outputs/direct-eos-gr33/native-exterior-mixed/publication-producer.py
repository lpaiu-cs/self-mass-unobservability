from pathlib import Path
import hashlib,json,shutil,subprocess,sys

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-exterior-mixed164-work';out=root/'outputs/direct-eos-gr33/native-exterior-mixed'
manifest=out.parent/'native-exterior-mixed-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)


def package():
    assert not manifest.exists();final=read(work/'result.json');audit=read(work/'audit.json')
    assert final['passed'] and audit['passed'] and not final['full_goal_complete']
    prefixes=read(root/'.phase164-doc-prefixes.json');previous=read(out.parent/'native-matched-mass-manifest.json')
    preserved={k:h for k,h in previous['sha256'].items() if not k.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    module=root/'verification/bound_native_exterior_mixed.py';assert sha(module)==sha(runtime/'verification/bound_native_exterior_mixed.py')
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==3 and all(v['source_sha256']==sha(module) and v['error'] is None for v in receipts.values())
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),registered_total_seconds=135,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies only. Imports, derivation, inspection, writing, publication and Git are not included. No new EOS roots, fluid steps or rays.')
    assert cost['recorded_action_seconds']<135
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.iterdir():
        assert p.is_file();q=dest/p.name;shutil.copyfile(p,q);assert sha(p)==sha(q);copies[p.name]=sha(q)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost,audit_passed=True);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bind_current_producer=True))
    before=sha(master);data=read(master)
    tails={
        'model-definition':'분류: Proven. 직접 방출 photon의 혼합 제약에는 밀도 측도·에너지·이동 delta와 초기 scalar 기울기의 계량 곱이 함께 필요하다. 광선 시간 Jacobian과 기울기 껍질을 더하면 해당 zeta 진폭 항은 상쇄된다. 분류: Counterexample candidate. 원 양의 방출과 primary+유한 Born 입력의 이 부문을 상계에 연결했다. 생성 배경 scalar와 이후 signed 반환은 별도다.',
        'observable-targets':'분류: Counterexample candidate. 새 외부 광자 혼합 부문의 전하 상계1.36903e-57을 기존 외부 상계와 합하면 선택 끝점-2.53758910e-51의0.7139523% 이하다. 조건부 구간[-2.55570627e-51,-2.51947192e-51]은 음수다. 전체 물리 오차 구간이나 정적/관측 비흡수성 판정은 아니다.',
        'adiabatic-limit':'분류: Conjectural. 지정한 짧은 펄스에서 외부 광자 기여를 제한해도 단열·궤도 규모의 정적 비교가 성립하거나 깨졌다는 결론은 따르지 않는다. 기존 동일 재고 비교와 미분 nuisance 경계는 유지한다.',
        'nonadiabatic-regime':'분류: Proven. 추가 외부 질량 제약의 바깥 꼬리는 반경별 시간 길이T+x/c로 제한해야 한다. 분류: Counterexample candidate. 이 가중치를 포함한 정적 potential 수축 상계2.14309e-6을 사용해 지정한 새 원천의 반복을 제한했다. 원 양의64/128 방출과 저장 입력을 재사용했으며 새 광선·유체 경로는 없다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 단계159의 원 점 추정 실패는 보존한다. 새 bound는 직접 방출 광자의 혼합 제약·기울기·측도·이동 경로에 대한 별도 절대 기여 결론이다. 추가 질량을0으로 맞추지 않았고, 전체 배경 scalar·signed 반환·ADM/flux·물리 원천 오차는 수락하지 않았다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. direct_photon_mixed_sector_enclosed와selected_sign_survives는true다. additional_exterior_mixed_stress_closed,exact_ADM_conservation_verified,physical_final_charge_solved,uniform_error_enclosure,complete_static_comparison,observable_identified,full_goal_complete는false다. 다음에는 생성된 배경 scalar의 외부 혼합 응력·계량 작용과 실제 질량 포트를 연결한다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';key=p.relative_to(root).as_posix();v=prefixes[key]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],key
        with p.open('ab') as f:f.write(('\n\n## 단계164 — 외부 광자 혼합 원천의 전하 상계\n\n'+line+' [단계164 보고서](../notes/REQUEST164_EXTERIOR_PHOTON_MIXED_BOUND_KO.md).\n').encode())
    note=root/'notes/REQUEST164_EXTERIOR_PHOTON_MIXED_BOUND_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 기록: 세 action 본문 합계{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트다. 원135초 한도·1thread·3GiB 안에서 완료했다. import·유도·열람·작성·발행/Git 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase163_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary={k:final[k] for k in ['classification','passed','audit_passed','new_sector_charge_envelope','combined_selected_exterior_envelope','combined_envelope_over_selected','conditional_selected_interval','selected_sign_survives','direct_photon_mixed_sector_enclosed','additional_exterior_mixed_stress_closed','exact_ADM_conservation_verified','physical_final_charge_solved','full_goal_complete']}
    summary.update(accepted_scope=final['scope'],cost=cost,previous_master_commit='53a1718d7ec4b7c7bf4ff01152e83dab882d8ec9',previous_master_sha256=before,
        next_decisive_lever='Connect generated-background scalar mixed exterior stress and metric action, signed later returns and actual mass ports before an ADM/flux claim. Reuse saved source/operator data first.',
        preserved_phase163_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_exterior_mixed']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)


def check(mode):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],name
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-exterior-mixed-164-paths').write_bytes(b'\0'.join(s.encode() for s in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        ref=':' if mode=='staged' else 'HEAD:'
        for name in paths:assert hashlib.sha256(git('show',ref+name)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,original_completion_preserved=True,full_goal_complete=False)))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
