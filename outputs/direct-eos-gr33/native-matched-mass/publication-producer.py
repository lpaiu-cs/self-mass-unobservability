from pathlib import Path
import hashlib,json,shutil,subprocess,sys

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-matched-mass163-work';out=root/'outputs/direct-eos-gr33/native-matched-mass'
manifest=out.parent/'native-matched-mass-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())

def package():
    assert not manifest.exists();final=read(work/'final-result.json');assert final['passed'] and final['audit_passed']
    assert not final['exact_ADM_conservation_verified'] and not final['full_goal_complete']
    prefixes=read(root/'.phase163-doc-prefixes.json');previous=read(out.parent/'native-mixed-return-manifest.json')
    preserved={k:h for k,h in previous['sha256'].items() if not k.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    module=root/'verification/read_native_matched_mass.py';assert sha(module)==sha(runtime/'verification/read_native_matched_mass.py')
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')};versions={sha(module):'verification/read_native_matched_mass.py',sha(work/'registered-producer.py'):'registered-producer.py'}
    for v in receipts.values():assert v['source_sha256'] in versions
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),registered_total_seconds=180,CPU_threads=1,virtual_GiB=3,
        scope='Recorded action bodies only; excludes imports, preflight, inspection, derivation, writing, publication and Git.')
    assert cost['recorded_action_seconds']<180
    dest=out/'completed';dest.mkdir();copies={}
    for p in work.iterdir():
        assert p.is_file();q=dest/p.name;shutil.copyfile(p,q);assert sha(p)==sha(q);copies[p.name]=sha(q)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final['cost']=cost;write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,receipt_source_versions={k:versions[v['source_sha256']] for k,v in receipts.items()}))
    before=sha(master);data=read(master)
    tails={
        'model-definition':'분류: Counterexample candidate. lapse에 사용하던 기존 보존 원천의 질량 잔여와 새 mixed homogeneous 성분을 같은 배경의 전하 분모에 적용했다. 보조 J를 새 방출 에너지로 세지 않는다. 추가 외부 혼합 원천과 정확한 ADM/flux 보존은 미완료다.',
        'observable-targets':'분류: Counterexample candidate. 질량 정규화 보정+3.4829714644393073e-57을 실제 선택 전하에 반영해 끝점-2.5375890992998707e-51을 얻었다. 변화는 기존 크기의 약1.37255ppm이며 음의 부호가 유지된다. 전체 물리 전하의 오차 인증이나 정적/관측 비흡수성은 아니다.',
        'adiabatic-limit':'분류: Conjectural. 같은 보존 원천의 질량항을 분리해 공통 정규화하더라도, 짧은 펄스의 전하 변화가 같은 재고의 정적 비교를 벗어난다는 결론은 따라오지 않는다. 단열·궤도·정적 비교 경계는 유지한다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 기존64/128 incident 경로의 질량 보정 상대차0.00688011,4/8 구적1.00867e-6,mixed 원천17/33 표현0.0100427을 원 기준으로 수락했다. 마지막은 새 진화 수렴이 아닌 원천 표현 대조다. 신규 진화 단계는0이다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 큰 직접장이 지배하는 전체 합의 고정밀 대조만으로 작은 질량·분모 보정을 인증할 수 없었다. 저장값의 개별 보정을110자리로 계산해 독립 열로 보존하고 각 유리식 차분을 대조했다. 원 판독·프로그램은 보존하며 남은 비영 질량 잔여를 열이나 광자 에너지로 넘기지 않았다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. stored_mass_ports_applied,common_background_denominator,compensated_component_checks는true다. homogeneous_component_readout만 수락했으며 additional_exterior_mixed_stress_closed,exact_ADM_conservation_verified,physical_final_charge_solved,uniform_error_enclosure,complete_static_comparison,observable_identified,full_goal_complete는false다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';key=p.relative_to(root).as_posix();raw=p.read_bytes();v=prefixes[key]
        assert len(raw)==v['bytes'] and sha(p)==v['sha256'],key
        with p.open('ab') as f:f.write(('\n\n## 단계163 — 질량항의 공통 전하 정규화\n\n'+line+' [단계163 보고서](../notes/REQUEST163_MATCHED_MASS_READOUT_KO.md).\n').encode())
    note=root/'notes/REQUEST163_MATCHED_MASS_READOUT_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 기록: 네 action 본문 합계{cost['recorded_action_seconds']:.3f}초, CPU{cost['CPU_seconds']:.3f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트다. 원180초 합산 한도·1thread·3GiB 가상 메모리 안에서 저장값만 재사용했다. import·사전 검사·열람·유도·작성·발행/Git 시간은 이 합계에 포함하지 않는다.\n").encode())
    write(out/'preservation.json',dict(previous_phase162_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary={k:final[k] for k in ['classification','passed','matched_selected_endpoint','matched_mass_charge_endpoint','mass_effect_over_selected','actual_stored_mass_ports_applied','same_background_denominator','homogeneous_component_readout','additional_exterior_mixed_stress_closed','exact_ADM_conservation_verified','physical_final_charge_solved','full_goal_complete']}
    summary.update(accepted_scope=final['scope'],cost=cost,previous_master_commit='b1e26c922e93f7068fa704ff49f57fe22efaa614',previous_master_sha256=before,
        next_decisive_lever='Propagate the exterior packet measure/energy/path variation together with scalar mixed constraints, then compare with initial ADM and outgoing flux. Preserve full physical/static/observational obligations.',
        preserved_phase162_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_matched_mass']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)

def check(mode):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],name
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-matched-mass-163-paths').write_bytes(b'\0'.join(s.encode() for s in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,subprocess.check_output(['git','diff','--cached','--name-only','-z'],cwd=root).decode().split('\0')))==set(paths)
        ref=':' if mode=='staged' else 'HEAD:'
        for name in paths:assert hashlib.sha256(subprocess.check_output(['git','show',ref+name],cwd=root)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,full_goal_complete=False)))

if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
