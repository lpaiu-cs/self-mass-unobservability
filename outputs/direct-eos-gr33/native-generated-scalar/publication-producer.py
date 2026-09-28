from pathlib import Path
import hashlib,json,shutil,subprocess,sys

root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-generated-scalar165-work';out=root/'outputs/direct-eos-gr33/native-generated-scalar'
manifest=out.parent/'native-generated-scalar-manifest.json';master=root/'paper/revision-manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args):return subprocess.check_output(['rtk','proxy','git',*args],cwd=root)


def package():
    assert not manifest.exists();final=read(work/'result.json');audit=read(work/'audit.json')
    assert final['passed'] and audit['passed'] and not final['full_goal_complete']
    prefixes=read(root/'.phase165-doc-prefixes.json');previous=read(out.parent/'native-exterior-mixed-manifest.json')
    preserved={k:h for k,h in previous['sha256'].items() if not k.startswith('docs/')}
    for name,h in preserved.items():assert sha(root/name)==h,name
    module=root/'verification/bound_native_generated_scalar.py';assert sha(module)==sha(runtime/'verification/bound_native_generated_scalar.py')
    receipts={p.name:read(p) for p in work.glob('*-receipt.json')}
    assert len(receipts)==4 and all(v['source_sha256']==sha(module) and v['error'] is None for v in receipts.values())
    cost=dict(recorded_action_seconds=sum(v['seconds'] for v in receipts.values()),CPU_seconds=sum(v['CPU_seconds'] for v in receipts.values()),
        maximum_recorded_RSS_bytes=max(v['peak_RSS_bytes'] for v in receipts.values()),registered_total_seconds=165,CPU_threads=1,virtual_GiB=3,
        scope='Action bodies only. Imports, derivation, inspection, writing, publication and Git are not included. No new EOS roots, fluid steps or rays.')
    assert cost['recorded_action_seconds']<165
    dest=out/'completed';dest.mkdir(parents=True);copies={}
    for p in work.iterdir():
        assert p.is_file();q=dest/p.name;shutil.copyfile(p,q);assert sha(p)==sha(q);copies[p.name]=sha(q)
    shutil.copyfile(__file__,out/'publication-producer.py')
    final.update(cost=cost,audit_passed=True);write(out/'final-result.json',final)
    write(out/'publication.json',dict(copies=copies,all_raw_files_copied=True,all_receipts_bind_current_producer=True))
    before=sha(master);data=read(master)
    tails={
        'model-definition':'분류: Proven. 외부 scalar 혼합 질량은 공간·시간 미분의 에너지 교차항과 반경 경계항으로 분해된다. 두 outgoing 성분은2/c 적분 U_Bt U_It의 signed scalar 방출 에너지를 만든다. 분류: Counterexample candidate. 추가 점근 질량0인 외부 특해와 이 에너지를 판독 상계에 연결했다. 일반 해의 homogeneous 함수와 내부 접합을0으로 설정하지 않았다.',
        'observable-targets':'분류: Counterexample candidate. 생성 scalar의 새 선택 상계2.28147e-57과 이전 외부 상계의 합은 선택 전하의0.7140422% 이하이며 조건부 구간[-2.55570856e-51,-2.51946964e-51]은 음수다. 새 내부 질량 접합 보정과 전체 물리/관측 오차는 이 구간 밖이다.',
        'adiabatic-limit':'분류: Conjectural. 이 짧은 펄스의 scalar 교차 에너지와 선택 전하 상계는 같은 재고 정적 비교나 단열·궤도 관측량을 대신하지 않는다. 전체 질량 보존과 구동·비교 모형의 연결이 필요하다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 현재531셀·33시각 compact 원천과 원 양의 방출 및 실제 세 signed SDIRK 반환을 생성 scalar 원천의 절댓값 상계에 포함했다. 입력은 기존 primary+유한 Born이다. 새로운 유체·광선 이력 없이 scalar 질량/flux를 외부 판독에 연결했다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 종료T에서 기존 paired 질량에 이번 광자·scalar 허용량을 모두 더해도 구간[-3.53941e-50,-1.53983e-50]cm는0을 제외한다. 지정 외부 항만으로 질량 잔여가 사라질 것이라는 기대를 배제한다. 나머지 같은 차수 원천과 초기/에너지 기준을 검사하기 전 전체 ADM 불가능성으로 확대하지 않는다. 원 점 추정 실패도 유지한다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. generated_scalar_sector_enclosed,outgoing_scalar_energy_channel_included,selected_sign_survives는true다. actual_mass_flux_balance_verified,additional_exterior_mixed_stress_closed,exact_ADM_conservation_verified,physical_final_charge_solved,uniform_error_enclosure,complete_static_comparison,observable_identified,full_goal_complete는false다. 다음은 질량 경계·기준 에너지 변환·scalar flux·원천 일의 동일 재고 대조다.'}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';key=p.relative_to(root).as_posix();v=prefixes[key]
        assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],key
        with p.open('ab') as f:f.write(('\n\n## 단계165 — 생성 scalar의 외부 혼합 질량과 에너지 흐름\n\n'+line+' [단계165 보고서](../notes/REQUEST165_GENERATED_SCALAR_EXTERIOR_KO.md).\n').encode())
    note=root/'notes/REQUEST165_GENERATED_SCALAR_EXTERIOR_KO.md'
    with note.open('ab') as f:f.write((f"\n운영 기록: 네 action 본문 합계{cost['recorded_action_seconds']:.6f}초, CPU{cost['CPU_seconds']:.6f}초, 최대 RSS{cost['maximum_recorded_RSS_bytes']}바이트다. 원165초 한도·1thread·3GiB 안에서 완료했다. import·유도·열람·작성·발행/Git 시간은 별도다.\n").encode())
    write(out/'preservation.json',dict(previous_phase164_nondoc=preserved,document_prefixes=prefixes,previous_master_sha256=before))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    summary={k:final[k] for k in ['classification','passed','audit_passed','new_sector_charge_envelope','combined_selected_exterior_envelope','combined_envelope_over_selected','conditional_selected_interval','selected_sign_survives','generated_scalar_sector_enclosed','outgoing_scalar_energy_channel_included','finite_T_selected_mass_port_interval_cm','finite_T_selected_mass_port_excludes_zero','additional_exterior_mixed_stress_closed','exact_ADM_conservation_verified','physical_final_charge_solved','full_goal_complete']}
    summary.update(accepted_scope=final['scope'],cost=cost,previous_master_commit='a57c0d5da7c295b8859660b4ccf7e1e8fdb5e9a1',previous_master_sha256=before,
        next_decisive_lever='Compare the same-interface mass constraint, reference-to-Killing photon energy conversion, scalar flux and material source work. Resolve remaining same-order inputs and actual mass matching before refining small exterior terms or running new evolution.',
        preserved_phase164_nondoc_files=len(preserved),preserved_document_prefixes=prefixes,raw_files=len(copies),sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_generated_scalar']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)


def check(mode):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],name
    completion='docs/dynamic-charge-completion.md'
    assert (root/completion).read_bytes().splitlines()[:20]==git('show','d5744f766719f0a29c7acb913c7984bdfd3624e2:'+completion).splitlines()[:20]
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-generated-scalar-165-paths').write_bytes(b'\0'.join(s.encode() for s in paths)+b'\0')
    if mode in ['staged','head']:
        if mode=='staged':assert set(filter(None,git('diff','--cached','--name-only','-z').decode().split('\0')))==set(paths)
        ref=':' if mode=='staged' else 'HEAD:'
        for name in paths:assert hashlib.sha256(git('show',ref+name)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),preserved_prefixes=6,original_completion_preserved=True,full_goal_complete=False)))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1])
