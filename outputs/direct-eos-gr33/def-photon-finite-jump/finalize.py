"""One-shot Phase62 evidence publication, preserving earlier verdicts."""
from pathlib import Path
import json
import hashlib

root=Path('.')
out=root/'outputs/direct-eos-gr33/def-photon-finite-jump'
result=json.loads((out/'result.json').read_text())
assert result['passed'] and not result['full_dynamic_charge_solved']
assert json.loads((root/'outputs/direct-eos-gr33/def-photon-spatial-coupling/result.json').read_text())['passed'] is False
texts={
'model-definition':'분류: Proven. 양의 각도 구적에서 |alpha_l|<=alpha_0인 유한 주파수 이동 충돌 블록은 양의 준정부호이며, l=0 반동 항은 같은 물질과 광자의 총 에너지·유한 광자 수를 보존한다. 분류: Counterexample candidate. 자유 전자 적분 커널로 기존 탄성 산란과 Kompaneets를 함께 교체해 실제 국소 시간식에 연결했다. 원 흡수·EOS·상태를 유지했으며 native 매질 산란의 복원으로 부르지 않는다.',
'observable-targets':'분류: Counterexample candidate. 동일 자유 전자 조건에서 유한 이동 커널과 확산식의 응답 차이는 초기 노름의 1.51638%, 유한 이동 커널과 탄성식만의 차이는 0.03342%다. 단계 61의 큰 확산식 스펙트럼 변화를 검증된 물리 신호로 채택하지 않는다. 분류: Conjectural. 현재 성분 입력은 실제 궤도 구동·동적 전하 관측량·정적 비교 및 nuisance 제거를 대체하지 않는다.',
'adiabatic-limit':'분류: Proven. 비국소 주파수 쌍에서도 상세평형 반동 외적은 LTE 온도와 Bose 화학퍼텐셜 영모드를 보존한다. 작은 theta는 모든 좁은 스펙트럼에 대한 미분 전개의 균일 타당성을 뜻하지 않는다. 분류: Counterexample candidate. 새 커널의 유한시간 대조는 단계 61의 무한시간 평형 미달을 통과로 바꾸지 않는다.',
'nonadiabatic-regime':'분류: Counterexample candidate. 같은 kH=1, 2.316801밀리초의 자유 전자 유한 이동 결합 경로는 시간 차수 2.01367, 마지막 시간 차이 1.19624e-5, 각도·산란 구적·셀 내 구적 대조를 통과했다. 원 셀 1/2점 가열 모멘트 오차는 단조 감소하지 않으므로 이 결과를 연속 연산자의 균일 미분 오차 인증으로 확대하지 않는다. 분류: Conjectural. 매질 산란·실제 온도 흡수·대기·GR·전하 연결은 남는다.',
'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 단계 61의 좁은 스펙트럼에 적용한 확산식은 같은 자유 전자의 유한 이동 커널과 1.51638% 다르다. 각도·구적 차이보다 큰 근사 의존성을 분리하고 적분 커널을 실제 결합식에 넣었다. 원 단계 61 집계·온도 조회·평형 실패를 보존한다. 상세평형 대칭화의 점별 상대 차이 최댓값 1과 두 점 구적에서 커진 초기 가열 모멘트 차이도 숨기지 않는다. 작은 가중/끝점 차이는 전 물리 영역의 인증이 아니다.',
'dynamic-charge-completion':'분류: Counterexample candidate. 단계 62 현재: 작은 주파수 이동 전개를 제거하고 자유 전자의 각도·주파수 적분 커널을 같은 EOS·흡수·물질 반동·공간 시간식에 연결했으며 지정 대조를 통과했다. 분류: Conjectural. 다음 물리 입력 병목은 현재 EOS 점유수와 일치하는 집단·결합 전자 산란 및 실제 온도 흡수 보간이다. 비균일 반경 수송·대기·GR 운동량/계량 되먹임, 기존 속도·전하 공간 기준·실제 구동·비교·관측과 전체 완료 요구사항은 유지한다.'}
paths=[]
for stem,body in texts.items():
    p=root/f'docs/{stem}.md';previous=p.read_bytes()
    assert b'## Phase62 ' not in previous
    addition='\n\n## Phase62 — 유한 주파수 이동 광자 커널의 결합\n\n'+body+'\n\n분류: Imported from prior work. 식·근사 판정·수치 대조·보존·예산·남은 경계는 [단계 62 보고서](../notes/REQUEST62_FINITE_JUMP_PHOTONS_KO.md)에 둔다.\n'
    p.write_bytes(previous+addition.encode());assert p.read_bytes().startswith(previous);paths.append(p)
report=root/'notes/REQUEST62_FINITE_JUMP_PHOTONS_KO.md'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,x:p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
source=[root/'verification/def_photon_finite_jump.py',root/'verification/def_photon_finite_jump_validate.py']
files=sorted(p for p in out.iterdir() if p.is_file())+source
manifest=out/'manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='892ea0265',passed=True,
    progress_class='loophole progress; conditional theorem progress',decision='FINITE_JUMP_RESPONSE_ACCEPTED_DIFFUSION_EFFECT_NOT_PHYSICALLY_ADOPTED',
    full_dynamic_charge_solved=False,sha256={p.as_posix():sha(p) for p in files}))
global_path=root/'outputs/direct-eos-gr33/gr-photon-finite-jump-milestone-manifest.json'
write(global_path,dict(classification='Counterexample candidate',checkpoint='892ea0265',
    decision='FINITE_JUMP_RESPONSE_ACCEPTED_DIFFUSION_EFFECT_NOT_PHYSICALLY_ADOPTED',
    prior_milestone_manifest_sha256=sha(root/'outputs/direct-eos-gr33/gr-photon-spatial-coupling-milestone-manifest.json'),
    full_dynamic_charge_solved=False,sha256={p.as_posix():sha(p) for p in paths+[report,manifest]}))
paper=root/'paper/revision-manifest.json';old=json.loads(paper.read_text());assert len(old)==170
entry=dict(classification='Counterexample candidate',progress_class='loophole progress; conditional theorem progress',
    decision='FINITE_JUMP_RESPONSE_ACCEPTED_DIFFUSION_EFFECT_NOT_PHYSICALLY_ADOPTED',finite_jump_coupled_response_passed=True,
    free_electron_diffusion_response_adopted=False,native_medium_scattering_certified=False,actual_temperature_opacity_certified=False,
    physical_atmosphere_closed=False,full_GR_photon_feedback_evolved=False,heat_velocity_time_order_passed=False,
    full_dynamic_charge_solved=False,new_stellar_steps=0,new_native_EOS_calls=0,new_physical_queries=0,source_bindings_verified=11,
    next_bottleneck='EOS-consistent collective and bound-electron scattering and actual-temperature absorption; then stratified transport, physical atmosphere, moving matter/GR feedback, original velocity/charge convergence, physical drive/comparator and observations.',
    report=report.as_posix(),report_sha256=sha(report),evidence_manifest=global_path.as_posix(),evidence_manifest_sha256=sha(global_path))
new=dict(old,request62_photon_finite_jump=entry);write(paper,new)
assert all(json.loads(paper.read_text())[k]==v for k,v in old.items())
print('Finalized',len(paths),'maintained notes; preserved',len(old),'paper entries; evidence files',len(files))
for manifest_path in [manifest,global_path]:
    for key,value in json.loads(manifest_path.read_text())['sha256'].items():assert sha(root/key)==value,key
print('All manifest hashes read back and verified')
