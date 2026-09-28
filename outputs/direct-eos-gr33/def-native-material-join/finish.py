"""Bind the inspected Phase118 outcomes; retain the rejected first trajectory."""
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parents[3];out=root/'outputs/direct-eos-gr33/def-native-material-join';center=out/'centered'
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
prod=read(center/'production.json');result=read(center/'result.json');audit=read(center/'audit.json');join=read(center/'join-audit.json')
assert all(x['passed'] for x in [prod,result,audit,join])
original=read(out/'result.json');failure=read(out/'join-audit.json');diagnostic=read(out/'diagnostic.json')
assert original['passed'] and not failure['passed']
gr=read(center/'gr/result.json');bound=read(center/'gr/bound.json');ga=read(center/'gr/audit.json')
assert all(x['passed'] for x in [gr,bound,ga])
fine=prod['paths'][-1];coarse=prod['paths'][0]
docs=[root/'docs'/x for x in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
text=f'''

## 수정된 실제 결합 진화와 전하 판정

분류: Counterexample candidate. 중앙 유속의64/128 경로가 모두3.434ms를 완주했고, 실제 공유 면의 양방향 물질 교환을 유지했다. 아래 값은 이 선언 모형에 대한 결과이며 전체 공간·EOS·동적 GR 오차 구간이 아니다.

| 항목 | 결과 |
|---|---:|
|64/128 생산 경로 시간 |{coarse['seconds']:.2f}s / {fine['seconds']:.2f}s |
|fine 결합 에너지 잔차 |{fine['energy_relative']:.5e} |
|fine 바리온 잔차 / 초기 대기량 |{fine['joint_baryon_relative_to_initial_atmosphere']:.5e} |
|fine 바리온 잔차 / 실제 교환량 |{fine['joint_baryon_relative_to_actual_port']:.5e} |
|fine 최대 상대 밀도 변화 |{fine['maximum_density']:.5e} |
|직접 전하의64/128 차이 |{100*result['controls']['time_direct']:.5f}% |
|저장 간격 대조의 직접 전하 차이 |{100*result['controls']['cadence_direct']:.5f}% |
|native 끝점 압력/에너지 앵커 오차 |{audit['native_endpoint_pressure_energy_relative']:.5e} |
|직접 전하 변화 |{result['endpoint_direct_relative']:+.10e} |
|직접+선언된 광자 질량 정규화 합 |{result['endpoint_direct_plus_photon_mass']:+.10e} |
|추가 GR와 외부 광자/퍼텐셜 상계 후 직접 하한 |{ga['arbitrary_outward_angular_mass_conditional_lower']:+.10e} |

분류: Counterexample candidate. 공유 면의 실제 누적 질량은 `{join['integrated_shared_mass_g']:.8e}g`, 반경 운동량 유속 적분은 `{join['integrated_shared_momentum_g_cm_s']:.8e}g cm/s`, Killing 에너지 유속 적분은 `{join['integrated_shared_Killing_energy_erg']:.8e}erg`다. 부호는 바깥 방향 양수다. 운동량 유속은 압력 응력을 포함하므로 순 이류 운동량이나 분출량으로 해석하지 않는다. 끝점의 확장 정밀도 유속 대조는 원 `1e-8` 기준을 통과했지만 전체 중성수소 반응 누적 원장은 여전히 미완료다.

분류: Counterexample candidate. 수정된 새 저장 원천에 단계117의 GR 식과 외부 광자·퍼텐셜 상계를 다시 적용했다. 실제 물질/광자의 반경 응력과 안쪽 광자 에너지 차감을 사용하며, 이전 수치 구간을 복사하지 않았다. 새 추가 GR 두 사건 오차 상계는 `{ga['observer_difference_GR_error_bound']:.8e}`, 직접+광자 질량 합의 명목 중심은 `{ga['endpoint_nominal_direct_plus_photon_mass']:.8e}`이다. 이 구간은 해당 추가 GR 항만 포함하며 원천의 공간 오차나 초기 Einstein 제약 오차를 포함하지 않는다.

분류: Counterexample candidate. 이 단계는 실제 자유 물질 유속을 결합 진화에 적용하고, 그 과정에서 발견한 인공 음향 질량 확산을 수정한 loophole progress다. 첫 시간 수렴 음수 결과와 strict 중성수소 합산 실패는 보존한다. 시간 대조·보존의 통과만으로 작은 정지질량 전하가 물리적인 것은 아니라는 구체적 반례다. 고친 결과 역시 반경/주파수/내부 각도 및 초기 제약이 닫히기 전까지 최종 항성 전하로 승격하지 않는다.

분류: Conjectural. 다음 우선순위는 같은 수정 유속의 실제 원천 공간 오차, 특히 약34km 중심-경계 재구성이 잔여 전하를 얼마나 바꾸는지 결정하는 일이다. 작은 GR 구적의 추가 반복은 우선순위가 낮다. 새 공간 계산은 필요한 정확도·최소 변경·저장 재사용·실측 자원 예산을 먼저 정하고 진행해야 한다.
'''
report=root/'notes/REQUEST118_NATIVE_MATERIAL_JOIN_KO.md';report.write_text(report.read_text()+text)
common='\n\n## 단계118 — 실제 공유 물질 경계와 인공 음향 확산 수정\n\n'
values={
 'model-definition.md':'분류: Counterexample candidate. 내부 셀 운동량을 실제 진화시키고400m 내부 공유 면의 질량·운동량·Killing 에너지·중성수소 유속을 동일 SSP 단계에서 양쪽에 반대 부호로 적용했다. 실제 내부 상태를 저장 native 면 배경에 재구성하여 대기 HLL 해법에 입력한다. 거친 내부 upwind 압력 점프 유속이 작은 정지질량 전하를 지배하여 실패했으므로, 같은 짧은 평활 구간의 내부 중앙 음향 유속으로 수정했다. 광자·EOS·원 경계 HLL 유속과 물리 초기 힘은 유지한다. 반경 재구성 및 일차 재고 구성 관계의 물리 오차는 미완료다.',
 'observable-targets.md':f'분류: Counterexample candidate. 수정된 공유 물질 진화의 직접 전하 변화는{result["endpoint_direct_relative"]:.8e}, 실제 원천의 추가 GR·외부 광자/퍼텐셜 상계 후 조건부 직접 하한은{ga["arbitrary_outward_angular_mass_conditional_lower"]:.8e}다. 첫 upwind 시간 수렴 음수 결과는 압력 점프 질량 확산이 거의 전부를 설명하여 물리 결과로 채택하지 않았다. 새 양수도 원천 공간·초기 Einstein 제약·동적 GR·정적 비교 모형이 닫히지 않은 조건부 성분 결과다.',
 'adiabatic-limit.md':'분류: Proven. 같은 수치 면 유속의 반대 부호 소거와 평탄 상수 계수 중앙 음향 연산자의 반대칭성을 검증했다. RK2 모드 증폭률 제곱은1+x^4/4이므로 무조건 안정이라고 주장하지 않는다.\n\n분류: Counterexample candidate. 이 짧은 구간의 물질 접합·전하 시간 대조는 정적 EFT 흡수 경계를 바꾸거나 실제 궤도 완화 검출을 뜻하지 않는다.',
 'nonadiabatic-regime.md':'분류: Counterexample candidate. 지정된 물질 복사 조건 대신 실제 양쪽 상태가 결정하는 공유 물질 유속을 사용해 두 시간 경로를 완주했다. 최초 내부 음향 확산의 큰 음수 전하는 수정 뒤 제거되었다. 같은 초기 과도 스펙트럼과 고정 계량의 비선형 열/H/광자 진화이며 물리적 주기 구동·전역 비선형 GR 또는 관측 식별성 완료가 아니다.',
 'failure-ledger-dynamic-chi.md':f'분류: Counterexample candidate. 첫 공유 유속 진화는 보존과 시간 대조를 통과했으나 직접 성분{original["endpoint_direct_relative"]:.8e}의 약100%가 거친 음향 압력 점프 질량 확산으로 재구성되었다. 원 음수·원 소스·전체 이력은 보존하고, 공유 HLL 경계를 유지한 내부 중앙 유속으로 실제 결합 경로를 수정했다. binary64 내부 중성수소 발산 합산은1.7804e-7로 실패했고 같은 면 유속의 확장 정밀도 대조는2.3214e-11이었다. 원 실패를 지우거나 전체 반응 원장 통과로 바꾸지 않는다. 압력 점프 HLL 초기 질량 검사 오류와 super 클래스 폐쇄 오류도 보존했다.',
 'dynamic-charge-completion.md':'분류: Counterexample candidate. actual_shared_mass_momentum_energy_neutral_flux, corrected_coupled64_and128_completed, original_time_energy_baryon_gates_passed, actual_source_GR_response_and_conditional_bound는true다. original_upwind_physical_charge_accepted와original_binary64_join_audit_passed는false로 유지한다. complete_neutral_trajectory_audit, physical_interface_reconstruction_certified, source_spatial_frequency_inner_angle_error_certified, initial_full_Einstein_constraints_matched, full_dynamic_GR_feedback, final_charge_solved, full_goal_complete는false다. 다음은 원천 공간·초기 제약의 결정적 오차를 닫는 일이다.'}
for p in docs:p.write_bytes(p.read_bytes()+(common+values[p.name]+'\n\n상세: [단계118 보고](../notes/REQUEST118_NATIVE_MATERIAL_JOIN_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/x for x in ['def_native_material_join.py','verify_native_material_join.py','def_native_material_centered.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith(('-checkpoint.npz','.tmp')))
files={str(p.relative_to(root)):sha(p) for p in paths}
manifest=root/'outputs/direct-eos-gr33/native-material-join-manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='b024c099c',actual_shared_material_interface_evolved=True,
    original_upwind_physical_charge_rejected=True,original_strict_neutral_sum_failure_preserved=True,corrected_time_and_conservation_passed=True,
    actual_source_GR_conditional_bound_passed=True,full_goal_complete=False,final_charge_solved=False,files=files))
master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_material_join']=dict(classification='Counterexample candidate',report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)),
    actual_shared_material_evolution=True,acoustic_artifact_rejected_and_repaired=True,conditional_positive_remnant=True,
    source_continuum_certified=False,initial_Einstein_constraints_matched=False,final_charge_solved=False,full_goal_complete=False)
write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,r in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:r['bytes']]).hexdigest()==r['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),passed=True)))
