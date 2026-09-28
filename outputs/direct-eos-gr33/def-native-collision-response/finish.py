from pathlib import Path
import hashlib
import json

root=Path.cwd();out=root/'outputs/direct-eos-gr33/def-native-collision-response';mono=out/'monolithic'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

aliases={
out/'plan.json':{'def_native_collision_response.py':out/'inspect-producer.py'},
mono/'plan.json':{'def_native_monolithic_response.py':mono/'first-pilot-producer.py'},
mono/'source-audit-first-plan.json':{'verify_native_collision_response.py':mono/'source-audit-first-producer.py'}}
for p in [out/'plan.json',out/'execution-plan.json',mono/'plan.json',mono/'execution-plan.json',mono/'source-plan.json',mono/'source-audit-first-plan.json']:
    for f,h in read(p)['bindings'].items():
        target=aliases.get(p,{}).get(Path(f).name,Path(f));assert sha(target)==h,(str(p),str(target))
split=read(out/'result.json');result=read(mono/'result.json');audit=read(mono/'source-audit.json');bank=read(out/'bank-result.json');owner=read(out/'check.json')
assert not split['passed'] and result['passed'] and audit['passed'] and bank['passed'] and owner['passed']
fine=result['paths'][1];tc=result['comparisons']['time'];bc=result['comparisons']['background_time']
report=root/'notes/REQUEST124_MONOLITHIC_COLLISION_RESPONSE_KO.md'
report.write_text(f'''# 단계124 — GR 구동 광자와 물질 에너지·H 응답의 동시 풀이

분류: Counterexample candidate. 앞 단계의 실제 lapse·기하 수송 입력을 기존 이동 물질의 흡수·방출·Thomson 산란 및 수정 native EOS의 에너지·중성수소 응답에 연결했다. 원531셀·8각도·152주파수와3.434431ms 구간의 세 경로를 완주했다. 첫 수송/충돌 분할식은 시간 기준에 실패했고, 같은64/128 시간 간격에서 한 단계의 전체 선형 방정식을 동시에 푸는 방식으로 바꿔 원2percent 기준을 통과했다. 실패를 덮거나 격자를 늘리지 않았다. 이전 goal turn의 실제 기하 수송 연결에 이어지는 loophole progress다.

## 계산한 응답의 정확한 범위

분류: Counterexample candidate. 밀도·유속·다른 이온 재고는 수락된 수정 배경의 실제 이력을 따른다. 광자 packet 변화와 물질의 reference 에너지 및 중성수소 변화는 별도 변수로 동시에 적분한다. 전달 운동량을 같은 충돌에서 계상했지만, 이번에는 그 운동량이 속도·질량 유속·운동에너지를 바꾸는 추가 응답을 아직 되돌리지 않았다. 따라서 에너지→온도/압력 변환은 처방된 물질 운동에서의 구성관계 응답이다. 실제 자유 유체의 완전한 열/역학 분리나 닫힌 복사 유체 GR 해라고 부르지 않는다.

분류: Proven. 고정 reference photon energy 좌표의 보존 packet 변분 x와 국소 점유수 변분은 delta_F=x-s*F0, s=3*alpha*delta_phi+delta_lambda로 연결된다. 충돌의 packet 원천 변분에는 체적·좌표 시간 변화, 이 점유수 변화, 고정 재고 압축에 따른 밀도/단열 온도 변화 및 추가 물질 에너지/H 변화가 함께 들어간다. bound-free 반응에서 광자 수와 중성수소 수는 같은 방향으로 변하므로 photon-minus-neutral이 보존 조합이다. photon-plus-material reference 에너지도 같은 충돌 원천을 상쇄한다. 주파수 상자 밖 packet은 별도 포트이며 물질 열로 넘기지 않는다.

분류: Counterexample candidate. 위 보존 대수는 유한 차원 반응 부분계의 항등식이다. 시간에 따라 처방된 밀도·유속의 추가 이류·기계 응답까지 물리적으로 닫혔다는 증명이 아니다. 배경에 직접 더하면 소실될 변화량을1e-26으로 재척도화하여 별도로 풀었으며, 구동의 물리적 크기를 바꾸지 않았다.

## 원 구현 대조와 시간 실패의 해결

분류: Counterexample candidate. 현재 이동 산란의 두 주파수 노드 packet 재분배를 희소 행렬로 표현했다. 세 실제 시각의 산란/주파수 출구 대조와 전체 bound-free+산란 원 구현 대조를 통과했다. 전체 충돌 상대 차이 최대는{max(r['full_collision_owner_relative'] for r in owner['rows']):.3e}다. native 보간에서 온도·중성비·밀도의 상대 미분1e-5/5e-6 대조 최대는{max(v['relative'] for v in bank['derivative_checks']):.3e}였다. 이는 실제 표에 대한 샘플 미분 대조이며 전체 native EOS의 엄밀한 미분 상계가 아니다.

분류: Counterexample candidate. 첫 Strang 수송/충돌 분할식은 모든 경로의 에너지·H 수지를 통과해도 광자 에너지 시간 차이{100*split['comparisons']['time'][0]:.4f}percent, 운동량 충격 차이{100*split['comparisons']['time'][3]:.4f}percent로 실패했다. 그 결과와 생산 코드·계획은 보존했다. 새 방식은 반경/각도 수송, 실제 이동 충돌, 에너지/H를 한 SDIRK 단계에 넣고 전체 잔차1e-10을 요구한다. 기존 분할 역산은 반복 해법의 전처리에만 쓰며 최종 방정식을 분할하지 않는다. 국소 선형 계수는 각 전체 시간 간격의 중간 시각에서 평가하고 기하 구동의 두 단계 시각은 유지한다.

| 분류: Counterexample candidate — 동시 풀이 | 64/128 시간 차이 | 64/128 배경 이력 차이 |
|---|---:|---:|
| 광자 에너지 | {100*tc[0]:.6f}% | {100*bc[0]:.6f}% |
| 물질 에너지 | {100*tc[1]:.6f}% | {100*bc[1]:.6f}% |
| 중성수소 수 | {100*tc[2]:.6f}% | {100*bc[2]:.6f}% |
| 전달 운동량 | {100*tc[3]:.6f}% | {100*bc[3]:.6f}% |
| 광자 반경 압력 | {100*tc[4]:.6f}% | {100*bc[4]:.6f}% |
| 물질 압력 | {100*tc[5]:.6f}% | {100*bc[5]:.6f}% |

분류: Counterexample candidate. 세 경로의 최대 에너지 수지 오차는{max(r['energy_balance_relative'] for r in result['paths']):.3e}, 광자 수−중성수소 수 수지 오차는{max(r['species_balance_relative'] for r in result['paths']):.3e}였다. 모두 원1e-8 문턱 안이다. 전체 단계 잔차 최대는{max(r['linear_relative'] for r in result['paths']):.3e}다. 새 방식의 최대 시간 차이는{100*max(tc):.4f}percent, 최대 배경 차이는{100*max(bc):.4f}percent다. 이 통과는 지정 격자/물리 연산자의 수치 비교이며 반경·각도·주파수 연속 한계나 모든 미시 채널의 인증이 아니다.

## 실제 물질 응답과 다음 GR 입력

분류: Counterexample candidate. fine 경로 끝점의 광자 reference 에너지 변화는{fine['endpoint_photon_reference_energy_erg']:.9e}erg, 물질 reference 에너지 변화는{fine['endpoint_material_reference_energy_erg']:.9e}erg다. 물질 에너지의 셀별 절댓값 합은{fine['endpoint_material_energy_L1_erg']:.9e}erg다. 반응으로 온도·중성비가 변하지만 그 양 자체를 최종 전하 또는 전체 되먹임의 상계로 해석하지 않는다.

분류: Counterexample candidate. 저장된 에너지/H에서 독립적으로 구성관계 변분을 복원한 압력 차이는{audit['pressure_reconstruction_relative']:.3e}다. 최대 abs(delta_logT)는{audit['max_delta_logT']:.3e}, 최대 abs(delta_log_neutral)는{audit['max_delta_log_neutral']:.3e}다. 비활성 셀에 남은 물질 응답은0이다. 추가 gas trace의 최대 L1은{audit['maximum_material_trace_L1_erg']:.9e}erg다. 이미 GR에 포함한 canonical 광자 변화를 빼고 남긴 광자 에너지/반경 압력 최대 L1은{audit['maximum_additional_photon_energy_L1_erg']:.9e}/{audit['maximum_additional_photon_radial_pressure_L1_erg']:.9e}erg다.

분류: Counterexample candidate. monolithic/additional-stress-128.npz는17개 시각의 gas 에너지·반경 압력·trace, 중복 제거한 광자 에너지·반경 압력, 전달 운동량을 포함한다. 동시 풀이 이력에는 광자 반경 압력도 저장했다. 첫 분할 이력은 해당 항을 저장하지 않아 완전한 GR 입력으로 사용할 수 없으며 반복 생산하지 않았다. 이 원천의 추가 GR 진화는 아직 실행하지 않았다.

## 실행 예산과 원 실패 보존

새 native 표를 만들거나 기존 비선형 배경을 반복하지 않았다. 첫 단순 속도 예측은 고정된 계수 행렬 준비 비용을 매 단계에 곱해 과대 추정했다. 준비와 반복 비용을 분리하고 같은 Schur 역산을 두 단계에서 재사용했다. 원240s/480s 계획의 dispatch 미수락 기록은 남아 있다. 각각 실측 예상169s/249s와 유지한 두 배 여유337s/498s를 검토하여 같은 세 경로에만360s/510s 상한을 먼저 재등록했다. 수락 기준이나 물리 격자는 바꾸지 않았다. 첫 분할 생산은{split['seconds']:.2f}s, 수정 동시 생산은{result['seconds']:.2f}s였으며 실패 계산 비용도 보존한다. 원본 및 warm prefix에서 이어 계산했다.

끝점 source audit은 배열 목록에 Python abs를 사용한 요약 코드 오류로 한 번 실패했다. 원 생산자·계획·오류를 남기고 NumPy 절댓값으로 고쳐 저장 결과 판독만 다시 실행했다. 물리 진화는 반복하지 않았다. 독립 끝점 에너지/H 수지 오차는{audit['independent_endpoint_energy_relative']:.3e}/{audit['independent_endpoint_photon_minus_neutral_relative']:.3e}로 통과했다.

분류: Conjectural. 다음 핵심 연결은 전달 운동량·변화한 압력이 만드는 추가 물질 속도와 보존 질량/에너지/재고 유속을 진화시키고, 그 에너지·응력·trace를 계량과 scalar에 되돌리는 일이다. 외부 생성/산란 scalar와 깊은 원천의 인과 경계도 남는다. full_coupled_metric_feedback, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 전체 연구 목표를 축소하거나 완료 처리하지 않는다.
''',encoding='utf-8')

docs=[root/'docs'/name for name in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
append={
'model-definition.md':'분류: Counterexample candidate. 실제 이동 흡수·방출·산란과 보상 photon packet, 물질 reference 에너지/H 변분을 같은 시간식에서 적분했다. 밀도·유속·다른 재고는 저장 배경을 따른다. 운동량 충격은 별도 저장하며 아직 유체 속도·운동에너지로 되돌리지 않았다. 따라서 열/압력 변분은 처방된 물질 운동에서의 구성관계 응답이다.',
'observable-targets.md':'분류: Counterexample candidate. 단순 기하 광자 응답을 실제 물질 에너지/H 응답에 연결하고 지정 시간/배경 대조를 통과했다. 다음 GR용 추가 에너지·반경 압력·trace를 저장했지만 아직 적용하지 않았다. 작은 성분 응답을 전체 물리 전하의 상계나 관측 검출로 사용하지 않는다.',
'adiabatic-limit.md':'분류: Proven. 같은 bound-free 원천에서 photon-plus-material reference 에너지와 photon-minus-neutral 수 조합이 상쇄된다. 주파수 출구는 별도 포트다. 이 부분계의 보존이 추가 유체 운동과 전체 단열/비단열 구분을 닫지는 않는다.',
'nonadiabatic-regime.md':f'분류: Counterexample candidate. 수송/충돌 분할 시간 기준 실패를 보존하고 전체 단계 동시 풀이로 최대 시간 차이{max(tc):.6e}, 배경 차이{max(bc):.6e}를 얻었다. 같은64/128 단계와 물리 격자에서 판정했으며 추가 밀도·속도·GR 되먹임은 미완료다.',
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 첫 Strang 분할은 에너지/H 보존을 통과해도 광자 에너지3.7903percent와 운동량4.6305percent의 시간 차이로 실패했다. 원 실패를 보존하고, 같은 시간 격자의 한 단계 전체 방정식을 동시 풀이해 해결했다. 실패한 자원 예측은 실행 전에 재평가했으며 더 촘촘한 시간 경로로 덮지 않았다. 끝점 감사의 목록 abs 오류는 판독만 수정했다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. actual_monolithic_radiation_thermal_H_response_evolved와additional_GR_source_exported는true다. 추가 물질 밀도·속도/운동에너지·공유 유속, additional_GR_source_applied, full_exterior_scalar, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 다음 병목은 운동량/압력 응답의 실제 물질 이동 및 GR 원천 양방향 연결이다.'}
for p in docs:p.write_bytes(p.read_bytes()+('\n\n## 단계124 — 광자·물질 에너지/H의 동시 응답\n\n'+append[p.name]+'\n\n상세: [단계124 보고](../notes/REQUEST124_MONOLITHIC_COLLISION_RESPONSE_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/n for n in ['def_native_collision_response.py','def_native_monolithic_response.py','verify_native_collision_response.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc','.tmp'] and not p.name.endswith(('-checkpoint.npz','-progress.json')))
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-monolithic-collision-response-manifest.json'
flags=dict(classification='Counterexample candidate',checkpoint='e71ed7eba',split_time_verdict=False,monolithic_response_verdict=True,prescribed_material_motion=True,additional_material_motion_evolved=False,additional_GR_source_exported=True,additional_GR_source_applied=False,final_charge_solved=False,full_goal_complete=False)
write(manifest,dict(flags,files=files));master=root/'paper/revision-manifest.json';data=read(master)
data['sha256'].update(files);data['sha256'][str(manifest.relative_to(root))]=sha(manifest)
data['native_monolithic_collision_response']=dict(flags,report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)));write(master,data)
for p,h in data['sha256'].items():assert sha(root/p)==h,p
for p,v in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
print(json.dumps(dict(passed=True,manifest_files=len(files),master_files=len(data['sha256']),full_goal_complete=False)))
