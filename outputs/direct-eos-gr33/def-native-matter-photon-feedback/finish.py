from pathlib import Path
import hashlib,json
import numpy as np

root=Path.cwd();out=root/'outputs/direct-eos-gr33/def-native-matter-photon-feedback';stage=out/'stage-time';ret=stage/'return';gr=ret/'gr'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

aliases={out/'plan.json':{'def_native_matter_photon_feedback.py':out/'first-check-producer.py'},
    out/'pre-history-execution-plan.json':{'def_native_matter_photon_feedback.py':out/'pilot-producer.py'},
    out/'return/plan.json':{'def_native_feedback_return.py':out/'return/planned-producer.py'},
    stage/'plan.json':{'def_native_stage_feedback.py':stage/'producer.py'},
    stage/'execution-plan.json':{'def_native_stage_feedback.py':stage/'producer.py'},
    ret/'plan.json':{'def_native_stage_feedback.py':stage/'producer.py','def_native_feedback_return.py':ret/'planned-producer.py'},
    ret/'execution-plan.json':{'def_native_feedback_return.py':ret/'production-producer.py'},
    ret/'probe-scaled/plan.json':{'def_native_feedback_return.py':ret/'source-producer.py'},
    ret/'probe-scaled/execution-plan.json':{'def_native_feedback_return.py':ret/'source-producer.py'}}
plans=list(out.rglob('*plan.json'))
for p in plans:
    for f,h in read(p).get('bindings',{}).items():
        target=aliases.get(p,{}).get(Path(f).name,Path(f));assert sha(target)==h,(str(p),str(target),h,sha(target))
failed=read(out/'result.json');photons=read(stage/'audit.json');material=read(ret/'material-audit.json');sources=read(ret/'sources.json');fields=read(gr/'result.json')
assert not failed['passed'] and all(x['passed'] for x in [photons,material,sources,fields])
fine=fields['paths'][0];old=read(root/'outputs/direct-eos-gr33/def-native-anisotropic-gr/audit.json')
ratio=fine['endpoint_compact_with_metric']/old['endpoint_compact_with_metric']
summary=dict(classification='Counterexample candidate',midpoint_frozen_time_verdict=False,stage_time_verdict=True,
    original_stage_time_conservation_verdict=False,refined_linear_solve_verdict=True,
    original_material_return_directional_verdict=False,material_probe_precision_verdict=True,
    actual_material_response_returned_to_photons=True,returned_collision_transfer_applied_to_material=True,returned_source_applied_to_represented_GR=True,
    additional_compact_charge=fine['endpoint_compact_with_metric'],additional_over_prior_compact_charge=ratio,
    material_energy_H_waveform_residual=sources['paths'][1]['energy_H_waveform_residual'],
    coupled_fixed_point_verified=False,updated_GR_reapplied_to_transport=False,full_exterior_scalar=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
write(out/'closure-summary.json',summary)
pc=photons['comparisons'];mc=sources['comparisons'];gc=fields['comparisons'];res=sources['paths'][1]['energy_H_waveform_residual']
report=root/'notes/REQUEST126_MATTER_PHOTON_GR_RETURN_KO.md'
report.write_text(f'''# 단계126 — 물질 운동을 광자 교환과 GR에 되돌린 실제 응답

분류: Counterexample candidate. 단계125의 실제 바리온·운동량·재고 이동을 움직이는 광자 충돌에 적용했다. 광자와 물질의 에너지·중성수소를 동시에 진화한 다음, 새 충돌 에너지·중성수소·운동량 전달을 실제 공유 물질 유속에 다시 적용하고 그 추가 에너지·압력·trace를 GR/scalar 방정식에 넣었다. 동일531셀·8각도·152주파수·3.434431ms와 원64/128 경로를 사용했다. 이 작업은 입력을 내보내는 단계에서 실제 귀환 응답을 계산하는 단계로 전진했다.

분류: Counterexample candidate. 이것은 저장 GR 입력에서 수행한 한 번의 waveform 귀환이다. 광자 단계의 바리온/운동량과 비충돌 수송은 이전 물질 이력에서 공급되며, 에너지/H만 광자와 동시 미지수다. 다음 물질 단계가 그 이력을 갱신한다. 새 GR를 다시 수송에 반영한 고정점이나 완전한 비선형 항성 진화를 완료했다고 해석하지 않는다.

분류: Counterexample candidate. 이전 물질 응답과 이번 귀환 이력의 차이는, 새 응답의 최대 공간L1으로 나눴을 때 기준 에너지{100*sources['paths'][1]['material_sweep_change'][2]:.4f}%, 중성수소 변화량{100*sources['paths'][1]['material_sweep_change'][3]:.4f}%다. 이는 미소 응답끼리의 상대 차이이며 항성 전체의 에너지나 수소량이 그 비율로 변했다는 뜻이 아니다. 실제 재교환을 생략한 응답과 구별되는 결과다.

## 실패를 해결한 실제 적분 변경

분류: Counterexample candidate. 첫 보존 변수 역산에서는 큰 reference rest-energy 이동과 작은 열에너지를 함께 다루며 소실이 발생했다. Etilde=Eref-(a_ref-a_surface)*cx*c²*B로 변수를 바꾸고 최종 판독에는 그 항을 복원했다. 밀도·속도·온도·중성비와 실제 재고 이동의 native 미분을 사용했다. 첫 실패 검사, 생산자, 직렬화 실패를 보존한다.

분류: Proven. 이 변수 변환은 가역적이다. 상대론적 비정지 에너지 K=ρcx·c²(W²-W)+ρuW²+p(W²-1)와 속도 미분을 기호 계산으로 확인했다. 여기서ρ는 보존 바리온 질량 밀도이며 cx는 선언 EOS의 정지에너지 계수다. 에너지 기준 변경으로 질량 수송이나 물리적 일을 제거하지 않는다.

분류: Counterexample candidate. 첫 전체 광자 귀환은 수지는 보존했지만 중성수소64/128 시간 차이{100*failed['comparisons']['time'][2]:.4f}%로 원2% 기준에 실패했다. 선형 연산자와 입력을 전체 시간 단계의 중간값으로 고정한 방식을, 각 SDIRK 단계의 실제 시각에서 평가하도록 수정했다. 원 실패 궤적과 판정을 보존하며 시간 간격·공간 격자·수락 기준을 바꾸지 않았다.

분류: Proven. y'=L(t-y)에서 L→∞일 때 중간값 고정 SDIRK는 y_(n+1)→t_n+h/2, 실제 단계 시각 SDIRK는 y_(n+1)→t_n+h다. 이 시험은 빠른 반응에서 첫 방식의 반 단계 평형 지연을 입증한다. 실제 다변수계의 모든 오차를 이 한 식으로 보증하지 않는다.

분류: Counterexample candidate. 수정 단계 시각의 마지막 배경64 경로는 전체 선형 잔차1e-10을 통과했지만 중성수소 수지2.804e-8로 원1e-8 기준에 실패했다. 통과한 두 배경128 경로는 그대로 재사용하고, 실패 경로만 전체 단계 잔차와 GMRES 기준을1e-12로 강화해 다시 풀었다. 원 단계 시각 result.json은 실패로 보존하며, 통과 경로를 결속한 audit.json만 수락 판정이다. 다른 가중치의 보존 수지가 전체 유클리드 잔차에서 자동 보장되지는 않는다.

분류: Counterexample candidate. 실제 물질 귀환의 첫 배경64 경로도 끝점 운동량 미분 대조0.305percent로 실패했다. 저장 상태에서 오차가 깊은0–8셀에 집중되고 탐침을 줄일수록 커짐을 확인했다. probe4/8 차이는0.0341percent,8/16 차이는0.0135percent였다. 큰 평형 힘의 차분 소실을 줄이도록 Richardson 탐침을8배로 잡고 실패 경로만 재진화했다. 원 물리 진폭과 donor 분기 및0.2percent 기준은 유지했고, 새 진화 끝점에서4/8과8/16을 모두 확인했다. 원 production.json은 실패로 남기고 material-audit.json으로 수락 경로를 결속한다.

| 분류: Counterexample candidate — 전체 이력의 L1 대조 | 시간 차이 | 배경 이력 차이 |
|---|---:|---:|
| 광자 에너지 | {100*pc['time'][0]:.6f}% | {100*pc['background'][0]:.6f}% |
| 동시 물질 기준 에너지 | {100*pc['time'][1]:.6f}% | {100*pc['background'][1]:.6f}% |
| 동시 중성수소 | {100*pc['time'][2]:.6f}% | {100*pc['background'][2]:.6f}% |
| 실제 귀환 물질 에너지 | {100*mc['time'][2]:.6f}% | {100*mc['background'][2]:.6f}% |
| 실제 귀환 물질 중성수소 | {100*mc['time'][3]:.6f}% | {100*mc['background'][3]:.6f}% |
| GR/scalar field 최대 norm | {100*gc['time']:.6f}% | {100*gc['background']:.6f}% |

분류: Counterexample candidate. 광자 경로의 최대 에너지 수지 오차는{max(r['energy_balance_relative'] for r in photons['paths']):.3e}, 종 수지 오차는{max(r['species_balance_relative'] for r in photons['paths']):.3e}다. 실제 물질 귀환의 보존/압력 판독과 원 방향 미분 기준도 통과했다. GR 적분4/8차 대조는{gc['quadrature']:.3e}, 독립 Jordan 반경 직접 적분과의 차이는{fields['independent_direct_relative']:.3e}다. 이는 선언된 유한 모형의 수치 대조이며 연속 EOS·미분의 엄밀한 전체 오차 보장이 아니다.

## 전하에 적용한 결과와 남은 오차

분류: Counterexample candidate. 추가 compact scalar 전하 성분은{fine['endpoint_compact_with_metric']:.12e}다. 이전 표현 compact 성분 대비 signed 비는{ratio:.9e}다. 이 값은 이번 귀환 원천의 추가 응답이며 최종 물리 전하나 기존 양의 하한의 승계가 아니다. 실제 광자 안쪽 포트 debit와 바깥쪽 유출을 함께 저장했고, 초기 GR 연산자의 canonical 응답과 시간에 따라 바뀐 배경의 체적 항을 구별해 source를 구성했다.

분류: Counterexample candidate. 동시 광자 단계의 물질 E/H와 그 전달을 실제 물질 유속에 다시 적용한 이력 사이의 상대 차이는 각각{res[0]:.6e}, {res[1]:.6e}다. 이는 다음 결합 반복에서 해결할 이력 불일치이며, 스스로 전체 해의 오차 상계를 주지 않는다. 같은 입력에서 원 두 시간 경로가 합치하는 것과 결합 고정점에 도달한 것은 별도 판정이다.

분류: Conjectural. 다음은 남은 물질 이력 차이와 추가 GR의 수송 귀환에 필요한 오차 한도를 최종 전하 규모와 연결하고, 아직 포함하지 않은 외부 생성/산란 scalar·깊은 원천을 닫는 일이다. 이번 추가 성분이 매우 작은 만큼 같은 국소 반복을 자동 확대하지 않고 저장 이력과 인과 경계에서 먼저 제한한다. 작은 물질 포트의 공간 수렴, 전체 EOS 미분 인증과 관측 연결도 남아 있다.

첫 중간값 고정 생산은{failed['seconds']:.2f}s, 단계 시각 수정 생산은{photons['seconds']:.2f}s, 실제 물질 귀환은{material['seconds']:.2f}s, source 구성은{sources['seconds']:.2f}s, GR 적용은{fields['seconds']:.2f}s였다. 새 EOS 표나 비선형 배경을 재계산하지 않았고, 각 실행 전 저장 결과와 prefix 실측으로 예산을 등록했다.

분류: Counterexample candidate. actual_material_response_returned_to_photons, returned_collision_transfer_applied_to_material, returned_source_applied_to_represented_GR는true다. coupled_fixed_point_verified, updated_GR_reapplied_to_transport, full_exterior_scalar, nonlinear_GR, final_charge_solved, full_goal_complete는false다.
''',encoding='utf-8')
docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in docs})
notes={
'model-definition.md':'분류: Counterexample candidate. 보존 물질의 밀도·속도·재고 이동을 실제 움직이는 광자 충돌에 적용하고, 새 교환량을 다시 실제 물질 유속과 GR/scalar에 적용했다. 광자/에너지/H는 단계 시각의 동시 SDIRK 미지수이며 바리온·운동량은 이전 물질 이력에서 공급된다. 한 귀환 sweep을 고정점으로 부르지 않는다.',
'observable-targets.md':f'분류: Counterexample candidate. 새 물질·광자 원천의 추가 compact 전하 성분은{fine["endpoint_compact_with_metric"]:.9e}다. 이는 이전 표현 성분에 대한 추가 응답이며 최종 물리 전하·관측 검출·양의 하한 승계가 아니다.',
'adiabatic-limit.md':'분류: Proven. y_prime=L*(t-y)의 무한히 빠른 완화 극한에서 중간값 고정 SDIRK는 t+h/2로 가지만 실제 단계 시각 SDIRK는 t+h로 간다. 빠른 반응의 평형 추종을 보존하려면 시간 의존 계수/입력의 평가 시각을 구별해야 한다.',
'nonadiabatic-regime.md':f'분류: Counterexample candidate. 실제 광자 귀환과 물질 재진화 및 compact GR 적용이 원 시간/배경 기준을 통과했다. 물질E/H 이력 불일치{res[0]:.6e}/{res[1]:.6e}를 남겼다. 새로운 GR의 수송 재반영과 전체 결합 고정점은 미완료다.',
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 첫 물질→광자 귀환의H 시간 대조6.6447percent 실패는 실제 SDIRK 단계 시각에 계수/입력을 평가해 수정했다. 다음 배경64 경로의 종 수지 실패는 선형 풀이 기준 강화로, 물질 귀환의 깊은 힘 차분 소실은 양쪽 크기 대조를 거친 탐침으로 실패 경로만 재계산해 수정했다. 가역적 비정지 에너지 좌표와 판독 배열 수정도 적용했다. 원 실패 판정·궤적·생산자는 보존한다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. 실제 물질→광자→물질→추가compact GR 귀환을 적용했다. coupled_fixed_point_verified, updated_GR_reapplied_to_transport, full_exterior_scalar, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 다음은 저장 이력의 물질/광자 불일치와 새로운 GR의 실제 수송 되먹임을 해결하는 일이다.'}
for p in docs:p.write_bytes(p.read_bytes()+('\n\n## 단계126 — 실제 물질·광자·GR 귀환\n\n'+notes[p.name]+'\n\n상세: [단계126 보고](../notes/REQUEST126_MATTER_PHOTON_GR_RETURN_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/n for n in ['def_native_matter_photon_feedback.py','def_native_stage_feedback.py','def_native_feedback_return.py','def_native_material_return_precision.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc','.tmp'] and not p.name.endswith(('-progress.json','-checkpoint.npz')))
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-matter-photon-gr-return-manifest.json'
write(manifest,dict(summary,checkpoint='d6ccddff9',files=files));master=root/'paper/revision-manifest.json';data=read(master);data['sha256'].update(files);data['sha256'][str(manifest.relative_to(root))]=sha(manifest)
data['native_matter_photon_gr_return']=dict(summary,report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)));write(master,data)
for p,h in data['sha256'].items():assert sha(root/p)==h,p
for p,v in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
print(json.dumps(dict(passed=True,manifest_files=len(files),master_files=len(data['sha256']),full_goal_complete=False)))
