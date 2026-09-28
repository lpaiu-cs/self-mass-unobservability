from pathlib import Path
import hashlib,json

root=Path.cwd();out=root/'outputs/direct-eos-gr33/def-native-material-response';good=out/'branch-aware'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

# Recover the exact initial dispatch producer from the saved execution source.
# These changes only added the later resource-reassessment entry point.
src=(out/'execution-producer.py').read_bytes()
src=src[:src.index(b'def reassess():')]+src[src.index(b"if __name__",src.index(b'def reassess():')):]
src=src.replace(b"plan=json.loads((OUT/'execution-plan.json').read_text());assert plan['eligible'];",b"plan=json.loads((OUT/'dispatch-plan.json').read_text());assert json.loads((OUT/'pilot.json').read_text())['eligible'];")
src=src.replace(b"signal.alarm(plan['production_seconds']);rows=[]",b"signal.alarm(160);rows=[]")
expected=next(v for p,v in read(out/'dispatch-plan.json')['binds'].items() if p.endswith('def_native_material_response.py'))
assert hashlib.sha256(src).hexdigest()==expected
(out/'dispatch-producer.py').write_bytes(src)
aliases={out/'plan.json':{'def_native_material_response.py':out/'producer-introspection-failure.py'},
    out/'dispatch-plan.json':{'def_native_material_response.py':out/'dispatch-producer.py'},
    good/'plan.json':{'def_native_material_branch_response.py':good/'pilot-producer.py'},
    good/'first-audit-plan.json':{'verify_native_material_response.py':good/'first-audit-producer.py'},
    good/'stable-audit-plan.json':{'verify_native_material_response.py':good/'stable-audit-producer.py'},
    good/'audit-plan.json':{'verify_native_material_response.py':good/'stable-audit-producer.py'}}
plans=[out/'plan.json',out/'dispatch-plan.json',out/'execution-plan.json',good/'plan.json',good/'execution-plan.json',good/'first-audit-plan.json',good/'stable-audit-plan.json',good/'audit-plan.json',good/'source-plan.json']
for p in plans:
    d=read(p)
    for f,h in d.get('bindings',d.get('binds',{})).items():
        target=aliases.get(p,{}).get(Path(f).name,Path(f));assert sha(target)==h,(str(p),str(target))
production=read(good/'production.json');result=read(good/'source-audit.json');failed=read(out/'production.json');first=read(good/'first-audit-result.json')
assert production['passed'] and result['passed'] and not failed['passed'] and not first['passed']
c=result['comparisons'];fine=production['paths'][1];branch=result['maximum_physical_branch_ratio']
report=root/'notes/REQUEST125_CONSERVED_MATERIAL_RESPONSE_KO.md'
report.write_text(f'''# 단계125 — GR·광자 입력으로 구동한 실제 보존 물질 응답

분류: Counterexample candidate. 앞 단계에서 저장만 했던 광자 운동량·에너지·중성수소 교환과 GR 기하 변화를 실제 내부–대기 공유 유속에 넣어, 추가 질량 이동·운동량·보존 에너지·중성수소를 원531셀과3.434431ms 구간 끝까지 진화했다. 원64/128 시간 경로와64/128 배경 대조를 통과했다. 이것은 실제 물질 이동의 loophole progress이며, 광자·계량이 그 이동에 다시 반응하는 전체 폐쇄는 아니다.

분류: Counterexample candidate. 새로운 비선형 항성 배경이나 EOS 표를 만들지 않았다. 저장 배경의17개 시각에서 실제 연산자의 방향 미분을 보간하고, 기존 hydro CFL을 유지한 SSP2를 사용했다. fine에서는 macro 시간과 CFL 상한을 함께 절반으로 했다. 실제 substep 수는{[r['substeps'] for r in production['paths']]}다. 물리 진폭1e-26은 바꾸지 않고 변화량을 별도 변수로 보관했다.

## 보존 변수와 계량 연결

분류: Counterexample candidate. 변수는 셀 바리온 질량B, 물리 반경 운동량 곱하기c, 기준 lapse 에너지Eref, 중성수소 수NH다. 같은 HLL 면 유량을 내부와 대기에 반대 부호로 적용한다. 진공 floor가 제거한 변화량은 독립 signed discard로 남긴다. 기존 내부 중앙 음향식, 실제 EOS 역산, limiter와 shared face를 재사용했다. 내부 밀도/물질 재고의 첫 차수 근사와 제한된 시간·물리 영역을 유지한다.

분류: Proven. shift가 없는 균질 비등방 셀의 물리 운동량 적분P는 dP/dt=-(h_t/h)P를, 기준 에너지는 -a_ref*V*[Pr*h_t/h+2*Pt*R_t/R]의 기하 일을 받는다. 공통 면 유속은 합산 시 상쇄된다. 이는 해당 국소 대수의 증명이며 비선형 GR 수치해 전체의 증명이 아니다. [3+1 보존식의 원문](https://arxiv.org/abs/gr-qc/0201064).

분류: Counterexample candidate. 원 압력 기준 유속을 source로 내보낼 때 누락했던 항과 rest-subtracted 에너지의 소실을 고쳤다. 기존 물질 rate 재현 오차는 최대3.796e-14다. 실제 반경·체적·lapse·반경 길이·lapse 기울기를 변화시키되, 원 압력 불균형을 정역학으로 재설정하지 않았다. 보존 역산에는 운동에너지와 재고 이동이 포함된다. 내부 배경 운동에너지의 기하 일에서 생략된 항의 저장 시각 기반 보수적 크기 추정은{result['omitted_deep_kinetic_metric_work_sampled_envelope_erg']:.6g}erg이며, 엄밀한 연속시간 상계로 부르지 않는다.

## 실제 실패와 수정

분류: Counterexample candidate. 처음 완주한64 경로는 보존 수지를 통과했지만 중간 시각의H 방향 미분 대조가5.54195percent로 실패했다. 확대된 수치 탐침이 작은 기존 유량 부호를 뒤집어 donor 셀을 바꾼 것이 원인이었다. 원 궤적·판정·생산자를 보존했다. 실제 비영 배경 유량의 donor를 유지하고, 정확히0인 면에서는 방향에 따른 donor를 사용해 미분했다. 실제 물리 진폭에서 유량 변화/기존 유량의 저장17시각 최대 비는{branch:.3e}로, 분기 유지 조건 안이다. 이를 모든 limiter나 EOS의 보편 미분 인증으로 승격하지 않는다.

분류: Proven. 비영 유량F의 donor는 abs(deltaF)<abs(F)에서 유지되며, 그 분기에서 delta(F*y)=deltaF*y+F*delta_y다. F=0에서는 deltaF의 방향으로 donor를 정한다. 이 대수와 실제 작은 진폭 검사가 확대 탐침의 잘못된 분기 변경을 분리한다.

분류: Counterexample candidate. 첫 압력 판독은 큰 배경을 빼는 방식 때문에 탐침 대조와 압력 시간 비교에 실패했다. 진화는 다시 돌리지 않았다. 바리온·운동량·에너지·NH 변화로부터 밀도·속도·온도·중성비의 선형 역산을 직접 풀고, 기존 native 밀도/열/화학 미분과 물질 재고 기울기를 적용했다. 개별 셀의 primitive 탐침으로 독립 대조한 압력 오차 최대는{max(v for row in result['pressure_probe_relative'] for v in row):.3e}다. 첫 실패 압력 source는 별도 보존한다.

| 분류: Counterexample candidate — 저장 이력 대조 | 시간 차이 | 배경 이력 차이 |
|---|---:|---:|
| 바리온 재배치 | {100*c['time'][0]:.6f}% | {100*c['background'][0]:.6f}% |
| 반경 운동량 | {100*c['time'][1]:.6f}% | {100*c['background'][1]:.6f}% |
| 기준 에너지 | {100*c['time'][2]:.6f}% | {100*c['background'][2]:.6f}% |
| 중성수소 | {100*c['time'][3]:.6f}% | {100*c['background'][3]:.6f}% |
| 추가 반경 압력 | {100*c['stress_time'][1]:.6f}% | {100*c['stress_background'][1]:.6f}% |
| 추가 trace | {100*c['stress_time'][2]:.6f}% | {100*c['stress_background'][2]:.6f}% |

분류: Counterexample candidate. 원2percent 시간/배경 기준을 유지했다. 네 보존량의 독립 저장 수지 오차 최대는{result['independent_balance_relative']:.3e}, 실제 생산의 half-probe 방향 대조 최대는{max(r['directional_relative'] for r in production['paths']):.3e}다. 전 시각 forward-probe 지표는 영에 가까운 rate에서 큰 상대값을 보이므로 별도 진단으로 보존하고, 등록한 evolved-state Richardson 대조와 혼동하지 않는다.

## 얻은 GR 입력과 남은 경계

분류: Counterexample candidate. fine 끝점의 셀별 바리온 변화 절댓값 합은{fine['endpoint_L1_physical'][0]:.9e}g, 기준 에너지 절댓값 합은{fine['endpoint_L1_physical'][2]:.9e}erg다. 작은 이동이어도 rest mass를 포함한 추가 에너지의 최대 공간L1은{result['maximum_additional_source_L1'][0]:.9e}erg다. 따라서 전달 열만 보고 이 응답을 생략할 수 없다. 이 수치를 순 에너지 생성이나 관측 신호로 해석하지 않는다.

분류: Proven. 기준 에너지의 정의에 따라 delta(VE)=[deltaEref+a_surface*cx*c^2*deltaB]/a_ref다. 고정 반경 밀도 source에는 여기서 기존 시간 배경의E*deltaV를 빼야 한다. 이전 GR 연산자가 포함한 초기 canonical 기체 반응을 다시 더해 분리하면, 추가 강제원천은 eF=delta(VE)-E_bg*deltaV+(E_initial+P_initial)*deltaV이며 압력은 pF=delta(VP)-P_bg*deltaV+K_initial*deltaV다. 같은 체적 기준에서 적용한다.

분류: Counterexample candidate. 세 material-source 파일에 정지질량 이동·운동량·reference 에너지·NH와 추가 proper 에너지·반경/접선 압력·trace를 저장했다. 초기 canonical 반응을 중복 제거했고 시간 배경과 초기 배경의 차이를 남겼다. 실패한 압력 판독의 최대 branch 비를 재사용해 측정하지 않은0을 보고하지 않는다. 최종 입력은 source-audit.json이며 원 result.json과 first-audit-result.json도 보존한다.

첫 실패 생산은{failed['seconds']:.2f}s, 수정 세 경로 생산은{production['seconds']:.2f}s였다. 초기 짧은 구간만으로 계산한 비용이 후반 CFL을 과소평가하여 dispatch를 거절한 기록을 보존했다. 실제 실패 전체 경로의 호출 수와 수정 prefix의 단위 비용을 재사용해 예상225.85s, 두 배 여유461.70s, 상한470s를 실행 전에 등록했다. 물리 기간·공간·시간 경로 수·수락 기준은 늘리지 않았다.

분류: Conjectural. 다음 결정적 작업은 이 새 물질 운동을 광자 충돌/수송과 GR/scalar에 되돌려 양방향 응답을 닫고, 최종 전하에서의 잔여 또는 제한을 다시 계산하는 것이다. 미표현 외부·깊은 원천, 작은 물질 포트의 공간 수렴, 전체 EOS 미분 인증과 관측 연결도 남는다. additional_material_motion_evolved와additional_GR_source_exported는true지만 full_photon_material_feedback, additional_GR_source_applied, nonlinear_GR, final_charge_solved, full_goal_complete는false다.
''',encoding='utf-8')
docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in docs})
notes={
'model-definition.md':'분류: Counterexample candidate. 실제 GR 기하와 광자 전달 에너지/H/운동량을 공유 물질 유속에 적용해 추가 바리온·운동량·기준 에너지·중성수소를 원 전체 구간 진화했다. native 보존 역산에 운동에너지·재고 이동을 포함한다. 광자와 GR에 되돌리는 다음 sweep은 미완료다.',
'observable-targets.md':'분류: Counterexample candidate. 추가 물질 이동의 실제 시간/배경 비교를 통과했다. 작은 질량 이동도 정지질량 에너지 source를 만들므로 전달 열만으로 생략하지 않는다. 새 원천의 GR 재적용과 최종 물리 전하는 아직 판정하지 않는다.',
'adiabatic-limit.md':'분류: Proven. 비영 배경 유량의 donor 분기를 유지할 때 delta(F*y)=deltaF*y+F*delta_y다. 정확히 영인 유량은 방향에 따라 donor를 선택한다. 기준 에너지에서 정지질량 이동을 복원하고 초기 canonical 압축 반응을 분리해야 GR source를 이중 계산하지 않는다.',
'nonadiabatic-regime.md':f"분류: Counterexample candidate. 같은531셀·3.434431ms의 물질 응답에서 시간 차이 최대{max(c['time']):.6e}, 추가 stress까지 배경 차이 최대{max(c['stress_background']):.6e}를 얻었다. 최대 독립 수지 오차{result['independent_balance_relative']:.3e}다. 이는 처방된 GR/광자 입력에 대한 실제 물질 sweep이며 양방향 전체 폐쇄가 아니다.",
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 처음 완주한 물질 응답은 확대 미분 탐침의 donor 부호 전환 때문에H 대조5.54percent로 실패했다. 원 donor 분기를 유지하는 미분으로 해결했다. 첫 압력 판독의 소실은 보존 변분의 직접 primitive 역산과 native 미분을 사용해 저장 해만 다시 읽어 해결했다. 두 원 실패·생산자·수치 판정은 보존한다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. additional_material_motion_evolved와additional_GR_source_exported는true다. actual_two_way_photon_material_metric_response, additional_GR_source_applied, full_exterior_scalar, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 다음 병목은 새 속도·밀도·압력·재고의 광자/GR 재반영이다.'}
for p in docs:p.write_bytes(p.read_bytes()+('\n\n## 단계125 — 실제 보존 물질 응답\n\n'+notes[p.name]+'\n\n상세: [단계125 보고](../notes/REQUEST125_CONSERVED_MATERIAL_RESPONSE_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/n for n in ['def_native_material_response.py','def_native_material_branch_response.py','verify_native_material_response.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc','.tmp'] and not p.name.endswith('-progress.json'))
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-conserved-material-response-manifest.json'
flags=dict(classification='Counterexample candidate',checkpoint='5acae2934',original_directional_verdict=False,original_pressure_readout_verdict=False,corrected_material_response_verdict=True,additional_material_motion_evolved=True,additional_GR_source_exported=True,full_photon_material_feedback=False,additional_GR_source_applied=False,final_charge_solved=False,full_goal_complete=False)
write(manifest,dict(flags,files=files));master=root/'paper/revision-manifest.json';data=read(master);data['sha256'].update(files);data['sha256'][str(manifest.relative_to(root))]=sha(manifest)
data['native_conserved_material_response']=dict(flags,report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)));write(master,data)
for p,h in data['sha256'].items():assert sha(root/p)==h,p
for p,v in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/p).read_bytes()[:v['bytes']]).hexdigest()==v['sha256'],p
print(json.dumps(dict(passed=True,manifest_files=len(files),master_files=len(data['sha256']),full_goal_complete=False)))
