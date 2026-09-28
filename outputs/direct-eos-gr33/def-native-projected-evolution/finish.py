"""Preserve actual projected evolution, its readout and its remaining limits."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

production=read(out/'production.json');result=read(out/'result.json');audit=read(out/'audit.json')
comparison=read(out/'comparison.json');check=read(out/'check.json');installed=read(out/'installed-source-audit.json')
ports=read(out/'material-ports.json')['ports']
assert production['passed'] and audit['passed'] and check['passed']
paths=production['paths'];fine=paths[-1];controls=result['controls']
assert len(paths)==2 and [p['completed_steps'] for p in paths]==[64,128]
plan=read(out/'plan.json');pilot_hash=sha(out/'pilot-producer.py')
assert next(v for k,v in plan['bindings'].items() if k.endswith('def_native_projected_evolution.py'))==pilot_hash
assert sha(root/'verification/def_native_projected_evolution.py')==read(out/'restart-receipt.json')['source_sha256']
for p,h in read(out/'readout-plan.json')['bindings'].items():
    assert sha(out/'readout-producer.py' if p.endswith('verify_native_projected_evolution.py') else p)==h,p

report=root/'notes/REQUEST121_PROJECTED_COUPLED_EVOLUTION_KO.md'
report.write_text(f'''# 단계121 — 수정 초기 GR 상태의 실제 물질·광자 결합 진화

분류: Counterexample candidate. 단계120의 실제 셀 재고를 보존하는 초기 질량·lapse·순간 스칼라 균형 해를 실제 진화기에 설치했다. 새 초기값에서 원19개 내부 셀+512개 대기 셀,152주파수,8각도,3.434431ms의64/128 두 경로를 모두 완주했다. 생산 벽시간은 {production['seconds']:.2f}s로 사전500s 상한 안이다. 새 초기값의 진화 연결은 완료했으나, 진화 중 계량·스칼라는 고정되어 완전한 GR 결합 진화는 아니다.

분류: Counterexample candidate. 바꾼 항은 Einstein 고정 반경의 Jordan 반경 대응, 실제 proper 셀 체적, 모든 수송 면적·행렬·원천 계수, 물질 중력, native 기준 밀도·온도, 같은 공유 물질/광자 면, 읽어내는 지연·스칼라 가중치와 외부 진공이다. 외부 진공은 원 가스 표면이 아니라 초기 광자가 끝나는 실제 바깥 면에서 시작한다. 내부 EOS 표의 원 기준은 유지하고 검증된 밀도·주파수 일차 미분으로 새 기준을 연결한다. 내부 온도 변수는 새 native 등엔트로피 온도에 대한 변화다. 초기 가스 압력 구배는 원 물리적 구배와 새 압력 변화의 구배를 합한다. 새 중력에 맞춰 정수압 균형을 정의상 강제하지 않았다.

분류: Proven. 같은 국소 광자 에너지의 Killing 주파수는 a_new/a_old만큼 옮겨진다. 점유수 표현의 체적·적색편이 인자는 (V_old/V_new)*(a_new/a_old)^3이며, 이 변환은 셀 광자 개수를 보존한다. 두 주파수 노드 사이 양의 보간은 개수와 일차 에너지 모멘트를 보존한다. 이 대수는 유한 주파수 표의 물리적 연속 정확도를 증명하지 않는다.

분류: Counterexample candidate. 유한 주파수 상자의 끝에서 광자를 버리지 않고, 끝점 절단에 따른 극소 에너지 차이는 양의 스펙트럼 기울기로 보정했다. 모든 셀·각도에서 광자 개수/에너지 대조 오차는 {max(max(x['number'],x['energy']) for x in check['remap']):.6e}, 초기 복사 반경 압력 대조는 {check['initial_photon_pressure']:.6e}다. 끝점의 최대 상대 기울기는 {max(x['endpoint_tilt'] for x in check['remap']):.6e}다. 새 입력의 스펙트럼 재구성 가정이며 관측된 유일한 분포가 아니다.

분류: Counterexample candidate. 실제 설치된 물질과 단계120 초기 제약 원천의 최대 상대 차이는 밀도 {installed['errors']['density']:.6e}, 압력 {installed['errors']['pressure']:.6e}, 정지에너지 포함 에너지 {installed['errors']['energy_including_rest']:.6e}다. 초기 native 압력/내부에너지 대조는 {max(check['native_pressure'],check['native_energy']):.6e}, 실제 fine 끝점 대조는 {audit['native_endpoint_pressure_energy_relative']:.6e}로 원0.002 기준 안이다. 유한 native 점 검사는 전역 EOS·미분 인증이 아니다.

분류: Counterexample candidate. fine 에너지 수지는 {fine['energy_relative']:.6e}, 바리온 수지는 초기 대기량 대비 {fine['joint_baryon_relative_to_initial_atmosphere']:.6e}, 실제 작은 공유 포트 대비 {fine['joint_baryon_relative_to_actual_port']:.6e}다. 최대 내부 상대 밀도 변화는 {fine['maximum_density']:.6e}, 최대 속도/c는 {fine['maximum_frame_velocity']:.6e}다. 수송·반응과 두 SSP 단계의 공유 물질 유속을 실제 적용했으며, 초기 입력 진단만 완료한 상태가 아니다.

## 새 전하 판독

분류: Counterexample candidate. 각 셀의 적분 원천은 좌표 반경에 균일하게 나누지 않고 정규화한 Jordan proper 체적에 따라 분포시켰다. 구적마다 셀 내용물 합은 보존한다. 이 규칙을 저장된 단계119 원천에도 적용해 초기 배경 변경과 원천 판독 변경을 분리했다. 원 단계119 궤적은 재계산하지 않았다. 원 판독·조건부 GR 하한은 이전 기록으로 보존한다.

| 분류: Counterexample candidate — 원2percent 기준의 관심량 대조 | 상대 차이 |
|---|---:|
|64/128 직접 전하 |{controls['time_direct']:.8%}|
|64/128 직접+광자 질량 항 |{controls['time_total']:.8%}|
|선형/자연 cubic 저장 이력 |{controls['history_direct']:.8%}|
|17/9 저장 시점 |{controls['cadence_direct']:.8%}|
|8/4 반경 구적 |{controls['quadrature_direct']:.8%}|

분류: Counterexample candidate. 새 원천의 직접 정규화 전하 변화는 {result['endpoint_direct_relative']:+.12e}, 지정 외부 광자 질량 항을 합친 명목값은 {result['endpoint_direct_plus_photon_mass']:+.12e}다. 원 사전 readout 전체 판정은 passed={result['passed']}다. 이 값은 직접 원천과 선언한 광자 정규화의 조건부 결과이며, 새 배경에서의 추가 GR 원천·퍼텐셜 상계를 아직 포함하지 않는다.

분류: Counterexample candidate. 같은 proper 체적 판독의 기존 배경 직접 변화는 {comparison['old_same_volume_rule_endpoint']:+.12e}다. 새 초기 배경 적용 전후의 최대 직접 전하 차이는 새 진폭 대비 {comparison['projected_vs_old_same_rule']:.8%}, 기존 저장 원천의 체적 판독 규칙 변경만의 차이는 {comparison['old_volume_rule_change']:.8%}다. 초기값 변경에 따른 차이는 이번 시간 대조 차이보다 작다. 시간 대조 자체도 엄밀한 오차 상계는 아니다. 이 유한 모형 비교를 전체 원천·주파수·반경 오차의 상계로 바꾸지 않는다.

분류: Counterexample candidate. 작은 공유 물질 포트는 별도다. 누적 공유 질량의64/128 값은 {ports['join_mass']['coarse']:.9e}/{ports['join_mass']['fine']:.9e}g이며 fine 대비 차이는 {ports['join_mass']['relative_to_fine']:.4%}다. 관심 전하가 원 시간 기준을 통과했다고 작은 물질 포트도 수렴했다고 부르지 않는다.

분류: Counterexample candidate. 생산 전 재시작 경로를 검사하면서, history의 기준 압력을 현재 밀도에서 다시 평가하는 문제를 발견했다. 초기 두 단계 시험은 영 변위에서 시작해 영향이 없었다. 저장된 새 초기 압력을 읽도록 수정하고 원 시험 코드를 보존한 뒤 동일 접두 구간을 이어갔다. RHS·초기값·수락 기준은 바꾸지 않았다. 바인딩과 변경 사유는 restart-receipt.json에 둔다.

끝점 audit의 첫 실행은 NumPy bool의 JSON 직렬화에서 종료되어 수치 보고와 호출 수가 저장되지 않았다. 해당 코드와 실패를 보존하고 bool 변환만 수정해 남은15s 안에서 audit만 재실행했다. 재실행은90 native호출·{audit['seconds']:.2f}s였다. 원 실행의 실제 호출 수는 미기록이며 원160호출 한도와 재실행90호출을 구분한다. 실제 궤적과 전하 판독은 반복하지 않았다.

## 판정 경계와 다음 병목

분류: Counterexample candidate. corrected_initial_state_installed, actual_projected_material_photon_trajectories_completed, new_direct_charge_readout는true다. 관심량 readout의 수락 여부는 위 원 판정 그대로이며, old_GR_lower_bound_transferred, dynamic_metric_scalar_evolved, final_charge_solved, full_goal_complete는false다. 이 단계는 loophole progress다.

분류: Conjectural. 다음 병목은 실제 물질·광자 에너지/반경 압력의 변화를 질량·lapse·스칼라 진화에 되돌리는 일이다. 초기 비영 외재곡률을 저장한 것만으로 시간에 따라 Einstein 제약을 유지했다고 할 수 없다. 먼저 저장한 새 궤적을 재사용해 동적 항과 필요한 오차를 정하고, 기존 등방·정적 배경용 상계를 새 비등방 광자 배경에 그대로 이식하지 않아야 한다. 이후 실제 동적 결합에 반영해야 최종 전하를 판정할 수 있다. 전역 반경/주파수/각도, 모든 반응·중성수소 누적 원장과 관측 식별성도 여전히 별도다.

근거: outputs/direct-eos-gr33/def-native-projected-evolution의 plan/check/pilot/production/result/comparison/audit 및 source/wave/coupled 파일. 원 실패와 계산 예산은 보존한다. 추가 경로·해상도·기간을 자동 실행하지 않았다.
''',encoding='utf-8')

docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
paragraphs={
'model-definition.md':'분류: Counterexample candidate. 재고 일치 초기 GR/scalar 상태를 실제19+512셀 물질·광자 진화의 기하·체적·EOS 기준·스펙트럼·중력·공유 면에 설치했다. 원64/128 경로를 모두 완주했다. 새 계량은 진화 중 고정이므로 전체 동적 GR는 아니다.',
'observable-targets.md':f"분류: Counterexample candidate. 새 배경의 직접 전하 변화는{result['endpoint_direct_relative']:+.8e}, 지정 광자 질량 항과 합친 명목값은{result['endpoint_direct_plus_photon_mass']:+.8e}다. 원 readout passed={result['passed']}. 기존 양수 GR 하한은 승계하지 않는다. 최종 물리 전하·정적 EFT 이탈·관측 검출은 미확정이다.",
'adiabatic-limit.md':'분류: Proven. 광자 셀 개수는 체적·적색편이 변환 인자(V_old/V_new)*(a_new/a_old)^3로 보존된다.\n\n분류: Counterexample candidate. 새 초기값에서 원 짧은 비평형 응답을 실제 진화시켰지만 이 결과는 단열 한계와 자유 정적 계수의 기존 흡수 경계를 변경하지 않는다.',
'nonadiabatic-regime.md':'분류: Counterexample candidate. 새로운 순간 scalar 균형 초기값에서 실제 물질·광자 교환과 공유 유속의 비단열 진화를 완주했다. 고정된 새 계량에서의 원천 응답이며, 초기 Krr의 비영 값은 시간 진화에 아직 적용하지 않았다.',
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 재시작 history 기준 압력을 현재 밀도에서 재평가할 위험을 생산 전에 고쳐 설치된 초기 기준을 사용했다. 기존 좌표 균일 셀 판독 대신 proper 체적 원천 판독을 명시하고 저장 기존 원천에도 같은 규칙을 적용했다. 원 기록·문턱은 보존하며 이전 GR 하한을 새 배경에 재사용하지 않는다.',
'dynamic-charge-completion.md':f"분류: Counterexample candidate. new_initial_state_installed_in_transport, new_coupled_trajectory_completed, new_direct_charge_readout는true다. new_readout_accepted={result['passed']}. full_dynamic_GR_feedback, full_native_continuum_EOS, complete_reaction_ledger, final_charge_solved, full_goal_complete는false다. 초기 제약 입력에서 실제 결합 진화까지의 연결 병목은 닫았고, 다음은 새 원천을 동적 계량/scalar에 되돌리는 일이다."}
for p in docs:
    added='\n\n## 단계121 — 수정 초기값의 실제 결합 진화\n\n'+paragraphs[p.name]+'\n\n상세: [단계121 보고](../notes/REQUEST121_PROJECTED_COUPLED_EVOLUTION_KO.md).\n'
    p.write_bytes(p.read_bytes()+added.encode('utf-8'))
paths=docs+[report,root/'verification/def_native_projected_evolution.py',root/'verification/verify_native_projected_evolution.py']
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and not p.name.endswith('-checkpoint.npz'))
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-projected-evolution-manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='f9987fc1a',actual_projected_coupled_evolution=True,
    readout_accepted=result['passed'],old_GR_lower_bound_transferred=False,full_dynamic_GR_feedback=False,final_charge_solved=False,full_goal_complete=False,files=files))
master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_projected_evolution']=dict(classification='Counterexample candidate',report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)),
    new_initial_state_installed_in_transport=True,new_coupled_trajectory_completed=True,new_readout_accepted=result['passed'],
    full_dynamic_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,p in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:p['bytes']]).hexdigest()==p['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),readout_passed=result['passed'],passed=True)))
