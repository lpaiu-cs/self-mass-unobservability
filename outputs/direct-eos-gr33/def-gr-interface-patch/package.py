"""Close only the verified four-component convergence objective."""
from pathlib import Path
import hashlib
import json
import subprocess

base=Path('outputs/direct-eos-gr33');out=base/'def-gr-interface-patch'
manifest=base/'gr-interface-convergence-manifest.json';assert not manifest.exists()
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
replay=read(out/'replay.json');spatial=read(out/'spatial-result.json');contrast=read(out/'contrasts-result.json')
operator=read(out/'operator-check.json');time=read(out/'p4-result.json');pilot=read(out/'cached-pilot.json')
assert replay['saved_replay_passed'] and replay['all_four_component_gates_passed'] and replay['histories']==21
assert spatial['spatial_passed'] and contrast['all_original_contrasts_passed'] and operator['passed']
assert read(out/'control.json')['passed']
assert contrast['parent_accounted_compute_seconds']+replay['seconds']<1800
assert all(replay['propagation'].values()) and all(replay['contrasts'].values())
fields=list(time['comparisons']);names=['전체 속도','scalar','기존 절단면 속도','새 지지 경계 속도']
requirements={
    'same_input_horizon_and_readouts':['plan.json','pilot-budget.json','operator-check.json'],
    'all_four_time_gates_on_degrees1_2_4':['p1-result.json','p2-result.json','p4-result.json','replay.json'],
    'all_four_spatial_relative_and_decrease':['spatial-result.json','replay.json'],
    'original_coefficient_gate':['coefficient-contrast.json','replay.json'],
    'original_outer_gate':['outer-contrast.json','replay.json'],
    'original_quadrature_gate':['quadrature-contrast.json','replay.json'],
    'original_abscissa_gate':['abscissa-contrast.json','replay.json'],
    'linear_residual_heat_balance_initial_state':['replay.json'],
    'symbolic_source_identity':['control.json'],
    'unmodified_equation_cache_check':['cached-plan.json','cached-pilot.json'],
    'computation_budget':['cached-execution-plan.json','contrasts-result.json'],
}
write(out/'completion-audit.json',dict(classification='Counterexample candidate',
    objective='전체 속도, scalar, 기존 절단면 속도, 새 지지 경계 속도의 동시 수렴',
    scope='Original declared numerical gates on the frozen-source linear GR model, same0.23080495568542375s horizon and65 native readouts; no full nonlinear/EOS/observational closure.',
    every_declared_numerical_gate_verified=True,
    requirements={k:dict(passed=True,evidence=[(out/p).as_posix() for p in ps]) for k,ps in requirements.items()},
    original_failed_grid_preserved=(base/'def-gr-full-weeks/beta1024/spatial-result.json').as_posix(),
    original_thresholds_preserved=True,component_splicing=False,new_EOS_calls=0))
result=dict(classification='Counterexample candidate',checkpoint='959bc940f',
    decision='ACCEPT_ALL_FOUR_COMMON_PROPAGATION_SPATIAL_AND_ORIGINAL_CONTRAST_GATES',
    progress_class='loophole progress; the requested common four-component numerical convergence is resolved for the declared frozen-source GR model',
    actual_same_input_applied=True,propagation_passed=True,spatial_passed=True,all_original_contrasts_passed=True,
    original_failure_resolved=True,goal_complete=True,full_nonlinear_star_solved=False,full_dynamic_charge_solved=False,
    time_comparisons=time['comparisons'],spatial_comparisons=spatial['comparisons'],contrasts=contrast['cases'],
    new_EOS_calls=0,next_physical_boundary='Full outer heat/photon support, nonlinear temperature feedback, physical EOS uncertainty, dynamical charge and observational inference remain separate incomplete tasks.')
write(out/'result.json',result)
table='\n'.join(f'| {name} | {100*time["comparisons"][key]["last"]:.8g}% | {100*spatial["comparisons"][key]["previous"]:.8g}% | {100*spatial["comparisons"][key]["last"]:.8g}% | 통과 |' for name,key in zip(names,fields))
ctable='\n'.join(f'| {name} | '+ ' | '.join(f'{100*contrast["cases"][label]["relative"][key]:.8g}%' for label in ['coefficient','outer','quadrature','abscissa'])+' |' for name,key in zip(names,fields))
cost=contrast['parent_accounted_compute_seconds']+replay['seconds']
report=Path('notes/REQUEST80_GR_COMMON_CONVERGENCE_KO.md')
text=f'''# 단계80 — 네 성분의 공통 수렴 기준 통과

분류: Counterexample candidate. **요청한 전체 속도·scalar·기존 절단면 속도·새 지지 경계 속도의 공통 수렴 기준을 모두 통과했다.** 같은 실제 결합 GR 경로에서 원 시간·공간·계수·외곽·구적·역변환 위치 기준을 확인했다. 단계79에서 남았던 새 지지 경계의 공간 차이는6.09016%에서1.67109%로 줄었다. 통과한 성분을 서로 다른 해에서 합치지 않았다.

분류: Imported from prior work. 이 완료는 같은4,012개 면 열 pole, 영 초기 상태, 물질 에너지 debit·운동량 lift, 고정 GR 배경,0.23080495568542375초와65개 native 읽기의 선형 모형에 대한 수치 수락이다. 상대2%·관측 차수1.5, 공간2%·차이 감소, 계수2%, 외곽·구적0.2%, 역변환 위치0.02% 기준을 유지했다. 전체 EOS 인증·비선형 항성·광자·전하·관측 추론 완료를 뜻하지 않는다.

## 공간 병목의 실제 수정

분류: Counterexample candidate. 원 고정 격자의1·2·4차 공간 대조는 실패 상태로 보존했다. 저장된 파동 이동 거리는 새 경계에서 원 셀 약1–3개였으며, 속도 제곱차이는 열 입력 종료 면 주변에 집중됐다. 직접 열 lift 뺄셈 차이는 그 차이의약7e-13으로 작았다. 이를 근거로 별도 계획을 세워 원 새 경계 mask 전체를 포함하는 셀과 양옆 한 셀, 총25개만4분할했다. 이 후보가 실패하면 추가 분할하지 않는 조건으로 실행했다. 단순히 실패할 때마다 전체 격자를 확대하는 반복은 하지 않았다.

분류: Counterexample candidate. 수정 범위는0.51486456–0.53188198R, 추가 격자점75개다. 4차 공간의 자유도는47,417에서48,017로1.265% 늘었다. 기계적 trial 공간을 바꿨으며 새 EOS 상태나 열 셀을 계산하지 않았다. 원 격자점을 모두 보존했고, 열 입력 종료 면과 격자점의 거리는 계산상0이었다. 원 열 이력·기간·native 반경·질량 가중치·두 국소 mask가 동일함을 확인했다. 같은 수정 격자에서 원1·2·4차 공간 대조를 모두 실행했다.

분류: Proven. 일정한 원천 g에 대해 셀을 나눈 약한 원천 적분은 `sum integral(test_prime*g)=g*(test(right)-test(left))`를 유지하며 내부 면 항은 상쇄된다. 같은 차수의 다항식은 각 하위 셀에도 같은 차수로 제한되므로, 이 분할은 원 trial 공간을 포함한다. 기호 검산은 원천의 이 항등식에 관한 것이며 실제 GR 수렴을 대신하지 않는다.

분류: Counterexample candidate. 시간 전파는 단계79에서 수락한 전체 행렬의 확장 정밀도 잔차와 Weeks/Laguerre 표현을 유지했다. beta1024,sigma12,4,096개 contour 표본,512/1024/2048개 계수와 세 공간 차수를 늘리지 않았다. 다음 표의 마지막 시간 차이는4차 공간 결과이며,1차와2차 공간도 네 시간·contour 기준을 각각 통과했다.

| 성분 | 마지막 시간 전개 차이 | 공간1→2차 차이 | 공간2→4차 차이 | 판정 |
|---|---:|---:|---:|---|
{table}

분류: Counterexample candidate. 거친512계수의 이력은 여전히 충분히 해상되지 않아 약64%의 첫 차이를 보인다. 그로부터 얻은 큰 log2 비율을 엄밀한 방법 차수로 부르지 않는다. 마지막 계수 차이와 contour 대조, 후속 sigma14 대조를 함께 사용했다. 이는 선언된 수치 수락 기준의 통과이며 연속 시간·연속 공간 참해의 엄밀한 오차 상계는 아니다.

## 원 조건부 대조의 완료

분류: Counterexample candidate. 다음은 같은 최종4차 공간의 이력에 대한 최대 상대 차이다. 열 계수는 원 coarse-bank, 외곽은2R→3R, 구적은6→8점, 역변환 위치는sigma12→14의 사전 선언한 대조다. 각 대조 경로 자체의 네 시간·contour 기준도 통과했다.

| 성분 | 계수: 기준2% | 외곽: 기준0.2% | 구적: 기준0.2% | 역변환 위치: 기준0.02% |
|---|---:|---:|---:|---:|
{ctable}

분류: Imported from prior work. coarse-bank는 이전28-anchor 입력의 민감도 대조이며 독립적인 실제 물리 오차 포락선이 아니다. 이번 통과가 미계산 외층 열유속을 물리적으로0이라고 인증하지도 않는다.

분류: Counterexample candidate. 실제 수정 행렬의 `K+144M`에 대해 부동소수점 Cholesky 잔차{operator['scaled_factor_residual']:.3g}를 얻었다. native 끝점 재구성과 저장 이력의 네 상대 차이는최대{max(operator['endpoint_readout_relative'].values()):.3g}였다. 이는 유한 행렬의 검사이며 구간 SPD 인증이나 연속계 안정성 정리가 아니다.

## 비용·독립 재생과 완료 경계

분류: Counterexample candidate. 4차 공간 실제 경로는{time['seconds']:.2f}초,1·2차 두 경로는합{spatial['seconds']:.2f}초였다. 최초 네 조건부 경로의 예측1,495.92초가 남은1,285.38초를 넘어 본 계산을 시작하지 않았다. 프로파일링에서 반복되는 희소행렬 절댓값·대각 스케일링·배열 변환 비용을 확인했다. 일정 배열과 희소 구조를 저장하고 같은 원천 식·세 잔차 보정을 유지한 구현은32개 선언 전달 대조에서 네 성분의 차이가 계산상0이었다. 재측정 예측990.40초가 남은1,188.58초 안에 들어서만 실행했다. 이 최적화는 수락 기준·물리 입력·표본 수를 바꾸지 않았다.

분류: Counterexample candidate. 네 후속 대조 실제 합계는{contrast['seconds']:.2f}초다. 파일럿·검사와 실패한 프로파일 import의60초 상한 전체를 포함해 보수적으로 계상한 계산 시간은{cost:.2f}초로 원1,800초 예산 안이다. 원 예산 초과 예측과 import 실패를 보존했다. 한 CPU/BLAS 스레드,4GB 상한, 새 EOS 호출0을 유지했다. 이 시간은 대화·편집 대기까지 포함한 세션 총 경과시간이 아니다.

분류: Counterexample candidate. 독립 재생은21개 이력의65개 시각, 영 초기값, 네 native RMS, 열 수지, 원 시간·공간·네 대조 판정과{replay['source_bindings']}개 입력·실행 결속을 확인했다. 최대 읽기 재생 상대 차이는{replay['maximum_readout_relative']:.3g}다. `completion-audit.json`은 각 요구사항을 실제 근거 파일에 대응시킨다. 과거 실패는 고치거나 삭제하지 않았다.

분류: Conjectural. 다음 물리 연구는 현재 수렴한 고정 입력 선형 응답을 출발점으로 전체 외층 열·광자, 비선형 온도 되먹임, 물리 EOS 불확실성, 동적 전하와 관측 추론을 연결하는 것이다. 현재의 네 성분 수렴 목표와 이 미완료 물리 목표를 혼동하지 않는다.

분류: Imported from prior work. 이 단계는 **loophole progress: 선언 모형에서 요청한 네 성분의 공통 수렴 해결**로 분류한다. 코드·원 실패·계획·실측 이력·완료 감사는 `gr-interface-convergence-manifest.json`에 SHA256으로 묶는다.
'''
report.write_bytes(text.encode('utf-8'))
details={
    'model-definition':'원 새 경계 창의25개 셀만 고정4분할하여 같은 열 입력·기간·초기값·양방향 결합을 실제 GR 경로에 적용했다. 입력과 읽기를 바꾸지 않았고 새 EOS 호출은0이다.',
    'observable-targets':'전체 속도·scalar·기존 절단면·새 지지 경계의 네 native RMS가 한 경로에서 시간·공간 및 원 조건부 대조를 함께 통과했다. 이 수치 통과를 새 관측 신호나 실제 물리 오차 인증으로 확대하지 않는다.',
    'adiabatic-limit':'이번 수치 수렴 해결은 정적 흡수와 단열 no-go 경계를 변경하지 않는다. 선언된 비영 입력의 실제 선형 응답을 수락한 것이며 물리적 관측 식별성은 별도다.',
    'nonadiabatic-regime':'네 마지막 공간 차이는 전체0.4630%,scalar0.0006105%,기존 절단면0.6588%,새 경계1.6711%다. 모두 감소하며 원2% 기준과 모든 시간·계수·외곽·구적·sigma 기준을 통과했다.',
    'failure-ledger-dynamic-chi':'단계79의 새 경계6.09% 공간 실패를 보존하고, 별도 계획한 국소 보정의 실제 경로에서1.67%로 낮췄다. 조건부 계산의 최초 예산 초과는 실행 보류한 뒤 동일 연산의 캐시 동등성 검증으로 수정했다. 문턱을 완화하지 않았다.',
    'dynamic-charge-completion':'요청한 네 성분의 공통 수렴 목표를 선언된 고정 입력 선형 GR 모형에서 완료했다. 전체 비선형 항성·물리 EOS·광자·동적 전하·관측 폐쇄는 미완료이며 이번 완료 판정에 포함하지 않는다.',
}
docs=[]
for stem,value in details.items():
    p=Path('docs')/(stem+'.md');old=p.read_bytes();assert b'Phase80' not in old
    addition=('\r\n## Phase80 — 네 성분의 공통 수렴 기준 통과\r\n\r\n분류: Counterexample candidate. '+value+
        '\r\n\r\n분류: Imported from prior work. 원 기준·실제 결과·예산·완료 경계는 [단계80 보고서](../notes/REQUEST80_GR_COMMON_CONVERGENCE_KO.md)에 둔다.\r\n')
    p.write_bytes(old+addition.encode('utf-8'));assert p.read_bytes().startswith(old);docs.append(p)
paths=[*docs,report,Path('verification/def_gr_interface_patch.py'),Path('verification/def_gr_cached_resolvent.py')]
audit=read(out/'completion-audit.json')
audit['requirements']['maintained_notes_and_classification']=dict(passed=True,evidence=[p.as_posix() for p in [*docs,report]],progress_class=result['progress_class'])
write(out/'completion-audit.json',audit)
paths.extend(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
hashes={p.as_posix():sha(p) for p in sorted(set(paths))}
write(manifest,dict(classification='Counterexample candidate',checkpoint='959bc940f',progress_class=result['progress_class'],
    decision=result['decision'],prior_milestone_manifest_sha256=sha(base/'gr-joint-propagation-manifest.json'),
    propagation_passed=True,spatial_passed=True,all_original_contrasts_passed=True,original_failure_resolved=True,
    goal_complete=True,full_nonlinear_star_solved=False,full_dynamic_charge_solved=False,sha256=hashes))
paper=Path('paper/revision-manifest.json');old=read(paper);head=json.loads(subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json']))
assert old==head and len(old)==188;key='request80_common_gr_convergence';assert key not in old
new=dict(old);new[key]=dict(classification='Counterexample candidate',progress_class=result['progress_class'],passed=True,
    propagation_passed=True,spatial_passed=True,all_original_contrasts_passed=True,original_failure_resolved=True,goal_complete=True,
    full_nonlinear_star_solved=False,full_dynamic_charge_solved=False,new_EOS_calls=0,scope=read(out/'completion-audit.json')['scope'],
    report=report.as_posix(),report_sha256=sha(report),evidence_manifest=manifest.as_posix(),evidence_manifest_sha256=sha(manifest))
write(paper,new);assert all(read(paper)[k]==v for k,v in head.items())
assert all(sha(Path(p))==h for p,h in hashes.items())
print(json.dumps(dict(files_bound=len(hashes),prior_paper_entries_preserved=len(head),paper_entries=len(new),document_prefixes_preserved=len(docs),numerical_goal_complete=True)))
