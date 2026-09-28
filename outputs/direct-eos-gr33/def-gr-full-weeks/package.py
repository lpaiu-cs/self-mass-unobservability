"""Bind the common propagation pass and the unchanged spatial failure."""
from pathlib import Path
import hashlib
import json
import subprocess

base=Path('outputs/direct-eos-gr33');out=base/'def-gr-full-weeks';b=out/'beta1024'
manifest=base/'gr-joint-propagation-manifest.json';assert not manifest.exists()
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
replay=read(out/'replay.json');spatial=read(b/'spatial-result.json')
assert replay['saved_replay_passed'] and replay['histories']==16 and replay['source_bindings']==64
assert not spatial['spatial_passed'] and all(read(b/f'p{p}-result.json')['propagation_passed'] for p in [1,2,4])
result=dict(classification='Counterexample candidate',checkpoint='45afcf302',
    decision='ACCEPT_COMMON_TIME_PROPAGATION_STOP_ON_NEW_BOUNDARY_SPATIAL_FAILURE',
    progress_class='loophole progress; common four-field propagation accepted on all original polynomial spaces, new-boundary spatial discrepancy remains',
    actual_same_input_applied=True,propagation_passed=True,propagation_passed_all_three_spaces=True,
    spatial_comparison_started=True,spatial_passed=False,original_failure_resolved=False,goal_complete=False,full_dynamic_charge_solved=False,
    time_comparisons=read(b/'p4-result.json')['comparisons'],spatial_comparisons=spatial['comparisons'],
    conditional_checks_not_run=['coefficient','outer3R','quadrature8','sigma14'],
    new_EOS_calls=0,next_bottleneck='Resolve the source-support interface spatial wave representation on the same coupled system; retain the accepted propagator and original thresholds. No automatic degree/grid expansion. Then complete coefficient, outer, quadrature and sigma contrasts.')
write(out/'result.json',result)
report=Path('notes/REQUEST79_GR_JOINT_PROPAGATION_KO.md')
details={
    'model-definition':'원 고차 결합 행렬을 축소하지 않고 확장 정밀도 잔차·FFT·Laguerre 전개로 실제 전파했다. 원 열 pole·기간·영 초기값·물질 debit·운동량 lift·양방향 결합을 유지했다.',
    'observable-targets':'같은 네 native RMS의 시간 전파 기준을 함께 통과했다. 작은 차이를 참해 오차 상계 또는 관측 가능한 신호로 부르지 않는다. 가장 거친512계수의 약64% 차이에서 얻은 큰 비율을 엄밀한 방법 차수로 해석하지 않는다.',
    'adiabatic-limit':'같은 물리 원천의 수치 전파 통과는 기존 단열 흡수 경계를 바꾸지 않는다. 공간·계수·외곽·구적·sigma 대조가 끝나지 않은 응답을 새로운 비단열 관측 신호로 승격하지 않는다.',
    'nonadiabatic-regime':'고정1·2·4차 공간 모두 네 시간 전파 기준을 통과했다. 마지막 공간 차이는 전체0.5096%,scalar0.0006105%,기존 절단면0.6588%,새 지지 경계6.0902%다. 새 경계만2% 기준에 미달했다.',
    'failure-ledger-dynamic-chi':'행 합 스케일링 기준해 실패·응답 기저의 기존 경계20.76% 실패·beta512 전체 속도 차수1.442 실패를 보존했다. beta1024는 공통 전파를 통과했으나 새 경계 공간 차이6.09%로 전체 수락은 실패했다. 직접 열 lift 뺄셈 차이는 그 공간 차이를 설명하지 못했다.',
    'dynamic-charge-completion':'네 성분의 공통 시간 전파 병목을 실제 같은 결합 경로에서 넘었지만 새 경계 공간 수렴과 후속 조건부 대조는 미완료다. 전체 비선형 항성·광자·전하·관측 목표로 확대하지 않는다.',
}
docs=[]
for stem,text in details.items():
    p=Path('docs')/(stem+'.md');old=p.read_bytes();assert b'Phase79' not in old
    addition=('\r\n## Phase79 — 공통 시간 전파 통과와 새 경계 공간 실패\r\n\r\n분류: Counterexample candidate. '+text+
        ' 전체 목표와 original_failure_resolved는 미완료로 유지한다.\r\n\r\n'
        '분류: Conjectural. 다음 수정 대상은 열 입력 종료 면 부근의 공간 파동 표현과 원천 점프 처리다. 통과한 전파를 유지하며 자동 격자 확대 없이 계획·비용·중단 기준을 재평가한다.\r\n\r\n'
        '분류: Imported from prior work. 실제 경로·원 실패·비용·독립 재생은 [단계79 보고서](../notes/REQUEST79_GR_JOINT_PROPAGATION_KO.md)에 둔다.\r\n')
    p.write_bytes(old+addition.encode('utf-8'));assert p.read_bytes().startswith(old);docs.append(p)
paths=[*docs,report,*[Path('verification')/name for name in ['def_gr_transfer_repair.py','def_gr_response_basis.py','def_gr_full_weeks.py']]]
for name in ['def-gr-transfer-repair','def-gr-response-basis','def-gr-full-weeks']:
    paths.extend(p for p in (base/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
paths=sorted(set(paths));hashes={p.as_posix():sha(p) for p in paths}
write(manifest,dict(classification='Counterexample candidate',checkpoint=result['checkpoint'],progress_class=result['progress_class'],
    decision=result['decision'],prior_milestone_manifest_sha256=sha(base/'gr-four-component-convergence-manifest.json'),
    actual_same_input_applied=True,propagation_passed=True,spatial_passed=False,original_failure_resolved=False,
    full_dynamic_charge_solved=False,goal_complete=False,sha256=hashes))
paper=Path('paper/revision-manifest.json');old=read(paper);head=json.loads(subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json']))
assert old==head and len(old)==187
key='request79_joint_gr_propagation';assert key not in old
new=dict(old);new[key]=dict(classification='Counterexample candidate',progress_class=result['progress_class'],passed=False,
    propagation_passed=True,spatial_passed=False,actual_same_input_applied=True,original_failure_resolved=False,
    full_dynamic_charge_solved=False,goal_complete=False,new_EOS_calls=0,next_bottleneck=result['next_bottleneck'],
    report=report.as_posix(),report_sha256=sha(report),evidence_manifest=manifest.as_posix(),evidence_manifest_sha256=sha(manifest))
write(paper,new);assert all(read(paper)[k]==v for k,v in head.items())
assert all(sha(Path(p))==h for p,h in hashes.items())
print(json.dumps(dict(files_bound=len(hashes),prior_paper_entries_preserved=len(head),paper_entries=len(new),document_prefixes_preserved=len(docs))))
