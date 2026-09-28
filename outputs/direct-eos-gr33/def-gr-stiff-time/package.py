"""Bind phase78 results without changing prior verdicts or source bytes."""
from pathlib import Path
import json
import hashlib
import subprocess

root=Path('.')
base=root/'outputs/direct-eos-gr33'
manifest=base/'gr-four-component-convergence-manifest.json'
assert not manifest.exists()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
stiff=read(base/'def-gr-stiff-time/stage-result.json')
precision=read(base/'def-gr-gauss-refined/arithmetic-result.json')
allocated=read(base/'def-gr-allocated-projection/stage-result.json')
replay=read(base/'def-gr-stiff-time/replay.json')
assert replay['saved_replay_passed'] and replay['histories']==8
assert not any(r['propagation_passed'] for r in [stiff,precision,allocated])
result=dict(classification='Counterexample candidate',checkpoint='caf5c5b',
    decision='STOP_ALL_THREE_CANDIDATES_WITHOUT_COMMON_ACCEPTANCE',
    progress_class='loophole progress; weak-readout arithmetic separated from unresolved global and old-interface velocity propagation',
    actual_same_input_applied=True,propagation_passed=False,spatial_comparison_started=False,
    original_failure_resolved=False,full_dynamic_charge_solved=False,goal_complete=False,
    sdirk_comparisons=stiff['comparisons'],same_step_arithmetic=precision['comparisons'],allocated_comparisons=allocated['comparisons'],
    new_EOS_calls=0,next_bottleneck='Compare full resolvent, projected direct resolvent and projected eigendecomposition on the same source/readouts before another time evolution; fix the dominant transfer error without component splicing or automatic size expansion.')
write(base/'def-gr-stiff-time/result.json',result)
report=Path('notes/REQUEST78_FOUR_COMPONENT_CONVERGENCE_KO.md')
details={
    'model-definition':'동일 물리 생성자·입력·기간의 전체 행렬 SDIRK2 및 원 Gauss 단계 산술 수정과 총차원 재배분을 적용했다. 물리 감쇠·원천 재적합·일방향 결합을 추가하지 않았다.',
    'observable-targets':'네 native RMS 정의를 유지했다. 약한 읽기의 기존 단계 차이는 같은 단계 산술 변경에 민감했고 두 속도의 차이는 그 변경보다 훨씬 컸다. 작은 상대 차이를 참해 오차 상계로 부르지 않는다.',
    'adiabatic-limit':'같은 단계 정밀도 대조와 시간 전파 검증은 단열 흡수 경계를 바꾸지 않는다. 다른 수치법에서 각각 통과한 성분을 합쳐 관측 후보로 승격하지 않는다.',
    'nonadiabatic-regime':'SDIRK2에서 scalar·새 경계는 약2차였지만 전체·기존 절단면 속도는1.17/0.99차였다. 같은2048단계 Gauss의 산술 변화는 두 약한 읽기의 기존 차이 약1.3배, 두 속도 차이의0.003–0.005%였다.',
    'failure-ledger-dynamic-chi':'SDIRK2는 두 속도 차수에서, 총차원7:1 재배분은 전체 속도 차수 및 기존 절단면3.42% 차이에서 실패했다. 정밀도 Gauss는 예산상 같은2048단계 한 경로만 실행하여 수렴을 판정하지 않았다. 세 후보 모두 미수락이며 원 실패를 보존한다.',
    'dynamic-charge-completion':'네 성분 동시 전파와 후속1/2/4차 공간 대조는 미완료다. 원 방정식·EOS·광자·비선형·전하·관측 폐쇄와 구분하며 사용자 목표는 활성 상태로 남긴다.',
}
docs=[]
for stem,text in details.items():
    p=Path('docs')/(stem+'.md');old=p.read_bytes();assert b'Phase78' not in old
    addition=('\r\n## Phase78 — 네 성분의 산술·전파 병목 분리\r\n\r\n분류: Counterexample candidate. '+text+
        ' 네 성분의 공통 수렴은 미통과이며 공간·계수·외곽·구적 대조는 시작하지 않았다.\r\n\r\n'
        '분류: Conjectural. 다음 수정 대상은 같은 원천의 전체 행렬 해·투영 직접 해·투영 고유분해 해의 전달 오차로 분리한다. 단계 수·기저 수를 먼저 확대하지 않는다.\r\n\r\n'
        '분류: Imported from prior work. 실제 코드·예산·실패 판정·저장 재생은 [단계78 보고서](../notes/REQUEST78_FOUR_COMPONENT_CONVERGENCE_KO.md)에 둔다.\r\n')
    p.write_bytes(old+addition.encode('utf-8'));assert p.read_bytes().startswith(old);docs.append(p)
paths=[*docs,report,Path('verification/def_gr_stiff_time.py'),Path('verification/def_gr_gauss_refined.py'),Path('verification/def_gr_allocated_projection.py')]
for name in ['def-gr-stiff-time','def-gr-gauss-refined','def-gr-allocated-projection']:
    paths.extend(p for p in (base/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
paths=sorted(set(paths));hashes={p.as_posix():sha(p) for p in paths}
write(manifest,dict(classification='Counterexample candidate',checkpoint='caf5c5b',
    progress_class=result['progress_class'],decision=result['decision'],
    prior_milestone_manifest_sha256=sha(base/'gr-coupled-repair-milestone-manifest.json'),
    actual_same_input_applied=True,propagation_passed=False,spatial_comparison_started=False,
    original_failure_resolved=False,full_dynamic_charge_solved=False,goal_complete=False,sha256=hashes))
paper=Path('paper/revision-manifest.json');old=read(paper);head=json.loads(subprocess.check_output(['git','show','HEAD:paper/revision-manifest.json']))
assert old==head and len(old)==186
key='request78_four_component_convergence';assert key not in old
new=dict(old);new[key]=dict(classification='Counterexample candidate',progress_class=result['progress_class'],passed=False,
    propagation_passed=False,spatial_comparison_started=False,actual_same_input_applied=True,original_failure_resolved=False,
    full_dynamic_charge_solved=False,goal_complete=False,new_EOS_calls=0,next_bottleneck=result['next_bottleneck'],
    report=report.as_posix(),report_sha256=sha(report),evidence_manifest=manifest.as_posix(),evidence_manifest_sha256=sha(manifest))
write(paper,new);assert all(read(paper)[k]==v for k,v in head.items())
assert all(sha(Path(p))==h for p,h in hashes.items())
print(json.dumps(dict(files_bound=len(hashes),prior_paper_entries_preserved=len(head),paper_entries=len(new),document_prefixes_preserved=len(docs))))
