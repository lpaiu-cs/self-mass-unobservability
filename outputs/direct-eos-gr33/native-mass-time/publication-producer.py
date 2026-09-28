"""Preserve the actual knot-side repair, its causal test and pending verdict."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase257-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-mass-time'
manifest=out.parent/'native-mass-time-manifest.json'
note=root/'notes/REQUEST258_MASS_WORK_ENDPOINT_KO.md'


def package():
    work=runtime/'native-mass-time258-work';actual=runtime/'native-left-metric258-work'
    affine=read(work/'canonical-affine-result.json');pressure=read(work/'pressure-work.json')
    assert affine['passed'] and all(v['max_saved_native_relative']==0 for v in affine['rows'])
    assert read(actual/'endpoint-check.json')['passed']
    assert read(actual/'check-receipt.json')['error'] is None
    module=root/'verification/repair_returned_metric_endpoint.py'
    assert sha(module)==sha(runtime/'verification'/module.name)
    for p,h in read(actual/'plan.json')['bindings'].items():assert sha(runtime/p)==h,p
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for p in work.rglob('*'):
        if p.is_file() and p.suffix in ['.json','.npz','.py','.log'] and not any(x in p.relative_to(work).parts for x in ['sweep-0','sweep-1','metric','gr','__pycache__']):
            copy(p,out/'diagnosis'/p.relative_to(work))
    for name in ['phase258-energy-owner.json']:
        copy(runtime/name,out/name)
    for name in ['.phase258-locate.py','.phase258-energy-owner.py','.phase258-native-projection.py','.phase258-affine-owner.py','.phase258-pressure-work.py','.phase258-followthrough.py']:
        copy(root/name,out/name.lstrip('.'))
    for p in runtime.glob('phase258-*.log'):copy(p,out/'logs'/p.name)
    for name in ['plan.json','prepare-receipt.json','check-receipt.json','endpoint-check.json','metric-result.json','controller-start.json','controller-status.json','pipeline-status.json','charge-reader.py','exterior-reader.py']:
        copy(actual/name,out/'actual-start'/name)
    capture=actual/'capture-64.json';count=read(capture)['actual_steps'] if capture.exists() else 0
    if capture.exists():copy(capture,out/'actual-start'/capture.name)
    copy(Path(__file__),out/'publication-producer.py')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 같은 해의 질량 시간 오차와 구간 끝 미분 수정

분류: Counterexample candidate. **최종 물리 전하는 미판정이다.** 257의 수정된 전체 해에서는 compact 음의 부호가 유지됐지만, 고정 외부·동질 질량 시간 대조는 2.005740886%로 원 2%를 넘었다. 이번에는 지배적인 압력일 원천을 확인하고 구간 끝 미분 선택을 수정한 입력으로 실제 원 119/231단계 반환을 시작했다. 이 문서는 실행 중 스냅샷이며 새 해의 질량·전하 통과를 주장하지 않는다. 게시 시각={now}, 저장된 coarse 수락 단계={count}.

사용자는 06:05 KST WSL 부팅 변경이 직접 실행한 재시작임을 확인했다. 해당 중단은 수치 수락 실패와 구분하며, 당시 완료 결과 재사용과 중단 증거는 단계251 기록을 유지한다.

분류: Counterexample candidate. 같은 저장 해의 기체 에너지 끝점 차이를 native 수송·충돌·floor로 재구성했다. native 차이는 2.30662396710335e−19erg, 충돌 차이는 −5.53568203077140e−21erg였다. 세 내부 셀의 기여 합은 순 기체 차이의 약 99.2%이며, 상쇄가 있는 부호 있는 합이다. 이를 확률이나 독립 오차 상계로 읽지 않는다.

분류: Counterexample candidate. 저장 coarse/fine의 700개 native 수송률을 모두 정확히 재현했고 같은 실제 분기를 재생했다. 각 분기에서 고정 계량 원천과 상태 의존 항으로 나눈 결과, 순 에너지 차이는 각각 {affine['thermal_difference_erg'][1]:.12e} / {affine['thermal_difference_erg'][2]:.12e}erg였다. 분해 항등식 잔차={affine['identity_relative']:.12e}. 이 분해는 각 저장 분기에 조건부이며 전체 native 선형성이나 유일 원인 정리가 아니다.

분류: Counterexample candidate. 최초 교차 경로 투영과 분기 재생 진단은 원 native anchor 1e−12에서 멈췄다. 진단이 전 구간에 한 모델을 재사용한 것이 문제였으며, 원214/226과 같이 canonical 구간마다 새 모델을 사용하자 저장률 전부가 정확히 재현됐다. 교차 경로의 다른 방향 때문이라는 앞선 추측은 확정 원인으로 채택하지 않는다. 초기 15분 예측 제한이 1037.92초 예상치를 거절한 기록, 재시작 디렉터리 오류 및 실패 소스·영수증·로그도 보존했다. 자원 예산을 충분히 늘리되 정밀도 기준은 유지했고, 성공 재생은 {read(work/'canonical-affine-receipt.json')['seconds']:.2f}초였다.

분류: Proven. 연속 조각별 선형 함수가 경계에서 미분 1에서 3으로 바뀌는 예에서, 이전 구간 Radau 적분은 왼쪽 미분을 사용해야 1이다. 다음 구간 미분을 닫는 단계에 넣으면 1.5가 된다. 별도로 부분적분 항등식을 기호 검산했다. 이 두 항등식은 실제 항성·수치 오차 보장이 아니다.

분류: Counterexample candidate. 실제 적용된 질량 제약 미분은 source 다항식의 오른쪽 미분을 사용했다. 독립 재계산에서 원 적용값을 최대 1.85e−14로 재현했다. 같은 원천·장·계량 값은 유지하고 knot에서만 왼쪽 미분 차이를 적용한 세 입력을 작성했다. 입력 시간 대조 최댓값 0.0440340%, 구적 차이 2.911e−14로 기존 기준을 통과했다. 원 미분 자료와 원 실패는 그대로 보존한다.

분류: Counterexample candidate. 명시적 압력일의 원 순 시계 차이는 {pressure['coarse_minus_fine_erg'][0]:.12e}erg였다. 왼쪽 미분만 적용한 원천 비교는 native 차이를 {pressure['estimated_native_difference_after_left_derivative']:.12e}erg로 줄였다. 약 3분의 1 감소이며, 나머지 coarse 시간 구적 오차까지 해결됐다고 주장하지 않는다. 부분적분 진단은 더 작은 차이를 보였지만, 이 값을 전하에 더하거나 원시함수 기반 적분기를 적용하지 않았다. 원시함수 제거가 stiff 유한 단계 정확도를 보장하지 않는 기존 실패 경계도 유지한다.

분류: Conjectural. 실제 원 광자·물질 결합 방정식에는 수정된 delta_lambda_rate를 적용한다. 원 primary 이력·세 GR장·변하지 않은 계량 값은 재사용한다. 수정 입력은 첫 구간부터 다르므로 기존 반환 상태를 그대로 이어 붙일 수 없으며, 원 두 반환 경로만 다시 푼다. 새로운 격자·기간·시계는 추가하지 않는다. 같은 257 오른쪽 전처리와 강화된 물질 잔차 기준을 사용하며, 원 단계·native·보존·구성·시간·dense·전하·질량 기준을 유지한다. 허용 벽시간은 coarse 6시간/fine 8시간, 각 16GiB/CPU3다. 과거 실측에 따른 가정은 coarse 40~100분/fine 90~240분이며 변경 입력의 후반 속도는 미측정이다. 기존 판독 약 50~90분은 별도 가정이다. 실패 시 그 결과를 보존하고 자동으로 더 촘촘한 경로를 추가하지 않는다.

분류: Conjectural. 실제 두 경로가 통과하면 같은 해의 광자·물질·에너지·경계 이력으로 원 연속 원천, compact 전하, 고정 외부·동질 질량의 독립 감사까지 자동 실행한다. 최종 연구 가치는 그 동일 해에서 기존 결론이 유지되는지로 판정한다. 물리 경계 에너지·외부 전파 중 일·동적 광선·질량/스칼라 접합과 같은 해의 재적용, 자기GR, EOS 인증·균일 오차·비선형·정적 비교·관측 연결은 여전히 열린다. 이번은 실제 병목에 대한 loophole progress이며 전체 완료가 아니다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계258 — 질량 압력일의 구간 끝 미분 수정\n\n분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 실제 저장 수송률700개를 정확히 재현하여 시간 차이가 명시적 계량 압력일에 집중됨을 확인했다. 닫는 Radau단계에서 다음 원천 구간의 미분을 읽는 문제를 수정한 입력으로 원119/231실제 결합 반환을 시작했다. 원257질량2.005740886%실패와 모든 기준을 보존하며 새 해의 질량·전하 판정은 아직 없다. 다른 해의 진단값이나 부분적분 보정값을 전하에 가산하지 않는다. 같은 새 해의 전하·질량 감사까지 후속 실행을 연결했고 물리 외부·EOS/균일 오차·자기GR·비선형/정적/관측 범위도 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=False,verdict='ACTUAL_CORRECTED_ENDPOINT_EVOLUTION_PENDING',
        actual_same_joint_solution_GR_computed=True,prior257_compact_sign_survives=True,
        new_actual_return_completed=False,new_final_charge_read=False,original257_mass_time_failure_preserved=True,
        current_coarse_saved_steps=count,original_scientific_gates_preserved=True,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes))
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_mass_time']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
