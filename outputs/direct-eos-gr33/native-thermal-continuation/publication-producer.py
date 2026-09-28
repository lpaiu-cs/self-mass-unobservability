"""Preserve the actual application of the thermal-coordinate precision fix."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,shutil,sys

helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,sha=b.read,b.sha
out=root/'outputs/direct-eos-gr33/native-thermal-continuation'
manifest=out.parent/'native-thermal-continuation-manifest.json'
note=root/'notes/REQUEST233_235_THERMAL_NATIVE_KO.md'


def write(p,value):
    p.parent.mkdir(parents=True,exist_ok=True);b.write(p,value)


def copy(src,dst):
    h=sha(src);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==h==sha(src)


def snapshot(base,number,names):
    for name in names:
        p=base/name
        if p.exists():write(out/str(number)/('snapshot-'+Path(name).name),read(p))


def package():
    assert not manifest.exists() and not note.exists();old=read(out.parent/'native-actual-stage-gr-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    out.mkdir(exist_ok=True);reused={}
    for number,name in [(233,'native-branch-consistency233-work'),(234,'native-thermal-precision234-work')]:
        base=runtime/name;assert read(base/'check-receipt.json')['error'] is None
        for src in base.iterdir():
            if src.is_file() and src.name not in ['normalization.json','photon-conservation-plan.json','check-result.json']:
                copy(src,out/str(number)/src.name)
    w=runtime/'native-thermal-continuation235-work'
    state=read(w/'controller-status.json');progress=read(w/'stage-progress-64.json')
    assert progress['solved_actual_steps']>=118 and read(w/'restart-check.json')['passed']
    for name in ['plan.json','prepare-receipt.json','check-receipt.json','restart-check.json','symbolic.json','controller-start.json','resource-affinity.json','controller-affinity.txt']:
        copy(w/name,out/'235'/name)
    snapshot(w,235,['controller-status.json','stage-progress-64.json'])
    if state['state'] in ['completed','failed']:
        for src in w.iterdir():
            if src.is_file() and src.suffix in ['.json','.log','.py']:copy(src,out/'235'/src.name)
        for name in ['rejected-joint-stage.npz','last-accepted-64.npz']:
            p=w/name
            if p.exists():reused['235/'+name]=dict(runtime_path=str(p),sha256=sha(p))
    # These two photons/ports belong to the accepted repaired step118.
    for i in [234,235]:copy(w/f'captured-64-{i:03d}.npz',out/'235'/f'captured-64-{i:03d}.npz')
    for number,name in [(232,'native-independent-fine232-work'),(224,'native-complete-radau224-work')]:
        base=runtime/name
        status=runtime/'native-complete-radau224-control/status.json' if number==224 else base/'controller-status.json'
        current=read(status);write(out/str(number)/'status-snapshot.json',current)
        snapshot(base,number,['stage-progress-128.json','result.json','fields-receipt.json','fine-receipt.json'])
        if current['state'] in ['completed','failed']:
            for src in base.iterdir():
                if src.is_file() and src.suffix in ['.json','.log','.py']:copy(src,out/str(number)/src.name)
            if number==224:
                for src in (base/'gr').iterdir():
                    if src.is_file():copy(src,out/'224/gr'/src.name)
    for name in ['.phase235-followthrough.py','.phase235-affinity.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    r=read(out/'234/result.json');original=read(out/'233/result.json')
    gr=read(out/'224/result.json');assert gr['GR_return_admitted']
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        native_thermal_arithmetic_repair=r,original_rejected_step=original,
        accepted_step118=True,coarse_progress_snapshot=progress,coarse_controller_snapshot=state,
        snapshot_KST=now,common_history_GR_result=gr,original_failures_preserved=True,paired_time_admitted=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    note.write_text(f'''# 내부 열 좌표 정밀도 수정의 실제 결합 진화 적용

분류: Counterexample candidate. **최종 전하의 기존 결론은 아직 판정 불가다.** 막혀 있던 실제118단계를 같은 결합 해에서 원 비선형·물리 모멘트 기준으로 통과했다. 이는 진단만의 개선이 아니라 실제 수락 이력의 전진이다. 전 기간·짝 시간 수렴·자기GR·최종 전하 수락과 구분한다.

분류: Counterexample candidate.233은231의 마지막 실패 제안 전체 비선형 잔차를 비트 단위로 재현했다. 선형 잔차{original['linear_relative']:.12e}는 통과했지만 실제 비선형 잔차{original['actual_relative']:.12e}는 실패했다. 실제 native 바리온 증가분과 affine 예측의 차이가 이를 설명하며 분해 나머지는{original['decomposition_remainder']:.12e}였다.40/80자리 재평가도 일치하여 자릿수 설정만 올리는 해결책은 지지되지 않았다.

분류: Counterexample candidate. 원인은 고정밀 native 평가 내부에서 보존량을 long double로 다시 낮춘 뒤 열 좌표의 상쇄 뺄셈을 하던 경로였다.234는 기존에 반올림된 고정 배경 계수 자체를 그대로 보존하고, 변하는 보존량의 곱셈·뺄셈을 고정밀로 유지했다. 같은 제안에서 native 증가분과 affine 예측의 차이는{r['affine_native_B_mismatch']:.12e}, 반/두 배 증분의 비선형성은{max(v['norm_stage_effect'] for v in r['increment_linearity']):.12e}로 줄었다.40/80자리 및 원0.5/2 구성식 대조도 통과했다. 그 제안은 새 우변에서 여전히 실패했으므로 물리 상태로 수락하지 않았다.

분류: Counterexample candidate.235는117수락 상태와 모든 이력의 정확한 재시작을 확인한 뒤, 수정 평가를 실제 Newton 우변과 독립 비선형 잔차에 함께 적용했다. 기존 제안은 초기 안내값으로만 사용하고 새 우변을 다시 풀었다. 실제118단계의 원 잔차는1.1053158577736875e-15, 물리 모멘트 최대1.2901845629864537e-17이었다. EOS·격자·원 시계·원 수락 기준은 바꾸지 않았다. 두 실제 Radau 순간의 광자 모멘트·출구 반환값도 수락 순간에 함께 저장했다.

분류: Counterexample candidate. 게시 시각 {now}, 거친 경로 수락 단계={progress['solved_actual_steps']}, controller={state['state']}. 실행 스냅샷은 전체 경로 완료 증거가 아니다. 별도로 진행하는232미세 경로는 고정된 이전 산술을 사용하므로, 끝나더라도 수정 산술의 명시적 교차 평가 없이235와 짝 수렴을 선언하지 않는다. 저장prefix 수지 대조 역시 이전 모든 벡터 방정식의 균일 오차 보증은 아니다.

계산 자원은 거친 속행2시간·저장prefix30분, 독립 미세 속행3시간, 공통 이력GR장2시간의 상한을 유지했다. 실제 거친·미세 작업이 모두CPU0에 묶였던 실행 배치를 확인하여235의 전 스레드와 controller를CPU2로 옮겼다.232는CPU0에서 계속한다. 각 적분은CPU1·8GiB이며 별도GR작업은12GiB다. 실행 중 소스·계획·수락 기준을 바꾸지 않았고, 중단이나 수락된 앞부분 재계산도 하지 않았다. 자원 배치 변경은 boot/PID/startticks와 함께 기록했다. 이 상한은 완료 예정시각이 아니다.

분류: Proven. 비선형 잔차와 affine 잔차의 차이가 같은 native 증가분 오차의 Radau 적분이라는 항등식을 기호 검사했다. 물리 수렴·EOS 인증·균일 미분 오차 정리가 아니다.

분류: Conjectural. 다음 결정은 마지막 실제 단계를 끝내고 수정 산술에서 두 경로·보존 이력이 양립하는지 확인한 뒤, 같은 해의 에너지·경계·GR를 최종 전하 판독에 연결하는 것이다. 앞서 남은 초기 국소 시간 차이, 공간·경계·EOS/미분·완전 비선형·정적 비교·관측·무한대 조건은 유지한다. 분류상 이번 진전은 loophole progress이며 최종 관측 가능성의 증명이 아니다.
''',encoding='utf-8')
    with note.open('a',encoding='utf-8') as f:
        f.write(f'\n분류: Counterexample candidate.224는 공통15/16기간3.219779ms의 같은 물질·광자·경계 이력에서 원천과 GR장을 완료했다. 원천 시간·표현 및 GR장 시간 기준을 통과했다. GR장 최대 시간 차이={max(gr["time"].values()):.12%},4/8적분 대조={gr["controls"]["quadrature"]:.12e}, 독립GR대조={gr["controls"]["independent_GR"]:.12e}다. 장 계산은2783.01초, 최대RSS6.08GB였다. 이 결과는 이 긴 구간의 실제GR반환을 시작할 수 있다는 판정이며, 반환 진화나 자기GR·최종 전하 완료가 아니다. 마지막1/16은 여전히 제외된다.\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계233–235 — 열 좌표 산술 수정의 실제118단계 적용\n\n분류: Counterexample candidate. 최종 전하 미판정. 이전 비선형 실패를 정확히 재현한 뒤, 고정밀 native 내부의 열 좌표 반올림을 제거하고 실제 같은 결합 해에 적용했다. 막혀 있던118단계가 원 비선형 잔차1.10532e-15와 물리 모멘트 기준을 통과했다. 원 실패·기준을 보존한다. 남은 단계와 다른 산술의 미세 경로 교차 평가, 같은 해의 전 기간·자기GR·최종 전하 수락은 별도다. [근거](../notes/'+note.name+').\n').encode())
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    modules=[root/f'verification/{name}.py' for name in ['align_native_branch_evaluation','repair_native_thermal_precision','continue_thermal_native']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_thermal_continuation']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
