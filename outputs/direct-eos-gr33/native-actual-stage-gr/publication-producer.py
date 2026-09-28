"""Freeze the paired GR verdict and actual native continuation evidence."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-actual-stage-gr';manifest=out.parent/'native-actual-stage-gr-manifest.json'
note=root/'notes/REQUEST223_232_ACTUAL_STAGE_GR_KO.md'


def copy(src,dst):
    h=sha(src);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==h==sha(src)


def package():
    assert not out.exists();old=read(out.parent/'native-dense-gr-return-manifest.json')
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    w=runtime/'native-stage-metric227-work';r=read(w/'result.json');assert read(w/'controller-status.json')['state'] in ['completed','failed']
    assert all(read(w/f'run-{n}.json')['passed'] for n in [64,128])
    out.mkdir();reused={}
    names={218:'native-rounding-continuation218-work',223:'native-captured-photon223-work',226:'native-interval-return226-work',227:'native-stage-metric227-work',228:'native-integer228-work',229:'native-integer-continuation229-work',230:'native-conserved-precision230-work'}
    for number,name in names.items():
        base=runtime/name
        for src in base.rglob('*'):
            if not src.is_file():continue
            rel=src.relative_to(base)
            # Preserve terminal records and unique numerical outputs. Large
            # already-saved input/checkpoint archives remain SHA-bound in place.
            if rel.parts[0] in ['sweep-0','sweep-1','clock-64','clock-128','checkpoints','snapshots']:continue
            if src.suffix=='.npz':
                reused[f'{number}/'+rel.as_posix()]=dict(runtime_path=str(src),sha256=sha(src));continue
            copy(src,out/str(number)/rel)
    for number,paths in {223:['recovered-64.npz','recovered-128.npz'],226:['sweep-1/photons/return-128.npz'],227:['recovered-64.npz','recovered-128.npz','sweep-1/photons/return-64.npz','sweep-1/photons/return-128.npz'],228:['proposal.npz'],229:['rejected-joint-stage.npz'],230:['comparison.npz']}.items():
        for path in paths:copy(runtime/names[number]/path,out/str(number)/path)
    # All GR inputs and outputs that were actually applied to the paired solve.
    for part in ['gr','metric']:
        for src in (w/part).glob('*.npz'):copy(src,out/'227'/part/src.name)
    live=runtime/'native-conserved-continuation231-work'
    for name in ['plan.json','restart-check.json','check-receipt.json','prepare-receipt.json','symbolic.json','controller-start.json']:
        copy(live/name,out/'231'/name)
    state=read(live/'controller-status.json');write(out/'231/status-snapshot.json',state)
    if state['state'] in ['completed','failed']:
        for src in live.glob('*.json'):copy(src,out/'231'/src.name)
        for src in live.glob('*.log'):copy(src,out/'231'/src.name)
        if (live/'rejected-joint-stage.npz').exists():copy(live/'rejected-joint-stage.npz',out/'231/rejected-joint-stage.npz')
    independent=runtime/'native-independent-fine232-work'
    for name in ['plan.json','restart-check.json','check-receipt.json','prepare-receipt.json','controller-start.json']:
        copy(independent/name,out/'232'/name)
    write(out/'232/status-snapshot.json',read(independent/'controller-status.json'))
    common=runtime/'native-complete-radau224-work'
    for name in ['plan.json','prepare-receipt.json','endpoints-receipt.json','source-receipt.json','endpoint-sources.json','sources.json','driver-polynomial-audit.json']:
        copy(common/name,out/'224'/name)
    write(out/'224/status-snapshot.json',read(runtime/'native-complete-radau224-control/status.json'))
    for number in [227,228,229,231,232]:
        source=root/(f'.phase{number}-run.py' if number==228 else f'.phase{number}-followthrough.py')
        copy(source,out/f'controller-{number}.py')
    copy(root/'.phase227-rate-alias.py',out/'rate-alias-producer.py');copy(root/'.phase230-probe.py',out/'defect-location-producer.py');copy(Path(__file__),out/'publication-producer.py')
    original=read(runtime/names[226]/'result.json');metric=read(w/'metric-result.json');integer=read(runtime/names[228]/'result.json');native=read(runtime/names[230]/'result.json')
    full=read(runtime/names[223]/'result.json');assert full['same_solution_accepted_history_recovered']
    note.write_text(f'''# 실제 단계 GR 반환과 전 기간 속행

분류: Counterexample candidate. **최종 전하의 기존 결론은 아직 판정 불가다.** 같은T/8 물질·광자·GR 반환의 두 실제 경로를 완료했다. 수정 후 원10채널 시간 대조 통과={r['passed']}, 최대 차이={max(r['time_relative']):.9%}. 자기GR 고정점·전체 기간·완전 비선형·무한대 최종 전하의 수락은 별도다.

분류: Counterexample candidate.226은 거친15/미세29실제 단계를 원 단계·물리·native·분기·보존 기준으로 완료했지만, 짝 시간 최대{max(original['time_relative']):.9%}>2%로 실패했다.188의 이전T/64 실패와226의 전역 실패를 모두 보존한다. 입력 미분을 원 Radau 시각에 합산한 대조에서 거친 시계의u 누적량5.05290%,lambda3.54655% 차이가 나타났고 미세 시계는2.36e-16이내였다. 이는 구체적인 입력 적분 오류이며 전체 응답 오차의 유일 원인 증명이 아니다.

분류: Counterexample candidate.227은 같은224밀집 원천과 지연 GR 전파를 기존 두 시계의 실제 단계·끝점73시각에 직접 계산했다. U_t를 전파하고 동일 원천 다항식의 미분을 중심 질량 제약에 넣어lambda_t를 구성했다. 원 끝점의 실제 바닥 투영은 보존했다. 최종 다항식의 투영 전 값과 저장된 투영 후 값은 서로 다르며 이를 매끄럽게 덮지 않았다. 외부 광자는 같은 두 Radau 각도 광도의 선형 collocation 표현으로 적분하여 같은 전체 에너지 측도를 유지했다. 입력의 최대 시간 차이는{max(metric['controls']['time'].values()):.9%}, 적분 차수 차이는{max(metric['controls']['quadrature'].values()):.9e}로 원2%/0.2% 기준을 통과했다. 필드의 기존 출력 재현, 같은 해의 에너지와 광선 대조도 통과했다.

분류: Counterexample candidate. 수정 입력을 실제 같은 물질·광자 단계 방정식에 적용했고, 두 경로 모두 실제 광자 모멘트와 원 경계 함수 반환값을 수락 순간에 저장했다. 다른 해의 전하를 사후 가산하지 않았다. T/8전역10채널은 통과했지만 최초T/64 끝점의 B/S/Eref/H 상대 차이는 각각0.402171%/0.121388%/2.421122%/18.803858%다. 특히 초기 Eref/H 차이는2%를 넘는다. 전역 성공으로 이 국소 차이를 지우지 않으며, 그 차이가 최종 전하 오차에 미치는 영향은 미판정이다. 입력 통과와 실제 단계 통과를 최종 전하 결론으로 대신하지 않는다.

분류: Counterexample candidate.223은 공통15/16구간3.219779ms의111/215단계 광자 이력 복원을 모두 끝냈고 원 끝점·물질 보존·native·반경 및 각도 출구 대조를 통과했다. 전체 주기의 마지막1/16은 여전히 제외된다. 이 실제 공통 이력은224전체 구간 소비자에 전달됐으며 원천 표현 검사를 완료하고 지연 GR을 계산 중이다. 이 실행의 스냅샷은 완료 판정이 아니다.

분류: Counterexample candidate.218은117단계 수락 뒤118선형 풀이에서7202.25초 예산으로 종료됐다. 저장된 정확한 우변·잔차를 재사용한228은 공동 정수 ULP 보정으로 선형 잔차{integer['linear_relative']:.9e}<1e-14를 통과했다. 이 제안을 실제 적용한229는8개의 선형 풀이를 통과했지만 실제 비선형 잔차8.97695e-10>1e-12로 실패했다. 새 물리 단계는 수락되지 않았고 원 실패를 유지한다.

분류: Counterexample candidate.230은229의 실제 전체 잔차를 정확히 재현했다. 고정밀 유속에 들어가기 직전의 보존량 변환을 동일 이진 입력으로 고정밀 계산하면 잔차가{native['arithmetic_effect']:.9e}만큼 바뀌었다.40/80자리 계산과 원0.5/2배 구성 대조는 통과했다. 기존 제안은{native['promoted_relative'][-1]:.9e}로 계속 실패하므로 이를 수락하지 않았다.231은 이 산술을 우변과 실제 비선형 평가에 함께 적용하여117상태와 모든 이력의 정확한 재시작 후 실제8번 풀이했지만1.21893e-6>1e-12로 다시 실패했다.263.89초이며 시간 예산 중단이 아니다. 보존량 변환 정밀도만으로118단계가 해결된다는 가설을 기각한다. 미세 경로와 저장prefix 재검사는 시작하지 않았고 새 물리 단계는 수락하지 않았다. 독립 실제 native 평가와 선형 분기/상태 표현 사이의 남은 오차를 해결해야 한다. 게시 시231 상태={state['state']}, 동작={state.get('action')}.

분류: Proven. 같은 Radau 선형 광도의 적분 가중치3/4,1/4와 보존 에너지 좌표 항등식을 기호 검사했다. 균일 미분 오차, EOS 인증, 물리 폐쇄의 일반 정리는 아니다.

분류: Counterexample candidate.232의 미세 경로 재시작은215수락 단계의 모든 저장 배열을 정확히 재현했다. 이미 승인된 원128시계의 마지막16단계만 독립 속행한다. 거친 경로 통과를 요구하던 실행 순서 조건만 제거했으며 어떤 거친 통과 기록도 만들지 않았다.231방정식과 경로별 원 수락 기준을 유지한다. 미세 경로가 끝나도 거친 실패·prefix 재검증·짝 시간 수락은 별도로 남는다.

계산 자원은 중간 중단의 재설정 비용을 줄이도록 여유 있게 유지했다.227은CPU1·12GiB, 각GR장30분·metric20분·실제 거친40분/미세80분이다.231은CPU1·8GiB, 거친2시간/미세3시간·저장 이력 재검사30분이다.232의 독립 미세 경로도CPU1·8GiB·3시간을 허용한다. 원 정확도·물리 기준과 실패는 유지하며 격자·물리 기간·경로 수를 자동 확대하지 않는다. namespace, 출력 비교 키, 원각도 배열의 prefix 범위 및 보고 normalizer 오류는 원 코드·계획·receipt와 함께 보존했다. 이미 계산된 GR 필드는 재사용했다.

분류: Conjectural. 남은 판단은 지배 오차를 수정한 동일 해의 전 기간·자기GR·경계/공간·EOS/미분·완전 비선형·정적 비교·관측 및 무한대 전하다. 이번 제한 구간의 수락만으로 기존 최종 전하의 부호나 결론을 계승하지 않는다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write((f'\n\n## 단계223–231 — 실제 단계 GR 및 마지막 결합 구간\n\n분류: Counterexample candidate. 최종 전하 미판정.226의 실제 두 경로는 완료했으나 시간 최대5.03484%로 실패했다. 동일 GR의 단계 시각 및 미분을 수정하여227의 실제 두 경로를 다시 완료했다. 원 시간 판정={r["passed"]}, 최대{max(r["time_relative"]):.9%}. 이전 실패·기준을 보존하며 자기GR 고정점이나 최종 전하로 확대하지 않는다. 공통15/16광자 이력은111/215단계 복원을 완료했다. 전 기간118단계는 선형 보정과 보존량 변환의 고정밀화를 각각 실제 풀이에 적용했으나 비선형 기준에서 계속 실패했다. 변환 정밀도만으로 해결됐다는 결론을 기각하고117수락 상태를 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,actual_stage_GR_return_paired_result=r,metric_result=metric,accepted_common_history_recovered=True,native_continuation_snapshot=state,original_failures_preserved=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    modules=[root/f'verification/{name}.py' for name in ['apply_actual_stage_gr','solve_saved_integer_native','continue_integer_native','repair_native_conserved_precision','continue_precise_conserved_native','continue_independent_fine_native']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_actual_stage_GR_return']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
