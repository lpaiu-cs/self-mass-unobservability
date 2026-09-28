"""Preserve the failed final stage, its actual repair, and parallel GR dispatch."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,shutil,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase235-publish.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-parallel-gr-momentum'
manifest=out.parent/'native-parallel-gr-momentum-manifest.json'
note=root/'notes/REQUEST236_238_PARALLEL_GR_MOMENTUM_KO.md'


def package():
    assert not out.exists();out.mkdir();old=read(b.manifest)
    preserved={p:h for p,h in old['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    reused={};snapshots={}
    for number,name in [(235,'native-thermal-continuation235-work'),(237,'native-true-momentum237-work')]:
        base=runtime/name
        assert read(base/('coarse-receipt.json' if number==235 else 'check-receipt.json'))['source_sha256']
        for src in base.iterdir():
            if src.is_file() and src.suffix in ['.json','.log','.py']:copy(src,out/str(number)/src.name)
        for name in (['rejected-joint-stage.npz','last-accepted-64.npz'] if number==235 else ['comparison.npz']):
            p=base/name
            if name=='last-accepted-64.npz':reused['235/'+name]=dict(runtime_path=str(p),sha256=sha(p))
            else:copy(p,out/str(number)/name)
    for number,name in [(236,'native-complete-return236-work'),(238,'native-true-momentum238-work')]:
        base=runtime/name
        for name in ['plan.json','prepare-receipt.json','check-receipt.json','restart-check.json','source-endpoint-check.json','symbolic.json','controller-start.json']:
            if (base/name).exists():copy(base/name,out/str(number)/name)
        for name in ['controller-status.json','stage-progress-64.json','coarse-receipt.json','prefix-receipt.json','prefix-result.json','path-64.json','failure-64.json']:
            if (base/name).exists():
                value=read(base/name);write(out/str(number)/('snapshot-'+name),value);snapshots[f'{number}/{name}']=value
        state=read(base/'controller-status.json')
        if state['state'] in ['completed','failed']:
            for src in base.iterdir():
                if src.is_file() and src.suffix in ['.json','.log','.py']:copy(src,out/str(number)/src.name)
    fine=runtime/'native-independent-fine232-work'
    for name in ['stage-progress-128.json','controller-status.json']:write(out/'232'/('snapshot-'+name),read(fine/name))
    full=runtime/'native-true-momentum238-work';coarse=read(full/'path-64.json')
    assert coarse['passed'] and coarse['actual_completed_steps']==119
    for name in ['path-64.json','coarse-receipt.json','stage-progress-64.json']:
        copy(full/name,out/'238'/name)
    for p in full.glob('captured-64-*.npz'):copy(p,out/'238'/p.name)
    for p in (full/'sweep-1/photons').glob('complete-64.*'):copy(p,out/'238/sweep-1/photons'/p.name)
    for n in [236,238]:copy(root/f'.phase{n}-followthrough.py',out/f'controller-{n}.py')
    copy(Path(__file__),out/'publication-producer.py')
    failed=read(out/'235/coarse-receipt.json');r=read(out/'237/result.json')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 긴 구간 실제 GR 반환과 마지막 운동량 단계

분류: Counterexample candidate. **최종 전하의 기존 결론은 아직 판정 불가다.** 열 좌표 수정으로118단계가 실제 통과했지만,235의 마지막119단계는 원 비선형 잔차 기준에서 실패했다. 시간 상한 중단이 아니며118수락 상태와 모든 이력·실패 제안을 보존했다. 긴 공통 구간GR의 수락 결과는 실제 반환 계산236으로 연결했고, 운동량 산술 수정은 실제 마지막 단계238로 연결했다.

분류: Counterexample candidate.235의 벽시간은{failed['seconds']:.3f}초로2시간 상한 이내였다.119의8개 Newton 제안 중 마지막 잔차는1.15246014647e-11>1e-12였다. 저장된 전체 잔차를237에서 정확히 재현했다. 잔차의 대부분은 반경 운동량S였으며, 기존 native 값을 고정하고 잔차 합산만 고정밀로 바꾼 결과{r['accumulation_only_relative']:.12e}는 실패했다. 동일한 S 면 유속·중력 원천을 기존 고정밀 tangent와 보존량 역변환으로 평가하면{r['precise_native_S_relative'][-1]:.12e}로 바뀌었다.40/80자리 차이={r['precision_change']:.3e}, 원0.5/2 구성 대조 최대={max(r['constitutive_S']):.12e}다. 저장 제안은 여전히 실패하므로 수락하지 않았다.

분류: Counterexample candidate.238은118상태와 모든 물리·광자·각도·원천·경계·floor 이력의 재시작을 새 물리 단계 없이 정확히 확인했다. 수정한 실제 운동량 평가를 Newton 우변, 독립 비선형 잔차, 면 유속, 중력 및 운동량 수지에 함께 적용했다. 기존 제안은 안내값일 뿐이며 새 우변을 다시 푼다. 이전 열 좌표·바리온 수정과 정확한 선형 운동량 연산자는 유지했다. 원 기준을 낮추거나 앞의 수락 이력을 재적분하지 않았다. 게시 시각 {now}의238상태는{snapshots['238/controller-status.json']['state']}이며 동작은{snapshots['238/controller-status.json'].get('action')}이다. 실행 상태만으로 단계 통과를 선언하지 않는다.

분류: Counterexample candidate.224의 공통15/16물질·광자·경계 이력에서 원천·GR 시간 대조가 통과했으므로,236은 그 동일 원천을 원111/215실제 단계와 끝점의 합집합535시각에 평가한다. 원 끝점의 투영값과 실제 원천 미분을 보존한다. 독립 GR장 세 개가 모두 원 출력 재현을 통과한 뒤 동일 에너지·각도 출구로 metric을 만들고, 실제 물질·광자 결합 방정식의 거친·미세 경로 및 원2%대조까지 자동 연결한다. high/low 표현은 같은 방정식의 보상 성분이며 별도 해의 진단 전하를 더하지 않는다. 이는 한 번의GR반환으로 자기GR 고정점·완전 비선형·무한대 전하 완료가 아니다. 마지막1/16과 이전 초기 국소 시간 실패도 남는다.

계산 배치를 넓혔다.236의 독립 장 계산은CPU4/6/8에서 병행하고 각16GiB·3시간을 허용한다. 초기화 생성 파일도 각각 소유하게 하여 충돌을 피했다. 이어 metric1시간·실제 반환 거친3시간/미세6시간이며 기존 물리 시계와 수락 기준을 유지한다.238은CPU2·12GiB에서 마지막 원 단계3시간과 저장prefix수지1시간,232는CPU0·8GiB에서 기존3시간 속행이다. 중간 중단과 재설정 비용을 줄이는 여유 상한이며 완료시각 보장은 아니다. 기존224실측으로 병행 장의 가장 긴 작업은약32–65분, 실제 반환은 기존 구간당 비용이 유지되면약41/79분으로 추산했지만 후반 분기 비용은 미측정이다.

분류: Proven. 원Radau source/collision 묶음의 대수 항등식을 기호 검사했다. 이 항등식이나 독립 정밀도 일치는 균일 미분·EOS·전체 물리 오차의 보증이 아니다.

분류: Conjectural. 다음 수락은 마지막 단계, 동일 산술의 경로 교차 평가와 같은 해의 수지·시간 대조, 실제GR반환으로 결정한다. 공간·경계·EOS/미분·자기GR·완전 비선형·정적 비교·관측·무한대 최종 전하 조건은 축소하지 않는다. 이번 작업은 loophole progress이며 최종 전하 판정을 대신하지 않는다.
''',encoding='utf-8')
    with note.open('a',encoding='utf-8') as f:
        f.write('\n분류: Counterexample candidate.238의 실제 마지막119단계가 원 비선형 잔차9.167083482563306e-15, 물리 모멘트 최대1.464944057759037e-17로 통과했다. 거친64시계의 원 전체 기간3.434431ms 적분을119실제 단계로 완료했고, 새 단계 풀이를 포함한 coarse 작업은116.38초였다. 같은 이력의 각도 출구 차이1.32040e-16, 국소 물질 수지 최대1.33915e-13도 원 기준을 통과했다. 이는235의 실제 마지막 실패를 수정 산술의 실제 풀이로 해소한 결과다. 마지막 수락 순간의 광자·출구도 저장했다. 다만 수정 산술에서의 저장prefix 수지 재검증, 미세 경로와 공통 산술의 교차 평가 및 짝 시간 판정은 별도로 남아 있으므로 전체 물리 폐쇄나 최종 전하 완료가 아니다.\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계236–238 — 긴 구간 실제 GR 반환 및 마지막 운동량 단계\n\n분류: Counterexample candidate. 최종 전하 미판정.235의 마지막119단계 실패를 정확히 재현한 뒤 운동량 수정과 수지를 실제238풀이에 적용했다. 마지막 단계가 원 잔차9.16708e-15로 통과하여 거친 경로의 원 전체 기간3.434431ms를119단계로 완료했다. 수정 산술의prefix 수지·미세 경로 교차 평가·짝 시간 판정은 별도다. 공통15/16의 원천·GR 시간 통과는 원111/215단계 실제GR반환에 연결했고 독립 장 계산 세 개를 별도 코어에서 병행한다. 원 실패·기준과 자기GR·최종 전하 미해결 조건은 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        last_coarse_failure_preserved=True,momentum_controls=r,actual_coarse_period_completed=True,
        coarse_path_result=coarse,live_snapshots=snapshots,snapshot_KST=now,
        paired_time_admitted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    modules=[root/f'verification/{n}.py' for n in ['return_complete_history_gr','repair_true_momentum_stage','continue_true_momentum']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_parallel_GR_momentum']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
