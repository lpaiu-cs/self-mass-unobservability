"""Record a checked full-history consumer without promoting its waiting state."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase243-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-full-captured-source'
manifest=out.parent/'native-full-captured-source-manifest.json'
note=root/'notes/REQUEST244_FULL_CAPTURED_SOURCE_KO.md'


def package():
    w=runtime/'native-full-captured244-work';assert read(w/'regression.json')['passed']
    assert not out.exists();out.mkdir();preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    for name in ['regression.json','check-receipt.json','initial-check-receipt.json','initial-check-producer.py','controller-start.json']:
        copy(w/name,out/name)
    snapshots={}
    for number,folder in [(236,'native-complete-return236-work'),(239,'native-common-arithmetic239-work'),(244,'native-full-captured244-work')]:
        for name in ['controller-status.json','stage-progress-128.json','capture-64.json','capture-128.json','result.json']:
            p=runtime/folder/name
            if p.exists():
                value=read(p);snapshots[f'{number}/{name}']=value;write(out/str(number)/('snapshot-'+name),value)
    for name in ['.phase244-followthrough.py','.phase244-launch.ps1']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    fine=snapshots['239/stage-progress-128.json'];ret=snapshots.get('236/capture-64.json',{})
    note.write_text(f'''# 전체 기간의 실제 캡처를 원천 계산으로 연결

분류: Counterexample candidate. **최종 전하의 기존 결론은 아직 미판정이다.** {now} 기준 수정 미세 경로는 {fine['solved_actual_steps']}/231단계, 실제 GR 반환의 거친 경로는 {ret.get('actual_steps')}/111단계다. 새244는239의 실제 완료와 원 짝 시간 수락을 기다리는 소비자이며, 대기 상태를 전체 기간 수락이나 최종 전하 성과로 세지 않는다.

분류: Counterexample candidate.243의 완성된119단계 거친 광자 이력과223의215단계 미세 광자 이력을 재사용한다. 미세의 남은16단계는239가 실제 수락 순간에 저장하는32개 캡처를 붙인다. 앞의 물질 단계·광자 조건부 풀이를 반복하지 않는다. 앞215단계 보존 상태/원천/시간의 정확 일치, 새 캡처의 시간·가중치·충돌률 정확 일치, 전체 물질 수지1e-8 및 각도/반경 출구1e-12를 요구한다.

분류: Counterexample candidate. 새 조립 함수는 이미 완료된243의117단계 조건부 prefix와 실제118/119캡처로 모든 출력 배열을 비트 단위로 재현했다. 시각을 한 표현값만 옮긴 입력은 거부했다. 첫 검사는 조건부로 재계산된 과거 충돌률까지 저장 물질률과 비트 일치를 요구하여 실패했다. 이는 기존243의 원 전체 방정식/수지 기준에 없던 잘못 추가한 조건이다. 실패 소스·receipt를 보존하고, 과거 배열 자체 보존과 원 전체 수지, 새 실제 캡처의 정확 일치로 원 기준을 유지했다. 물리 허용오차는 바꾸지 않았다.

분류: Proven. 기존 소비자의 원 Radau 모멘트 및 같은 affine 배경 분할 항등식을 재사용한다. 분류: Counterexample candidate. 기존 원천 다항식의 저장 예제 회귀는 동일 결과를 냈다. 이 대조는 전체 물리 오차나 실행 대기 중인 미세 경로의 수락 증거가 아니다.

분류: Conjectural.239가 통과하면 조립 → 끝점 원천 → 실제 Radau/구동 다항식 원천을 자동 실행한다.224의 기존 소비자를 그대로 재사용하여 마지막1/16을 포함한119/231단계의 같은 에너지·압력·광자·경계를 읽는다. 원천 시간2%, 압력0.2%, mapping/dense/polynomial1e-12, 수지1e-8을 유지한다. 실제 반환236의 기하를244의 일차 구동 모드로 대체하지 않으며,236반환 해의 최종 전하는 별도로 같은 해에서 판독해야 한다.

CPU2·16GiB를 배정했다.224실측은 끝점185.92초, 원천570.36초였다. 단계 수는326에서350으로7.4% 늘며 단위 비용이 유지될 때 끝점4–8분, 원천11–20분으로 계획하되 후반 비용은 미측정이다. 불필요한 중단을 줄이도록 각각45/90분을 허용한다. 기다림은6시간 상한이며 의존 실행의PID·시작 시각·WSL부팅ID를 확인한다. 실패하면 자동 장 계산이나 격자 확대를 하지 않는다. 기존236/239의 소스·계획·실행은 변경하지 않았다.

분류: Conjectural. 연구 가치 기준은 지배 오차를 해결한 동일 결합 해의 최종 전하 결론이다. 이번 연결은 그 판독을 위한 작업이며 EOS/미분·공간/경계·자기GR 고정점·완전 비선형·정적 비교·관측·무한대 정규화의 미해결 범위를 줄여 선언하지 않는다. 이번 작업은 loophole progress다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계244 — 실제 캡처로 전체 기간 원천 연결\n\n분류: Counterexample candidate. 최종 전하 미판정.239의 실제 미세 완료와 원 짝 수락 이후, 기존 거친 광자 이력·미세215단계와 실제32개 새 캡처를 재사용하는 소비자를 연결했다. 모든 배열이 같은 기존 거친 완료를 재현하고 시각 불일치를 거부하는 대조를 통과했다. 대기 상태는 실제 전체 조립이나 물리 수락이 아니다. 원 기준·이전 실패·236실제GR반환과 전체 완료 조건을 유지한다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=False,
        checked_existing_history_assembly=True,full_period_source_started=False,live_snapshots=snapshots,snapshot_KST=now,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused={}))
    module=root/'verification/read_full_captured_history.py';assert sha(module)==sha(runtime/'verification'/module.name)
    files=[module,note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_full_captured_source']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
