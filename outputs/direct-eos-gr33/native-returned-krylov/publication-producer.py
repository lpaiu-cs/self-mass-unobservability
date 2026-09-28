"""Preserve the actual late-stage failure and its direct coupled continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase252-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-returned-krylov'
manifest=out.parent/'native-returned-krylov-manifest.json'
note=root/'notes/REQUEST254_ACTUAL_RETURN_KRYLOV_KO.md'


def package():
    current=runtime/'native-returned-krylov254-work';failed=runtime/'native-returned-refinement253-work'
    assert read(failed/'controller-status.json')['state']=='failed'
    assert read(failed/'replay-111.json')['passed'] and read(failed/'replay-112.json')['passed']
    assert read(current/'replay-111.json')['passed'] and read(current/'replay-112.json')['passed']
    modules=[root/'verification'/name for name in ['continue_returned_linear_refinement.py','read_refined_return_charge.py',
        'continue_returned_krylov.py','read_krylov_return_charge.py','read_krylov_return_exterior.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    plan=read(current/'plan.json')
    for p,h in plan['bindings'].items():assert sha(runtime/p)==h,p
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for p in failed.iterdir():
        if p.is_file() and (p.suffix in ['.json','.log'] or p.name.startswith('original-limit-')):copy(p,out/'failed253'/p.name)
    old=runtime/'native-full-return249-work'
    for name in ['controller-status.json','controller-start.json','coarse-receipt.json','coarse.stderr.log','metric-receipt.json','full/metric-result.json','full/capture-64.json']:
        copy(old/name,out/'original249'/name)
    copy(runtime/'native-refined-charge253-work/full/controller-status.json',out/'failed253-charge.json')
    copy(runtime/'native-full-charge251-work/full/controller-status.json',out/'failed251-charge.json')
    copy(runtime/'native-returned-exterior252-work/full-controller-status.json',out/'failed252-exterior.json')
    for folder,label in [(current,'current254'),(runtime/'native-krylov-charge254-work/full','charge254'),(runtime/'native-krylov-exterior254-work','exterior254')]:
        for p in folder.iterdir():
            if p.is_file() and p.suffix in ['.json','.log']:copy(p,out/label/p.name)
    for n in [64,128]:
        for name in [f'recovered-{n}.npz',f'sweep-1/photons/return-{n}.npz']:
            if (current/name).exists():copy(current/name,out/'current254'/name)
    for name in ['.phase253-followthrough.py','.phase253-charge-followthrough.py','.phase254-followthrough.py','.phase254-charge-followthrough.py','.phase254-exterior-followthrough.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(runtime/'phase254-live-verification.json',out/'live-verification.json')
    copy(runtime/'phase254-symbolic-check.json',out/'symbolic-check.json')
    copy(Path(__file__),out/'publication-producer.py')
    state=read(out/'current254/controller-status.json');counts={str(n):read(out/f'current254/capture-{n}.json')['actual_steps'] for n in [64,128] if (out/f'current254/capture-{n}.json').exists()}
    crossed=counts.get('64',0)>=113
    expanded=read(out/'current254/krylov-expansion.json')['events']
    completed=[v for v in expanded if v['state']=='completed']
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 실제 반환의 후반 선형 풀이 수정과 속행

분류: Counterexample candidate. **최종 물리 전하는 미판정이다.** 같은 긴15/16해에서252의 고정 외부·질량 판독까지 음의 조건부 부호가 유지된 결과는 보존한다. 전체 기간 반환은 실제113단계에서 실패했다. 이번에는 그 실패에 사용된 실제 결합 풀이를 직접 수정해 원119/231단계 끝까지 이어가는 중이다. 게시={now}. 게시 시 수락 단계={counts}, 실행 상태={state['state']}, 기존113단계 통과={crossed}. 이 스냅샷을 전체 경로의 완료나 최종 전하 판정으로 대체하지 않는다.

분류: Counterexample candidate.249의 전체 계량은 원 시간·구적·광선 및 과거 구간 접합 기준을 통과해 완료됐다. 그 계량을 적용한 coarse는111단계의 canonical15구간을 완료하고112단계 쌍을 저장했다. 실제113단계 첫 선형 풀이가 원4회 보정 제한에서 벡터1.881796669776e−6,광자 수8.591209781709e−6으로 실패했다. 벽시간·메모리 종료가 아니다. 원 상태·로그·거절 벡터를 보존한다.

분류: Counterexample candidate.253에서 방정식·전처리·Krylov공간20·각 호출5주기·수락 기준을 유지하고 잔차 보정 횟수를12까지 늘렸다. 기존111단계 저장 배열 전부와112단계 초기상태·단계 쌍이 정확히 재현됐다.113단계 첫 Newton선형 풀이는12회로 벡터2.202494133788e−19,최대 물리 잔차2.719491840941e−19까지 내려갔다. 그러나 같은 단계 다음 Newton선형 풀이가 벡터1.267402725341e−5로 실패했다. 후반 보정의 자체 잔차 비율이0.92~0.99에 머무는 호출들이 있어, 같은 작은 Krylov공간의 반복 횟수만 늘리는 것으로 병목이 해소됐다고 할 수 없다.

분류: Conjectural.254에서는 원GMRES(restart20,maxiter5)를 먼저 실행한다. 수렴 실패(info>0)한 호출만 같은 미지수·연산자·RHS·전처리와 허용오차에서restart80,maxiter20으로 이어간다.12회 extended-residual보정과3회 Newton제안,선형1e−14/물리1e−13,실제 비선형1e−12/물리1e−13,constitutive0.2%,보존1e−8,port1e−12,시간2%는 그대로다. 큰 탐색 공간의 효과는 원 실제 단계 수락으로 판정하며 격자·물리 기간을 늘리지 않는다.

분류: Counterexample candidate.254에서도 원111/112재현은 정확히 통과했다. 전체249계량과248GR,전체 주 해239를 다시 계산하지 않았다. 저장 누락된 transient branch/anchor진단을 복원하기 위한 coarse의canonical15구간만 재연산한다. 높은/낮은 성분과 동일 해의 광자·물질·실제 출구·적용 계량을 유지한다. 새254전하 후속은 이미229배열을 정확히 재현한251판독 연산자에서 입력 경로만 바꿨으며,254실제 물리 쌍과 감사가 통과해야 출발한다. 같은 해의 외부·질량 판독 및 독립 감사도 그 뒤에 연결했다. 이전251/252/253후속의 의존성 실패 기록은 남는다.

분류: Counterexample candidate. 수정은 진단에만 적용한 것이 아니다. 실제 원 실패 단계113을 원 비선형·물리 잔차 및 constitutive·보존 검사까지 통과해 다음 물리 단계로 진행했다. 게시 시 끝난 큰 공간의 풀이 수={len(completed)}, 첫 추가 풀이={completed[0]['seconds']:.3f}초/{completed[0]['iterations']}회다. 이 중간 수락을 전체 coarse/fine 시간 수렴이나 최종 전하의 유지로 바꾸어 주장하지 않는다.

실행 예산은 coarse2시간·fine4시간·각16GiB·CPU3을 유지했다. 원20공간100반복은 실측10~12초다.80공간1600반복의3~5분 추정은 가정이며, 후반 수렴 속도는 아직 미측정이다. 전하 판독은 물리 재적분 없이 endpoint90분·기하60분·원천120분·GR세 경로 각3시간과16GiB,외부30분이다. 시간 상한의 여유는 정확도 완화가 아니다. 통과한 물리 구간을 이유 없이 반복하거나 자동 새 해상도를 추가하지 않는다.

분류: Conjectural. 물리 에너지 변환·시간 의존 외부 전파·배경 질량 재정규화·자기GR고정점과 전체 EOS/균일 미분·시간/공간/경계/비선형/정적 EFT·관측 연결은 미완료다.252의 고정 외부 부호와253의 선형 한 개 통과는 이 조건들을 대신하지 않는다.188의 초기 국소 실패와 모든 원 수락 기준을 보존한다. 분류상 실제 병목에 대한 loophole progress이며 전체 목표는 진행 중이다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계254 — 실제 반환의 후반 선형 풀이 속행\n\n분류: Counterexample candidate. 최종 물리 전하 미판정.249의 완료 계량을 실제 반환에 적용했으나113단계의 선형 풀이가 실패했다.253은12회 보정으로 첫 선형계를 통과했으나 다음 계에서 잔차가 정체됐다.254는 실패한 호출만 더 큰 Krylov공간으로 이어 풀며 원 물리 방정식·수락 기준을 유지한다.111/112저장 상태 재현은 정확히 통과했다. 전체 기간 쌍이 통과하면 동일 해의 전하 및 고정 외부·질량 판독을 자동 수행하도록 연결했다. 전체 물리 전하·자기GR와EOS/미분/공간/비선형/정적·관측 조건 및 원 실패는 남는다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,
        original_249_linear_failure_preserved=True,twelve_solve_253_failure_preserved=True,
        every_accepted_111_112_replay_exact=True,original113stage_passed=crossed,actual_step_snapshot=counts,
        same_full_solution_charge_and_exterior_followers_connected=True,scientific_gates_changed=False,
        live_snapshot=state,snapshot_KST=now,full_period_actual_return_completed=state['state']=='completed',
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=plan['bindings']))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_returned_krylov']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
