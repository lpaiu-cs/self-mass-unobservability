"""Publish the completed actual returned pair, preserving failures and scope."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase254-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-returned-right'
manifest=out.parent/'native-returned-right-manifest.json'
note=root/'notes/REQUEST256_FULL_RETURN_RESIDUAL_KO.md'


def package():
    actual=runtime/'native-returned-right256-work';old=runtime/'native-returned-krylov254-work'
    assert read(actual/'controller-status.json')['state']=='completed'
    result=read(actual/'result.json');assert result['passed']
    assert read(actual/'replay-228.json')['passed'] and read(actual/'actual-cache-reuse.json')['rhs_guess_exact']
    assert read(actual/'prefix-native.json')['passed']
    for p,h in read(actual/'plan.json')['bindings'].items():assert sha(runtime/p)==h,p
    assert sha(actual/'sweep-1/photons/return-64.npz')==sha(old/'sweep-1/photons/return-64.npz')
    assert sha(old/'recovered-128.npz')==read(actual/'private-output-copy.json')['original_hash']
    modules=[root/'verification'/n for n in ['finish_returned_right_residual.py','read_right_return_charge.py','read_right_return_exterior.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for p in actual.rglob('*'):
        if p.is_file() and p.suffix in ['.json','.log'] and 'sweep-0' not in p.parts and 'metric' not in p.parts and 'gr' not in p.parts:
            copy(p,out/'actual256'/p.relative_to(actual))
    for n in [64,128]:
        for name in [f'recovered-{n}.npz',f'sweep-1/photons/return-{n}.npz']:copy(actual/name,out/'actual256'/name)
    for name in ['accepted-linear.npz','initial-closure-failure/producer.py','initial-json-failure/producer.py']:copy(actual/name,out/'actual256'/name)
    for name in ['fine-receipt.json','fine.stderr.log','controller-status.json','linear-232.json','krylov-expansion.json','coarse-krylov-expansion.json','coarse-receipt.json','run-64.json','original-limit-232-30.npz','original-limit-232-30.json']:
        copy(old/name,out/'failed254'/name)
    for folder,label,status in [('native-right-charge256-work/full','charge256','controller-status.json'),('native-right-exterior256-work','exterior256','full-controller-status.json')]:
        p=runtime/folder;copy(p/status,out/label/status);copy(p/status.replace('status','start'),out/label/status.replace('status','start'))
    physical=runtime/'native-physical-port255-work'
    for p in physical.iterdir():
        if p.is_file():copy(p,out/'physical-port255'/p.name)
    for name in ['.phase255-port-input.py','.phase256-followthrough.py','.phase256-charge-followthrough.py','.phase256-exterior-followthrough.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    rows={str(n):read(actual/f'run-{n}.json') for n in [64,128]}
    linear=read(actual/'right-result.json');cost=read(actual/'right-calls.json')['calls']
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    note.write_text(f'''# 전체 기간 실제 반환과 같은 해의 전하 판독

분류: Counterexample candidate. **최종 물리 전하는 아직 미판정이다.** 이번 성과는 원 전체 기간의 실제 광자·물질·GR 반환 쌍119/231단계를 완료하고 원 시간 대조까지 통과한 것이다. 고정 초기 외부의 조건부 전하와 완전한 물리 전하는 구별한다. 게시={now}.

분류: Counterexample candidate.254의 coarse119는 전체 기간을 통과했다. fine는228까지 수락한 뒤229의 선형계가 벡터9.500747175713496e−11,최대 물리2.4558336450717364e−11로 실패했다. 원 제한1e−14/1e−13을 넘었으며 시간·메모리 중단이 아니다. 마지막 왼쪽 전처리 호출은 내부 callback0에서도 실제 보정 잔차1.91676을 남겼다. 실패 로그·거절 상태·기존 수락 기준은 그대로 보존한다.

분류: Proven. 오른쪽 전처리의 보정식 A P y=b−A x0, x=x0+P y는 원 선형계와 같은 잔차를 갖는다. 저장한 symbolic 검사는 이 항등식만 증명하며 수치 오차나 물리 결론을 인증하지 않는다.

분류: Counterexample candidate. 기존191의 오른쪽 전처리 구현을 재사용했다. 실패229의RHS·초기 추정·저장 잔차를 정확히 재현한 뒤 같은 연산자와 전처리로 풀었다. 최종 선형 잔차={linear['relative']:.9e},최대 물리 잔차={max(linear['moments']):.9e}; 첫 통과는10반복/1.524초였다. closure 전달 누락과numpy불리언의JSON직렬화 오류 두 건은 별도 보존했다. 첫 수치 통과 뒤 저장된 해를 다시 사용했고, 출력용recovered128은 독립 파일로 분리하여 원254바이트가 변하지 않았음을 검사했다.

분류: Counterexample candidate. 이 수정을 실제 결합 진화에 적용했다. coarse119와 fine215의 전체 저장 상태를 재사용하고, 기록되지 않았던215까지의 후반 native/branch진단은 같은 저장 단계에서32회 평가하여 원 native율과 차이0을 확인했다.228까지의 transient상태·전체 이력은 저장되지 않아216~228을 재생했고,228의 초기상태와 단계 쌍은 모두 정확히 일치했다.229에서는RHS와추정이 정확히 일치한 저장 수치해를 실제 단계에 투입하고 원 비선형·물리·분기·구성·보존 기준을 적용했다. coarse최대 비선형={rows['64']['maximum_true_stage']:.9e},fine={rows['128']['maximum_true_stage']:.9e}; 최대 물리 잔차는각각{rows['64']['maximum_true_physical_stage']:.9e},{rows['128']['maximum_true_physical_stage']:.9e}. 전체 기간 두 경로와 시간 감사가 통과했다.

분류: Conjectural. 실행 전 예산은 저장 선형계1시간,실제 fine6시간,16GiB,CPU3이었다. 직전 fine가2906.52초 걸린 측정에 근거해 마지막13단계 복원25~45분을 가정했으며 뒤의 오른쪽 전처리 비용은 미측정으로 적었다. 같은 물리 기간·격자·시계만 실행했고 정확도는 완화하지 않았다. 실제 종료 시간·메모리는 receipt에 있다. 완료된 high원천·계량·coarse·fine215를 재적분하지 않았다.

분류: Counterexample candidate. 같은256해의endpoint→dense기하·원천→3경로GR전하→고정 외부·질량 판독이 자동 연결돼 있다. 원251의229배열 재현과252의 독립 판독 감사 코드를 그대로 사용하며 입력 경로만 바꾼다. 이 게시 시점의 후속 상태는 별도 스냅샷이며 실행 중이라는 이유로 전하 결론을 추정하지 않는다.

분류: Proven. 같은 경계의 일차 물리 에너지 변환은Lphys=Lref+(u+nu)_face L0이다. 구동 계량이 시간 의존적이면 경계에서 변환한 에너지를 무한대까지 보존된Killing에너지로 취급할 수 없다. 전파 일·광선/각도·도착시간 및 질량·스칼라 접합이 별도 필요하다.

분류: Counterexample candidate.255에서 실제119/231Radau시계와 이번에 적용된 같은128-g8반환 계량으로 이 경계 입력을 구성했다. q8누적 변환 high64/128은−0.008695808385643418/−0.008695842724288512erg,low는−7.579858102076161e−28/−7.5798795698511635e−28erg이다. 시간·구적 입력 대조는 통과했지만 이 값을 옛 전하에 가산하지 않았다. 실제 계량에 적용한 물리 외부 연결도 아직 완료하지 않았다.

분류: Conjectural. 자기GR고정점,물리 외부 전파·에너지/질량 접합,배경 재정규화,EOS인증과 균일 미분·시간/공간/경계 오차,완전한 비선형·정적EFT비교·관측 연결은 남는다.188의 초기 국소 실패도 보존한다. 분류상 실제 적분 병목을 해소한loophole progress이며 전체 연구 완료나 최종 전하의 유지 판정은 아니다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계256 — 전체 기간 실제 반환 완료\n\n분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 원119/231실제 반환을 완료하고 원 시간 대조를 통과했다.254fine229의 선형 실패는 기존 오른쪽 전처리와 저장 해 재사용으로 실제 단계에서 해소했다. 완료coarse와fine215를 재사용하고228단계 재생의 정확한 일치를 확인했다. 같은 해의 전하와 고정 외부·질량 판독을 연결했다.255물리 경계 에너지 입력은 실제 시계로 준비했지만 전파 일·물리 외부/질량 접합 및 실제 계량 적용은 미완료이며,옛 진단을 전하에 가산하지 않았다. 자기GR·EOS/미분/공간/비선형/정적·관측 조건과 원 실패는 그대로 남는다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',passed=True,actual_same_joint_solution_GR_computed=True,actual_returned_steps=[119,231],full_period_actual_return_completed=True,
        original254fine_failure_preserved=True,coarse_reused=True,fine215_reused=True,actual228replay_exact=True,
        saved229linear_applied_to_actual_stages=True,same_solution_charge_followers_connected=True,scientific_gates_changed=False,
        physical_launch_input_prepared=True,physical_launch_input_applied_to_return=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,snapshot_KST=now)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_returned_right']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
