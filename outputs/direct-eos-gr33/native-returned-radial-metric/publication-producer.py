"""Archive the matched radial returned metric and the registered live continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('previous',Path('.phase261-publish.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-returned-radial-metric'
manifest=out.parent/'native-returned-radial-metric-manifest.json'
note=root/'notes/REQUEST262_RETURNED_RADIAL_METRIC_KO.md'


def package():
    work=runtime/'native-returned-moments262-work';old=runtime/'native-returned-exterior262-work'
    check=read(work/'check.json');arithmetic=read(work/'arithmetic-check.json')
    assert check['passed'] and arithmetic['passed'] and read(work/'check-repair-receipt.json')['error'] is None
    assert read(old/'symbolic.json')['passed'] and read(work/'symbolic.json')['passed']
    bindings=read(work/'execution-plan.json')['bindings']
    for p,h in bindings.items():assert sha((root if p.startswith('outputs/') else runtime)/p)==h,p
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for folder,label,names in [
        (old,'original',['plan.json','prepare-receipt.json','symbolic.json','check.json','check-receipt.json','cost-stop.json']),
        (work,'continued',['plan.json','prepare-receipt.json','symbolic.json','arithmetic-check.json','check.json','check-receipt.json','check-repair-receipt.json','execution-plan.json','controller-start.json'])]:
        for name in names:copy(folder/name,out/label/name)
    write(out/'pipeline-snapshot.json',read(work/'pipeline-status.json'))
    for name in ['.phase262-check-repair.py','.phase262-production.py','.phase262-controller.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    final=dict(classification='Counterexample candidate',passed=True,
        verdict='RETURNED_RADIAL_METRIC_MATCHED_PHOTON_PROPAGATION_PENDING',
        actual_same_joint_solution_GR_computed=False,same_accepted_returned_boundary_reproduced=True,
        returned_metric_photon_propagation_complete=False,physical_GR_boundary_applied=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        snapshot_KST=now,controls=check,arithmetic=arithmetic)
    note.write_text(f'''# 반환 계량의 실제 외부 반경 연결

분류: Counterexample candidate. **기존 조건부 음의 전하는 보존했으며, 새로운 물리 외부를 포함한 최종 전하는 아직 미판정이다.** 수락한 같은 high/low 이력에서 실제 물질에 적용했던 반환 경계를 외부의 각 반경까지 연장했다. 반환 계량으로 배경 광자를 전파하는 실행도 연결했다. 이번의 수락 범위는 외부 반경 계량과 기존 경계·미분의 일치이며 광자 전 기간의 완료나 새 물질 결합 해의 판정은 아니다. 기록 시각={now}.

분류: Counterexample candidate. 같은 575개 저장 시각의 출구 U, 같은 signed Radau 방출, 같은 homogeneous 질량 잔차를 사용했다. 출구 U의 outgoing 연장은 기존 적용 경계의 가정과 동일하다. 새로운 외부 산란 스칼라장을 구한 것으로 해석하지 않는다. 원 방출 다항식과 각도 flux 가중치를 유지하고 광자 질량·압력 및 안쪽 질량 debit을 같은 반경의 제약에서 합쳤다. 기존 출구 lapse의 상대 차이는 {check['applied_boundary_relative']:.12e}, 독립 광행시간 역변환 차이는 {check['inverse_ray_age_fraction']:.12e}다. 시간·반경 미분의 독립 차분 대조 최대 차이는 {max(check['independent_metric_derivative_relative'].values()):.12e}였다. 이 대조는 유한 표본 검사이며 균일 미분 인증이 아니다.

분류: Proven. C(r,t)=C∞(t)−G/c⁴ Σw F(t−τ(r))를 쓰면 반환 질량 부분은 λ=C/(ra), ν=−C/(ra)−G/c⁴ Σ외부패킷 E(1+μ²)/(rₚaₚ)다. 움직이는 외부 적분 하한을 미분하면 ν′의 광자 압력 항이 나온다. signed 선형 방출과 다항 진공 커널의 모멘트 적분도 독립 Gauss 적분 및 기호 항등식으로 검산했다. 임의의 완전 비선형 배경에서 성립한다고 확대하지 않는다.

분류: Counterexample candidate. 정적 진공 커널의 값과 시간 미분은 독립 실제 ray와 대조했다. 입출력 시각·원천을 평활화하지 않고 정적 커널만 다항식으로 표현했다. 최초 대표 경로는 확장정밀도 거듭제곱 계산 비용이 지배적이었다. 모멘트 누적은 longdouble로 유지하고 별도로 보존된 작은 반환 성분의 조회만 float64로 평가하자, 같은256조회가 {arithmetic['old_seconds']:.6f}초에서 {arithmetic['new_seconds']:.6f}초로 줄었다. 같은 계량·미분 최대 상대 차이는 {max(arithmetic['relative'].values()):.12e}다. 큰 배경 차분이나 high/low 합산을 float64로 바꾼 것이 아니다. 원 대표 실행은 완료 cohort 없이 식별자를 확인하고 비용 사유로 종료했으며 계획·소스·중단 기록을 보존했다.

분류: Counterexample candidate. 빠른 조회 대조 뒤 같은 프로세스에서 legacy 초기화를 두 번 호출해 검사가 중단된 원 실패도 보존했다. 초기화가 생성 함수를 변형하는 구조 때문에 생긴 실행 오류였다. 독립 경계 검사는 새 프로세스에서 한 번만 초기화하도록 실행해 통과했다. 실제 광자 대표 구간 및 생산 경로도 프로세스당 한 번만 초기화한다. 이 오류를 물리 수락 실패로 해석하지 않는다.

분류: Conjectural. 현재 실행은 원 기간의 대표 방출 셀 세 개를 먼저 측정한다. 항등식 통과 및 경로당 보수적 비용 예상4시간 미만일 때만 네 기존 각도·방출 구적·기하 대조 경로를 최대4 CPU로 진행한다. 원 에너지 항등식0.2%, 각운동량1e−10, 구적0.2% 기준을 유지한다. 메모리 상한은 경로당16GiB다. 원 EOS·유체 해를 재계산하지 않는다. 대표 측정 전 전체 종료시각은 확정하지 않는다. 이 기록은 실행 예약을 완료된 물리 결과로 세지 않는다.

분류: Conjectural. 다음 판정에는 반환 계량 광자 응력의 전 기간 수락, 별도 반환 원천 시계 대조, 외부 생성·산란 스칼라와 reciprocal stress 및 배경 연산자 항, 물질 질량·일 접합이 필요하다. 이들을 합친 새 경계를 실제 결합 해에 적용하고 같은 해의 전하를 읽어야 한다. 과거 진단 전하나 이번 특수해를 기존 전하에 사후 가산하지 않는다. 이번은 누락돼 있던 반환 반경 입력을 실제 특성식에 연결한 loophole progress이며 최종 목표 완료가 아니다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as stream:stream.write(f'\n\n## 단계262 — 반환 계량의 외부 반경 연결\n\n분류: Counterexample candidate. 같은 수락 이력의 반환 질량·광자 압력·outgoing 스칼라 경계를 실제 외부 반경으로 연장했고 기존 출구 lapse와 독립 미분 대조를 통과했다. 이 계량을 배경 광자 전파에 연결했으며 대표 비용 검사 후에만 원 네 구적 경로를 실행한다. 최초 비용 중단과 중복 초기화 오류를 보존했다. 새 외부 광자 전 기간, reciprocal 스칼라·배경 연산자·접합·새 물질 결합 및 최종 물리 전하는 아직 미판정이다. 기존 조건부 음의 전하 수락은 유지한다. [근거](../notes/{note.name}).\n'.encode())
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    modules=[root/'verification'/n for n in ['propagate_returned_exterior.py','propagate_returned_moments.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,final);m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_returned_radial_metric']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
