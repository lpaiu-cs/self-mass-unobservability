"""Preserve exact characteristic transformations and measured continuation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('prior',Path('.phase262-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-primitive-photon'
manifest=out.parent/'native-primitive-photon-manifest.json'
note=root/'notes/REQUEST263_PRIMITIVE_PHOTON_CONTINUATION_KO.md'


def package():
    work=runtime/'native-primitive-photon263-work';comparison=read(work/'method-comparison.json');pilot=read(work/'pilot.json')
    check=read(work/'check.json');plan=read(work/'execution-plan.json')
    assert comparison['passed'] and check['passed'] and read(work/'pilot-receipt.json')['error'] is None
    assert not pilot['cost_eligible'] and pilot['forecast_upper_seconds_per_full_path']<plan['maximum_path_seconds']
    bindings=plan['bindings']
    for p,h in bindings.items():assert sha((root if p.startswith('outputs/') else runtime)/p)==h,p
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for name,label in [('native-smooth-photon263-work','mass-transform'),('native-outgoing-photon263-work','outgoing-transform'),('native-primitive-photon263-work','primitive')]:
        folder=runtime/name
        assert read(folder/'pilot-receipt.json')['error'] is None and read(folder/'symbolic.json')['passed']
        for f in folder.iterdir():
            if f.is_file() and f.suffix in ['.json','.npz'] and f.name not in ['pipeline-status.json']:
                copy(f,out/label/f.name)
    copy(work/'direct-reference/plan.json',out/'direct-reference/plan.json')
    for name in ['.phase263-direct-reference.py','.phase263-controller.py','.phase263-production.py','.phase262-production.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    write(out/'pipeline-snapshot.json',read(work/'pipeline-status.json'))
    old=runtime/'native-returned-moments262-work'
    for name in ['pipeline-status.json','pilot-receipt.json']:
        if (old/name).exists():copy(old/name,out/'original262'/name)
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    maximum=max(v for group in comparison['relative'].values() for row in group.values() for v in row.values())
    stress=max(row['stress_weak_moments_erg'] for group in comparison['relative'].values() for row in group.values())
    final=dict(classification='Counterexample candidate',passed=True,
        verdict='EXACT_PRIMITIVE_CHARACTERISTICS_AND_ACTUAL_COHORTS_ACCEPTED_DIRECT_REFERENCE_PENDING',
        actual_same_joint_solution_GR_computed=False,actual_returned_metric_photon_cohorts_computed=3,
        unchanged_ODE_tolerances=True,independent_direct_trajectory_acceptance_pending=True,
        whole_period_production_complete=False,physical_GR_boundary_applied=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        snapshot_KST=now,actual_transformed_rhs_relative=check['independent_actual_transformed_rhs_relative'],
        method_comparison_maximum=maximum,stress_comparison_maximum=stress,
        forecast_upper_seconds_per_path=pilot['forecast_upper_seconds_per_full_path'],maximum_path_seconds=plan['maximum_path_seconds'])
    note.write_text(f'''# 동일 광자 특성식의 누적량 적분과 실제 대표 전파

분류: Counterexample candidate. **기존 조건부 음의 전하 결론은 보존했고, 새 물리 경계에서의 최종 전하는 미판정이다.** 반환 계량에서 실제 광자를 전파하는 비용 병목을 수정했다. 같은 배경·반환 계량·방출 이력을 보존하면서 세 가지 정확한 변수변환으로 각각 대표 방출 구간0·7·15를 실제 적분했다. 누적량 방식은 세 구간을 모두 완료했으며 원 시간 미분식으로 적분하는 독립 두 광선 대조 뒤 전체 원 기간의 네 구적 경로를 실행하도록 연결했다. 기록={now}.

분류: Proven. 원 특성식은 ḣ=−vν′−μ²λ_t, δṫ=−ζ−(1−μ²)h/μ²다. k=h+μ²λ로 옮기면 λ_t가 소거된다. 같은 outgoing U(t−x/c)에 대해서는 그 시간 미분도 끝점 항으로 옮길 수 있다. 더 나아가 원 방출 primitive F를 쓰면 광선 교차의 dF/dt=(1−μ/μ_source)L과 압력 항의 같은 인자가 정확히 상쇄된다. q=h−μΦU+μ²K C∞+(G/c⁴)KμΣw μ_source F, K=1/(rN√b)로 두어 시간 미분과 압력 L을 모두 누적량·장 값으로 대체했다. 비영 초기 q와 끝점 복원 항을 모두 유지했다. 이는 현재의 선형 특성식과 outgoing 표현에서 성립하는 대수 항등식이다. 완전 비선형 별 또는 다른 외부 스칼라 표현의 정리로 확대하지 않는다.

분류: Counterexample candidate. 실제 계량 값 ν·λ·ζ와 배경 값은 기존 계산과 표본에서 동일했다. 독립 차분을 사용한 실제 변환 방정식 대조 차이는 {check['independent_actual_transformed_rhs_relative']:.12e}다. 질량만 변환한 방법, outgoing까지 변환한 방법, 모든 광자 누적량을 사용한 방법의 실제 세 cohort를 비교했다. 반경·방향·에너지·도착·일·응력 중 최대 차이는 {maximum:.12e}, 응력 모멘트 최대 차이는 {stress:.12e}로 원0.2%기준을 통과했다. 각 cohort는 원32개 구적 방향×8개 방출 구적, 총256개 실제 패킷이다. 이 대표 비교는 연속 전체 구간 오차 보증이 아니다.

분류: Counterexample candidate. 변환 후의 계량 일은 끝점 항등식으로 복원된다. 따라서 그 자체가 매우 작은 에너지 항등식 잔차를 보이더라도 **독립적인 일 적분 검증으로 세지 않는다.** 원 ν_t·λ_t를 그대로 쓰는 별도 실제 광선 두 개(같은 가장 이른 cohort의 근접 접선·근접 방사 방향)를 원 허용오차로 전파해 비교한다. 이 대조가 통과하기 전에는 전체 생산을 시작하지 않는다. 별도 대조의 유효 범위와 전체 시간 시계 대조 미완료도 유지한다.

분류: Counterexample candidate. 질량만 변환한 첫 대표 구간은352.34초, 전체 누적량 방식은{pilot['rows'][0]['seconds']:.2f}초였다. 누적량 방식의 세 구간 실측은 {[round(v['seconds'],2) for v in pilot['rows']]}초, 총 실행은{read(work/'pilot-receipt.json')['seconds']:.2f}초다. DOP853의 rtol=2e−8, atol=2e−11은 유지했다. 시간 미분식을 직접 쓰는 원 단계262대표 실행의 중간 계산은 완료로 세지 않았으며 그 실패·종료 기록을 보존한다.

분류: Conjectural. 원4시간 경로 예산의 비용 부적합은 그대로 남긴다. 가장 느린 대표 구간의 패킷당 비용을 전체 원17개 출력에 적용하고2배 여유를 둔 예상은 경로당{pilot['forecast_upper_seconds_per_full_path']/3600:.3f}시간이다. 후반 대규모 묶음에서도 패킷당 비용이 이 범위에 든다는 가정이며 보장된 종료시간이 아니다. 사용자의 반복 중단 방지 지시에 따라 정확도 기준을 유지하고 경로당8시간, 최대4개 CPU 경로로 재등록했다. 원 물리 격자·기간·경로 수는 늘리지 않는다. 완료 출력은 각각 저장하며 새 EOS·물질 재계산은 없다. 본문의 예상은 이 외부 전파 작업의 예산으로, 연구 전체 완료 시각이 아니다.

분류: Conjectural. 다음 실질 판정은 반환 광자 응력과 질량·lapse 원천의 원 전체 기간 구적 통과다. 이후 반환 원천의 별도 시간 대조, reciprocal 스칼라 응력·배경 연산자·같은 물질 질량과의 접합을 합친 경계를 실제 결합 해에 적용하고 전하를 읽어야 한다. 이 누적량 변환이나 대표 통과를 새 결합 해·최종 전하의 완료로 바꾸지 않는다. 이번은 실제 광자 전파의 계산 병목을 줄인 loophole progress이며 목표는 계속 진행 중이다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as stream:stream.write(f'\n\n## 단계263 — 같은 특성식의 누적량 전파\n\n분류: Counterexample candidate. 시간 미분과 광선 교차 압력 항을 정확한 변수변환으로 누적 방출량·장 값에 옮기고, 원 허용오차의 실제 대표 cohort 세 개를 완료했다. 서로 다른 변환의 궤적·응력 대조가 원0.2%기준을 통과했다. 복원된 계량 일 항등식은 독립 일 적분 검증으로 세지 않으며 원 시간 미분식의 별도 실제 광선 대조를 전체 생산의 전제조건으로 둔다. 실측 보수적 경로 예상6.09시간을 근거로 원4시간 비용 부적합을 보존하고8시간 예산으로 원 네 경로를 연결했다. 물질·EOS를 다시 계산하지 않았다. 전체 물리 전하는 미판정이다. [근거](../notes/{note.name}).\n'.encode())
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    modules=[root/'verification'/n for n in ['propagate_smooth_characteristics.py','propagate_outgoing_characteristics.py','propagate_primitive_characteristics.py']]
    assert all(sha(f)==sha(runtime/'verification'/f.name) for f in modules)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_primitive_photon']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
