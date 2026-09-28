"""Preserve the one-hour reference stop, its single-ray completion and the full-period propagation."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('prior',Path('.phase263-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-primitive-photon-full'
manifest=out.parent/'native-primitive-photon-full-manifest.json'
note=root/'notes/REQUEST264_DIRECT_REFERENCE_RETRY_KO.md'
ACTIONS=['fine','angular','temporal','geometry']


def package():
    work=runtime/'native-primitive-photon263-work';old=work/'direct-reference';reference=work/'direct-reference-264'
    plan=read(work/'execution-plan-264.json');status=read(work/'pipeline-status-264.json');failed=read(work/'pipeline-status.json')
    ref=read(reference/'result.json');collected=read(work/'result.json');stopped=read(old/'receipt.json')
    assert status['state']=='completed' and ref['passed'] and collected['passed'] and not collected['physical_final_charge_solved']
    assert failed['state']=='failed' and failed['completed']==[] and 'TimeoutError' in stopped['error']
    receipts={a:read(work/f'{a}-receipt.json') for a in ACTIONS+['collect']};receipts['reference']=read(reference/'receipt.json')
    assert all(v['error'] is None for v in receipts.values())
    bindings=plan['bindings']
    for p,h in bindings.items():assert sha((root if p.startswith('outputs/') else runtime)/p)==h,p
    rows={a:read(work/a/'result.json')['rows'] for a in ACTIONS};assert all(len(v)==16 for v in rows.values())
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for name in ['execution-plan.json','controller-start.json','pipeline-status.json']:copy(work/name,out/'phase263'/name)
    for name in ['receipt.json','progress.json','packet-0.json','packet-0.npz']:copy(old/name,out/'phase263/direct-reference'/name)
    for name in ['.phase263-reference.stderr.log','.phase263-controller.stderr.log']:copy(root/name,out/'phase263'/name.lstrip('.'))
    for f in sorted(reference.iterdir()):
        if f.is_file():copy(f,out/'direct-reference-264'/f.name)
    for name in ['execution-plan-264.json','controller-start-264.json','pipeline-status-264.json','stress-history.npz']+[f'{a}-receipt.json' for a in ACTIONS+['collect']]:
        copy(work/name,out/name)
    copy(work/'result.json',out/'collect-result.json')
    for f in sorted(work.glob('*.log')):copy(f,out/'logs'/f.name)
    for a in ACTIONS:
        for f in sorted((work/a).iterdir()):copy(f,out/a/f.name)
    for name in ['.phase264-direct-reference.py','.phase264-controller.py','.phase264-production.py']:copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
    ray=ref['rows'][1];first=ref['rows'][0];rel=ref['relative']
    e0=max(v for k,v in rel.items() if k.startswith('0:'));e31=max(v for k,v in rel.items() if k.startswith('31:'))
    stress=max(max(max(row) for row in v) for v in collected['controls'].values())
    boundary=max(max(v) for v in collected['boundary_controls'].values())
    seconds={a:receipts[a]['seconds'] for a in ACTIONS};terminal=collected['terminal']
    identity=max(r['energy_identity_relative'] for v in rows.values() for r in v)
    angular=max(r['angular_invariant_relative'] for v in rows.values() for r in v)
    final=dict(classification='Counterexample candidate',passed=True,
        verdict='INDEPENDENT_REFERENCE_AND_FULL_PERIOD_RETURNED_PHOTON_STRESS_ACCEPTED_BOUNDARY_FEEDBACK_PENDING',
        actual_same_joint_solution_GR_computed=False,phase263_reference_stopped_by_budget=True,reused_packet_indices=[0],independent_direct_reference_passed=True,
        direct_reference_maximum_relative={'0':e0,'31':e31},whole_period_production_complete=True,
        stress_control_maximum=stress,boundary_control_maximum=boundary,path_seconds=seconds,
        maximum_energy_identity_relative=identity,maximum_angular_invariant_relative=angular,
        returned_source_time_control_complete=False,reciprocal_scalar_source_complete=False,physical_GR_boundary_applied=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',previous_conditional_charge_sign='negative_retained',
        full_goal_complete=False,snapshot_KST=now)
    note.write_text(f'''# 독립 광선 대조 완료와 반환 계량 광자 응력의 전 기간 전파

분류: Counterexample candidate. **최종 전하는 아직 판정 불가다.** 기존 조건부 음의 전하 결론(단계260, 고정 외부 모형·한 차례 GR 피드백)은 유지한다. 이번에 완료한 반환 계량 광자 응력은 아직 결합 해에 되먹이지 않았으므로 그 결론을 바꾸거나 확장하지 않는다. 단계263에서 1시간 예산으로 멈춘 원 시간 미분식 독립 대조를 광선 하나만 다시 계산해 통과시켰고, 등록된 네 구적 경로와 collect 대조까지 원 기간 전체에서 완료했다. 기록={now}.

분류: Counterexample candidate. 단계263의 원 대조는 packet 0(μ=0.019842)을 {first['seconds']:.3f}초에 끝낸 뒤 packet 31(μ=0.995180) 도중 3600초 상한에서 `TimeoutError`로 멈췄고 생산 경로는 시작하지 않았다. 원 실행 계획·실패 상태·receipt·traceback과 저장된 packet 0를 변경 없이 보존·게시했다. 등록 바인딩이 모두 그대로임을 확인한 뒤 packet 0는 저장본을 재사용하고 packet 31만 같은 코드·cohort·DOP853 rtol=2e−8·atol=2e−11로 {ray['seconds']:.3f}초({ray['rhs_evaluations']}회 RHS)에 완료했다. 두 광선의 합은 {first['seconds']+ray['seconds']:.1f}초로, 원 1시간 관문 예산이 부족했다는 판단이 실측으로 확인됐다. 생산 경로 상한을 8시간으로 늘리면서 관문 예산은 늘리지 않은 실행 계획 오류였고, 물리·수치 기준의 실패가 아니다.

분류: Counterexample candidate. 누적량 방식 pilot-0와의 원 대조식 최대 상대 차이는 packet 0 {e0:.6e}, packet 31 {e31:.6e}로 원 0.2% 기준을 통과했다. packet 31의 자체 에너지 항등식 차이는 {ray['energy_identity_relative']:.3e}, 각 불변량은 {ray['angular_invariant_relative']:.3e}다. 적분 진행률은 같은 RHS 값을 같은 적분기에 넘기는 관측 전용 래퍼로 기록했다. 두 대표 광선의 궤적 대조이며 연속 전체 구간의 오차 보증은 아니다.

분류: Counterexample candidate. 네 경로(fine·angular·temporal·geometry)가 원 16개 출력 시각을 모두 완료했다. 경로별 실행 시간은 {', '.join(f"{a} {seconds[a]:.0f}초" for a in ACTIONS)}로, 최악 대표 cohort 비용을 모든 패킷에 적용하고 2배 여유를 둔 보수적 예상 6.09시간보다 훨씬 짧았다. 전 출력의 에너지 항등식 최대 {identity:.3e}, 각 불변량 최대 {angular:.3e}다. collect의 응력 모멘트 대조 최대 {stress:.6e}, 광자 질량 원천·lapse 특수해 경계 대조 최대 {boundary:.6e}로 원 0.2% 기준을 통과했다. 마지막 시각 fine 경로의 광자 질량 원천은 {terminal['photon_J_source_cm']:.12e} cm, lapse 특수해는 {terminal['photon_lapse_particular']:.12e}다. 격자·기간·경로 수·허용오차·기준은 바꾸지 않았다. 경로 상한 12시간은 부분 재개가 불가능한 점을 이유로 둔 여유이며 실제로 쓰이지 않았다.

분류: Conjectural. 이번 통과는 반환 계량에서 전파한 외부 광자 응력과 질량·lapse 원천이 원 기간 전체의 구적·기하 대조를 통과했다는 뜻이다. 반환 원천의 별도 64/128 시계 대조는 이번 실행에 없다. reciprocal 스칼라 응력·배경 연산자·같은 물질 질량과의 접합을 합친 경계를 실제 결합 해에 적용하고, 같은 해의 에너지·경계 이력으로 전하를 읽어야 최종 전하를 판정할 수 있다. 과거 다른 해의 진단값을 사후 가산하지 않는다. 이번은 외부 광자 전파의 실행 병목을 닫은 loophole progress이며 최종 목표 완료가 아니다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        f=root/f'docs/{name}.md';prefixes[f.relative_to(root).as_posix()]=dict(bytes=f.stat().st_size,sha256=sha(f))
        with f.open('ab') as stream:stream.write(f'\n\n## 단계264 — 독립 광선 대조 완료와 전 기간 광자 전파\n\n분류: Counterexample candidate. 단계263의 원 시간 미분식 대조가 1시간 관문 예산에서 멈춘 실패와 원 기록을 보존했다. 등록 바인딩 확인 뒤 저장된 packet 0를 재사용하고 packet 31만 같은 코드·허용오차로 완료했으며, 두 광선 대조 최대 {max(e0,e31):.3e}로 원0.2%기준을 통과했다. 이어 등록된 네 구적 경로가 원 16개 출력 시각을 모두 완료했고 응력·광자 질량 원천·lapse 경계 대조가 원0.2%기준을 통과했다. 격자·기간·경로 수·기준은 바꾸지 않았다. 반환 원천 64/128 시계 대조와 reciprocal 스칼라·배경 연산자·물질 질량 접합, 결합 해 적용이 남아 최종 전하는 미판정이며 단계260의 조건부 음의 전하 결론을 유지한다. [근거](../notes/{note.name}).\n'.encode())
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    files=[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={f.relative_to(root).as_posix():sha(f) for f in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_primitive_photon_full']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
