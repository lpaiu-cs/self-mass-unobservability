"""Preserve actual exterior photon propagation and its GR stress source."""
from pathlib import Path
from types import FunctionType
import datetime,importlib.util,sys
spec=importlib.util.spec_from_file_location('b',Path('.phase260-publish.py'));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha,copy=b.read,b.write,b.sha,b.copy
out=root/'outputs/direct-eos-gr33/native-dynamic-photon'
manifest=out.parent/'native-dynamic-photon-manifest.json'
note=root/'notes/REQUEST261_PHYSICAL_EXTERIOR_PHOTONS_KO.md'


def package():
    work=runtime/'native-dynamic-vacuum261-work';old=runtime/'native-dynamic-photon261-work'
    result=read(work/'result.json');boundary=read(work/'photon-boundary/result.json')
    assert result['passed'] and boundary['passed'] and read(work/'pipeline-status.json')['state']=='completed'
    assert read(old/'pipeline-status.json')['completed'][0]['returncode']==-15
    assert read(work/'symbolic.json')['passed'] and read(work/'photon-boundary/symbolic.json')['passed']
    bindings={}
    for plan in [work/'plan.json',work/'controller-start.json',work/'result.json',work/'photon-boundary/plan.json']:
        for p,h in read(plan)['bindings'].items():
            # The runtime outputs link resolves to this same repository on /mnt/e.
            source=(root if p.startswith('outputs/') else runtime)/p
            assert sha(source)==h,(str(plan),p)
            assert p not in bindings or bindings[p]==h,p
            bindings[p]=h
    for p in work.rglob('*receipt.json'):assert read(p).get('error') is None,p
    assert read(work/'photon-boundary/receipt.json')['error'] is None
    modules=[root/'verification'/n for n in ['propagate_physical_exterior.py','propagate_exterior_vacuum.py']]
    assert all(sha(p)==sha(runtime/'verification'/p.name) for p in modules)
    preserved={p:h for p,h in read(b.manifest)['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir()
    for folder,label in [(work,'actual'),(old,'cost-stop')]:
        for p in folder.rglob('*'):
            if p.is_file() and p.suffix in ['.json','.npz','.py','.log'] and 'initialization' not in p.parts:
                copy(p,out/label/p.relative_to(folder))
    for name in ['.phase261-controller.py','.phase261-collect.py','.phase261-vacuum-controller.py','.phase261-vacuum-collect.py','.phase261-boundary.py']:
        copy(root/name,out/name.lstrip('.'))
    copy(Path(__file__),out/'publication-producer.py')
    opt=read(work/'optical-check.json');check=read(work/'check.json');terminal=result['terminal']
    final=dict(classification='Counterexample candidate',passed=True,
        verdict='ACTUAL_DYNAMIC_PHOTON_STRESS_AND_GR_PARTICULAR_SOURCE_ACCEPTED_COMPLETE_BOUNDARY_PENDING',
        actual_same_joint_solution_GR_computed=False,
        same_accepted_reference_emission_applied=True,actual_photon_worldlines_propagated=True,
        actual_photon_stress_inserted_in_GR_constraint=True,
        previous_conditional_negative_charge_preserved=True,new_dynamic_final_charge_read=False,
        physical_GR_boundary_applied=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        source_controls=result['controls'],source_time_control=result['time_control'],boundary_controls=boundary['controls'],
        snapshot_KST=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat())
    note.write_text(f'''# 실제 외부 광자 전파와 GR 제약 원천

분류: Counterexample candidate. **단계260의 조건부 음의 전하 수락은 보존했다. 새 동적 외부를 포함한 최종 전하는 아직 판정하지 않았다.** 이번에는 같은 수정 EOS 배경과 수락한 high/low 방출 이력을 사용해, 입사 계량이 배경 광자의 에너지·반경·방향에 미치는 실제 1차 전파를 전체 원 기간에 계산했다. 네 각도·방출 구적·기하 대조를 완료했고 기존 0.2% 구적 기준 및 2% high 방출 시간 기준을 통과했다. 이 결과의 광자 응력 분포를 기존 반경의 GR 질량·lapse 제약에 직접 투영했다. 물리 광자 원천을 실제로 생성한 loophole progress이며, 완전한 새 경계를 물질에 재적용한 결과는 아니다.

분류: Counterexample candidate. 원 네 방출 각도 구간, 17개 배경 이력 시각, 전체 기간 0.0034344311179287023초, 기존 high/low Radau 방출을 재사용했다. 새로운 유체·EOS 계산이나 물리 격자·기간 확대는 없다. 보존된 참조 에너지 방출과 배경 광자의 기하 변화는 별도로 전파했다. 완성된 high 응력에는 두 항을 함께 넣었고, low의 참조 방출도 보존했으나 반환 계량에 의한 low 기하 항은 아직 적용하지 않았다.

분류: Proven. 구면 대칭의 각운동량 보존으로 h=δlnH−δν−δlnL, dh/dr=−δν′−μδλ_t/(cN√b), δt′=[−ζ−(1−μ²)h/μ²]/(cN√bμ)를 얻는다. 고정 시각에서 δr=−cN√bμδt와 δμ=(1−μ²)[h+(1/r−ν₀′)δr]/μ이며, δlnH=δν+u_launch+h다. 코드의 각운동량 항등식, 패킷 응력 측도, 평탄한 시간 의존 lapse 해를 검산했다. 이는 선형 특성식이며 완전 비선형 GR 정리가 아니다.

분류: Counterexample candidate. 기존 입사장 U 값은 동일하게 재현했다. 저장된 bilinear Born U를 일관되게 미분했을 때 기존 미분 표와의 차이는 {max(check['incident_jet_representation_relative']):.12e}, 입사 경계 lapse 값의 차이는 {check['incident_boundary_value_relative']:.12e}였다. 이 수치는 지정된 표현의 대조이며 연속 EOS/균일 미분 오차의 인증이 아니다.

분류: Counterexample candidate. 최초 대표 실행은 일반 내부 반경 변환을 외부 진공에도 반복해 비용이 컸다. 정확한 PID·시작 tick·boot·명령을 확인하고 225초 무렵 의도적으로 종료했으며 소스·계획·확인 결과와 -15 종료를 보존했다. 같은 진공의 dr/dx=N√b와 dx/dr=1/(N√b)를 한 번 풀어 재사용했다. 독립 원 변환 차이는 {opt['independent_original_optical_fraction']:.12e}, 계량 값·미분의 최대 차이는 {max(opt['identical_metric_relative'].values()):.12e}였다. 32패킷의 동일 조회는 {opt['old_seconds']:.6f}초에서 {opt['new_seconds']:.6f}초로 줄었다. 이 비용 비율을 전체 계산이나 물리 정확도 증명으로 확대하지 않는다. 새 대표 구간 세 개의 실측으로 전체 경로 예산을 확인한 뒤 네 독립 경로를 병렬 실행했으며 원 수락 기준은 유지했다.

분류: Counterexample candidate. 미세 경로 종단에서 기하에 의한 물리 출구 에너지 증분은 {terminal['physical_launch_energy_increment_erg']:.12e}erg, 현재 패킷의 좌표 에너지 증분은 {terminal['instantaneous_packet_energy_increment_erg']:.12e}erg, 전파 중 계량 일은 {terminal['propagated_work_erg']:.12e}erg였다. 종단의 에너지 항등식 상대 차이는 {terminal['energy_identity_relative']:.12e}, 각운동량 항등식 차이는 {terminal['angular_invariant_relative']:.12e}였다. 최대 반경 이동={terminal['maximum_radius_shift_cm']:.12e}cm, 최대 방향 변화={terminal['maximum_direction_change']:.12e}. 이동과 에너지를 큰 배경에 더해 반올림으로 소실시키지 않고 별도 상태로 보존했다. 출구 debit과 이후 전파 일을 혼동하지 않았다.

분류: Proven. 같은 스칼라 진공에서는 K(r)=∫ᵣ∞ds/(Nb^(3/2)s²)=1/(rN√b)다. 원 진공 선형 제약에 넣는 패킷 광자 원천의 lapse 특수해는 −G/c⁴ ΣE₀K[(1+μ²)(δlnH−δν−δλ−δr/(rb))+2μδμ]다. 에너지뿐 아니라 응력의 위치와 방향 변화도 포함한다. 이 공식과 질량 측도의 미분을 기호 검산했다. 스칼라 상호작용과 배경 연산자 변화가 자동으로 포함된 공식이 아니다.

분류: Counterexample candidate. 실제 패킷 이력으로 위 GR 광자 원천을 계산했다. 미세 경로 종단의 photon J 원천은 {boundary['rows'][0]['terminal_photon_J_source_cm']:.12e}cm, lapse 특수해는 {boundary['rows'][0]['terminal_photon_lapse_source']:.12e}다. 같은 기존 적용 계량에서 독립적으로 읽은 질량 커널과의 차이는 {boundary['independent_saved_mass_kernel_relative']:.12e}였다. 새로운 photon lapse 특수해 최댓값과 기존 적용 전체 lapse 최댓값의 비는 {boundary['maximum_particular_lapse_over_existing_applied_lapse']:.12e}다. 이 비는 원천 크기 비교이며 최종 전하 오차의 상한이 아니다.

분류: Conjectural. 다음 실제 연결에는 반환 계량이 배경 광자에 주는 기하 응답, 스칼라의 상호 에너지 교환 및 배경 연산자 변화, 같은 물질 질량 경계와의 접합이 필요하다. 그 합으로 만든 경계를 결합 해에 적용하고 전하를 다시 읽어야 한다. 이 photon 특수해나 과거 다른 해의 진단값을 기존 전하에 사후 가산하지 않는다. 배경 질량 정규화, 자기GR, 원 EOS/균일 오차, 완전 비선형·정적 EFT·관측 연결도 열린다.
''',encoding='utf-8',newline='\n')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        text=f'\n\n## 단계261 — 같은 이력의 동적 외부 광자 원천\n\n분류: Counterexample candidate. 기존 수락 해의 방출과 같은 EOS 진공·입사장을 사용해 광자의 에너지·반경·방향 변화를 원 전체 기간에서 실제 전파했고 네 원 구적 대조와 high 방출 시간 기준을 통과했다. 실제 패킷 응력을 원 반경의 GR 질량·lapse 제약 특수해에 연결했다. 단계260의 조건부 음의 전하 수락은 보존하되, 새 원천에 대한 완전한 경계·물질 재적용과 최종 전하는 미판정이다. 반환 계량 기하 응답, 상호 스칼라 일·배경 연산자·질량 접합을 누락한 채 전하에 보정값을 가산하지 않는다. [근거](../notes/{note.name}).\n'
        with p.open('ab') as f:f.write(text.encode())
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,verified_runtime_bindings=bindings))
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    m['native_dynamic_photon']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')
