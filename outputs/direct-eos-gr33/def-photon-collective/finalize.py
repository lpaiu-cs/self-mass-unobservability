"""Publish the recovered collective response without changing any failed verdict."""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp
from scipy.integrate import quad
import def_photon_collective as m

root=Path('.');out=root/'outputs/direct-eos-gr33/def-photon-collective';recovery=out/'continuum'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
write=lambda p,x:p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=read(recovery/'result.json');old=read(out/'result.json')
assert r['passed'] and not old['passed'] and not r['full_dynamic_charge_solved']
assert not read(root/'outputs/direct-eos-gr33/def-photon-spatial-coupling/result.json')['passed']
t,f0,slope,v0=sp.symbols('t f0 slope v0',real=True)
M=f0*t+slope*t*t/2;J=v0*M+f0*t*t/2+slope*t**3/3
assert sp.simplify(sp.diff(M,t)-(f0+slope*t))==0
assert sp.simplify(sp.diff(J,t)-(v0+t)*(f0+slope*t))==0
write(out/'continuum-symbolic.json',dict(classification='Proven',passed=True,
    claim='Exact zeroth and first antiderivatives for each linear nonnegative density interval; adding positive residue atoms and mixing adjacent q measures preserves positivity. The previously checked detailed-balance recoil/Legendre assembly then retains finite energy and number.',
    limits='A theorem about the declared finite representation, not about its distance from exact nonideal plasma physics.'))

checks=[]
for plan_path in [out/'plan.json',out/'coupling-plan.json',recovery/'pilot-plan.json',recovery/'coupling-plan.json']:
    for path,value in read(plan_path)['bindings'].items():
        target=root/path
        if plan_path.name=='pilot-plan.json' and path=='verification/def_photon_collective_continuum.py':
            target=recovery/'pilot-initial-source.py'
        assert sha(target)==value,(str(plan_path),path)
        checks.append(dict(plan=plan_path.as_posix(),source=target.as_posix(),sha256=value))
s=np.load(out/'inventory.npz');b=np.load(m.old.OUT/'bank.npz')
up=float(s['plasma_u']);Z=m.inventory_reader.g.d.CHARGES
bound=float((s['ni']*(Z[:,None]-s['charges'])).sum()/s['ne'])
fraction=float(quad(lambda u:u*u*(u/(2*np.sinh(u/2)))**2,0,up)[0]/(4*np.pi**4/15))
write(out/'source-audit.json',dict(classification='Counterexample candidate',passed=True,bindings_verified=len(checks),checks=checks,
    bound_electrons_per_free_electron=bound,Planck_capacity_below_plasma=fraction,
    first_cell_edges=b['edges_u'][:2].tolist(),first_cell_representative=float(b['u'][0]),
    first_cell_capacity_fraction=float(b['Ci'][0]/b['Ci'].sum()),
    warning='Zero count of representatives below the plasma frequency is not zero sub-plasma support; the first source cell begins at zero.',
    preserved_failures=dict(discrete_response=False,initial_arctan_pilot=False,
        fixed_count_log_grid_max_moment_error=.0003875915143740638),
    source_snapshot_note='The initial pilot plan source is preserved as continuum/pilot-initial-source.py. The fixed-count-log trial is pilot-log-source.py. Later uniform-log controls and production are bound by the recovery plan.'))

report=root/'notes/REQUEST63_COLLECTIVE_PHOTONS_KO.md';before=report.read_bytes()
assert '## 최종 결합 판정'.encode() not in before
rows=r['rows'];base=rows[0];comp=r['comparisons']
addition=f'''
## 최종 결합 판정

분류: Counterexample candidate. 수정 생산의 `continuum/result.json`은 **passed=true**다. 원 이산 속도 `result.json`의 **passed=false**를 변경하지 않았다. 아래 응답 차이는 같은 초기 물질 섭동 1K의 에너지 노름으로 정규화했다.

| 항목 | 결과 | 사전 기준 |
|---|---:|---:|
| 16/32/64 시간 차수 | {base['time_order']:.9f} | ≥1.8 |
| 마지막 시간 차이 | {base['time_differences_initial'][-1]:.8e} | <1e−3 |
| 연속 밀도 128/256 대조 | {comp['density-mesh']:.8e} | <1e−6 |
| q 절점 129/257 대조 | {comp['q-grid']:.8e} | <1e−6 |
| 각도 24/48·모드 8/12 대조 | {comp['angular']:.8e} | <1e−6 |
| 집단−독립 전자 효과의 32/64 시간 차이 | {r['paired_time_difference_initial']:.8e} | <1e−6 |
| 최대 에너지 수지 결함 | {max(z['balance'] for z in rows):.8e} | <1e−9 |
| 최대 원 에너지식 잔차 | {max(z['energy_residual'] for z in rows):.8e} | <1e−9 |
| 최대 선형 방정식 잔차 | {max(z['solver_residual'] for z in rows):.8e} | <1e−11 |
| 에너지·광자 수 영모드 잔차 | {max(max(z['algebra']['energy_number_null_relative']) for z in rows):.8e} | <1e−10 |
| 에너지 노름 증가 | {max(z['entropy_growth'] for z in rows):.1f} | <1e−10 |

분류: Counterexample candidate. 채택한 q 세분 경로와 같은 셀 적분의 독립 전자 경로 차이는 초기 노름의 **{r['collective_effect_initial']:.9e}**다. 물질 온도 섭동의 차이는 **{r['collective_temperature_effect_K']:.9e}K**, 끝점 온도 섭동은 **{r['final_temperature_perturbation_K']:.12f}K**다. 원 미수락 이산 속도 경로와의 차이는 {r['rejected_discrete_difference_initial']:.9e}로 저장했다. 이 결과는 지정한 충돌 모형의 효과이며 실제 항성 관측 신호나 새 자유낙하 전하 검출이 아니다.

분류: Counterexample candidate. 수정 생산은 **{r['seconds']:.3f}초, 최대 {r['peak_memory_GB']:.3f}GB**로 사전 상한 안에서 끝났다. 같은 방정식의 저장 독립 전자 경로를 재사용했고, 최초 입력·초기 시험·두 생산 계획의 {len(checks)}개 SHA 결속을 독립적으로 재확인했다. 이번 단계는 **loophole progress 및 선언 유한 모형의 conditional theorem progress**다. 전체 동적 전하 완료로 표시하지 않는다.
'''
report.write_bytes(before+addition.encode('utf-8'));assert report.read_bytes().startswith(before)

texts={
'model-definition':'분류: Proven. 공통 온도 Maxwell/RPA 전자·다종 이온 밀도 스펙트럼은 양수이며, 같은 양의 셀 적분과 상세평형 반동 블록에 연결하면 유한 광자 수·에너지 보존을 유지한다. 분류: Counterexample candidate. 현재 분자 기체 EOS의 점유수와 같은 물질/흡수/공간 시간식으로 집단 산란을 연결했다. 중성 원자 질량 공통 규약·비충돌·진공 광자 기하와 극점 근사는 명시적으로 남긴다.',
'observable-targets':f"분류: Counterexample candidate. 같은 셀 적분의 독립 전자 경로와 집단 산란 경로의 지정 응답 차이는 초기 노름의 {r['collective_effect_initial']:.7e}다. 원 이산 속도 구적 차이가 효과 크기와 비슷한 경로는 미수락으로 보존했고 연속 유전 스펙트럼의 결합 대조로 수정했다. 분류: Conjectural. 이 성분 효과는 실제 궤도 구동·동적 전하·정적 비교 및 관측 신호의 완료를 뜻하지 않는다.",
'adiabatic-limit':'분류: Proven. 선언 공통 온도 RPA의 정적 전자 구조인자는 (x²+Σz)/(x²+1+Σz)이며, 지정 이산 상세평형 반동은 LTE·광자 수 영모드를 보존한다. 분류: Counterexample candidate. 정적 합·이차 모멘트·허수 주파수 응답이 맞아도 좁은 광자 셀의 유한시간 응답 수렴은 별도로 실패할 수 있었다. 정적 수락을 동적 수락으로 바꾸지 않는다.',
'nonadiabatic-regime':f"분류: Counterexample candidate. 같은 2.316801밀리초 결합에서 연속 RPA 밀도의 시간 차수는 {base['time_order']:.6f}이고 밀도·파수·각도 대조 및 짝지은 집단 효과의 시간 차이가 모두 원 1e-6 응답 기준 안이다. 분류: Conjectural. 유한 셀 기하·극점·물리 매질 오차의 엄밀한 전 구간 상계, 대기와 전체 GR·전하 연결은 별도 요구사항이다.",
'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 단계 63 원 이산 속도 생산은 속도 구적 1.02750e-5·각도 3.39369e-6로 1e-6 기준에 미달했다. 양의 SVD 모멘트 통과로 이를 구제하지 않았다. 연속 유전 스펙트럼으로 표현을 바꾸고 공명 날개의 arctan/고정 점 수 보간 오류를 로그 간격 제한으로 수정해 같은 결합 기준을 통과했다. 원 실패·문턱·소스는 보존했다. 낮은 plasma 주파수의 대표점 수 0은 원 첫 셀에 해당 영역이 없다는 뜻이 아니다.',
'dynamic-charge-completion':'분류: Counterexample candidate. 단계 63 현재: 같은 EOS 전자·다종 이온 점유수의 집단 광자 산란을 물질·흡수·공간 시간식에 연결하고, 원 구적 실패를 연속 유전 스펙트럼 표현으로 수정해 결합 대조를 통과했다. 분류: Conjectural. 결합 전자 산란·실제 온도 흡수와 비충돌 등 물리 근사 오차, 비균일 반경·대기·전체 GR 되먹임 및 기존 속도/전하 수렴·실제 구동·비교·관측은 미완료다.'}
paths=[]
for stem,body in texts.items():
    p=root/f'docs/{stem}.md';data=p.read_bytes();assert b'## Phase63 ' not in data
    text='\n\n## Phase63 — EOS 점유수와 일치하는 집단 광자 결합\n\n'+body+'\n\n분류: Imported from prior work. 식·원 실패·표현 수정·수치 판정·예산·남은 경계는 [단계 63 보고서](../notes/REQUEST63_COLLECTIVE_PHOTONS_KO.md)에 둔다.\n'
    p.write_bytes(data+text.encode());assert p.read_bytes().startswith(data);paths.append(p)

source=[root/f'verification/{name}.py' for name in ['def_photon_collective','def_photon_collective_validate','def_photon_collective_continuum','def_photon_collective_continuum_validate']]
files=sorted(p for p in out.rglob('*') if p.is_file() and p.name!='manifest.json')+source
decision='CONTINUOUS_EOS_MATCHED_RPA_RESPONSE_ACCEPTED_DISCRETE_ALIASING_REJECTED'
manifest=out/'manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='9b1fce346',passed=True,
    original_discrete_passed=False,progress_class='loophole progress; conditional theorem progress',
    decision=decision,full_dynamic_charge_solved=False,sha256={p.as_posix():sha(p) for p in files}))
global_path=root/'outputs/direct-eos-gr33/gr-photon-collective-milestone-manifest.json'
write(global_path,dict(classification='Counterexample candidate',checkpoint='9b1fce346',decision=decision,
    prior_milestone_manifest_sha256=sha(root/'outputs/direct-eos-gr33/gr-photon-finite-jump-milestone-manifest.json'),
    full_dynamic_charge_solved=False,sha256={p.as_posix():sha(p) for p in paths+[report,manifest]}))
paper=root/'paper/revision-manifest.json';raw=paper.read_bytes();previous=read(paper);assert len(previous)==171
entry=dict(classification='Counterexample candidate',progress_class='loophole progress; conditional theorem progress',decision=decision,
    EOS_population_matched_collective_coupled_response_passed=True,original_discrete_response_passed=False,
    bound_electron_kernel_complete=False,actual_temperature_opacity_certified=False,
    nonideal_collisional_quantum_error_certified=False,native_scalar_kernel_identified=False,
    physical_atmosphere_closed=False,full_GR_photon_feedback_evolved=False,heat_velocity_time_order_passed=False,
    full_dynamic_charge_solved=False,new_stellar_steps=0,new_native_EOS_calls=4,new_physical_queries=0,
    source_bindings_verified=len(checks),next_bottleneck='Bound-electron redistribution and actual-temperature absorption, remaining kinetic approximations; then stratified transport, physical atmosphere, full GR feedback, original velocity/charge convergence, drive/comparator and observations.',
    report=report.as_posix(),report_sha256=sha(report),evidence_manifest=global_path.as_posix(),evidence_manifest_sha256=sha(global_path))
tail=json.dumps(dict(request63_photon_collective=entry),ensure_ascii=False,indent=2).encode()
paper.write_bytes(raw[:raw.rfind(b'}')].rstrip()+b',\n'+tail[2:]+b'\n')
assert all(read(paper)[k]==value for k,value in previous.items())
for p in [manifest,global_path]:
    for key,value in read(p)['sha256'].items():assert sha(root/key)==value,key
print('Verified source bindings',len(checks),'; evidence files',len(files),'; preserved paper entries',len(previous),'; passed',r['passed'])
