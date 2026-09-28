from verify_direct_eos_gr import append_scoped_material_progress
def verify():
 import json
 from pathlib import Path
 import numpy as np
 import gr_conservative_composition_tangent as c
 plan=c.bindings()
 star=c.initialize(None)
 folder=c.OUT/'path-1'
 times=c.parent.wall.prior.time_nodes(plan,1)
 previous=older=None
 energy_flux=baryon_flux=np.zeros(star.n+1,dtype=c.ld)
 older_energy,older_baryon=energy_flux.copy(),baryon_flux.copy()
 rows=[]
 for step in range(8):
  cp=np.load(folder/f'step-{step:04d}.npz')
  assert cp['time_seconds']==times[step]
  d=cp['delta'].copy(); y=star.base+d
  star.material_cache={c.e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
  z=star.evaluate(d)
  for key in ['m','mf','a','N','Q','aux','dU']:
   assert np.array_equal(z[key],cp[key]),(step,key)
  value,energy_defect,baryon_defect=c.budget(star,z,cp['integrated_energy_flux'],cp['integrated_baryon_flux'])
  assert np.array_equal(energy_defect,cp['normalized_energy_defect'])
  assert np.array_equal(baryon_defect,cp['normalized_baryon_defect'])
  assert c.budget_passed(value,plan)
  cone=c.parent.cones(z); assert cone['sampled_cone_inside_light_cone']
  norm=0.
  if step:
   h=times[step]-times[step-1]
   coeff=c.weights(h,None if step==1 else times[step-1]-times[step-2])
   native,_=c.residual(star,d,previous,older,h,coeff)
   norm=float(np.max(abs(native)/c.ATOL))
   assert norm<=1,(step,norm)
   c0,c1,c2=coeff; area=4*np.pi*star.rf**2
   energy_flux,older_energy=(-c1*energy_flux-c2*older_energy+h*c.e.C*area*z['fluxes'][1])/c0,energy_flux
   baryon_flux,older_baryon=(-c1*baryon_flux-c2*older_baryon+h*c.e.C*area*z['fluxes'][0])/c0,baryon_flux
   older,previous=previous,(d,z)
  else:
   assert np.all(d==0)
   previous=older=(d,z)
  assert np.array_equal(energy_flux,cp['integrated_energy_flux'])
  assert np.array_equal(baryon_flux,cp['integrated_baryon_flux'])
  assert np.all(energy_flux[[0,-1]]==0) and np.all(baryon_flux[[0,-1]]==0)
  rows.append(dict(step=step,time_seconds=float(times[step]),maximum_native_residual=norm,**value,**cone))
 print('ACTUAL COMPOSITION COUPLED STEP SEVEN',json.dumps({k:v for k,v in rows[-1].items() if k!='isotope_inventory_defects'}),flush=True)
 
 oldplan=c.previous.bindings()
 for key in ['method','nonlinear_absolute_tolerances','nonlinear_relative_tolerance','maximum_stage_iterations','finite_conservation_gates','time_refinement_gate','conduction_time','duration_seconds','refinements','coordinate_edges_seconds']:
  assert plan[key]==oldplan[key],key
 op=plan['operator_check']
 assert op['passed'] and op['corrected_native_residuals']['without_composition']>1
 assert op['corrected_native_residuals']['native_composition']<=1
 for parent_path,score in [(c.HISTORY.OUT,72.5041790539357),(c.previous.OUT,3066.268985402333)]:
  assert json.loads((parent_path/'path-1/failure.json').read_text())['step']==6
  oldlog=[json.loads(s) for s in (parent_path/'path-1/iterations.jsonl').read_text().splitlines()]
  assert oldlog[-1]['iteration']==23 and oldlog[-1]['residual_norm']==score>1
 current=[json.loads(s) for s in (folder/'iterations.jsonl').read_text().splitlines()]
 for step,count in [(6,8),(7,11)]:
  logs=[r for r in current if r['step']==step]
  assert [r['iteration'] for r in logs]==list(range(count))
  assert logs[-1]['residual_norm']==rows[step]['maximum_native_residual']<=1
 return rows
rows=verify()
r,s=rows[6],rows[7]
body=f'''분류: Counterexample candidate. 실제 보존형 결합 진화가 이전 두 풀이의 여섯 번째 단계 실패를 넘어 일곱 번째 단계까지 진행했다. 여섯 번째 단계는 8개 잔차 평가 뒤 {r['time_seconds']:.17g}초에서 수락됐고, 31개 native 방정식의 최대 잔차/허용오차는 {r['maximum_native_residual']:.12g}로 기준 1을 통과했다. 최대 국소 에너지 수지 결함/초기 열용량은 {r['maximum_local_energy_defect']:.12g}, 최대 국소 바리온 상대 결함은 {r['maximum_local_baryon_defect']:.12g}, 26종 재고 결함/전체 초기 바리온은 {r['maximum_isotope_inventory_defect']:.12g}이다. 최대 국소 정지계 특성속도는 {r['maximum_local_rest_characteristic_speed_over_c']:.12g}c이며 모든 표본 특성값이 실수다. 일곱 번째 단계도 11개 평가 뒤 {s['time_seconds']:.17g}초, 잔차 점수 {s['maximum_native_residual']:.12g}로 수락됐다. 일곱 단계의 상태·EOS 반환값·계량·보존량·BDF 잔차와 공유 면 교환의 시간 누적을 저장 자료에서 재계산해 대조했다.

분류: Counterexample candidate. 앞 절의 매 반복 접선 갱신만 사용한 경로는 다섯 번째 단계 뒤 여섯 번째 단계의 24회 반복 상한에서 종료했다. 마지막 잔차 점수 72.5041790539357과 실패 상태를 보존한다. 물질 면 흐름의 유한차분이 donor 선택 경계를 가로지르는 별도 결함도 확인했다. 실제 기준 상태에서 선택한 donor를 반복 행렬의 미분 동안 유지하자, 지정 방향의 에너지 행렬 작용 상대 차이가 0.0016594838460034416에서 2.0202545076012503e-11로 줄었다. 그러나 그 수정만 적용한 별도 실제 경로도 여섯 번째 단계에서 24회 뒤 점수 3066.268985402333으로 실패했다. 이 결과는 방향 미분 대조의 통과가 실제 비선형 진화의 완료를 대신할 수 없음을 보여 준다. 두 원 실패·소스·계획·실패 반복값을 모두 보존한다.

분류: Counterexample candidate. 실제 병목을 넘긴 수정은 수송에 따른 조성 변화가 native EOS의 압력·내부에너지와 복사/전도 불투명도에 미치는 응답을 비선형 반복 행렬에 포함한 것이다. 각 단계의 시작 상태에서 좌우 이웃 조성으로 향하는 두 방향을 사용하며, 같은 EOS/불투명도 평가기로 구한 응답을 그 단계 안에서 유지한다. 최대 절대 핵종 변화 1e-5, 혼합 비율 1 이하의 볼록 혼합을 사용했다. 보존된 실패 상태에서 9,454개 native 조성 평가를 수행한 같은 한 번의 보정 대조에서는, 조성 응답을 생략하면 잔차 점수 21.299866784224427, 포함하면 0.005876933788246076이었다. 후자는 실제 31개 잔차·BDF 에너지/바리온/26종 수지·특성속도 기준도 통과했다. 이 한 번의 보정은 별도의 구현 대조이며, 실제 시간 경로를 그 보정 상태로 시작하지 않았다.

분류: Counterexample candidate. 새 동결 소스 gr_conservative_composition_tangent.py의 실제 경로는 재검증한 수락 상태 0–5를 이어받고, 여섯 번째 단계를 수락된 다섯 번째 상태에서 시작했다. 두 방향의 국소 조성 기저와 고정한 계량 미분·폭 제한 희소 행렬은 반복 풀이의 근사일 뿐, 실제 26종 상태나 EOS를 축소하지 않는다. 실제 수락 잔차에서는 전체 EOS, 물리적 donor 선택, 모든 조성, 계량 및 31개 방정식을 다시 계산한다. 같은 물리식·시간 절점·24회 반복 상한·여덟 선 탐색 보폭·Anderson 이력·잔차 허용오차·보존 문턱을 유지했다. 수락 상태를 자르거나 재규격화하지 않았다. 실제 시간 적분은 첫 후방 Euler와 가변 간격 BDF2이다. 계획의 상속된 SDIRK 기호 검사 항목은 이전 방법의 검사 기록이며, 현재 방법의 근거는 별도 BDF2 일관성 검사와 실제 BDF 잔차다.

분류: Conjectural. 전체 0.42117120910640804초와 1/2/4 시간 세분화 대조는 아직 완료되지 않았다. 이후 구간의 실제 진화를 같은 연속 실행으로 이어간다. 현재 조건부 공통 LTE 두 열수송 성분·지정 전도 시간·반사 벽은 실제 복사/대기/외부를 완성한 모형이 아니다. 이번 유한 방향 대조와 native 단계 통과를 물리 EOS 인증, 연속/엄밀 미분 오차 보증, 핵반응 포함 진화나 실제 구동·관측 폐쇄로 확장하지 않는다. 이 사안들은 계속 남는 연구 작업이다.

재현: verification/gr_conservative_composition_tangent.py의 동결 계획으로 run --refinement 1/2/4 --workers 15, 전체 경로 완료 뒤 compare. 원 실패는 gr-conservative-refreshed-tangent 및 gr-conservative-branch-tangent의 path-1/failure.json이다. 새 실제 상태는 gr-conservative-composition-tangent/path-1/step-0006.npz 및 step-0007.npz이다.
'''
brief='분류: Counterexample candidate. 접선 갱신과 donor 경계 수정만으로는 여섯 번째 실제 단계를 풀지 못한 원 실패를 보존했다. native EOS의 수송 조성 응답을 반복 행렬에 연결한 새 경로가 같은 시간·물리식·문턱에서 여섯 번째와 일곱 번째 단계를 수락했다. 실제 31개 잔차 점수는 각각 0.263285와 0.120397이며 보존·특성속도 기준을 통과했다. 두 방향 근사는 반복 행렬에만 적용되고 실제 26종과 EOS는 유지된다. 분류: Conjectural. 전체 시간 수렴·연속/미분 오차·물리 EOS·복사/외부/반응·관측 폐쇄는 미완료다.'
append_scoped_material_progress('conservative-composition-tangent','실제 EOS 조성 응답을 포함한 보존형 GR 비선형 진화',body,brief,
 'native 조성 응답 연결로 기존 여섯 번째 단계 실패를 넘어 일곱 번째까지 실제 31개 잔차·보존·특성속도 문턱 통과. 전체 시간 수렴·물리 EOS/외부/반응/관측 미완료.',verify)
