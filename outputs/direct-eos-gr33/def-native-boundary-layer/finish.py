"""Record the local radial verdict without claiming global convergence."""
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
result=read(out/'result.json');radial=read(out/'radial.json');prod=read(out/'production.json');audit=read(out/'audit.json');local=read(out/'local-audit.json')
gr=read(out/'gr/result.json');bound=read(out/'gr/bound.json');ga=read(out/'gr/audit.json');mixed=read(out/'contrast.json')
assert all(x['passed'] for x in [result,radial,prod,audit,local,gr,bound,ga])
fine=prod['paths'][-1];coarse=prod['paths'][0];ports=radial['shared_material_ports']
docs=[root/'docs'/name for name in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
report=root/'notes/REQUEST119_NATIVE_BOUNDARY_LAYER_KO.md'
text=f'''

## 실제 결합 결과와 판정 범위

분류: Counterexample candidate. 새19셀64/128 경로는 모두 원 구간을 완주했다. 생산은 합계{prod['seconds']:.2f}s로500s 상한 안이었다. 실제 native 끝점 압력/에너지 오차는{audit['native_endpoint_pressure_energy_relative']:.5e}, fine 에너지 잔차는{fine['energy_relative']:.5e}, 바리온 잔차는 초기 대기량 대비{fine['joint_baryon_relative_to_initial_atmosphere']:.5e}였다. 처음15개 위치의 기하·thermal/rate·밀도/주파수 미분 표가 비트 단위로 유지되었고 공유 광자 면 면적도 일치했다.

| 관심량의 대조 | 상대 차이 |
|---|---:|
|새19셀 직접 전하의64/128 시간 대조 |{100*result['controls']['time_direct']:.6f}% |
|16→19셀 직접 전하 공간 대조 |{100*radial['controls']['space_direct']:.6f}% |
|직접+광자 질량 합의 공간 대조 |{100*radial['controls']['space_total']:.6f}% |
|전체 출사 주파수 스펙트럼의 공간 대조 |{100*radial['controls']['space_surface_spectrum']:.6f}% |
|64/128 두 시계에서 공간 차이의 혼합 대조 / fine 직접 전하 |{100*mixed['mixed_difference_over_fine_direct']:.8f}% |

분류: Counterexample candidate. 원2% 기준을 적용한 직접 전하·총합·출사 스펙트럼의 국소 공간 대조가 통과했다. 추가 혼합 차이는 두 시계에서 공간 대비가 크게 바뀌는지 확인한 기술 통계이며, 전역 시간/공간 오차 상계가 아니다. 관측 시계는 이전16셀 결과와 비트 단위로 같고, 새 계산만 다른 시각에 평가해 차이를 줄이지 않았다.

분류: Counterexample candidate. 새 정규화 직접 전하 변화는 `{result['endpoint_direct_relative']:+.10e}`, 광자 질량 항과 합친 명목값은 `{result['endpoint_direct_plus_photon_mass']:+.10e}`다. 같은 새 원천의 실제 GR 응력·질량 기여와 외부 광자/퍼텐셜 상계 뒤에도 조건부 직접 하한은 `{ga['arbitrary_outward_angular_mass_conditional_lower']:+.10e}`로 양수다. 추가 GR 항만 포함한 명목 합 구간은 `[{ga['nominal_GR_interval'][0]:.8e},{ga['nominal_GR_interval'][1]:.8e}]`다. 이 좁은 구간에 원천의 공간·EOS·초기 제약 오차가 포함된다고 표시하지 않는다.

분류: Counterexample candidate. 작은 공유 물질 포트는 위 전하와 다르게 아직 수렴했다고 할 수 없다. fine 누적 질량은 기존 `{ports['join_mass']['original']:.8e}g`에서 새 `{ports['join_mass']['refined']:.8e}g`로 부호까지 바뀌었다. 새64/128 경로끼리도 `{mixed['local_ports'][0]['shared_mass_g']:.8e}/{mixed['local_ports'][1]['shared_mass_g']:.8e}g`여서 상대 시간 차이는{100*mixed['shared_mass_time_relative']:.2f}%다. 바깥 방향이 양수다. 누적 Killing 에너지·중성수소 포트도 큰 상대 차이를 보였다. 이를 전체 물질 유출 예측의 통과로 세지 않는다. 관심 전하의 원 기준 통과와 작은 포트 자체의 미수렴을 함께 보존한다.

분류: Counterexample candidate. 이번 성과는 기존34km 면 재구성을 실제 native 경계층으로 바꾸어 결합 진화를 완주하고, 그 국소 세분이 직접 전하를 거의 바꾸지 않는다는 유한 대조를 얻은 loophole progress다. 같은 국소 격자를 자동으로 더 세분할 근거는 약해졌다. 전체 방사·물질 원천이 연속 한계에서 인증되거나 최종 물리 전하가 정해진 것은 아니다.

분류: Conjectural. 다음 우선순위는 실제 초기 물질·광자 상태가 저장 GR 배경의 Einstein 제약과 맞는지 확인하고 불일치가 있으면 초기 상태에 반영하는 일이다. 나머지15개 내부 셀의 공간 오차도 남는다. 작은 추가 GR 구적이나 이번 국소 경계를 계속 자동 세분하는 대신, 초기 제약과 전체 원천 오차 중 잔여 전하를 실제로 흔들 수 있는 항을 먼저 닫아야 한다.
'''
report.write_bytes(report.read_bytes()+text.encode('utf-8'))
common='\n\n## 단계119 — 실제 native 경계층과 국소 공간 대조\n\n'
paragraphs={
 'model-definition.md':'분류: Counterexample candidate. 원16개 내부 셀의 마지막 한 셀을 실제 native 네 셀로 바꿔 총19개를 진화했다. 처음15개 입력은 그대로 유지했고, 최근접 점–공유 면 거리는33,975m에서300m로 줄었다. 새 네 위치의 초기 재고·열/H/복사 계수는 직접 native 평가했다. 기존 중앙 물질 유속·공유 HLL·양방향 광자·고정 계량·일차 재고 모형은 유지한다.',
 'observable-targets.md':f'분류: Counterexample candidate. 같은 관측 시계에서 국소16→19셀 차이는 직접 전하{100*radial["controls"]["space_direct"]:.6f}%, 직접+광자 질량 합{100*radial["controls"]["space_total"]:.6f}%다. 새 실제 GR 원천·외부 광자/퍼텐셜 상계 후 조건부 직접 하한은{ga["arbitrary_outward_angular_mass_conditional_lower"]:.8e}로 양수다. 전체 공간 오차막대나 최종 항성 전하로 승격하지 않는다.',
 'adiabatic-limit.md':'분류: Proven. 반경 셀을 나누면서 추가한 동일한 내부 면 유속은 인접 체적 합에서 소거되며 원 두 외부 포트는 유지된다.\n\n분류: Counterexample candidate. 이 항등식과 국소 공간 대조의 통과는 정적 유한 차수 흡수 정리를 변경하거나 실제 궤도 완화 검출을 증명하지 않는다.',
 'nonadiabatic-regime.md':'분류: Counterexample candidate. 실제19셀 열·수소·광자·물질 결합의 원64/128 경로가 완주했고 전하·출사 스펙트럼의 국소 공간 대조는 통과했다. 초기 과도 Cauchy 스펙트럼과 고정 계량이라는 전제를 유지한다. 경계 질량 교환량 자체는 큰 상대 시간/공간 차이가 남아 수렴 판정하지 않는다.',
 'failure-ledger-dynamic-chi.md':f'분류: Counterexample candidate. 작은 공유 질량 포트의 부호는16→19셀에서 바뀌었고 새64/128 시간 차이는{100*mixed["shared_mass_time_relative"]:.2f}%다. 관심 전하의 공간 차이{100*radial["controls"]["space_direct"]:.6f}%가 작다는 이유로 포트나 전체 상태 수렴으로 표시하지 않는다. 원 upwind 인공 전하와 이전 strict 중성수소 합산 실패도 그대로 보존한다. 이번 공간 대조는 마지막 셀만 대상으로 하므로 원15셀 전체 반경/초기 제약의 인증이 아니다.',
 'dynamic-charge-completion.md':'분류: Counterexample candidate. actual_native_boundary_layer_evolved, retained15_banks_bitwise, common_observer_clock, original_two_time_paths_completed, local_direct_total_spectrum_radial_gates_passed, updated_GR_conditional_lower_positive는true다. small_material_port_time_space_converged, global_radial_error_certified, full_neutral_trajectory_audit, initial_Einstein_constraints_matched, full_dynamic_GR_feedback, final_charge_solved, full_goal_complete는false다. 다음은 실제 초기 Cauchy 상태의 제약과 나머지 원천 공간 오차다.'}
for p in docs:p.write_bytes(p.read_bytes()+(common+paragraphs[p.name]+'\n\n상세: [단계119 보고](../notes/REQUEST119_NATIVE_BOUNDARY_LAYER_KO.md).\n').encode('utf-8'))
paths=docs+[report,root/'verification/def_native_boundary_layer.py',root/'verification/verify_native_boundary_layer.py']
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith(('-checkpoint.npz','.tmp')))
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-boundary-layer-manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='bf060cff9',actual_native_boundary_layer_evolved=True,
    local_charge_radial_gates_passed=True,small_material_port_converged=False,global_radial_error_certified=False,
    final_charge_solved=False,full_goal_complete=False,files=files))
master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_boundary_layer']=dict(classification='Counterexample candidate',report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)),
    last_cell_replaced_by_four_actual_native_cells=True,local_charge_and_spectrum_comparison_passed=True,small_material_port_converged=False,
    global_radial_error_certified=False,initial_Einstein_constraints_matched=False,full_goal_complete=False,final_charge_solved=False)
write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,r in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:r['bytes']]).hexdigest()==r['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),passed=True)))
