"""Bind the actual radial coupling and preserve the failed first mesh."""
from pathlib import Path
import hashlib
import json
import subprocess

OUT=Path('outputs/direct-eos-gr33/def-photon-radial-gr-resolved')
OLD=OUT.parent/'def-photon-radial-gr'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
r=read(OUT/'result.json');a=read(OUT/'replay.json');s=read(OUT/'state-check.json');i=read(OUT/'input.json')
assert r['passed'] and a['numerical_target_passed'] and a['saved_input_bindings_passed']
assert not read(OLD/'result.json')['passed']
rows=[read(OUT/(name+'.json')) for name in ['fine-32','fine-64','fine-128','spatial-p2','one-way']]
balance=max(x['exchange_balance'] for x in rows);gr=max(x['GR_residual'] for x in rows);ph=max(x['photon_residual'] for x in rows)
history=rows[2]['history'];final=history[-1];thermal=a['endpoint']['fine-128']
cost=sum([read(OLD/'setup-failed.json')['wall_seconds_upper_bound'],read(OLD/'pilot.json')['seconds'],
    read(OLD/'result.json')['seconds'],read(OUT/'pilot.json')['seconds'],r['seconds'],
    read(OLD/'replay.json')['seconds'],a['seconds'],s['seconds']])
assert cost+60<900
names=['질량 가중 속도 RMS','scalar RMS','내부 면 속도','내부 면 광도']
table='\n'.join(f'| {name} | {100*v["last"]:.6g}% | {v["order"]:.6g} | {100*v["spatial"]:.6g}% |' for name,v in zip(names,r['comparisons'].values()))
report=Path('notes/REQUEST82_PHOTON_MATTER_RADIAL_GR_KO.md')
text=f'''# 단계82 — 광자·물질 교환의 실제 반경 GR 연결

분류: Counterexample candidate. **실제 외곽 반경 r/R={i['coordinate_radius_R']:.10f}에서 광자·물질 에너지와 운동량을 기존 전체 선형 GR 행렬에 넣고, GR 압축·속도를 다시 광자·물질에 돌려주는 경로를 계산했다.** {i['coordinate_duration_seconds']*1e6:.6f}μs 구간에서 원 시간·공간2% 기준을 통과했다. 별도 광자 계산의 출력을 표시하는 데 그치지 않고 매 적분 단계에서 함께 풀었다. 완료한 범위는 **두 물질 부피·한 내부 광자 면의 국소 P1 조건부 모형**이다. 전 반경 광자 수송·실제 대기 경계·완전 비선형 항성 진화는 완료하지 않았다.

## 실제로 연결한 상태와 식

분류: Imported from prior work. 단계70의24,134개 광자 주파수 셀, 각 물질 부피의18개 유한 이온 점유, 원 충돌·흡수·방출·유한 주파수 이동 표를 재사용했다. 내부 준위와 미보유 화학 좌표의 조건부 LTE, 빠진 충돌 반응 등 이전 물리 한계는 유지한다. 단계80의 전 반경 GR 공간 연산자와 압력·질량 제약을 사용한다. 기존 전도 응답과 새 광자 결합을 동시에 모두 진화한 것은 아니다.

분류: Counterexample candidate. 국소 구간은 실제 최외곽 native 물질 셀 안에 있으며 전체 폭은 약1.389km다. 물질 부피2개와 가운데 광자 유속1개에 P1 에너지·유속 모멘트를 두었다. 초기 구동은 저장된 이웃 상태의 A*N*T 기울기를 국소 Taylor 근사한 것으로, 임의의 새 완화시간이나 속도 배율을 넣지 않았다. 양 끝의 생략된 광자 유속은 이 부분 계산에서 닫았다. 실제 주변 물질·대기가 단열 또는 반사 경계라는 물리적 주장은 아니다. 두 초기 온도 편차는 {i['initial_temperature_offsets_K'][0]:.6f}, {i['initial_temperature_offsets_K'][1]:.6f}K다.

분류: Counterexample candidate. 최종 GR 자유도는{i['GR_degrees_of_freedom']:,}개이며 기계 모드 축약 없이 적분했다. 두 부피의 비단열 에너지 S, 압력 보정 P와 공유 면 유속 F가5개 연결 변수다. 각 SDIRK2 단계에서 전체 GR 선형 풀이, 전체 주파수 광자·유한 이온 풀이,5×5 연결식을 함께 푼다. 아래 q는 물리적 변위/scalar 좌표이며 w는 광자 운동량을 포함한 속도 변수다.

```text
V_j = 4*pi*(100R)^3 integral_j N*a*A^4*r^2 dr
L_infinity = 4*pi*(100R*r_face)^2*N_face^2*A_face^4*F

q_dot = w - H_F*F
M*w_dot = -K*q + B_S*S + B_P*P
y_dot + C*y = H_rho*rho_dot + C*E_rho*rho - B_v*v_dot
S,P,F = O*(y-y_initial) + D_rho*rho
rho = D_q*q + D_S*S
```

분류: Counterexample candidate. y에는 물질 온도, 유한 이온 점유, 광자 에너지·유속이 포함된다. 광자 수송의 배경 계량 인자·면 위치·계수는 고정했고, 압축과 속도의 선형 역응답을 연결했다. 전체 계량 변화에 따른 광자 궤적·불투명도 변화나 모든 고차 각 모멘트의 되먹임을 연결한 것은 아니다.

분류: Proven. 선언 모형에서 내부 면 광도는 두 부피의 에너지 장부에 반대 부호로 들어가 합이 상쇄된다. 기존 LTE EOS의 광자를 그대로 둔 채 독립 광자를 더하지 않고, 물질 항과 독립 광자 응력을 한 번씩 센다. 등방 복사의 관성계 변환에는 (E+P)v=4Ev/3가 필요하다. 이 에너지·운동량 구분은 [Park의 GR 복사 유체 방정식](https://academic.oup.com/mnras/article/367/4/1739/1747567)의 공변 보존식과 관측자 변환을 참고했다. 해당 문헌이 이번 이산화나 P1 정확도를 인증하는 것은 아니다.

분류: Counterexample candidate. 새 MHD 물질 EOS와 저장 native EOS의 같은 상태 압력·열용량 차이는 각각 {i['pressure_native_difference']:.8g}, {i['heat_capacity_native_difference']:.8g}다. 국소 자유에너지 밀도에 T와 무관한 a+b*rho를 더해 정적 값을 맞추는 접선 모형을 명시했다. a={i['free_energy_density_anchor']['constant_erg_cm3']:.8g}erg/cm³, b={i['free_energy_density_anchor']['coefficient_erg_g']:.8g}erg/g이다. 새로운 물리 EOS나 전 반경 평형의 인증으로 읽지 않는다. 정지 에너지에는 기존 동위원소 질량 기준 CX={s['isotope_rest_CX']:.12g}를 유지한다. 이를1로 바꿔 생기는 약0.7% 차이는 모형의 에너지 불일치가 아니다. 같은 CX에서 실제 GR 압력·에너지와 native 값의 차이는 각각 {s['pressure_match_relative']:.6g}, {s['energy_match_relative']:.6g}다.

분류: Proven. 위 자유에너지 보정은 delta p=-a, delta(rho*u_lnrho)=-a이므로 압축 열 계수 p-rho*u_lnrho를 바꾸지 않는다. 압력의 밀도 접선과 열용량도 바꾸지 않는다. 이 항등식과 면 에너지 상쇄·복사 운동량 변환을 기호 검사했다. 실제 자유에너지의 전 상태 범위 오차 보증은 아니다.

## 최초 실패와 직접 수정

분류: Counterexample candidate. 최초 균일 패치의 시간 대조는 통과했지만 GR 속도·scalar 공간 차이는65.0–65.7%여서 실패했다. 원 코드·계획·모든 경로·실패 판정은 def-photon-radial-gr에 보존했다. 이 구간의 음향 이동 폭은39.484cm인데 격자 간격은2170.49cm였다. 전체 부피에12점만 둔 압축 구적도 좁은 음향 층을 놓쳤다. 따라서 단순히 적분 단계 수를 늘려서는 해결되지 않는 공간 문제였다.

분류: Counterexample candidate. 별도 수정 계획에서 세 압력 경계의 ±1.25 음향 이동 폭에32구간씩, 바깥 빛의 전파 영역에256구간을 한 번 배치했다. 압축 평균은 실제 유한요소 구적점을 모두 사용하고, 운동량 유속 보간도 반경 선형 대신 동일한 누적 redshift 부피로 맞췄다. 광자 격자·물질 부피·적분 기간·수락 기준은 유지했다. 최초 입력 파일의 mass_scale 연결 오류도 계산 전 수정했고 원 계획·소스를 보존했다. 기준 미달 후 격자를 반복해서 늘리지 않았다.

## 실제 경로와 독립 재계산

분류: Counterexample candidate. 아래 시간 오차는32/64/128단계 중64–128의 공통 시각 차이를128경로의 최대 절댓값으로 나눈 값이고, 차수는 앞뒤 차이의 비다. 공간은 같은128단계의2차–4차 GR 대조다. 모든 행에 시간·공간2%, 시간 차수1.5 이상 기준을 동일하게 적용했다.

| 읽어낸 양 | 시간 차이 | 시간 차수 | 공간 차이 |
|---|---:|---:|---:|
{table}

분류: Counterexample candidate. 별도로 물질 온도의 끝점 시간·공간 차이는 {s['material_temperature_time_relative']:.8g}, {s['material_temperature_spatial_relative']:.8g}, 시간 차수는{s['material_temperature_time_order']:.6g}였다. 전체 경로의 최대 GR 선형 잔차는{gr:.6g}, 광자 선형 잔차는{ph:.6g}, 두 부피의 비단열 에너지 수지 오차는{balance:.6g}다. 이는 지정 이산 모형과 유한 대조의 수치 수락이며 엄밀한 연속체 미분 오차 상계가 아니다.

분류: Counterexample candidate. 끝점의 물질 온도 변화는 {thermal['material_temperature_increment_K'][0]:.8f}K와 {thermal['material_temperature_increment_K'][1]:.8f}K였다. 안쪽으로 전달된 비단열 에너지는 약{abs(thermal['nonadiabatic_energy_erg'][0]):.8g}erg다. 같은 저장 상태에서 물질 EOS와 이온 Hessian으로 압력·에너지 연결항을 별도 재구성한 상대 차이는 최대{max(a['independent_thermodynamic_ports_relative']):.6g}였다. 끝점 광자–물질 충돌 교환률은 부피당 약{abs(a['material_collision_power_erg_s'][0]):.8g}erg/s이고, 광자의 손실과 물질의 획득이 {a['collision_power_cancellation_relative']:.6g} 상대 오차로 일치했다.

분류: Counterexample candidate. GR 역응답을 끈 수치 대조와 비교하면 끝점 물질 온도는 최대{max(abs(x) for x in a['GR_return_temperature_increment_K']):.8g}K 달라진다. 따라서 실제 GR→물질 되먹임이 계산에 들어갔다. 이 대조는 압축·속도와 기준 단열 차감을 함께 끈 코드 대조이며 다른 보존 항성의 예측으로 해석하지 않는다. 내부 면 속도의 최대는{max(abs(h['face_velocity_m_s']) for h in history):.8g}m/s다. scalar 가중 RMS 최대는{max(abs(h['scalar_mass_RMS']) for h in history):.6g}로 매우 작다. 그 상대 대조율을 동적 전하 검출·관측 가능성으로 승격하지 않는다.

## 비용·판정 경계

분류: Counterexample candidate. 최초 실패 설정 비용을120초 상계로 포함한 계산 타이머 합은{cost:.2f}초이며 실행 시작 비용60초 여유를 더해도 원 총900초 안이다. 수정 경로의 실측 기반 예상은{read(OUT/'pilot.json')['forecast_seconds']:.2f}초, 실제 생산 계산은{r['seconds']:.2f}초였다. CPU/BLAS1스레드·메모리6GB 제한, pilot 관측 메모리는{read(OUT/'pilot.json')['memory_GB']:.3f}GB다. 새 EOS 호출과 주파수 격자 확대는0이다. 단계80·81의 이전 성공과 실패 판정은 바꾸지 않았다.

분류: Counterexample candidate. 이번 성과는 **loophole progress: 광자·유한 이온 물질·전체 선형 GR의 실제 국소 상호 결합과 지정 구간 수치 수렴**이다. 분류: Conjectural. 다음 물리 병목은 닫아 둔 두 바깥 광자 면을 이웃의 인증된 상태 및 실제 표면·대기와 연결하는 것이다. 한 외곽 상태의 표를 전 반경에 그대로 적용하거나 이번 짧은 구간을 전 항성 장기 안정성·물리 EOS 인증·정적 EFT 비흡수성·관측 추론 완료로 해석할 수 없다. 전체 비선형·계수/계량 되먹임·동적 전하·관측은 미완료다.
'''
report.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
paragraphs={
 'model-definition.md':'분류: Counterexample candidate. 실제 외곽 r/R=0.994284의 두 물질 부피·한 광자 면을 기존 전체 선형 GR에 연결했다.24,134개 주파수와 부피당18개 유한 이온 좌표, 비단열 에너지·압력·광자 운동량, GR 압축·속도 역응답을 매 단계 함께 풀었다. P1·고정 수송 계량/계수·국소 접선 매칭·닫힌 패치 경계라는 범위를 유지한다.',
 'observable-targets.md':'분류: Counterexample candidate.18.5344μs 실제 결합에서 물질 온도 변화는 약±0.25352K, GR 역응답을 끈 대조와의 끝점 차이는 최대0.0007521K다. 이는 국소 상태 응답이며 자유 계수의 정적 모형을 벗어나는 관측량이나 동적 전하 검출이 아니다. 전파된 scalar 가중 RMS는 약6e−52로 상대율만으로 의미를 부여하지 않는다.',
 'adiabatic-limit.md':'분류: Proven. LTE 복사 응력을 한 번만 계수하고 등방 복사 운동량(E+P)v를 보존해야 한다. 국소 T독립 자유에너지 밀도 a+b*rho의 정적 매칭은 p-rho*u_lnrho를 보존한다. 분류: Counterexample candidate. 실제 이온·광자 비평형 경로의 수치 통과는 전체 물리 단열 극한이나 정적 EFT 이탈의 증명이 아니다.',
 'nonadiabatic-regime.md':'분류: Counterexample candidate. 전체 선형 GR51353자유도와 국소 광자·물질을18.5344μs 결합했다.32/64/128단계 시간 대조 최대0.190%,2차–4차 공간 대조 최대0.191%로 원2% 문턱을 통과했다. 에너지 수지 오차 약1.1e−12와 별도 EOS 압력·에너지 재구성을 확인했다. P1 각도/두 부피의 수송 공간 자체나 연속체 오차 인증은 아니다.',
 'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 최초 균일 패치는 시간 통과·GR 공간65% 실패를 보존한다.18.5μs 음향 이동39.5cm보다21.7m 격자가 컸고 전 부피12점 압축 구적도 경계 층을 놓쳤다. 세 음향 경계와 빛의 전파 영역을 한 번 해상하고 유한요소 구적·보존 부피 유속을 연결한 별도 경로가 원 기준을 통과했다. 닫힌 패치 바깥의 유속을 물리적으로 영이라고 인증하지 않는다.',
 'dynamic-charge-completion.md':'분류: Counterexample candidate. local_radial_photon_material_GR_coupling=True, local_numerical_target_passed=True다. 전 반경 광자·실제 대기·전체 비선형·동적 전하·관측은 false를 유지한다. 새 물리 연결은 실제 반경의 유한 광자·물질 교환을 GR에 넣고 GR 압축·속도를 되돌려 진화한 것이다.'}
prefixes={}
for name,paragraph in paragraphs.items():
    p=Path('docs')/name;before=p.read_bytes();prefixes[p.as_posix()]=dict(bytes=len(before),sha256=sha(p))
    addition='\r\n\r\n## Phase82 — 광자·물질의 실제 반경 GR 결합\r\n\r\n'+paragraph+'\r\n\r\n분류: Imported from prior work. 식·원 실패·보정·수치 판정·범위는 [단계82 보고서](../notes/REQUEST82_PHOTON_MATTER_RADIAL_GR_KO.md)에 둔다.\r\n'
    p.write_bytes(before+addition.encode('utf-8'));assert p.read_bytes()[:len(before)]==before
write(OUT/'documentation-prefixes.json',prefixes)
write(OUT/'completion.json',dict(classification='Counterexample candidate',local_coupling_passed=True,
    original_failure_preserved=True,accounted_compute_seconds_upper_bound=cost,startup_allowance_seconds=60,
    maximum_energy_balance=balance,independent_thermodynamic_ports_passed=True,
    photon_material_collision_balance_passed=True,actual_GR_state_match_passed=True,
    full_radial_photon_atmosphere=False,full_nonlinear_evolution=False,full_dynamic_charge=False))
manifest=OUT.parent/'photon-radial-gr-manifest.json';assert not manifest.exists()
files=[Path('verification')/name for name in ['def_photon_radial_gr.py','def_photon_radial_gr_resolved.py','def_photon_radial_gr_audit.py']]
files += [report,*[Path('docs')/name for name in paragraphs]]
files += [p for directory in [OUT,OLD] for p in directory.iterdir() if p.is_file()]
write(manifest,dict(classification='Counterexample candidate',checkpoint='f470378f3',
    progress_class='loophole progress; actual local photon/material/radial GR feedback and finite-interval convergence',
    numerical_target_passed=True,full_radial_photon_atmosphere=False,full_nonlinear_evolution=False,
    files={p.as_posix():sha(p) for p in sorted(files)}))
paper=Path('paper/revision-manifest.json');previous=read(paper);assert len(previous)==190 and 'request82_photon_radial_gr' not in previous
entry=dict(classification='Counterexample candidate',local_radial_photon_material_GR_coupling=True,
    numerical_target_passed=True,full_radial_photon_atmosphere=False,full_nonlinear_evolution=False,
    full_dynamic_charge_solved=False,new_EOS_calls=0,report=report.as_posix(),report_sha256=sha(report),
    evidence_manifest=manifest.as_posix(),evidence_manifest_sha256=sha(manifest))
write(paper,dict(previous,request82_photon_radial_gr=entry));now=read(paper);assert all(now[k]==v for k,v in previous.items())
stage=[p.as_posix() for p in files]+[manifest.as_posix(),paper.as_posix()]
subprocess.run(['git','-c','core.autocrlf=false','add','-f','--',*stage],check=True)
print(json.dumps(dict(manifest_files=len(files),staged_files=len(stage),paper_entries=len(now),report=report.as_posix()),indent=2))
