"""Preserve dynamic lapse and its actual compensated photon transport."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for part in iter(lambda:f.read(1<<20),b''):h.update(part)
    return h.hexdigest()

lapse=read(out/'result.json');repair=read(out/'corrected/audit.json');tr=read(out/'transport/result.json');ci=read(out/'transport/collision-input-audit.json')
assert lapse['passed'] and repair['passed'] and tr['passed'] and ci['passed']
# Frozen failed/pre-repair producers still explain their original bindings.
aliases={out/'first-plan.json':{'def_native_dynamic_lapse.py':out/'first-symbolic-producer.py'},
    out/'plan.json':{'def_native_dynamic_lapse.py':out/'lapse-producer-before-boundary-sign.py'},
    out/'transport/plan.json':{'def_native_compensated_transport.py':out/'transport/producer-before-time-alignment.py',
        'def_native_dynamic_lapse.py':out/'lapse-producer-before-boundary-sign.py'}}
for plan in [out/'first-plan.json',out/'plan.json',out/'corrected/plan.json',out/'transport/plan.json',out/'transport/execution-plan.json',out/'transport/collision-input-plan.json']:
    for p,h in read(plan)['bindings'].items():assert sha(aliases.get(plan,{}).get(Path(p).name,Path(p)))==h,(plan,p)
fine=tr['paths'][1];raw=lapse['paths'][0];fixed=repair['paths'][0]
report=root/'notes/REQUEST123_DYNAMIC_LAPSE_TRANSPORT_KO.md'
report.write_text(f'''# 단계123 — 실제 출사 광자 lapse와 보상 광자 수송

분류: Counterexample candidate. 수정된 초기 배경과 단계122의 비등방 GR 변화량에 실제 바깥 광자 패킷을 연결해 무한원 정규화 lapse를 구했다. 이 입력으로 원531셀·8각도·152주파수의 공유 광자 수송 연산자에서 반경·방향·주파수 변화량을 원3.434431ms 끝까지 실제 적분했다. 계량 입력만 저장한 상태에서 한 단계 나아갔다. 다만 이번에 진화한 것은 기하 수송 부분이며, 변화한 광자의 충돌·열·H 반응 및 물질 유체가 다시 변하는 전체 되먹임은 아직 미완료다. 전체 목표는 유지한다.

## 실제 외부 광자와 lapse

분류: Proven. shift가 없는 구면 계량에서 null Hamiltonian H=a*sqrt(p_r^2/h^2+L^2/R^2)를 미분하면, 국소 광자 에너지 epsilon과 방향 코사인 mu는 다음 식을 따른다. h는 반경 proper 길이 계수, R은 면적 반경이다.

```text
dr/dt       = c*a*mu/h
dmu/dt      = (1-mu^2)*[c*a/h*(R'/R-a'/a)+mu*(R_t/R-h_t/h)]
dln(eps)/dt = -c*mu*a'/h-(h_t/h)*mu^2-(R_t/R)*(1-mu^2)
```

분류: Proven. 고정 Einstein 반경에서 a=A*N,h=A/sqrt(b),R=A*r, u=alpha*delta_phi, zeta=delta_nu-delta_lambda로 두고, 원 주파수 좌표 Eref=a0(r)*epsilon을 유지하면 다음 계량 구동을 얻는다. 모든 delta는 배경과 별도로 저장한다.

```text
delta_vr       = c*N*sqrt(b)*mu*zeta
delta_mudot    = (1-mu^2)*[c*N*sqrt(b)*((1/r-nu')*zeta-delta_nu')-mu*delta_lambda_t]
delta_lnE_rate = -c*N*sqrt(b)*mu*(delta_nu'+u')-u_t-mu^2*delta_lambda_t
```

분류: Imported from prior work. 구면 GR 질량·압력·복사 수송의 연속 방정식은 [Salgado,2002](https://arxiv.org/abs/gr-qc/0201064)를 참조한다. 위 특성식은 현재 코드에서 Hamiltonian으로 별도 symbolic 검증했다.

분류: Counterexample candidate. 실제 바깥 면의 네 양의 각도 빈 안에서 점유수를 상수로 재구성하고, 실제 저장 스펙트럼·시간 이력의 패킷을 새 정확한 진공 계량의 null 궤적으로 전파했다. 초기 외부 광자는0이라는 원 조건을 유지했다. 모든 각도의 출력 에너지 적분은 원 포트 적분과{raw['emitted_energy_relative']:.3e}로 일치하며, 곡률 안의 impact 불변량 대조는{raw['ray_invariant']:.3e}다. 방출 에너지의 내부 질량 감소와 패킷 반경 압력을 함께 lapse에 넣었다. 수치 면에서 lapse=0을 강제하지 않는다.

분류: Proven. 한 패킷의 Killing 에너지 dE와 현재 위치/방향(rp,mup)의 lapse 기여는 (G*dE/c^4)*[integral(r0,rp)dr/(N*b^(3/2)*r^2)-mup^2/(rp*Np*sqrt(bp))]이다. 동일한 패킷의 내부 질량 감소가 첫 항에, 반경 압력이 둘째 항에 들어간다. scalar 진공 항은 delta_nu_scalar(r0)=+Phi0*U0-integral(r0,infinity)2*Phi*U/(r*b)dr다.

분류: Counterexample candidate. 외부 scalar 부분에는 현재 표현한 compact 파동을 바깥 radial 특성으로 연장한 항을 사용했다. 외부 광자가 새로 만드는 scalar와 전체 산란 되먹임까지 포함한 경계는 아니다. 실제 내부·출사 에너지의 잔차가 만드는 점근 질량 오차도 따로 저장했으며 보존을 맞추는 상수를 적합하지 않았다.

| 분류: Counterexample candidate — 실제 fine 경로 | 값 |
|---|---:|
| 최대 delta_nu | {fixed['maximum_delta_nu']:.9e} |
| 끝점 광자 lapse | {raw['endpoint_outer_photon_lapse']:.9e} |
| 부호 수정 후 끝점 scalar lapse | {fixed['endpoint_outer_scalar_lapse']:.9e} |
| 최대 점근 질량 잔차 | {raw['maximum_asymptotic_mass_residual_cm']:.9e} cm |
| 위 잔차에 의한 최대 lapse | {raw['maximum_ADM_residual_lapse']:.9e} |
| lapse4/8 구적 차이 | {lapse['controls']['quadrature']:.9e} |
| lapse64/128 원천 차이 | {lapse['controls']['time']:.9e} |

## 실제로 적용한 수송과 다음 입력

분류: Counterexample candidate. 배경의 물질·광자 비선형 이력은 재사용하고, 그 이력의 실제 각도·주파수 분포로 계량 구동항을 만들었다. conserved phase-space packet 수의 변화량 delta_N에 대해 기존 공유 반경·각도 수송 행렬과 두 단계 L-stable SDIRK를 적용했다. 변화량이 음수인 것은 기존 양의 분포의 감소를 뜻하므로 이를0으로 자르지 않았다. 깊은 면의 주어진 입사 packet flux 변화는0, 바깥 입사는0으로 지정했다. 이 조건이 더 깊은 별의 물리적 폐쇄는 아니다.

분류: Proven. Eref 노드 사이의 주파수 구동은 패킷 수를 옆 노드로 옮긴다. 이동률을 abs(omega)*E/abs(E_next-E)로 두면 이산 에너지 변화는 정확히 omega*E가 된다. 유한 상자 양 끝의 ghost 패킷은 자기 수·에너지로 별도 계상한다. 이를 물질 열로 넘기면 안 된다. 측정한 metric frequency work와 복사 충돌 열은 서로 다른 항이다.

| 분류: Counterexample candidate — 보상 수송 | 값 |
|---|---:|
|64/128 수송 시간 대조 | {tr['comparisons']['time']:.7%} |
|64/128 비선형 원천 이력 대조 | {tr['comparisons']['background_time']:.7%} |
| 독립 반경 포트·주파수 일 수지의 최대 오차 | {max(v['conservation_relative'] for v in tr['paths']):.9e} |
| 주파수 패킷 수·에너지 moment 오차 | {max(v['frequency_moment_relative'] for v in tr['paths']):.9e} |
| fine 끝점 Eref 에너지 변화 | {fine['endpoint_reference_energy_erg']:.9e} erg |
| fine 끝점 절대 에너지 변화 norm | {fine['endpoint_absolute_energy_erg']:.9e} erg |
| 누적 계량 주파수 일 | {fine['integrated_metric_frequency_work_erg']:.9e} erg |
| 주파수 ghost 에너지 | {fine['endpoint_spectral_energy_erg']:.9e} erg |

분류: Counterexample candidate. 이 수송의 reference-frequency 에너지는 최종 관측 전하가 아니다. 현재 L0에 추가 충돌·열·H·유체 응답이 들어가 있지 않으므로, 위 작은 값을 전체 되먹임의 오차 상계로도 사용하지 않는다. 현재 수치 판정은 세 지정 수송 경로에 대해서만passed=true다.

분류: Proven. 충돌 연산자에 넣을 실제 점유수 변분은 저장된 packet 변수 x에 대해 delta_F=x-(3*u+delta_lambda)*F0다. Phase122가 이미 포함한 canonical 체적 응답의 적분 에너지는 -u*E0-delta_lambda*Pr0, 적분 반경 압력은 -u*Pr0+delta_lambda*(R4_0-2*Pr0)다. 전체 수송 moment에서 이 항을 빼야 다음 GR 원천에 같은 변화가 두 번 들어가지 않는다.

분류: Counterexample candidate. collision-input-128.npz와additional-photon-sources-128.npz에 위 입력을 저장했다. 점유수 복원 대조는{ci['occupation_reconstruction_relative']:.3e}, canonical 항을 뺀 추가 광자 에너지와 압력의 최대 L1 norm은 각각{ci['maximum_additional_photon_energy_L1_erg']:.9e}/{ci['maximum_additional_radial_pressure_L1_erg']:.9e}erg다. 이것은 충돌·유체 되먹임을 이미 진화했다는 뜻이 아니다.

## 실패·수정과 실행 범위

분류: Counterexample candidate. 첫 Hamiltonian 검사는 시간 의존 각도 항의 부호 오류를 실제 계산 전에 잡았다. 올바른 항은mu*(R_t/R-h_t/h)다. 원 코드·계획과 실패를 보존했다. 이후 기존 lapse 소유자와 대조하면서 scalar 진공 적분 경계의 부호도 잘못 옮긴 것을 확인했다. 반경 미분 항등식은 맞았지만 경계식의 부호가 틀렸다. 독립 직접 적분으로 이를 확인하고 모든 lapse에2*Phi0*U0만 보정했으며, 이미 완료한 광자 궤적은 반복하지 않았다. 최대 수정은{fixed['maximum_boundary_sign_correction']:.3e}, 독립 부호 검사 오차는{repair['independent_boundary_sign_absolute_error']:.3e}다. 과거 수치 비교 통과가 이 부호를 증명했다고 부르지 않는다.

새 수송 파일 이름이 기존 모듈과 충돌해 첫 import가 실패했다. 기존 추적 파일의 바이트를 즉시 복원하고 고유 파일명을 사용했다. 물리적 단계는 실행되지 않았다. 그 뒤에는 전체 저장 배열에 pilot/restart 시각도 들어 있어 단순 시간 배열 비교가 실패했다. 기존 source export와 동일하게17개 기준 시각의 유일한 저장 snapshot을 선택하고1e-18 일치 기준을 적용했다. 실패한 binding 갱신도 hash 검사에서 실행 전에 차단됐으며 원 계획·코드는 보존한다. 이 과정에서 물리 이력이나 정확도 문턱을 변경하지 않았다.

lapse 세 경로{lapse['seconds']:.2f}s, 실제 보상 수송 세 경로{tr['seconds']:.2f}s, 충돌 입력 변환{ci['seconds']:.2f}s였다. 각각 사전45s, 실패 전 준비 비용을 차감한40s,15s 한도 안이다. 추가 native 상태 호출과 비선형 유체 재계산은0이며, import/준비 실패 비용을 성공 경로 시간에 숨기지 않는다.

분류: Conjectural. 다음 병목은 이번 delta_F와 보존 체적 변화가 실제 흡수·방출·산란, 물질 열/H 반응과 유체 운동을 바꾸고 다시 GR 원천으로 돌아오는 양방향 연결이다. 기존 공통 충돌·열·유체 소유자의 선형 응답 또는 보상 변수 진화를 사용해야 한다. 동시에 외부 생성 scalar와 깊은 원천의 인과적 경계를 닫아야 한다. 실제 기하 수송을 완료했지만 full_coupled_metric_feedback, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 이 단계는 loophole progress다.

근거: outputs/direct-eos-gr33/def-native-dynamic-lapse와corrected/transport 하위 계획·원 실패·수치 배열·독립 검사 및 네verification 생산/검사 파일.
''',encoding='utf-8')
docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
text={
'model-definition.md':'분류: Proven. null Hamiltonian에서 얻은 방향 변화의 시간 항은mu*(R_t/R-h_t/h)다. 고정 Einstein 반경과 원 주파수 좌표에서 반경·각도·주파수의 계량 구동을 별도로 계산한다. 충돌 점유수는packet 변분에서(3*alpha*delta_phi+delta_lambda)*F0를 빼야 한다.',
'observable-targets.md':'분류: Counterexample candidate. 실제 출사 광자의 질량 감소와 반경 압력을 함께 lapse에 연결하고 같은 내부·대기 광자 수송의 기하 응답을 완주했다. 수송의 reference-frequency 에너지 변화는 최종 전하가 아니다. 추가 충돌·유체·외부 scalar 폐쇄 뒤에 다시 전하를 읽어야 한다.',
'adiabatic-limit.md':'분류: Proven. 계량 주파수 일과 실제 흡수·방출에 의한 물질 열은 별도 항이다. canonical 체적 응답을 포함한 뒤 동일한 photon moment를 새 원천에 다시 더하면 이중 계산이다.\n\n분류: Counterexample candidate. 이 보존 분리가 기존 단열 흡수 경계의 이탈이나 전체 EOS 인증을 뜻하지 않는다.',
'nonadiabatic-regime.md':f"분류: Counterexample candidate. 별도 저장 계량 변화량을 원531셀·8각도·152주파수의 공유 광자 수송 연산자에 적용해 원 시간 구간을 완주했다. 시간 대조{tr['comparisons']['time']:.3e}, 원천 경로 대조{tr['comparisons']['background_time']:.3e}다. 추가 충돌·열/H·유체 되먹임은 아직 적용하지 않았다.",
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. Hamiltonian 검사로 동적 각도 항 부호를 실제 계산 전에 수정했다. scalar lapse 적분 경계의 부호 오류는 독립 적분과 기존 소유자 대조로 확인하고 완료한 광자 궤적을 반복하지 않고 수정했다. 원 생산 코드·실패·raw 수치 verdict는 보존한다. 수송 파일 이름 충돌은 기존 파일 복원 후 해소했으며, 저장 이력은 원 기준 시각 선택 방식으로 읽는다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. actual_emitted_photon_lapse, compensated_geodesic_transport_completed, collision_occupation_input_exported는true다. full_collisional_thermochemical_hydrodynamic_feedback, full_exterior_scalar, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 다음은 실제 충돌·물질 보존 상태로의 되먹임과 추가 GR 원천의 양방향 폐쇄다.'}
for p in docs:p.write_bytes(p.read_bytes()+('\n\n## 단계123 — 동적 lapse와 실제 보상 광자 수송\n\n'+text[p.name]+'\n\n상세: [단계123 보고](../notes/REQUEST123_DYNAMIC_LAPSE_TRANSPORT_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/n for n in ['def_native_dynamic_lapse.py','verify_native_dynamic_lapse.py','def_native_compensated_transport.py','verify_native_compensated_transport.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-dynamic-lapse-transport-manifest.json'
flags=dict(classification='Counterexample candidate',checkpoint='2fc1fb85f',actual_photon_lapse=True,actual_compensated_geodesic_transport=True,
    collision_input_exported=True,full_coupled_metric_feedback=False,full_exterior_scalar=False,final_charge_solved=False,full_goal_complete=False)
write(manifest,dict(flags,files=files));master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_dynamic_lapse_transport']=dict(flags,report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)));write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,p in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:p['bytes']]).hexdigest()==p['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),actual_transport=True,passed=True)))
