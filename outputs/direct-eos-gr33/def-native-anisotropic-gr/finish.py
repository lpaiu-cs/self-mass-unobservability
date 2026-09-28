"""Record the anisotropic field evolution and preserve both failed raw gates."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

raw=read(out/'fields.json');fixed=read(out/'characteristic/fields.json');audit=read(out/'audit.json');source=read(out/'sources.json')
fine=fixed['paths'][0];symbolic=read(out/'symbolic.json');box=read(out/'characteristic/check.json')
assert not raw['passed'] and not fixed['passed'] and audit['passed'] and symbolic['passed'] and box['passed']
for plan in [out/'plan.json',out/'characteristic/plan.json',out/'audit-plan.json']:
    for p,h in read(plan)['bindings'].items():assert sha(root/p)==h,p
ratio=fine['endpoint_compact_with_metric']/fine['endpoint_direct']-1
report=root/'notes/REQUEST122_ANISOTROPIC_DYNAMIC_GR_KO.md'
report.write_text(f'''# 단계122 — 실제 비등방 원천의 GR·스칼라 변화량 진화

분류: Counterexample candidate. 단계121의 수정 초기 배경에서 완주한 실제 물질·광자 원천을 사용해, 반경 압력과 각도4차 모멘트를 포함한 GR 질량 제약 및 스칼라 파동의 시간 응답을 계산했다. 배경에 더하면 사라지는 작은 변화량을 별도 변수로 보존한다. 물질·광자 궤적을 다시 계산하지 않았으며, 공간 수송이 변화한 계량에 반응하는 마지막 되먹임은 아직 적용하지 않았다. 이 단계는 loophole progress이며 전체 목표 완료가 아니다.

## 비등방 광자와 보존 체적의 연결

분류: Proven. 순간 scalar 균형, 초기 Pi_phi=0인 선언 배경의 일차 변분에서, 물질 좌표 바리온·엔트로피·이온 재고와 광자의 canonical momentum을 고정한 계량 부분 단계는 다음 식을 따른다. f=delta_phi, dl=delta_m/(r*b), Hg=Eg+Pg, Kg=(dPg/dlnrho)_s, R4=integral(E*mu^4)다.

```text
delta_Eg = eFg - Hg*(3*alpha*f + dl)
delta_Pg = pFg - Kg*(3*alpha*f + dl)
delta_Er = eFr - 4*alpha*Er*f - (Er+Pr)*dl
delta_Pr = pFr - 4*alpha*Pr*f - (3*Pr-R4)*dl
```

분류: Proven. 광자 수밀도 변화는 -3*alpha*f-dl, 각 패킷의 국소 에너지 변화는 -alpha*f-mu^2*dl, 방향 코사인 제곱 변화는 -2*mu^2*(1-mu^2)*dl이다. 이들을 함께 변분하면 위 반경 압력 식을 얻으며 광자 trace는0이다. 실제 비등방 광자를 등방 기체의 Gamma 하나로 대체하는 폐쇄는 이 식과 같지 않다.

분류: Proven. 이 폐쇄에서는 비등방 복사를 포함해도 delta_m=r^2*b*Phi*f+J와 J'+(nu'+lambda')J=4*pi*r^2*A^4*eF가 성립한다. U=r*f, dx=dr/(N*sqrt(b))로 두면 U_tt/c^2-U_xx+Veff*U=Keff*J-4*pi*r*N^2*A^4*[alpha*tF+r*Phi*(eF-pF)]다. 구현한 Veff/Keff의 완전 전개와 질량 항 상쇄, 기체만 남기는 극한을 symbolic.json에서 검증했다. 이 대수 증명은 전체 비선형 항성 진화나 미지의 경계 원천을 증명하지 않는다.

분류: Imported from prior work. 극면적 좌표의 Einstein 질량·반경 압력 식과 scalar 파동식은 [Novak,1997](https://arxiv.org/abs/gr-qc/9707041)의2.23–2.26 및 [Salgado,2002](https://arxiv.org/abs/gr-qc/0201064)의233–237에 대응한다. 현재의 비등방 보존 변분은 위 식을 사용해 별도로 유도했다.

## 실제 저장 원천을 사용한 결과

분류: Counterexample candidate. 원64/128 시간 경로,531개 물질/광자 셀,17개 저장 시점을 그대로 사용했다. 모든 양의 바깥 방향 광자와 실제 안쪽 유입의 Killing 에너지를 적분했다. 포트 이력 대조는 coarse/fine에서{source['paths'][0]['port_history_relative']:.7%}/{source['paths'][1]['port_history_relative']:.7%}로 원2percent 기준 안이다. 추가 EOS 호출·유체 단계는0이다.

| 분류: Counterexample candidate — 수정 특성 적분의 수치 결과 | 값 |
|---|---:|
| 직접 정규화 전하 끝점 | {fine['endpoint_direct']:+.12e} |
| 표현한 영역의 GR 질량·응력·퍼텐셜을 포함한 끝점 | {fine['endpoint_compact_with_metric']:+.12e} |
| 위 GR 항의 직접값 대비 상대 변화 | {ratio:+.8%} |
| 최대 delta_phi | {fine['maximum_delta_phi']:.9e} |
| 셀 중심 직접 제약 적분의 최대 delta_m | {audit['maximum_delta_mass_cm']:.9e} cm |
| 최대 delta_lambda | {audit['maximum_delta_lambda']:.9e} |
| 최대 delta_nu_prime | {audit['maximum_delta_nu_prime_per_cm']:.9e} /cm |
| 최대 proper 중력 가속도 변분 | {audit['maximum_delta_proper_acceleration_cm_s2']:.9e} cm/s^2 |
| 전체 저장 scalar 장의64/128 상대 차이 | {fixed['controls']['time']:.9e} |
| 수정 반경 구적4/8 상대 차이 | {fixed['controls']['quadrature']:.9e} |

분류: Counterexample candidate. 질량 변화는 배경 binary64 ulp의 약{fine['maximum_mass_increment_over_background_ulp']:.6e}, 스칼라 변화는 약{fine['maximum_scalar_increment_over_background_ulp']:.6e}다. 배경 배열에 직접 더해 진화기를 재실행하면 이 신호를 표현하지 못한다. 별도 저장한 field/center-constraints 배열이 다음 되먹임의 입력이다. 작다는 사실만으로 물리적 영향이0이거나 최종 전하가 닫혔다고 판정하지 않는다.

분류: Counterexample candidate. 표현한 compact 영역의 퍼텐셜 연산자 norm 추정은{fine['compact_potential_contraction']:.9e}, 한 번의 되먹임에 따른 상대 변화는{fine['potential_iteration_relative'][0]:.9e}다. 공간의 셀별 상수 장·저장 시각의 선형 이력에 대한 수치다. 이 작은 값은 바깥 진공/깊은 층의 전체 연산자 상계나 물질·광자 수송의 안정성 보증을 대신하지 않는다.

## 원 실패와 적분 수정

분류: Counterexample candidate. 첫 fields.json은 반경 구적 차이{raw['controls']['quadrature']:.9e}가 사전0.002를 넘어 passed=false다. 원 코드·계획·결과를 수정하지 않았다. 광원뿔의 절댓값 거리 cusp와 저장 시각의 retarded knot를 기존 셀 안에서 분할해 적분했다. 셀·시간·주파수·각도 수를 늘리지 않았고 원4/8 공간 보간 차수를 유지했다. 각 매끈한 조각의 다항식 차수에 필요한 Gauss 적분만 수행한다.

분류: Proven. 상자 안의 선형 시간 원천에 대해 관측점과 파면이 셀 내부에 있는 해석해를 사용한 최소 실행 검사의 절대 차이는{max(box['exact_box_source_absolute_errors']):.9e}다. 다항식 특성 적분의 확인이며 GR 연속해 전체의 오차 정리가 아니다.

분류: Counterexample candidate. 수정 characteristic/fields.json 역시 passed=false다. 구적 문제는 해결됐지만 기존 미분할 직접 판독과의 재현 차이{fixed['controls']['direct_reproduction']:.9e}가 원1e-9 호환 기준을 넘기 때문이다. 기존 숫자와의 호환 실패를 삭제하거나 기준을 완화하지 않았다. 이후 audit-plan.json을 별도로 등록해 광학 좌표의 원천 보간을 사용하지 않는 독립 Jordan 반경 적분으로 같은 물리량을 계산했다. 원1e-9 동등성 기준에서 실제 차이{audit['independent_equivalence_relative']:.9e}로 audit.json은passed=true다. 이는 독립 정확도 대조 통과이며 앞선 두 raw verdict의 사후 통과 전환이 아니다.

분류: Counterexample candidate. J를 Gauss 지점들 사이에서 선형 보간했던 최초 center 질량은 전체 J norm 대비{audit['exact_center_J_vs_old_interpolation_relative']:.9e} 차이가 있었다. center-constraints.npz에서는 각 중심까지의 proper 체적·lapse 적분을 직접 수행하고 바깥 면에는 전체 적분값을 사용한다. 바깥 압력은 진공 쪽 값을 사용하며 마지막 셀의 압력을 복사하지 않는다. 이에 따른 질량·체적·압력·중력의 변화량을 저장했다.

원천 추출{source['seconds']:.2f}s, 최초 장 계산{raw['seconds']:.2f}s, 특성 수정 세 경로{fixed['seconds']:.2f}s, 독립 audit{audit['seconds']:.2f}s였다. 각각 사전25/90/90/30s 한도 안이다. 원 symbolic 및 해석해 확인의 준비 시간과 WSL 시작 시간은 이 계산 시간과 구분한다.

## 남은 실제 병목

분류: Counterexample candidate. represented_anisotropic_scalar_mass_response_evolved와independent_characteristic_accuracy_audit는true다. legacy_compatibility, full_spatial_material_photon_feedback, exterior_deep_source_enclosure, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 전하 끝점은 표현한 원천·선형 GR 폐쇄에서 양수이며, 외부 광자 질량 정규화까지 합친 최종 전하가 아니다. 이전 배경의 양수 GR 하한도 승계하지 않는다.

분류: Conjectural. 다음 구현은 이 별도 변화량을 lapse 경계까지 연결한 뒤 물질 유속·압력 및 광자 공간/주파수 수송에 일관되게 반영해야 한다. 기존 배경에 작은 수를 더하는 방식은 사용할 수 없다. 보상 변수의 실제 수송 되먹임 또는 그 전체 연산자에 대한 오차 상계가 필요하다. 또 표현하지 않은 외부·깊은 층의 인과적 원천과 퍼텐셜을 새 비등방 배경에서 닫아야 한다. 이번 장 진화를 전체 결합 완성으로 부르지 않는다.

근거: outputs/direct-eos-gr33/def-native-anisotropic-gr의 원본 및characteristic 하위 결과, audit와center-constraints, verification/def_native_anisotropic_gr.py 및def_native_characteristic_gr.py, verify_native_anisotropic_gr.py.
''',encoding='utf-8')

docs=[root/'docs'/n for n in ['model-definition.md','observable-targets.md','adiabatic-limit.md','nonadiabatic-regime.md','failure-ledger-dynamic-chi.md','dynamic-charge-completion.md']]
assert not (out/'documentation-prefixes.json').exists()
write(out/'documentation-prefixes.json',{str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in docs})
paragraphs={
'model-definition.md':'분류: Proven. 고정 canonical 광자 운동량에서 delta_Er=eFr-4*alpha*Er*f-(Er+Pr)*dl, delta_Pr=pFr-4*alpha*Pr*f-(3*Pr-R4)*dl이다. 실제 비등방 반경 압력의4차 모멘트를 포함해야 보존 체적과 GR scalar 변분을 연결한다. 이 폐쇄에서 delta_m=r^2*b*Phi*f+J의 질량 상쇄식이 유지된다.',
'observable-targets.md':f"분류: Counterexample candidate. 특성 경계를 분할한 새 직접 끝점은{fine['endpoint_direct']:+.8e}, 표현한 비등방 GR 항 포함 끝점은{fine['endpoint_compact_with_metric']:+.8e}다. 독립 Jordan 반경 적분과{audit['independent_equivalence_relative']:.3e}로 일치한다. 원 호환 실패는 보존하며 최종 전하·정적 EFT 이탈·관측 검출은 미확정이다.",
'adiabatic-limit.md':'분류: Proven. 비등방 광자의 canonical momentum 고정 변분은 반경·접선 응력을 함께 변화시키며 trace=0을 보존한다.\n\n분류: Counterexample candidate. 이 순간 계량 부분 단계의 보존 폐쇄는 모든 원천의 단열성이나 기존 정적 계수 흡수 경계의 이탈을 뜻하지 않는다.',
'nonadiabatic-regime.md':'분류: Counterexample candidate. 실제64/128 물질·광자 원천에서 retarded scalar 파동과 질량 제약의 별도 변화량을 계산했다. compact 퍼텐셜 되먹임은 포함했으나, 물질·광자 공간 수송의 동적 계량 되먹임과 미표현 영역의 원천은 아직 포함하지 않았다.',
'failure-ledger-dynamic-chi.md':'분류: Counterexample candidate. 첫 장 계산은 광원뿔/시간 knot를 셀 안에서 분할하지 않아0.002 구적 기준에 실패했다. 분할 수정은 구적 대조를 해결했지만 기존 미분할 전하와의1e-9 호환 기준은 실패했다. 두 원 verdict를 보존하고 별도 등록한 독립 Jordan 적분 정확도 audit만 통과로 기록한다. 배경 ulp 아래의 변화량을 배경에 단순 합산하는 경로는 사용하지 않는다.',
'dynamic-charge-completion.md':'분류: Counterexample candidate. represented_anisotropic_GR_scalar_response_evolved, independent_characteristic_accuracy_audit는true다. 원 raw 필드/호환 verdict는false다. full_spatial_material_photon_feedback, exterior_deep_source_enclosure, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 다음은 별도 변화량의 실제 수송 되먹임과 lapse/인과적 외부 연결이다.'}
for p in docs:p.write_bytes(p.read_bytes()+('\n\n## 단계122 — 비등방 GR·스칼라 변화량 진화\n\n'+paragraphs[p.name]+'\n\n상세: [단계122 보고](../notes/REQUEST122_ANISOTROPIC_DYNAMIC_GR_KO.md).\n').encode('utf-8'))
paths=docs+[report]+[root/'verification'/n for n in ['def_native_anisotropic_gr.py','def_native_characteristic_gr.py','verify_native_anisotropic_gr.py']]
paths+=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
files={str(p.relative_to(root)):sha(p) for p in paths};manifest=root/'outputs/direct-eos-gr33/native-anisotropic-gr-manifest.json'
write(manifest,dict(classification='Counterexample candidate',checkpoint='c57efa4b6',represented_anisotropic_fields_evolved=True,
    raw_fields_passed=False,legacy_compatibility_passed=False,independent_accuracy_audit_passed=True,
    spatial_material_photon_feedback=False,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False,files=files))
master=root/'paper/revision-manifest.json';d=read(master);d['sha256'].update(files);d['sha256'][str(manifest.relative_to(root))]=sha(manifest)
d['native_anisotropic_gr']=dict(classification='Counterexample candidate',report=str(report.relative_to(root)),manifest=str(manifest.relative_to(root)),
    represented_anisotropic_fields_evolved=True,raw_fields_passed=False,legacy_compatibility_passed=False,independent_accuracy_audit_passed=True,
    spatial_material_photon_feedback=False,full_GR_feedback=False,final_charge_solved=False,full_goal_complete=False)
write(master,d)
for path,h in d['sha256'].items():assert sha(root/path)==h,path
for path,p in read(out/'documentation-prefixes.json').items():assert hashlib.sha256((root/path).read_bytes()[:p['bytes']]).hexdigest()==p['sha256'],path
print(json.dumps(dict(manifest_files=len(files),master_files=len(d['sha256']),independent_audit_passed=True,preserved_failed_raw_verdicts=2,passed=True)))
