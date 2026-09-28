# 단계122 — 실제 비등방 원천의 GR·스칼라 변화량 진화

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

분류: Counterexample candidate. 원64/128 시간 경로,531개 물질/광자 셀,17개 저장 시점을 그대로 사용했다. 모든 양의 바깥 방향 광자와 실제 안쪽 유입의 Killing 에너지를 적분했다. 포트 이력 대조는 coarse/fine에서0.2189959%/0.1540199%로 원2percent 기준 안이다. 추가 EOS 호출·유체 단계는0이다.

| 분류: Counterexample candidate — 수정 특성 적분의 수치 결과 | 값 |
|---|---:|
| 직접 정규화 전하 끝점 | +1.021000070584e-27 |
| 표현한 영역의 GR 질량·응력·퍼텐셜을 포함한 끝점 | +1.020984763205e-27 |
| 위 GR 항의 직접값 대비 상대 변화 | -0.00149925% |
| 최대 delta_phi | 1.905957331e-31 |
| 셀 중심 직접 제약 적분의 최대 delta_m | 1.943689194e-16 cm |
| 최대 delta_lambda | 2.854127132e-26 |
| 최대 delta_nu_prime | 1.528336904e-34 /cm |
| 최대 proper 중력 가속도 변분 | 1.373597110e-13 cm/s^2 |
| 전체 저장 scalar 장의64/128 상대 차이 | 1.058750316e-03 |
| 수정 반경 구적4/8 상대 차이 | 1.447544486e-15 |

분류: Counterexample candidate. 질량 변화는 배경 binary64 ulp의 약5.342711e-05, 스칼라 변화는 약8.789677e-13다. 배경 배열에 직접 더해 진화기를 재실행하면 이 신호를 표현하지 못한다. 별도 저장한 field/center-constraints 배열이 다음 되먹임의 입력이다. 작다는 사실만으로 물리적 영향이0이거나 최종 전하가 닫혔다고 판정하지 않는다.

분류: Counterexample candidate. 표현한 compact 영역의 퍼텐셜 연산자 norm 추정은1.026861694e-09, 한 번의 되먹임에 따른 상대 변화는6.113168547e-11다. 공간의 셀별 상수 장·저장 시각의 선형 이력에 대한 수치다. 이 작은 값은 바깥 진공/깊은 층의 전체 연산자 상계나 물질·광자 수송의 안정성 보증을 대신하지 않는다.

## 원 실패와 적분 수정

분류: Counterexample candidate. 첫 fields.json은 반경 구적 차이4.966891360e-03가 사전0.002를 넘어 passed=false다. 원 코드·계획·결과를 수정하지 않았다. 광원뿔의 절댓값 거리 cusp와 저장 시각의 retarded knot를 기존 셀 안에서 분할해 적분했다. 셀·시간·주파수·각도 수를 늘리지 않았고 원4/8 공간 보간 차수를 유지했다. 각 매끈한 조각의 다항식 차수에 필요한 Gauss 적분만 수행한다.

분류: Proven. 상자 안의 선형 시간 원천에 대해 관측점과 파면이 셀 내부에 있는 해석해를 사용한 최소 실행 검사의 절대 차이는8.326672685e-17다. 다항식 특성 적분의 확인이며 GR 연속해 전체의 오차 정리가 아니다.

분류: Counterexample candidate. 수정 characteristic/fields.json 역시 passed=false다. 구적 문제는 해결됐지만 기존 미분할 직접 판독과의 재현 차이3.322330808e-07가 원1e-9 호환 기준을 넘기 때문이다. 기존 숫자와의 호환 실패를 삭제하거나 기준을 완화하지 않았다. 이후 audit-plan.json을 별도로 등록해 광학 좌표의 원천 보간을 사용하지 않는 독립 Jordan 반경 적분으로 같은 물리량을 계산했다. 원1e-9 동등성 기준에서 실제 차이1.809663530e-14로 audit.json은passed=true다. 이는 독립 정확도 대조 통과이며 앞선 두 raw verdict의 사후 통과 전환이 아니다.

분류: Counterexample candidate. J를 Gauss 지점들 사이에서 선형 보간했던 최초 center 질량은 전체 J norm 대비1.705938360e-05 차이가 있었다. center-constraints.npz에서는 각 중심까지의 proper 체적·lapse 적분을 직접 수행하고 바깥 면에는 전체 적분값을 사용한다. 바깥 압력은 진공 쪽 값을 사용하며 마지막 셀의 압력을 복사하지 않는다. 이에 따른 질량·체적·압력·중력의 변화량을 저장했다.

원천 추출4.98s, 최초 장 계산12.43s, 특성 수정 세 경로18.88s, 독립 audit5.32s였다. 각각 사전25/90/90/30s 한도 안이다. 원 symbolic 및 해석해 확인의 준비 시간과 WSL 시작 시간은 이 계산 시간과 구분한다.

## 남은 실제 병목

분류: Counterexample candidate. represented_anisotropic_scalar_mass_response_evolved와independent_characteristic_accuracy_audit는true다. legacy_compatibility, full_spatial_material_photon_feedback, exterior_deep_source_enclosure, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 전하 끝점은 표현한 원천·선형 GR 폐쇄에서 양수이며, 외부 광자 질량 정규화까지 합친 최종 전하가 아니다. 이전 배경의 양수 GR 하한도 승계하지 않는다.

분류: Conjectural. 다음 구현은 이 별도 변화량을 lapse 경계까지 연결한 뒤 물질 유속·압력 및 광자 공간/주파수 수송에 일관되게 반영해야 한다. 기존 배경에 작은 수를 더하는 방식은 사용할 수 없다. 보상 변수의 실제 수송 되먹임 또는 그 전체 연산자에 대한 오차 상계가 필요하다. 또 표현하지 않은 외부·깊은 층의 인과적 원천과 퍼텐셜을 새 비등방 배경에서 닫아야 한다. 이번 장 진화를 전체 결합 완성으로 부르지 않는다.

근거: outputs/direct-eos-gr33/def-native-anisotropic-gr의 원본 및characteristic 하위 결과, audit와center-constraints, verification/def_native_anisotropic_gr.py 및def_native_characteristic_gr.py, verify_native_anisotropic_gr.py.
