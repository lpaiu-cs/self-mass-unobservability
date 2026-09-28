# 단계117 — 실제 물질·광자 원천의 GR 기여와 조건부 잔여

분류: Counterexample candidate. 단계116에서 완주한 양방향 물질·광자 이력에 Einstein 질량 제약과 계량을 통한 scalar 원천을 적용했다. 광자를 trace0이라고 생략하지 않고 실제 에너지와 반경 압력을 사용했다. 추가 선형 GR 항은 현재 직접 전하의 양의 부호를 뒤집지 못한다. 이로써 새 저장 이력에 대한 해당 누락항의 크기를 제한했으며, **전체 물리 전하의 완료나 새로운 관측 검출을 뜻하지 않는다.**

분류: Counterexample candidate. 아래 구간은 명시한 추가 GR 항에 대한 조건부 구간이다. 원천의 공간·주파수·내부 각도 오차, EOS 구성 관계 절단, 초기 완전한 Einstein 제약 및 동적 계량에 따른 물질·광자의 재진화 오차를 포함하지 않는다.

| 같은 관측 끝점의 정규화 전하 성분 | 값 |
|---|---:|
| 단계116 직접 성분 재현 | +1.1439211186e-27 |
| 실제 물질+광자의 계량 응력 기여 | −1.46401978e-32 |
| 실제 에너지 질량 제약 기여 | −1.39861290e-34 |
| 위 세 성분의 합 | +1.1439063385e-27 |
| 외부 광자 응력 기여의 절대 상계 |1.03652051e-30 |
| 외부 광자 질량 제약 기여의 절대 상계 |9.83237704e-33 |
| scalar 퍼텐셜의 모든 반복 차수 상계 |5.52017150e-32 |
| 저장 경계 유속 이력의 영향 상계 |1.46585608e-35 |
| 단일 끝점의 미계산 GR 항 합 상계 |1.10156926e-30 |

분류: Counterexample candidate. 끝점 하나의 직접+GR 구간은 `[1.14280477,1.14500791]×10^-27`다. 실제 보고하는 관측 변화는 `u=0`과 끝점의 차이이므로 두 사건의 미계산 항을 각각 허용해 **두 배 상계**를 적용했다. 계산된 계량 항까지 합친 관측 변화의 GR 영향 상계는 직접 성분의 **0.19389%**다. 기존 빈 내부 상수 출사 광자 재구성에서 직접+광자 질량 합의 중심은 `5.29686546e-27`, 추가 GR 항만의 구간은 `[5.29466233,5.29906861]×10^-27`이다. 이를 전체 수치·물리 오차막대로 표시하지 않는다.

분류: Counterexample candidate. 고정된 실제 내부 이력과 총출사 에너지를 유지하면서 외부 각도 분포를 임의의 비음수 외향 분포로 허용해도, 외부 광자 계량 기여의 상계를 포함한 관측 변화의 직접 하한은 `1.14170320e-27`로 양수다. 배경 정규화 전하가 양수이므로 양의 광자 질량 손실 항은 이 하한을 지우지 못한다. 이것은 **외부 자유 전파의 각도 분포**에 대한 결론이며, 내부 광자 각도 해상도를 바꾸어 진화 원천 자체가 달라지는 경우까지 포함하지 않는다.

분류: Proven. 고정된 정적·등방 배경의 선형 polar-areal 제약에서 `f=delta phi`, `Phi=phi_prime`, `b=1−2m/r`라 쓰고, 추가 계량 변화 중 좌표 바리온 재고·엔트로피·조성을 고정하면

```text
delta m = r² b Phi f + J
delta lnV = (3 alpha+r Phi) f + J/(r b)
delta E = E_forced − (E+P) delta lnV
delta P = P_forced − Gamma1 P delta lnV
J_prime + (nu_prime+lambda_prime) J = 4 pi r² A⁴ E_forced
```

가 된다. 따라서 `J`는 실제 누적 Killing 에너지에 `sqrt(b)/N`을 곱해 복원된다. 본문의 식은 광자·물질을 계량 변화와 함께 다시 수송했다는 뜻이 아니다. 해당 보존 체적 tangent 폐쇄를 적용한 추가 응답 식이다. 유도와 기존 symbolic 검사를 재사용했다.

분류: Counterexample candidate. 질량 제약의 안쪽 경계는 실제 저장된 광자 유속의 시간 적분을 음의 질량 변화로 반영했다. 깊은 영역에서 나온 광자 에너지를 외부 공급 에너지로 남겨 놓지 않는다. fine 경로의 안쪽 누적 광자 에너지는 `7.80603e32erg`, 바깥으로 나간 에너지는 `2.75889e31erg`다. 실제 저장된 영역 에너지와 두 유속의 별도 대조 차이는0.15397%이며 원2% 기준 안이다. 이 차이를 맞추려고 보상 상수를 넣지 않았고, 차이의 GR 영향은 위 상계에 포함했다. 누락된 깊은 내부 경계의 직접 scalar 도착 시간은3.66924ms여서 이번3.43043ms 관측 끝점 뒤다.

분류: Proven. 진공 외향 null ray에서 `r nu_prime≤kappa<1`이면, 접선에 가까운 출사까지

```text
mu(r)² ≥ (1−kappa) [1−r0²/r²]
integral_r0^infinity dr / [r² sqrt(1−r0²/r²)] = pi/(2 r0)
```

로 제한된다. 실제 양의 누적 출사 에너지와 이 식으로 외부 광자의 `E−Pr` 응력 및 질량 제약 원천을 제한했다. 인위적인 외부 반경 절단이나 최소 출사 각도는 필요하지 않다. 배경 계수와 광자 보존 전제는 유지해야 한다.

분류: Proven. `U_tt/c²−U_xx+Veff U=S`의 지연 적분 연산자 노름이 `eta=(c T/2) integral |Veff| dx<1`이면, 절대 자유 원천 노름 `M0`에 대해 모든 퍼텐셜 반복의 차이는 `eta M0/(1−eta)` 이하이다. 이미 상쇄된 관측 전하를 `M0`로 사용하면 안 된다.

분류: Counterexample candidate. 이번 frozen 배경·저장 이력에서 계수 극값을 독립적으로 허용한 보수적 계산은 `eta≤3.30286e-8`, `M0≤4.87505e-20cm`를 주었다. 입력 binary64 합산에1e-9 여유를 둔 뒤 상계 산술은 바깥 방향 구간 연산을 사용했다. 샘플된 경계 에너지 차이를 연속 상계라고 오해하지 않도록 각 저장 구간의 이차식 내부 극점까지 포함했다. 원 source 함수나 EOS의 연속 오차를 인증한 것은 아니다.

분류: Counterexample candidate. 검증은 다음과 같다.

| 검사 | 결과 |
|---|---:|
| 이전 직접 전하 재현 상대차 |1.568e-15 |
|64/128 직접+영역 내 GR 성분 시간 차이 |0.94265% |
| 선형/3차 이력 보간 차이 |0.56366% |
| 영역 내4/8점 반경 구적 차이 |0.00479% |
| 해석적으로 알려진 지연 ramp 적분 차이 |5.552e-16 |
| 물질·질량·scalar 배경0일 때 GR 계수0 환원 |통과 |
| 접선 광선 적분의 해석식과 상계 대조 |통과 |
| 비영 퍼텐셜의 알려진 resolvent 나머지 대조 |통과 |
| 원 생산 소스·데이터 해시 및 경계 이력 상계 대조 |통과 |

실행: 원 유체·광자 경로, native 표, 반경/각도/주파수 해상도 및 시간 구간을 재계산하거나 늘리지 않았다. 새 원천 추출·GR 읽기·구간 계산은 각25/45/30초 예산 안에서 수행했고 GR 읽기는3.718초, 수정 상계 계산은2.509초, 독립 감사는1.196초였다. 구성자/import/WSL 시작·해시 등의 시간은 별도다. 첫 상계는 보고용 NumPy bool 직렬화에서 실패했으며 해당 소스를 보존했다. 보고 형 변환을 고치고 저장 원천을 재사용했다. 수락 기준은 완화하지 않았다.

분류: Conjectural. 이 정도의 선형 GR 누락항을 더 정밀하게 구적하는 일은 현재 우선순위가 낮다. 다음 물리 병목은 지정된 물질 경계를 실제 내부·대기 상태가 함께 결정하도록 연결하고, 초기 Einstein 제약 및 원16개 내부 셀의 공간 오차를 확인하는 것이다. 이번 frozen 원천의 양수 판정을 실제 항성의 최종 전하로 옮기려면 이 원천 자체의 변화가 제한돼야 한다. 전체 중성수소 누적 원장, 내부 각도·주파수 및 전체 반응 채널도 남는다.

분류: Counterexample candidate. `actual_matter_and_photon_metric_source_applied`, `actual_inner_energy_debit_applied`, `declared_all_orders_scalar_potential_bound`, `external_outward_angle_independent_GR_bound`는true다. `initial_full_Einstein_constraints_matched`, `free_mechanical_interface`, `source_continuum_error_certified`, `full_dynamic_GR_feedback`, `final_charge_solved`, `full_goal_complete`는false다. 이 단계는loophole progress와 명시한 조건의theorem progress이며 전체 목표는active다.

재현 파일은 `verification/def_native_feedback_gr.py`, `verification/verify_native_feedback_gr.py`다. 계획·입력·GR 성분·상계·감사는 `outputs/direct-eos-gr33/def-native-feedback-gr/`에 보존했다. 기존 출력에 덮어쓰려고 명령을 반복하지 않는다.
