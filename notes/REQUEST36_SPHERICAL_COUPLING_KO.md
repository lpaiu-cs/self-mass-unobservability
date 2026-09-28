# 단계 36 — 전 구역 scalar–물질–계량 동시 계산

분류: Counterexample candidate. 새 분자 EOS의 5,735셀, 26핵종 이류, 두 열수송 채널을 구면 DEF scalar 파동과 같은 비선형 잔차에 연결했다. scalar와 방사 방향 Einstein 제약은 물질 잔차를 평가할 때마다 함께 푼다. 단방향 후처리나 고정 계량의 scalar source 실험으로 대체하지 않았다. 이 문서의 수치 결과는 유한 격자 계산이며 연속 이론의 인증이 아니다.

## 고정한 실험

분류: Counterexample candidate. 저장된 분자 EOS 초기상태에서 `phi=Pi=Phi=0`으로 시작한다. 두 경로는 같은 초기상태, 같은 시간 간격 `h=3.2132813194763797e-6 s`를 쓴다. 물질 경계는 기존 반사벽이고, scalar 외곽 경계는 `phi_b(t)=amplitude*t/h`다. 진폭 0과 1e-6을 비교한다. 이는 외부 파동을 주는 수치 대조이며 동반천체·궤도·관측 자료에 맞춘 구동이 아니다. `c*h`는 약 9.63e4 cm로 외곽 셀 폭 3.69e6 cm보다 작다. 이 첫 단계의 경계 응답으로 전 항성 전달함수나 물리적 온도 변화를 추론하지 않는다.

분류: Counterexample candidate. 기존 물질의 31개 native 잔차, 바리온·핵종·국소 에너지 수지, 격자점 특성속도 문턱을 유지했다. scalar 잔차는 경계 진폭 정규화에서 2e-17, scalar 상대 에너지 수지는 2e-15다. 총 에너지 수지는 물질 열용량 기준 1e-9와 scalar 에너지 규모의 2e-15를 합친 허용량으로 따로 평가한다. scalar의 큰 에너지 때문에 직접 총량 차분만으로 작은 물질 열 오차를 보증하지 않는다.

분류: Counterexample candidate. 전 단계의 native 단일 GR 실행 156.5초를 예산 근거로 썼다. CPU 8개·BLAS 1 thread, 모든 본 계산의 합계 벽시간 상한 1,200초를 정했다. 동일 무구동 결과를 수정 경로에서 재계산하지 않으며, 각 수정의 원인·결정·중단 기준을 별도 계획으로 고정했다. 기간·격자·진폭은 자동 확대하지 않았다. 각 계획과 실행 소스는 실행 중 바꾸지 않았다.

## 방정식과 명시적인 이산화 변경

분류: Proven. `A=exp(-2 phi^2)`, `rho_J=rho_E/A^3`, `T_J=T_E/A`를 사용한다. 물질 native EOS와 불투명도는 Jordan 변수로 호출하고, `P_E=A^4 P_J`, `u_E=A u_J`, `rest_E=A rest_J`, `kappa_E=kappa_J/A^2`, `tau_cond,E=tau_cond,J/A`로 변환한다. 이 항목은 앞 단계에서 유도한 frame 관계를 구현한 것이다.

분류: Counterexample candidate. `Pi=a phi_t/(Nc)`는 셀 중심, `Phi=phi_r`는 셀 경계에 둔다. 경계 기울기 에너지는 내부 면의 두 셀에 절반씩, 외곽 면은 마지막 셀에 전부 배분한다. 물질에는 backward Euler, scalar에는 midpoint를 쓴다. 총 에너지로 질량 제약을 푼다. lapse는 `log(N*a)`의 구간별 상수 적분을 사용한다. 이는 앞 단계의 lapse 구적과 다른 공간 이산화다. 초기 무 scalar 상태의 압력/중력 보정만 고정하며, 새 scalar 힘을 평형 보정으로 제거하지 않는다.

분류: Proven. 연속 방정식의 scalar가 받는 물질 교환은 직접항 `alpha*T_trace*phi_t`와 기하 교차항 `4*pi*g*c*r*N*a*[2*E_phi*S_m-S_phi*(E_m+R_m)]`이며, `g=G/c^4`다. 물질은 반대 부호를 받는다. 진공에서는 두 물질 교환항 모두 0이다. 따라서 진공에서 남는 scalar 수치 에너지 오차를 물질 source로 옮기는 방법은 물리적 교환항의 구현이 아니다.

## 보존한 실패와 원인

분류: Counterexample candidate. 최초 구현은 scalar의 유한 에너지 곱셈식에서 남는 기하 항을 전부 물질로 넘겼다. 무구동은 3개 반복 기록, native 정규화 잔차 0.000441744로 통과했다. 구동은 잔차 7.41069e18에서 24회 한도를 소진했다. 소요 시간 393.54초다. 실패 상태와 source는 `def-spherical-coupled/` 및 `verification/def_spherical_coupled.py`에 보존한다.

분류: Counterexample candidate. 4셀 진공 대조에서 실제 물질 교환은 0인데 최초 방법은 최대 1.19514e12 erg/cm^3의 가짜 물질 에너지 증가량을 만들었다. 해당 방법 자체의 scalar 수지 잔차는 4.87e-20로 작았다. **수지의 상쇄만으로 물리적인 상호 교환이 증명되지 않는 반례**다. 전체 격자의 실패 상태에서는 이 가짜 항이 초기 열용량의 약 7.41069e9배였다.

분류: Counterexample candidate. 물리적인 교차항만 사용한 별도 수정에서 전 구역 native 잔차는 8개 반복 기록으로 0.000355360까지 내려갔다. 물질 에너지 잔차 3.554e-13, 바리온 잔차 3.929e-24, 최대 특성속도 0.407216c였다. 그러나 scalar 상대 에너지 수지는 8.36449e-10으로 문턱 2e-15에 미달했다. 이 경로의 최종 판정은 `passed=false`로 유지한다. 실행은 262.51초였고, 이 결과를 결합 진화 인증으로 승격하지 않는다.

## scalar 쪽에서 처리하는 수치 에너지 보정

분류: Proven. 유한 곱셈식의 scalar 기하 잔여항을 `q_disc`, 물리 기하 교환항을 `q_phys`, `b=1/a^2`, `k=N/a`, midpoint scalar 운동량을 `p`라 하자. scalar midpoint 방정식 우변에

```text
delta_rhs = h^2 * pi * g * c * k * (q_phys-q_disc)/(b*p)
```

를 더하면 그 추가 운동량 변화가 scalar 에너지식에 공급하는 항은 정확히 `q_phys-q_disc`다. 이 대수 항등식은 symbolic check로 검증했다. 물질 교환항은 바꾸지 않는다.

분류: Counterexample candidate. 이 방법은 **에너지 보정을 포함한 새로운 scalar 시간 갱신식**이다. 원 midpoint 방정식을 그대로 풀었다고 주장하지 않는다. 수정 식의 잔차와 함께 원 midpoint 식의 잔차 및 보정 크기를 별도로 저장한다. 진공 대조의 물질 교환은 정확히 0이고 상대 scalar 수지는 5.711e-20이었다. native 물질 반복의 접선에도 물리 교환항의 물질 의존성을 포함했다.

분류: Conjectural. 이 에너지 보정의 일반적인 연속 방정식 합치성·시간/공간 수렴은 아직 증명하지 않았다. 특히 `Pi_mid=0`인 비정상 회전점에서 나눗셈의 정칙성을 보장하지 않는다. 구현은 아주 작은 운동량에서 보정을 0으로 두므로 이 경계도 검증해야 한다. 단조 경계 ramp의 한 단계 성공을 반복 궤도 구동이나 일반 파형의 성공으로 확장하지 않는다.

## 최종 수치 판정

분류: Counterexample candidate. **수정된 이산 모형의 전 구역 동시 결합 첫 단계를 수락했다.** 아래 모든 문턱과 24개 반복 기록 한도를 유지했으며 실제 사용은 8개 기록이었다. 단순한 source 분할이나 한 방향 계산이 아니라 native 물질, scalar, 질량 및 lapse가 같은 수렴 상태를 이룬 결과다.

| 항목 | 측정값 | 문턱 |
|---|---:|---:|
| 물질 native 정규화 잔차 | 3.65633e-4 | 1 |
| 국소 물질 에너지 잔차 / 초기 열용량 | 3.65633e-13 | 1e-6 |
| 국소 바리온 잔차 | 4.79335e-24 | 1e-9 |
| 핵종 재고 잔차 / 초기 총 바리온 | 4.28199e-26 | 1e-9 |
| 수정 scalar 갱신식 잔차 / 경계 진폭 | 5.59436e-25 | 2e-17 |
| scalar 상대 에너지 수지 | 1.32535e-19 | 2e-15 |
| 총 에너지 오차 / 등록 허용량 | 3.65633e-4 | 1 |
| 최대 격자점 물질 특성속도 | 0.407216c | c 미만, 허수부 1e-10 미만 |

분류: Counterexample candidate. 위 보존 문턱은 상속한 계획의 별도 예산 검사다. native 방정식 검사는 에너지 성분 절대 허용오차 1e-9, 바리온 성분 1e-13 등 더 엄격한 원 ATOL로 정규화하며 함께 통과했다.

분류: Counterexample candidate. 원 midpoint 갱신식의 정규화 잔차와 수치 보정 크기는 각각 8.36637e-10이다. **원 midpoint 식의 통과가 아니다.** scalar 에너지 총량은 1.39381e50 erg, 물질이 받은 직접 교환의 적분은 -6.90273e20 erg, 기하 교환 적분은 -9.44377e26 erg였다. 이들은 선언한 외곽 수치 구동의 값이며 실제 계의 물리적 규모로 채택하지 않는다.

분류: Counterexample candidate. 무구동 끝점과 비교한 최대 변화는 `delta ln rho=7.41089e-10`, `delta ln T=0.0842616`, `delta(v/c)=2.90822e-11`이었다. 수정 전 물리 교환 경로와의 온도 차이는 최대값 비교에서 약 1.51e-10이다. 이 결과는 유한 이산 모형의 응답이며, 미해상 경계 파동·공간 이산화·수치 보정의 오차를 포함할 수 있어 실제 항성의 온도 응답으로 주장하지 않는다.

분류: Counterexample candidate. 무구동·수정 구동 두 저장 상태의 native 호출 인자 키를 정확히 재사용하여 모든 물질 및 scalar 잔차를 독립 재생했다. 밀도·온도·EOS·질량·계량·열유속·보존 증가량·scalar 변수와 에너지의 재구성 배열이 모두 저장 배열과 정확히 일치했다. 새 EOS 호출은 0회였다. 이는 저장 상태 재생이며 독립된 EOS 제공자의 인증은 아니다.

분류: Counterexample candidate. 최종 수정은 native material 호출 31,439회, 벽시간 330.23초였다. 최초 실패와 물리 교환 수정까지 합한 측정 벽시간은 **986.28초(약 16분 26초)**로 등록한 1,200초 안에 끝났다. 독립 재생은 약 6초였고 장기 작업은 남아 있지 않다.

## 완료 경계와 다음 목표

분류: Counterexample candidate. 이번 성과는 **native EOS–scalar–물질–열수송–계량을 같은 전 구역 이산 진화에 연결한 loophole progress**다. 최초 방법의 가짜 열원을 제거하고, scalar 에너지 보정을 scalar 갱신식 안으로 옮겨 유한 단계의 수락 기준과 저장 재생을 충족했다. 최초 반복 실패와 물리 교환만 적용한 에너지 실패는 그대로 보존했다.

분류: Conjectural. 다음 우선 병목은 **운동량 0의 회전점에서도 정칙하고 연속 방정식에 합치하는 scalar–계량 보존 이산화**다. 현재 `1/Pi_mid` 보정의 일반 정칙성, 시간·공간 수렴, Einstein 운동량/각방향 제약과의 일치가 먼저 필요하다. 이 조건을 작은 진공·매끄러운 파동 대조에서 해결하기 전에는 전체 native 장시간 적분을 확대하지 않는다. 그다음에 해상 가능한 외부 경계 및 동반천체 구동을 맞추고 정적 비교를 제거한 잔여 응답과 관측 추론을 연결한다.

실행 근거: `verification/def_spherical_coupled.py`, `def_spherical_physical_exchange.py`, `def_spherical_balanced_wave.py`, `def_spherical_replay.py`; `outputs/direct-eos-gr33/def-spherical-{coupled,physical-exchange,balanced-wave}/`의 계획·실패·결과·저장 상태와 최종 `replay.json`. 통합 결속은 `gr-spherical-coupling-milestone-manifest.json`에 둔다.
