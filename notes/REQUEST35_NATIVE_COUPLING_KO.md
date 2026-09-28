# 단계 35 — native EOS와 물질–scalar 결합

분류: Counterexample candidate. 실제 분자 EOS의 Einstein/Jordan frame 변환과 상호적 에너지 교환을 구현했다. 대표 5상태에서 16/32/64단계 구동 경로, 무구동 경로, 부호 반전 경로를 실행했다. 약 2차 시간 수렴, native 엔트로피 잔차, 총 에너지–외부 일 수지 및 부호 대칭 기준을 통과했다. 실행은 CPU 1개로 약 60초였다.

분류: Proven. scalar 응력과 직접적인 물질 교환항 외에 두 성분의 에너지 유속이 다를 때의 계량 교차항도 필요하다. 연속 방정식에서 두 성분의 직접 교환과 계량 교차항이 각각 상쇄됨을 symbolic check로 확인했다. 완전한 quadratic-entropy 열수송식의 conformal 변환에서는 `tau_E=tau_J/A`, `K_E=A^2 K_J`, `q_E=A^4 q_J`가 필요하다.

분류: Conjectural. 이 대표 상태 실험은 고정 부피·밀도·조성의 가역 source substep이다. 임의의 scalar spring과 시간 단위를 쓴다. 열수송·구면 공간 결합·계량 진화·실제 외부 구동은 포함하지 않으므로 실제 항성의 열 모드 또는 관측 신호로 해석하지 않는다.

분류: Counterexample candidate. 새 EOS 보존 초기화의 기존 고정밀 적분은 복구본에서 1,920/5,735개 셀까지 저장되어 있다. 저장된 독립적인 GR shell mass 증가량과 원 바리온 재고로 finite-volume 초기상태를 구성하는 별도 방법을 실행하여, **5,735개 셀 전체 초기화가 완료됐다.** 과거 8/16-node 적분과 실패·중단 기록은 보존한다. 기존 보존 및 native 역산 수락 기준을 바꾸지 않았다.

실행 근거: `verification/def_native_coupling.py`, `verification/def_native_source_step.py`, `outputs/direct-eos-gr33/def-native-source-step/plan.json`, `result.json`, `manifest.json`.

## 실제로 해소한 병목

분류: Proven. 고정된 finite-volume 바리온 및 좌표 에너지 모멘트 `B_i,U_i`가 주어지면, `E_i=U_i/V_i`로 질량 제약을 구성하고 `rho_i=B_i/(a_i V_i)`를 정한 뒤 native EOS의 `u(rho_i,T_i,X_i)=E_i/rho_i-C_X c^2`를 풀어 두 보존량에 맞는 원시 상태를 얻는다. 양의 열용량과 실제 역산 해의 존재가 필요하다. 이는 초기화 대수식이며 물리적 EOS 정확성의 정리가 아니다.

분류: Counterexample candidate. 이번 방법은 기존 RK4의 독립적으로 누적된 각 껍질 질량 증가량으로 `U_i=c^4 delta_m_i/G`를 정한다. 누적 질량의 근접 차분이나 목표 질량으로의 재조정은 하지 않았다. 기존 구적 모멘트 대신 **초기 finite-volume 모형의 정의를 명시적으로 바꾼 것**이다. 기존 8/16-node native 구적 통과, 셀 내부 엔트로피 보존 또는 연속 공간 정확도를 얻었다고 주장하지 않는다.

분류: Counterexample candidate. 아래는 같은 기존 문턱으로 실제 전체 native EOS 역산을 수락한 결과다.

| 항목 | 측정값 | 기존 문턱 |
|---|---:|---:|
| 최대 셀 바리온 상대 오차 | 1.247e-18 | 1e-14 |
| 최대 셀 에너지 상대 오차 | 4.377e-14 | 1e-13 |
| 최대 native 열 역산 잔차 | 9.957e-9 | 1e-8 |
| 전체 원 바리온 상대 오차 | 0 | 1e-11 |
| 전체 원 질량 상대 오차 | 1.959e-14 | 2e-12 |

분류: Counterexample candidate. 대표 7셀은 native 호출 2–3회로 역산됐고 전체에서는 최대 4회였다. 7셀 실측에서 추정한 CPU 4개 벽시간은 120초, 실제 전체 조립은 135.4초(프로세스 벽시간 137.4초)였다. 파일을 저장한 후 프로세스가 정상 종료했다. 재구적의 고정밀 엔트로피 역산을 반복하지 않고, 이미 보유한 구조의 보존량으로 실제 진화 입력을 만들었다.

초기화 근거: `verification/gr_molecular_shell_initial.py`, `outputs/direct-eos-gr33/gr-molecular-shell-initial/{plan,restriction,runtime,initial-manifest}.json`, `initial.npz`.

## 물질–scalar source의 의미와 한계

분류: Proven. `A=exp(-2 phi^2)`, `rho_J=rho_E/A^3`, `T_J=T_E/A`, `P_E=A^4 P_J`, `epsilon_E=rho_E A(C_X c^2+u_J)`다. 정지에너지를 고정해 두고 힘만 추가하면 같은 frame 변환이 아니다. 고정 `rho_E,s,X`에서 `d(epsilon_E/rho_E)/dphi=-alpha T_trace/rho_E`이며 이 도함수와 실제 native 에너지 차분을 source 양쪽에서 공유했다.

분류: Counterexample candidate. 대표 상태의 16/32/64단계 끝점 scalar 변수 차이에서 관측한 수렴 차수는 1.98386이었다. 외부 일 대비 등록 진폭 에너지로 정규화한 총 수지 오차는 최대 3.428e-12로 사전 문턱 1e-9 미만이었다. 무구동 경로의 scalar 및 온도 변화는 정확히 0, 배경과 구동의 동시 부호 반전에서 scalar는 홀수·온도는 짝수 대칭을 저장 정밀도에서 만족했다. 이 값은 실제 항성의 응답 크기가 아니다.

분류: Proven. `Pi=a phi_t/(N c)`, `Phi=phi_r`일 때 `E_phi=R_phi=(Pi^2+Phi^2)/(8 pi (G/c^4) a^2)`, `S_phi=-Pi Phi/(4 pi (G/c^4) a^2)`다. 물질 에너지 source에는 `-alpha T_trace phi_t`와 `-4 pi (G/c^4)c r N a[(E_phi+R_phi)S_m-S_phi(E_m+R_m)]`가 모두 들어간다. 후자는 계량 교차 교환이며 scalar 에너지식에서 반대 부호다. 질량 제약에는 총 에너지, lapse 식에는 총 방사 방향 응력을 넣어야 한다.

## 새 EOS의 실제 구면 GR 첫 단계

분류: Counterexample candidate. 새 초기상태를 기존의 수정된 압력/중력 공간 연산자, 상호 일관된 donor 반복, native 분자 EOS와 불투명도, 26핵종 이류에 연결해 **전 5,735셀의 첫 GR 시간 단계를 실제로 수락했다.** 초기 backward-Euler 단계의 길이는 3.213281319e-6초다. 이전 EOS의 시간 상태를 가져오지 않았다. 24개 반복 기록 한도와 모든 원 native/보존/특성속도 문턱을 유지했다.

분류: Counterexample candidate. 반복 잔차는 780920.978 → 391.477 → 0.289357로 내려가 3개 기록 안에서 기준 1을 통과했다. 최대 국소 에너지 수지 잔차는 초기 열용량 정규화에서 2.577e-10, 최대 바리온 잔차는 5.270e-24, 최대 핵종 재고 잔차는 2.786e-26이었다. 격자점 최대 특성속도는 초기 0.356563c에서 끝점 0.356293c로 바뀌었고 허수부는 0이었다. 최대 log T 변화는 7.818e-4로 실제 상태 변화가 있었다. 12,905회의 native material 평가와 초기 준비를 포함한 수락까지의 벽시간은 CPU 8개에서 156.5초였다. 저장 파일의 독립 재생도 같은 잔차 0.289357, 상태 배열·누적 수지 배열의 정확한 일치로 통과했다. 재생을 포함한 최종 프로세스 벽시간은 207.8초로 10분 예산 안에서 정상 종료했다.

실행 근거: `verification/gr_molecular_compatible_step.py`, `outputs/direct-eos-gr33/gr-molecular-compatible-step/{plan,result,replay,manifest}.json`, `step-0000.npz`, `step-0001.npz`.

## 완료 경계와 다음 결정

분류: Counterexample candidate. 이번 단계는 **새 EOS 초기화와 검증된 구면 GR 계산기의 실제 연결**이라는 병목을 해소한 loophole progress다. 단순히 작은 인증 기록을 추가한 것이 아니라, 미완료였던 전체 초기상태를 생성하고 그 상태에서 원 방정식의 첫 진화를 수행했다.

분류: Conjectural. 새 EOS의 전체 시간 수렴과 scalar–물질–계량의 전 구역 동시 진화는 아직 완료하지 않았다. 대표 상태의 scalar source 대조와 이번 무 scalar GR 단계를 합쳐 이미 동시 진화했다고 주장하지 않는다. 다음 목표는 두 구현을 동일한 전 구역 잔차와 물리적 외부 구동에 연결한 짝 비교다. 기존 장시간 자유 열 이력을 자동 재실행하거나 더 촘촘한 경로를 추가하지 않는다.
