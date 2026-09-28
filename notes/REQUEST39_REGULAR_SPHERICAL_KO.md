# Request 39 — 정칙 scalar와 native 구면 열유체의 실제 조립

작성: 2026-09-17. 선행 결과: [Request 38](REQUEST38_HEAT_COUPLING_KO.md).

## 결론

분류: Counterexample candidate. **새 정칙 scalar 갱신식을 실제 5,735구역의 물질·두 열유속·26종 조성·방사 계량과 같은 잔차에 조립하고, 구동·무구동 첫 단계를 완료했다.** 역운동량 보정이나 scalar 수치 결함을 물질 열로 넘기는 항을 사용하지 않는다. 원 native 수락 문턱을 유지했고 새 방법에 맞춰 사전 등록한 질량·경계 에너지 수지를 통과했다. 전 구역 계산은 227.69초였다. 이는 유한 첫 단계 조립이라는 **loophole progress**다.

분류: Proven. 같은 방사 수반, scalar 공간 가중치, 물질 계량 일과 공유 면 수송을 쓰면 총질량 증가에서 내부 교환이 소거된다. 물질 에너지 잔차를 남겨 둔 정확한 유한 항등식을 기호 검산했다. 양의 계량·수반 가지와 선언한 단계 방정식에 조건부인 **theorem progress**다.

분류: Counterexample candidate. **이번 경계 첫 단계는 직접적인 trace 물질 피드백의 분리 검출에는 부족하다.** 해당 항을 저장 상태의 물질 에너지 잔차에서 제거해도 수락 판정이 바뀌지 않는다. 결합식을 실제 실행했다는 사실과 그 물리 효과를 수치적으로 식별했다는 주장을 구분한다.

분류: Conjectural. 새 결합법의 다단계 시간 수렴, 물질을 포함한 공간 수렴, 전체 Einstein 제약의 수렴, 해상된 외부 파동, 실제 동반천체 구동과 관측 추론은 남는다. 이전 GR 전체 기간 결과나 진공 파동의 2차 차수를 이번 새 모형에 이전하지 않는다.

## 같은 수반으로 연결한 식

분류: Imported from prior work. 기존 native 분자 EOS·불투명도·정지에너지와 conformal 변환 `A=exp(-2 phi²)`, 바리온 및 26종 이류, 두 완전한 열수송식, 호환 압력·음향 응력·정적 평형 기준을 재사용했다. 반응 네트워크의 시간 source를 새로 결합한 것은 아니다. 이번의 조성 진화는 보존 이류다.

분류: Proven. 물질 부피 `V=4π/3 Delta(r_face³)`와 기존 셀 안 질량 분율 `f`는 유지한다. scalar의 별도 구적은 Phase37의 정칙 파동 구조를 따라

```text
V_scalar_i = 4π r_i² Delta(r_face)_i
A_scalar_(i+1/2) = 4π r_i r_(i+1)
A_scalar_outer = 4π r_last R_outer,    A_scalar_center=0
w_i = V_scalar_i/(4π)
K_i = w_i Pi_i²/2
gradient_face = A_scalar_face distance_face Phi_face²/(8π)
```

로 정의한다. 내부 gradient energy는 이웃 두 셀에 반씩, 외곽 것은 마지막 셀에 전부 배분한다. 이 둘을 `L_i,Rg_i`라 한다. 물질과 scalar 부피가 같다고 가정하지 않는다. 이 구현은 처음 `phi=Pi=0`인 한 단계를 대상으로 한다.

분류: Proven. `g=G/c⁴`, `b=1-2m/r`, `b_face=1-2m_face/r_face`에서 끝점 질량 증분은

```text
C_i = 1-2K_i(1-f_i)/r_i-2L_i/r_left
D_i = 1+2K_i f_i/r_i+2Rg_i/r_right
D_i Delta(m_right)-C_i Delta(m_left)
  = g V_i Delta(E_m)+K_i b0_i+L_i b0_left+Rg_i b0_right.
```

를 방사 순서로 푼다. 중심의 `L/r_left` 항은 0이다. 큰 절대 질량이나 정지에너지를 직접 빼지 않고 증분을 사용한다. 구한 질량에서 `a`, scalar 에너지와 물질 상태를 함께 갱신한다.

분류: Proven. 두 끝점 평균을 위줄로 쓰면, Phase38의 물질 일을

```text
Rstar = mean(aR) log(a1/a0)/(a1-a0)
gstar = mean(a alpha T_trace)/abar
chi = 4/(sqrt(b0)+sqrt(b1))²
Delta(E_m)+(Ebar_m+Rstar) Delta(a)/abar+gstar Delta(phi)=F+epsilon
```

로 쓴다. logarithmic divided difference는 대각선에서 `1/a`이고 작은 증분은 정칙 급수로 평가한다. 물질 metric work를 질량 차분식에 넣은 계수는

```text
W_i = g V_i (Ebar_m+Rstar) chi/r_i
Cstar=(1+C)/2-W(1-f),    Dstar=(1+D)/2+Wf
lambda_outer=1
H_i=lambda_(i+1)/Dstar_i,    lambda_i=H_i Cstar_i
```

다. 이번 실행은 `Cstar,Dstar,b>0`을 검사한다. scalar 면 수반은 내부에서 이웃 `H`의 산술 평균, 외곽에서 `H_last`다. 물질 면 수반은 `lambda_face`다. 두 면 가중치의 역할이 다르며 각각 같은 유속을 양쪽 셀에서 공유한다.

분류: Proven. `k_i=H_i bbar_i`, `kf_face=Hgrad_face bbar_face`를 쓰고 scalar 운동량의 직접 source를 `c H_i g V_i gstar_i/w_i`로 둔다. 물질 에너지 유속은 기존 유체 유속 `fE`에 `lambda_face`를 곱하여

```text
F_i = -h c Delta_face(area * lambda_face * fE)/(H_i V_i)
N_i=H_i/a_i,    a_t=(a1-a0)/h
```

로 조립한다. 따라서 이전 lapse와 이전 에너지식을 그대로 둔 채 scalar 항만 추가한 방법이 아니다. 공간 운동량·열수송식은 이 새 `N,a_t`를 사용한다. 나머지 실제 native 잔차 및 수락 문턱은 그대로 검사한다.

분류: Proven. 위 식들의 곱셈 항등식과 면별 합을 취하면

```text
Delta(M_outer)/g
  = W_scalar_boundary - (I_m_outer-I_m_center)
    + sum_i H_i V_i epsilon_i

W_scalar_boundary
  = H_last bbar_outer A_scalar_outer/(4πg)
    * Phibar_outer * Delta(phi_boundary)
I_m_face = h c area_face lambda_face fE_face
```

이다. scalar의 내부 gradient·운동 에너지 교환, reciprocal trace 일 및 물질의 metric work가 이 합에서 상쇄된다. 중심 regularity와 반사 물질 경계에서 중심 scalar port와 물질 경계 유속은 0이다. scalar 운동량이나 장의 변화량으로 나누지 않는다. 유한 반복과 부동소수점 오차는 아래의 별도 수치 문턱으로 검사한다.

## 고정한 실행과 결과

분류: Counterexample candidate. 먼저 원 초기 상태의 안쪽 16셀을 잘라 새 외곽에 반사 벽과 lapse 기준을 둔 작은 공간 조립 대조를 했다. 원 전체 항성을 축소해 물리적으로 동등한 모형으로 만든 것이 아니다. `h=3.2132813194763796874e-6 s`, 외곽 선형 scalar ramp 진폭 `0,1e-6`, 초기 `phi=Pi=0`의 두 경로를 90초 상한으로 등록했다. 2.506초, native 호출 30+32회로 통과했다.

분류: Counterexample candidate. 작은 대조 뒤 전체 5,735셀의 같은 첫 단계 두 경로를 등록했다. 기존 142–330초의 전 구역 단계와 작은 조립의 실측을 근거로 150–850초를 예상했고, CPU worker 8개·BLAS 1스레드·900초 상한·한 쌍으로 제한했다. 기존 무구동·보정 구동 끝점은 초기 추정값과 정확 native 입력 캐시로만 재사용했다. 새 방정식으로 모두 다시 판정했다.

분류: Counterexample candidate. 전체 구역 결과는 다음과 같다. native 정규화 잔차의 기준은 `<=1`이다.

| 항목 | 무구동 | 구동 |
|---|---:|---:|
| native 최대 정규화 잔차 | 0.000441744 | 0.0577209 |
| 기록 수 / Newton 보정 수 | 1 / 0 | 4 / 3 |
| 새 native 물질 호출 | 0 | 13,434 |
| scalar 방정식 잔차 | 0 | 4.186e-24 |
| 질량–경계–잔차 항등식 상대 오차 | 3.639e-20 | 6.028e-20 |
| 총 에너지 결함 / 사전 허용량 | 3.074e-7 | 7.307e-7 |
| 최대 국소 바리온 상대 결함 | 5.931e-24 | 5.764e-24 |
| 핵종 총재고 결함 / 초기 총바리온 | 2.788e-26 | 4.806e-26 |
| 최대 국소 정지계 특성속도 / c | 0.3562931 | 0.3562931 |

분류: Counterexample candidate. scalar 잔차 기준은 `2e-17`, 항등식 상대 오차 기준은 `2e-15`였다. 총 에너지 허용량은 사전 등록한 `sum(H V heat0)*1e-9 + 2e-15*(scalar_energy+abs(boundary_work))`다. 이는 새 물질 일 방정식의 native 허용량과 scalar 부동소수점 오차를 함께 계상한다. 기존의 서로 다른 에너지 방정식에 대한 수지 통과로 바꾸어 부르지 않는다.

분류: Counterexample candidate. 구동 경계의 scalar 일과 총질량 에너지 증가는 약 `1.39339e50 erg`, 그 차이는 `3.24519e32 erg`였다. 이 절대 에너지는 물리적으로 보정한 외부 천체 구동의 크기가 아니다. 최대 내부 scalar 값은 `3.05702e-10`이었다. 구동−무구동 끝점의 최대 `ln rho,ln T,v/c` 차이는 각각 `7.41730e-10,3.64408e-10,4.36836e-11`이다. 관측 신호나 오차가 보장된 응답 계수로 해석하지 않는다.

분류: Counterexample candidate. 전체 쌍의 내부 실측은 227.695초, 프로세스 시작을 포함한 `time` 벽시간은 229.35초였다. user/system CPU 시간은 688.80/7.18초, 보고된 최대 RSS는 407,896 KiB였다. GPU는 사용하지 않았다. 무구동 경로는 저장 상태가 새 잔차 문턱을 이미 만족해 EOS 재평가를 하지 않았지만 초기 provenance·binding 확인을 포함하여 90.62초가 걸렸다. 반복 적분 재사용과 검증 비용을 구분한다.

분류: Counterexample candidate. 저장된 네 끝점의 입력·실행 manifest SHA, 방사 질량 증분과 계량 secant를 재계산했다. 최대 방사 제약 오차는 행의 항들로 정규화하면 `4.686e-20`, 절대값은 `8.674e-19 cm`였다. 인접 질량 증분의 차가 거의 0인 행을 그 작은 차만으로 나누면 최대 `1.248e-8`이 된다. 이를 상대적인 전체 질량 오차로 해석하지 않는다. 이 저장 대조는 새 EOS 호출이나 전체 native 잔차의 독립 재생을 포함하지 않는다.

## 이번 통과로 아직 말할 수 없는 것

분류: Counterexample candidate. 전체 구동 경로의 수반 가중 직접 trace 일은 `6.89938e20 erg`로, 질량–경계 결함 `3.24519e32 erg`보다 작다. 셀별 직접 trace 항의 최댓값도 원 native 에너지 허용치의 `0.0219279`배다. 저장 상태에서 물질 trace 항만 제거한 잔차의 최대 노름은 여전히 `0.0577209`로 통과한다. 따라서 현재의 전체 항성 첫 단계는 이 항의 누락을 수치 수락 기준만으로 검출하지 못한다. 항의 필요성과 부호는 위 조건부 항등식 및 앞선 국소 대조로 뒷받침하며, 이번 결과를 직접 피드백의 수치적 분리 검출로 세지 않는다.

분류: Counterexample candidate. `c h`는 외곽 셀 폭의 약 `0.0261`배다. scalar는 이번 단계에서 경계 가까이에 주로 머물고 물질 내부와의 결합이 작다. 큰 경계 gradient 에너지가 총질량 수지를 지배하는 대조다. 전 구역을 실제 계산했다는 것과 항성 전체의 파동 응답을 해상했다는 것은 다르다.

분류: Conjectural. 다음 우선순위는 **물질 내부와 겹치는 해상된 scalar 상태·구동을 작은 보존 구면 모형에서 만들고, trace 교환을 제거한 대조와 분리되는지 확인하는 것**이다. 이어 비영 이전 scalar 상태를 받는 다단계 식과 시간 대조를 닫는다. 원 바리온·에너지 보존량을 재사용하는 작은 공간 모형, 선형 파동의 시간 척도 및 국소 물질 교환 척도를 먼저 검토한다. 이번 미해상 ramp를 장기간 반복하거나 native 문턱을 자동으로 조이는 방식은 택하지 않는다. 실제 동반천체 구동과 정적 비교 제거·관측 추론은 그 뒤다.

## 재현과 보존

분류: Imported from prior work. WSL Python의 기존 native 의존성 경로와 `OPENBLAS_NUM_THREADS=1`을 사용한다. 계산 커널은 `verification/def_spherical_regular.py`, 기호 항등식은 `verification/def_spherical_regular_identity.py`, 저장 제약·분리 가능성 대조는 `verification/def_spherical_regular_replay.py`다. 새 출력 폴더에서만 `prepare --cells 16`, `run --cells 16`, `prepare --cells 5735`, `run --cells 5735`를 순서대로 허용하며 이미 저장된 결과를 덮어쓰지 않는다.

분류: Counterexample candidate. 소스·사전 계획·네 끝점·반복 기록·수락 결과·기호 검산·저장 대조·자원 기록은 `outputs/direct-eos-gr33/def-spherical-regular/`와 이번 milestone manifest에 결속했다. 이전 실패와 수락 결과의 원문 및 SHA를 유지한다. 이번 네 경로에서 새로운 실행 실패는 없었다. 이전 변동 문서의 스냅샷은 checkpoint `9bb45bb2`에 있다.
