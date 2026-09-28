# Request 38 — 비평형 열 상태를 유지하는 정칙 scalar–물질 교환

작성: 2026-09-17. 선행 결과: [Request 37](REQUEST37_REGULAR_SCALAR_KO.md).

## 이번에 해결한 것

분류: Proven. **열유속과 평형 엔트로피를 고정해야 한다는 제약을 제거했다.** 실제 보존식과 두 열수송식을 함께 풀면 `v=0,Q≠0`에서도 계량·scalar의 일을 전달하는 정칙한 국소 상태 변화가 정의된다. 비평형 엔트로피 전류의 수지도 맞는다. 앞선 고정 상태의 반례는 보존하지만, 이를 모든 결합법에 대한 장애로 확대하지 않는다. 전역 정준 열 좌표를 먼저 찾는 것은 이 방법의 필수 조건이 아니다.

분류: Counterexample candidate. 기존 native 분자 EOS·불투명도와 26종 조성을 그대로 사용하여 온도·속도·두 열유속·scalar를 같은 비선형 단계에서 풀었다. 8/16/32단계와 무열유속·가역 대조의 5개 경로가 원 기준을 통과했다. 동결된 척도의 끝점 최대 노름 차수는 2.01053, 최대 에너지 수지 오차는 초기 열용량 척도 대비 5.660e-15였다. 이는 **국소 결합 적분을 실행한 loophole progress**다.

분류: Conjectural. 공간적으로 균일한 셀의 교환 대조이며, 계량 변형과 scalar 진동자의 시간 척도는 시험용이다. 새 방법으로 전 구역 구면 GR 항성을 진화시킨 것이 아니다. 실제 외부 구동·관측 폐쇄와 물리적 열수송 계수 보정도 남는다.

## 고정 엔트로피 가정에서 동적 열 상태로

분류: Imported from prior work. 기존 구현은 [Maartens, 식 2.17·2.20·2.22](https://arxiv.org/html/astro-ph/9609119v1)의 열유속 이차 엔트로피 모형을 두 독립 열 채널에 적용한다. 점성·교차계수·회전 항은 생략하고 전체 수송 미분항은 유지한다. 이 모형의 물리 유효성을 새로 인증하지 않는다.

분류: Proven. `c=1`에서 `W=(1-v²)^-1/2`, `Q=q_1+q_2`, `w=epsilon+P`, `ell=ln a`라 하자. 공간적으로 균일한 방사 방향 변형 `ds²=-dt²+a(t)²dx²+dy²+dz²`를 사용한다. 응력 `E,S,R`은 이전 단계와 같으며 보존량과 교환식은

```text
B=a rho W,             J=a² S
dB=0,                  dJ=0
dE=-(E+R) d ell-alpha T_trace dphi
d epsilon=w d ln rho+rho T ds-alpha T_trace dphi
```

이다. `d ln rho=-d ell-W²v dv`를 대입하면

```text
rho T ds = -v dQ-2Q W² dv-2Qv d ell
```

를 얻는다. 따라서 앞서 사라지지 않던 열유속 항은 평형 엔트로피의 변화로 처리된다. `ds=0,dQ=0`을 추가로 강제하지 않는다.

분류: Proven. 각 열 채널에 `C_j=K_j T/tau_j`, `beta_j=1/(C_j T)`를 두면, 완전한 공변 열법칙의 균일 셀 형태는

```text
dq_j + C_j [W² dv+v(d ell+d ln T)]
     + q_j/2 [d ln beta_j+d ln W+d ell]
  = -q_j dt/(tau_j W).
```

이다. 구현에서는 `q_j`가 열 에너지 유속을 `c`로 나눈 에너지 밀도 단위이므로 `C_j=K_j T/(c² tau_j)`를 사용한다. 밀도·온도·scalar에 따른 `beta_j` 변화를 실제 EOS와 불투명도 상태에서 모두 평가한다.

분류: Proven. 비평형 엔트로피의 좌표 부피 밀도는

```text
Sigma=a W [rho s+vQ/T-(beta_1 q_1²+beta_2 q_2²)/2]
dSigma/dt = a [beta_1 q_1²/tau_1+beta_2 q_2²/tau_2] >= 0
```

다. 이 항등식은 Gibbs 관계·보존식·전체 열법칙 아래에서 기호 검산했다. 완화 항을 제거한 가역 부분에서는 우변이 0이다. 평형 엔트로피 `s`를 고정하는 것과 전체 비평형 전류의 보존은 다르다. 연속 모형의 항등식이며 아래 유한 단계의 엔트로피를 정확히 보존한다는 주장은 아니다.

## 정지 상태에서의 정칙성

분류: Proven. `v=0`에서 바리온 보존을 소거한 뒤 미지수를 `(d ln T,dv,dq_1,dq_2)`로 잡자. `c_T=∂epsilon/∂ln T`, `b_Tj=∂ln beta_j/∂ln T`이면 시간 미분 행렬은

```text
[ c_T,           2Q,  0, 0 ]
[   0,            w,  1, 1 ]
[ q_1 b_T1/2,   C_1,  1, 0 ]
[ q_2 b_T2/2,   C_2,  0, 1 ]

det = c_T(w-C_1-C_2)+Q(q_1 b_T1+q_2 b_T2).
```

이다. 이 determinant가 0이 아니고 EOS가 매끄러우면 국소 상태 변화가 유일하다. `Q=0`에서는 `c_T>0`, `w>C_1+C_2`가 충분하다. 일반 비영 열유속에서 열유속 보정항을 버리면 안 된다. `v`, `Q`, scalar 운동량 또는 장의 변화량으로 나누지 않는다.

분류: Counterexample candidate. 저장된 native 초기 5,735셀의 실제 두 열유속·EOS 도함수·수송 계수에서 이 국소 행렬을 구성했다. `det/(c_T w)`의 최솟값은 0.99999995950, 연립 잔차는 1.046e-16, 계량 일의 상대 에너지 항등식 잔차는 1.083e-19였다. 새 EOS 호출 없이 저장 상태를 재사용했다. 이것은 초기 상태의 국소 변형 가능성 검사이며 전체 공간 진화·특성속도·전역 유일성의 증명이 아니다.

## 실제 native EOS로 수행한 유한 단계

분류: Counterexample candidate. 기존 shell 초기 자료의 셀 2762를 선택했다. `rho≈2.2629e3 g/cm³`, `T≈1.9479e7 K`와 26종 조성을 가져오고 `v=0`, scalar 운동량 `p=0`에서 시작했다. 두 열유속은 각각 `.01 sqrt(C_j c_T)`로 주어 각 채널의 초기 이차 엔트로피 감소량을 같은 작은 크기로 맞췄다. **실제 항성 열유속을 사용한 경로는 아니다.** 저장된 해당 셀의 열유속은 훨씬 작다.

분류: Counterexample candidate. 전도·복사의 초기 완화시간은 약 `4.113e-4 s`, `7.496e-15 s`로 크게 다르다. 같은 열유속 크기를 두 채널에 주면 작은 전도 관성 때문에 이차 엔트로피 감소가 과도해지므로, 각 채널의 `C_j`에 맞춘 위 초기값을 실행 전에 정했다. 초기 `q_cond/w=3.342e-16`, `q_rad/w=3.688e-10`이다.

분류: Counterexample candidate. 무차원 시간 단위는 초기 복사 완화시간, 종료 시간은 `2π`, 외부 계량은 `a(t)=exp(.001 sin t)`다. scalar 초기값은 `.01 sqrt(c_T/w)=1.38059e-5`, 인공 관성은 `I=16w`다. 이 약 `4.71e-14 s` 구간의 진동자는 시험 장치이며 궤도 구동이나 물리적 항성 모드의 시간이 아니다.

분류: Proven. 유한 단계는 바리온을 소거하고, 동일한 새 온도·속도·두 열유속·scalar 상태에서 다음을 함께 푼다. 위줄은 두 끝점 산술 평균이다.

```text
rho=B/(aW),                 Delta(a²S)=0
Delta(aE)=-mean(aR) Delta(ln a)-gbar Delta(phi)
gbar=mean(a alpha T_trace)
Delta(phi)=h pbar,          Delta(p)=h[-phibar+gbar/I]
```

분류: Proven. 두 열 채널에는 위 공변식의 대칭 midpoint 차분을 사용한다. 따라서 scalar 진동자 에너지 `I(p²+phi²)/2`의 증가는 `gbar Delta(phi)`이고, 물질의 직접 교환과 정확히 상쇄된다. 총 에너지의 유일한 외부 일은 등록된 계량 일이다. native EOS의 정지에너지는 직접 두 큰 수를 빼지 않고 고정 바리온과 `expm1`으로 차분했다. 유한 반복 오차는 별도로 검사한다.

분류: Counterexample candidate. 사전 기준과 최종 측정은 아래와 같다.

| 항목 | 결과 | 원 기준 |
|---|---:|---:|
| 끝점 고정 척도 최대 노름 수렴 차수 | 2.01053 | >=1.7 |
| 모든 경로의 최대 단계 잔차 | 4.10e-12 이하 | <2e-8 |
| 최대 총 에너지–계량 일 수지 / 초기 `c_T` | 5.660e-15 | <1e-8 |
| 32단계 엔트로피 증가, `c_T/T` 정규화 | 4.99998e-5 | 양수 |
| 32단계 엔트로피 수지 상대 차이 | 0.009529 | <=0.15 |
| 비영 열유속 세 경로의 scalar 운동량 부호 전환 | 각각 1회 | >0 |
| 무열유속 대조의 속도·두 열유속 | 0 유지 | 고정 척도 <1e-12 |

분류: Counterexample candidate. 끝점 척도는 `(Delta ln T/.001, v/v_scale, q_j/q_j_scale, phi/phi_scale,p/phi_scale)`다. 수락 차수 2.01053은 이 벡터의 최대 노름이며 각 성분 모두가 같은 차수라는 뜻이 아니다. 보조 개별 차수는 온도 1.8334, 속도 1.8968, 전도 1.8345, 복사 1.8968, scalar 위치 2.0105, scalar 운동량 1.7358이다. 전도 끝점 변화는 매우 작아 이 짧은 경로가 전도 완화시간을 해상한 것은 아니다.

분류: Counterexample candidate. 엔트로피 수지 상대 차이는 8/16/32단계에서 0.13343 → 0.03706 → 0.009529로 감소했다. 이는 엔트로피 변화와 양의 생산률 적분의 유한 시간 차분 오차를 포함한다. 가역 대조의 절대 엔트로피 변화는 `c_T/T` 척도로 -3.05e-14였다. 생산이 0인 이 대조에서 저장한 상대 비율은 의미가 없으며 절대 오차로 해석한다. 유한 단계의 무조건적인 엔트로피 증가 정리를 주장하지 않는다.

분류: Counterexample candidate. 예비 경로는 6.29초·native 호출 138회, 나머지 계산은 35.71초·823회였다. CPU 1개·BLAS 1스레드로 실행했고 본 실행의 프로세스 최대 RSS는 약 157 MiB였다. 소요시간은 실측 예비 경로로 추정한 52초 및 240초 상한 안에 있었다. 기존 장기 GR 이력을 재계산하지 않았다.

분류: Counterexample candidate. 모든 경로를 저장한 후 최종 JSON 요약의 NumPy 논리값 직렬화에서 실패했다. 실패와 소스를 보존하고 저장된 5개 경로만으로 같은 수락 기준을 재계산했다. 추가 적분은 0회다. 별도 native 호출 75회로 저장 경로를 재구성한 최대 에너지 수지 결함은 5.660e-15, 엔트로피 기록 차이는 6.005e-19였고 끝점 수렴 차수도 일치했다. 중간 상태는 binary64로 저장했으므로 bitwise 재생이나 모든 단계 잔차의 독립 재검사라고 하지 않는다.

## 구면 질량 제약과 연결할 정확한 유한 항등식

분류: Proven. 두 양의 계량 상태에서 다음 정칙 항등식이 성립한다.

```text
Delta(a)/abar = chi Delta(m_i)/r_i
chi = 4/(sqrt(b_old)+sqrt(b_new))²,  b=1/a²
L(a0,a1)=log(a1/a0)/(a1-a0),       L(a,a)=1/a
Rstar=mean(aR)*L(a0,a1)
gstar=mean(a alpha T_trace)/abar
Delta(E_m)=-(Ebar_m+Rstar) Delta(a)/abar-gstar Delta(phi)+F
```

분류: Proven. `L`은 양의 상태에서 매끄러운 logarithmic divided difference이며 작은 증분에서는 `log1p(x)/x`의 연속 연장으로 계산할 수 있다. 이것은 방금 시험한 로그 계량 일과 같은 유한 교환식이다. Phase37의 scalar 질량 제약에 대입하면

```text
Dstar=Dbar+Ggeom V_i(Ebar_m+Rstar) chi f_i/r_i
Cstar=Cbar-Ggeom V_i(Ebar_m+Rstar) chi(1-f_i)/r_i
```

를 얻는다. 이 계수로 같은 방사 방향 수반 `H_i`를 구성하면 물질의 trace source가 scalar 갱신에 `H_i Ggeom V_i gstar / w_scalar_i`로 들어가며 역운동량 보정이 필요 없다. 물질 수송 증가량을 `F_i=-Delta_face(I_m)/(H_i V_i)`, `I_m=h H_face area k S_m`으로 정의하면 수반 가중 총합은 정확히 경계 유속으로 망원 소거된다. 모든 식을 기호 검산했다.

분류: Conjectural. 이 마지막 구면 구성은 **조건부 대수적 연결식**이다. 정칙 계량 가지, 같은 수반·체적·공유 면 사용, 선언한 물질 일의 실제 이산 잔차 반영이 필요하다. native 구면 운동량·열수송·조성 이류와 함께 구현하거나 실행하지 않았다. 현재의 구면 lapse/에너지 연산자를 그대로 두고 source만 바꾸면 이 증명이 적용되지 않는다.

## 다음 목표와 근거

분류: Conjectural. 다음 목표는 **위 수반과 공유 면 에너지 유속을 실제 구면 잔차에 조립해, 정칙 scalar 파동과 native 물질·두 열유속을 동일 단계에서 푸는 것**이다. 기존 수락 초기 상태와 저장 이력을 재사용한다. 먼저 작은 공간 대조의 물질·scalar 교환과 경계 수지를 확인하고, 그 실측으로 전 구역 첫 단계의 예산을 정한다. 새 장기 이력이나 더 촘촘한 전체 기간 경로를 자동 실행하지 않는다.

근거: `verification/def_heat_coupling.py`, `def_heat_coupling_replay.py`, `def_heat_metric_secant.py`; `outputs/direct-eos-gr33/def-heat-coupling/`의 계획·기호식·5개 저장 경로·직렬화 실패·복구 판정·native 재구성·구면 항등식. 통합 결속: `gr-heat-coupling-milestone-manifest.json`. 이전 문서 snapshot은 commit `6cd55662`에 남긴다. 재검산은 `symbolic()`, `native_lift()`, `replay()`, `check()`를 호출하며 기존 결과 파일을 덮어쓰지 않는다.
