# Adiabatic limit and exact collapse

Status: Proven. There are two separate reductions. The charge oscillator I deltaQddot+Gamma deltaQdot+kappa deltaQ=delta phi approaches a one-pole response only if epsilon_I=I omega^2/|kappa+i Gamma omega| is small; its relative transfer error is bounded by epsilon_I/(1-epsilon_I). The resulting pole approaches a static derivative expansion under the different condition |omega Gamma/kappa|<<1. Neither reduction follows just from the presence of dissipation. See [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md).

Status: Proven. For the settled response `G(z)=c_Y+beta/(1+tau_chi z)`, the degree-N Taylor derivative comparator has exact residual

```math
G(z)-\left[c_Y+\beta\sum_{n=0}^{N}(-\tau_\chi z)^n\right]
=\frac{\beta(-\tau_\chi z)^{N+1}}{1+\tau_\chi z}.
```

Status: Proven. On `z=i omega` and `|omega tau_chi|<=rho<1`, its modulus is at most `|beta| rho^(N+1)`. A convergent infinite geometric expansion requires this frequency restriction; a formal inverse operator is not convergence for arbitrary forcing.

Status: Proven. At zero lag the readout is `(c_Y+beta)F`. At zero frequency only the settled solution becomes constant; an initial transient may still evolve. Zero beta removes the driven-state contribution but does not remove an independently initialized transient when c_chi is nonzero.

Status: Proven. At a single known frequency, `a0 F+a1 dot F` exactly matches the pole with `a0=c_Y+beta/(1+omega^2 tau_chi^2)` and `a1=-beta tau_chi/(1+omega^2 tau_chi^2)`. This is exact finite-sample collapse even away from small lag.

Status: Proven. With K distinct positive carriers and freely shared real coefficients, degree `N>=2K-1` interpolates the pole exactly. Thus low-frequency truncation is one collapse mechanism, not a necessary condition for all finite local-comparator collapse. See [the exact sampling proof](../paper/manuscript.md).

Status: Proven. Nuisance absorption and finite precision give additional boundaries independent of adiabaticity. A projected signal must survive the specified covariance-weighted nuisance span before any response distinction becomes measurable.

Status: Proven. A reciprocal positive relaxation spectrum with all rates >=Lambda supplies a physically specified fast comparator without selecting a Taylor order. Its two-frequency quadrature-ratio inequality is violated by a positive pole with tau_chi>1/Lambda. The gap is a separate premise; it is not inferred from an observed lag or introduced as a favorable prior.

Status: Proven. Settling to a fraction epsilon of an independent homogeneous amplitude takes tau_chi*log(1/epsilon). For tau_chi=500 days, one-per-mille settling takes 3453.88 days. Without an initial-amplitude bound, this relative decay supplies no uniform absolute signal bound.


## Request 12 follow-through

Status: Proven. Equilibrium susceptibility alone does not fix damping or inertia; Gamma=kappa*tau permits arbitrary relaxation at fixed equilibrium in the admitted EFT class. A specific microscopic theory can relate these coefficients and must be matched separately.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Imported from prior work. For the specified SLy/DEF candidate, static susceptibility/c is about 0.1387 ms and the computed scalar-pole decay time is about 0.1995 ms. This supplies a numerical fast-response example. It does not establish a uniform fast-rate gap over EOS, gravity coupling or stellar branches. See [matching equations and numerical checks](../notes/REQUEST13_STELLAR_DERIVATION.md).

## Request 14 validated flow

Status: Proven. The validated GR IVP has no dynamic-chi coupling. Its numerical error enclosures therefore neither constrain tau_chi nor change the analytic adiabatic-collapse boundary. A certificate for a nonzero chi signal must also enclose the matched force and readout. See [scope and missing links](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).


## Request 15 후속 검증

분류: Proven. 초기화 규약의 수정과 매개변수 재매핑은 새 완화 pole을 만들지 않는다. 이번 질량·Kepler·대수 관측식의 부분 인증은 기존 단일 pole의 단열 붕괴 조건을 변경하지 않는다. 분류: Conjectural. EOS 응답의 전체 이체 구동·힘·관측량 연결은 여전히 별도의 물리 조건이다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).


## Request 16 다체 관측식과 영 구동 경계

분류: Proven. 비균일 셀에서 끝점 오차 eps, 곡률 잔차 rho, 셀 폭 h를 알면 정확한 cubic에 대한 값·시간 미분·셀 적분 오차는 각각 eps+rho*h²/8, 2eps/h+rho*h/2, eps*h+rho*h³/12 이하이다. 자연 경계조건에 C4 전역 오차 공식을 강제로 적용하지 않는다. native 산술 반올림과 이동 격자 매개변수 미분은 별도 항이다. 영 구동 가지의 정확한 영 응답은 단열 근사를 요구하지 않는다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).


## Request 17 비영 구동과 동반성 응답

분류: Proven. 명시한 감수율·거리 구간에서 결합 행렬의 가역성과 정적 Hessian 양의 정부호, 전하의 역거리 일·이차 미분 구간을 확인했다. 관성을 생략한 선도 monopole 복사의 rank-one 감쇠 아래 총전하 S는 (Ceff/c) Sdot+S=Ceff*phi를 만족한다. 지정 구간의 완화시간은 0.150614–0.150631 ms이며, 평형 초기자료에서 각 전하의 연속 추적 오차는 1.37e−15 m 이하라는 조건부 상계를 얻었다. 관성·고차 복사·실제 궤도 오차는 이 상계에 포함되지 않는다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).


## Request 18 열 백색왜성 구조와 응답 경계

분류: Proven. 비음수·구대칭·고정 밀도에서 eta=4pi|beta|G/c² integral rho*r dr<1이면 scalar 적분 연산자가 축약 사상이어서 정적 해가 유일하다. Q0=|beta|GM/c²에 대해 Q0<=chi<=Q0/(1−eta), |f(k)−chi|/chi<=(kR)²/[3(1−eta)]+|k|Q0/(1−eta)²를 얻었다. 두 공개 열 구조에서 정의한 구각 모형은 eta<0.000166이다. 실제 시간 의존 항성 구조를 이 두 끝점으로 포괄한 것은 아니다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).


## Request 19 열 진화 재현과 질량 보정 대조

분류: Proven. 재현한 광학 모델 18969에서 정의한 고정 밀도 구각 scalar 모형은 eta=0.000165267742<1이고, 궤도 주파수의 상대 정적 차이 상계는 2.03086945e-10이다. 질량 보정 진화의 기록된 최적 구조에 동일한 조건부 식을 적용한 상계는 2.088720088e-10이다. 두 결과 모두 실제 시간 의존 열 별의 전체 미분 오차 보장이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).


## Request 20 열 구조 민감도와 질량 정의

분류: Proven. 두 고정 양의 구대칭 퍼텐셜의 eta_i<1 조건에서 정적 응답 차이의 해석 상계를 도출했다. 질량이 같으면 선도 Q 차이가 소거된다. 실제 변형 구조의 구각 합집합에서 65자리 구간 연산으로 상계를 계산했다. 이는 물리적 항성 오차 또는 전구간 미분 오차 보장이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).


## Request 21 열 EOS와 GR 질량 매칭

분류: Proven. 정지 에너지와 내부에너지의 반대 방향 기준 이동은 총 에너지 밀도를 보존한다. FreeEOS의 H2 영점 이동을 소스의 정확한 상수로 제거하고 중성 원자 질량을 더했다. 압력 좌표 TOV 변환과 외곽 폴리트로프의 제1법칙도 기호적으로 검증했다. 기존 평탄·고정 구각 scalar 오차 상계는 이 새 GR 유체 모형의 오차 보장이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).


## Request 22 GR 후보의 열수송과 물질 보존

분류: Proven. 양의 국소 전도계수와 준정적 구대칭 계량에서 열유속은 -K sqrt(f) (dT/dr + T dnu/dr)이며 영유속 조건은 T exp(nu)=상수다. 정확한 정적 대각 계량과 정지 물질, 다른 상쇄 흐름 부재에서는 G_tr=0이 순 열유속 0을 요구한다. 유한 광도를 쓰는 항성에는 준정적 근사 또는 시간 진화를 명시해야 한다. 이 경계는 궤도 완화시간이나 보편 동적 no-go를 증명하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).


## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성

분류: Proven. 고유 바리온 질량을 독립 변수로 한 TOV 식과 고정 조성에서의 엔트로피 정규화를 기호 검산했다. 조성과 단위 바리온 질량당 엔트로피를 유지한 균일 질량 배율은 비율을 보존하지만 총 원소 재고와 총 엔트로피는 그 배율만큼 바꾼다. 정확한 정적 계량의 순 열유속 no-go와 기존 평탄 구각 scalar 오차 상계의 적용 경계는 그대로다.

세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).


## Request 24 조성 변화의 반응 에너지 기준

분류: Proven. 조성이 변하면 du_B+P dv_B=T ds_B+Σ μ_i^th dY_i이다. 따라서 ds_B=0만으로 du_B+P dv_B=0을 결론낼 수 없다. 고정 조성의 Request23 엔트로피 재구성은 보존되지만, 반응이 켜진 후 동일한 엔트로피를 고정하는 것은 일반적인 열 진화 해가 아니다.

세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).
