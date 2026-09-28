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


## Request 25 새 GR 상태의 미세물리 평가

분류: Proven. 구역 내부에서 조성·엔트로피가 고정되면 TOV의 계량 퍼텐셜은 dν=−d ln h를 만족한다. 이를 같은 수학적 대기와 연결해 적색편이를 계산했다. 국소 에너지와 고유 시간의 변환은 무한원 광도 원천에 e^(2ν)를 준다. 반응이 켜진 후의 엔트로피 보존이나 정확한 정적 열유속의 허용을 뜻하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).


## Request 26 남은 폐쇄 조건의 실행 검산

분류: Proven. 명시한 계량·trace 값 영역에서 scalar 에너지의 음의 퍼텐셜 항은 양의 구배 항의 0.018670 미만으로 상계되어 음의 주파수 제곱 모드를 배제한다. 같은 고정 배경의 정적 beta 도함수 상계도 얻었다. 전체 EOS/TOV 연속 해의 영역 소속이나 유체·열 안정성 인증은 아니다. 영 scalar 가지에서는 열 교란이 선형 scalar 방정식에서 분리된다.

세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).


## Request 27 직접 반응 벡터와 명시적 PP 상태

분류: Proven. 안정한 실수 모드들의 응답이 D+Σaⱼ/(1+iωτⱼ)이면 빠른 모드를 정적 계수로 치환한 오차는 Σfast|aⱼ|ωτⱼ 이하이다. 느린 모드도 ωτ≪1에서는 선도 정적 계수로 붕괴하며, 그 구동 또는 출력 잔차가 0이면 관측되지 않는다. 영 scalar 배경의 선형 열 분리 경계는 새 핵 반응 상태의 존재로 해제되지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).


## Request 28 반응 에너지 보존과 scalar 미분 구간

분류: Proven. 내부 핵 전환은 닫힌 계의 총에너지에 독립적인 원천으로 더해지지 않는다. 반응열·정지에너지·조성별 내부에너지와 경계 복사 손실을 함께 계수해야 한다. 이전 PP 준정상 부분계의 가열은 벌크 연료를 필요로 하므로 네 중간 핵종의 상태를 닫힌 질량 저장고로 읽을 수 없다. 영 scalar 가지와 유한 carrier 정적 보간 경계는 유지한다.

세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).


## Request 29 공통 EOS와 GR 열 경로의 관측 연결

분류: Proven. 같은 핵 반응의 정지에너지 손실과 열량을 별도 총질량 원천으로 중복 계수할 수 없다. 고정 바리온 열 경로에는 화학 에너지, 적색편이된 중성미자·광자 손실과 표면 압력 일을 포함해야 한다.

분류: Conjectural. 표면 광도가 정해졌다고 내부 수송이 닫히거나 단일 orbital relaxation 상태가 식별되는 것은 아니다. 영 scalar 가지와 정적 carrier 보간 경계는 그대로 유지한다.

세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).


## Request 30 EOS 역산과 약반응 미분 수정

분류: Proven. 선형 혼합 반응률의 온도 미분에는 가중치 미분 항이 필요하다. 고정 표·혼합 분기에서 지수 일차식의 3차 미분으로 중심 차분 오차를 제한할 수 있다.

분류: Conjectural. 반응 미분을 수정해도 전체 항성의 빠른 모드 제거·orbital 완화시간·비영 읽기 잔여량이 자동으로 보장되지 않는다.

근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).


## Request 31 조성 미분과 보존형 GR 재적분

분류: Proven. 바리온 질량 좌표에서 구대칭 준정적 Fourier 식은 면적 제곱·적색편이·온도 기울기로 표현되며 지정 면 이산식은 일정한 적색편이 온도에서 영유속을 준다.

분류: Conjectural. 이 항등식이나 유한 시간 대조는 orbital 완화시간, 빠른 모드 제거와 물리 오차의 연속 보증을 대신하지 않는다.

근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).


## Request 32 직접 원소 EOS 중간 검증

분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.

분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.

근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 32 구조 피드백 엔탈피 시간 대조

분류: Proven. 양의 원소 치환 가중치가 이온 수·평균 전하·전하 제곱합을 모두 보존하면 치환 전하의 분산이 0이어야 한다. 따라서 지원하지 않는 전하를 다른 전하들로 치환하는 현재 사상은 완전 이온화 한계의 Coulomb 입력까지 정확히 보존할 수 없다. 조성만의 에너지 기준 이동으로 이 한계를 제거할 수 없다.

분류: Counterexample candidate. 직접 원소 후보의 CT 미분 대조는 실패다. 직접 Fermi 적분의 원래 정밀도 및 강화 정밀도 전구역 대조는 각각 실패, 통과다. 이 유한 대조는 물리·연속 EOS 오차 보증이 아니다. 교환 항의 CT 호출까지 직접 적분으로 바꾼 후보의 전구역 대조는 통과다.

분류: Proven. 보관 CT 근사식의 η=1 함수 및 일차 미분 점프가 0을 배제함을 계수 구간으로 보였다. 또한 서로 다른 동위원소 질량의 부분 이온화 평형은 조성만의 엔트로피 이동으로 원소 단위 EOS에서 정확히 복원되지 않는다. 새 공통 Fermi 평가의 η=1·4 지정 원시 값·미분 12개는 구간 구적으로 반폭 1e-11 이하에 감쌌고, 해당 정규화 출력의 절대오차 상계는 모두 1e-9 이하다. 이는 지정 원시 함수의 점별 보증이며 전체 EOS나 GR 오차 보증이 아니다.

근거: [한글 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 33 수정 EOS와 GR 연결 진행

분류: Counterexample candidate. 수정 EOS의 GR 진입 대조에서 미정 Coulomb 회절 출력 gamma_e를 발견했다. 비활성 회절 경로가 이 슬롯을 정의하지 않으므로 이를 반환 계약에서 제외하고, 나머지 21개 출력의 같은 역산·호출 이력 기준을 유지해 84개 지정 상태에서 통과했다. 원래 22출력 대조 실패는 보존한다. 이전 열역학 미분 검사는 이 미정 슬롯을 사용하지 않았다. 새 EOS의 GR 표와 고정 상태 반응 입력을 계산 중이며, 그 구조나 전체 시간 경로의 완료를 주장하지 않는다.

분류: Proven. 다른 항과 허용 집합을 고정한 동위원소 병진 질량 변경에 대해 모든 이온화 분율에 적용되는 조건부 최소 자유에너지 차이 상계를 얻었다. 압력·미분 또는 전체 물리 EOS 오차의 보증은 아니다.

분류: Conjectural. 물리 EOS·연속 전구간 오차, 자체 수송·유체·계량·대기, 비영 scalar 구동·전하와 완전한 관측 추론은 계속 진행한다.

근거: [단계 33 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 불투명도 입력 검증

분류: Counterexample candidate. 실제 호출 경로는 Type1이었다. Type2를 계측한 첫 시도의 무호출 실패를 보존한 뒤, 수정 경로에서 5,735구역의 밀도·온도·26핵종 조성과 반환 불투명도를 저장 프로파일에 대조했다. 반환 불투명도 차이는 0이었다. 동일 전자 입력을 다시 전달한 대조에서도 입력·조성·반환값이 모두 비트 단위로 같았다. 새 EOS 전자수·미분 입력의 실제 개입과 불투명도 유한 미분 대조는 별도로 실행 중이다.

분류: Proven. 동일한 양의 상태밀도를 갖는 이상 전자·양전자 기체에서 eta>=-50, 0<kBT/(m_e c^2)<=0.005이면 총 자유입자수와 순 자유입자수의 상대 차이는 1.03e-130보다 작다. 이는 선언한 이상기체 값의 조건부 상계이며 상호작용 플라스마, 미분 또는 항성 궤적의 인증이 아니다. 실제 EOS의 표본이 이 영역에 드는지는 따로 검사한다.

분류: Conjectural. 입력 계약 일치는 물리적 불투명도·자체 수송·전체 GR 시간 진화의 완료를 뜻하지 않는다. 근거: [단계 33 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 원천 연결과 불투명도 미분 검증

분류: Counterexample candidate. 새 EOS의 5,735구역 반응 입력을 실제 native 함수에 연결했다. 첫 계측은 180초 제한으로 실패했다. 불변 압축 입력을 한 번만 읽도록 바꾸고, 작성 중인 프로파일의 불완전 헤더를 기다리도록 수정한 별도 계측 코드에서 같은 제한을 유지했다. 이전 원천의 비트 일치를 확인한 뒤 새 입력 계측을 약 15.3초에 마쳤다. 새 전자 입력의 불투명도 개입은 반환값 변화 0이었다. 실제 온도 범위에서 이 입력을 사용하는 Compton 경로가 비활성이었다.

분류: Counterexample candidate. 원래 불투명도 미분은 사전 0.1% 기준을 실패했다. 실제 단조 3차 보간 19,242회를 읽어, 값에 적용된 선택 분기의 미분과 미분값을 별도로 보간한 결과의 불일치를 확인했다. 선택된 식을 직접 미분해 복사·전도 결합에 전달했다. 이어 원본 계수를 읽은 이중정밀도 계산에서 전도 혼합항 누락을 확인했다. 정수 z36th에 1/36을 대입하는 보관 코드가 그 원인이었다.

분류: Proven. 저장된 네 보간값의 선형 변화 경로 12,828개를 정확한 유리수 연산으로 판정했다. 12,826개는 단일 식으로 보증했고 나머지 두 경로는 분할했으며 두 비영 미분 점프를 보존했다. 누락된 전도 혼합항은 전체 29,400개 저장 격자 구간에서 값·두 미분의 조건부 상계를 얻었다. 실제 표 값의 물리 오차나 전체 EOS 미분의 보증으로 확대하지 않는다.

분류: Counterexample candidate. 혼합항을 복구한 계산의 유한 미분 최대 점수는 5.404e-4로 통과했지만, 원래 native 값과의 최대 차이 6.323e-4는 원래 값 대조 기준 1e-4를 실패했다. 혼합항을 0으로 둔 별도 재현은 native 값 차이 6.251e-6으로 값 대조를 통과했으나 미분 최대 점수 1.0624e-3으로 전체 판정은 실패했다. 이 원래 판정들을 유지한다. 아직 생성되지 않은 새 GR 상태에 적용할 별도 검증 기준과 모델 해시를 고정했다.

분류: Conjectural. 물리 EOS·표 데이터 오차·전체 연속 미분, 자체 수송·유체·계량·대기의 시간 진화, 실제 비영 scalar 구동·전하 및 완전한 관측 추론은 남는다. 동적 chi 후보의 계산 기반을 보완했으며 신규 관측량 판정과 최종 투고본 갱신은 보류한다. 근거: [단계 33 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 새 GR 상태의 수송·부력 검증 준비

분류: Counterexample candidate. 새 EOS의 5,735구역, 구역당 17점 GR 표를 완성했고 별도 직렬 재계산 아홉 점을 비트 단위로 재현했다. 같은 바리온·핵종 재고를 유지하는 GR 구조 계산을 진행한다. 아직 생성되지 않은 상태를 대상으로 반응·불투명도 보조 입력의 공통 EOS 평가와 수송·부력 검증 기준을 고정했다. 공통 입력 계산기의 제조 해 대조는 입력 순서와 밀도·온도 미분 방향을 검증했다. 이 준비를 새 GR 해나 시간 진화의 완료로 세지 않는다.

분류: Proven. 고정 GR 배경의 구역 경계에서 양의 계수 B와 양의 적색편이 온도 theta를 사용하는 L=B(theta_inner^4-theta_outer^4) 교환은 경계 밖 유속이 0일 때 총 적색편이 에너지를 보존한다. 구역쌍의 엔트로피 생성은 B(x-y)^2(x+y)(x^2+y^2)/(xy)>=0이다. B가 상태에 의존해도 이 순간 항등식은 유지된다. 일정 계수 선형 대조에서 비선형 양의 수송 법칙까지 확장한 조건부 정리이며, 유체·계량 진화와 복사 대기를 포함한 정리는 아니다.

분류: Proven. 매끄럽고 제1법칙을 만족하는 EOS, 순간 압력 평형, 보존된 유체 요소의 엔트로피·조성, 무시한 열교환·계량 섭동 아래 부력 가속도는 g*A*xi이다. A=(d epsilon/dell)/(epsilon+P)-(dP/dell)/(Gamma1*P), 따라서 고유시간의 국소 부력 주파수 제곱은 -g*A다. epsilon에는 실제 핵 정지 에너지를 포함한다. 저장 격자의 유한 차분으로 A를 계산하는 절차의 오차 보증은 별도다.

분류: Imported from prior work. 국소 고유 기준계의 혼합 길이 변수와 균일 조성 대류 식은 [Thorne의 원문](https://ntrs.nasa.gov/api/citations/19760017020/downloads/19760017020.pdf), 인쇄 쪽 6-7의 식 (8)-(10)과 대조했다. 새 EOS·조성에서의 난류 수송 해를 이 인용만으로 얻은 것으로 간주하지 않는다.

분류: Conjectural. 새 GR 상태의 실제 반응·수송 대조, 조성 혼합과 대류, 복사 대기, 유체·계량의 시간 진화, 전체 물리 EOS·연속 오차, 비영 scalar 구동 및 완전한 관측 추론은 계속 수행한다. 수송 모듈의 유효 불투명도에는 전도가 포함되므로 이를 광자 평균 자유행로나 복사 광학 깊이로 사용하지 않는다. 근거: [단계 33 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 EOS 오차 전달 정리

분류: Proven. 내부 평형이 동일하고 강한 볼록성을 유지해도 작은 자유에너지 값 오차가 임의로 큰 매개변수 미분 오차와 공존하는 반례를 검산했다. 별도로 제약을 제거한 내부 좌표의 강한 볼록성, 잔차, 모델 일차·이차 미분 및 블록 Lipschitz 상계가 주어지면 평형 위치·감소 자유에너지의 기울기·Hessian 오차를 제한하는 식을 얻었다. 압력·엔트로피·열용량으로의 오차 전달도 기호 검산했다.

분류: Conjectural. 실제 EOS의 전체 평형 분율, 투영 Hessian 하한과 물리적 미분 오차는 아직 인증해야 한다. 값 상계나 두 간격 차분으로 이를 대체하지 않는다. 정리의 가정·증명·정확한 대조와 남은 측정 항목은 [EOS 오차 전달 노트](../notes/REQUEST33_EOS_ERROR_TRANSFER_KO.md)에 보존한다. 전체 물리 EOS·GR 진화·비영 scalar 및 관측 추론 목표는 계속 열려 있다.


## Request 33 Fermi 연속 영역 인증

분류: Proven. eta∈[-17,24], beta∈[0,0.006]의 새 합성 Gauss 적분법에서 Fermi 값·혼합 미분 28성분의 절단·구적 오차를 1e-11 이하로 제한했다. 세 지정점 84성분의 기존 네이티브 비교도 원래 문턱을 통과했다. 소수 끝점 저장과 오차 필드의 불일치 실패를 보존한 뒤, 저장 끝점까지 포함하는 정확한 유리수 상계를 별도로 내보내고 감사했다. 최대 추가 저장 오차는 약 4.033e-55다.

분류: Conjectural. 이 결과는 기존 적응형 함수의 연속 오차나 암시적 평형·물리 EOS·GR 궤적 인증을 대신하지 않는다. 전자 밀도 역산과 열역학 오차 전달로 연결하며 전체 진화·scalar·관측 목표를 유지한다. 가정·증명·실패 및 재현 근거는 [Fermi 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 이상 전자 열역학·역산 인증

분류: Proven. Fermi 인증을 전자 밀도·압력·에너지의 27개 연속 값·미분 오차로 전달했다. 양의 밀도 기울기 하한으로 선언한 공통 구적법의 역산 오차를 약 3.50e-10으로 제한했다. 원시 구간에서 유리수로 독립 재계산한 30개 열역학 구간과 세 유일 역산 근의 bracket도 통과했다. 저밀도 네이티브 대조의 약 5.49e-5 반지름을 보존하며 작은 구적 오차 상계와 구별한다.

분류: Conjectural. 이 인증은 이상 전자 성분에 한정된다. beta=0 대조는 감소 함수의 우측 극한이고, 이온화·비이상 플라스마·전체 GR 진화·비영 scalar·관측 추론은 계속 열려 있다. 식·정량 상계·독립 감사와 재현 근거는 [Fermi 및 전자 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 일관된 전자 퍼텐셜·이온화 경계

분류: Proven. 동일 압력 구적 퍼텐셜을 미분해 밀도를 정의하고 Legendre 쌍의 자유에너지 값·기울기·곡률 오차를 제한했다. 밀도 목표가 같은 서로 다른 두 근의 역함수 도함수 상대 오차도 직접 밀도 구적에서 약 1.014e-9, 퍼텐셜 구적에서 약 3.575e-8로 제한했다. 지정 구간 내부의 정확한 구적 근이라는 조건을 유지한다. 이상 이온의 가중 Hessian 양의 하한과 단조 전하 중성 방정식의 유일성을 검산했다.

분류: Conjectural. 비이상 항·분자·압력 이온화의 곡률과 물리적 오차는 추가 인증 대상이다. 이상 모형의 전역 유일성이나 원래 EOS의 수렴으로 이를 대신하지 않는다. 전체 GR 진화·scalar·관측 목표도 남는다. 가정·증명·대조와 FreeEOS 원문 경계는 [전자 및 이온화 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 실제 EOS 이온화 재고 전수 추출

분류: Counterexample candidate. 고정 기준 상태 5,735구역의 실제 중성 원자·이온·H2·H2+ 재고를 추출했다. 기존 21출력은 독립 저장 EOS 표와 전 구역에서 비트 단위로 같았고, 직렬·호출 이력·병렬 대조도 통과했다. 독립 원소 재고·전하 잔차는 각각 약 4.85e-15, 9.06e-15 이하이다. 실제 전체 분율에 접근할 수 없다는 이전 한계는 해소했다.

분류: Proven. 원본 이온화 루틴은 작은 성분의 명시적 절단과 이전 0 마스크의 재사용을 포함한다. 수치 0을 참 물리적 0이나 모든 호출에 공통인 누락 상계로 해석하지 않는다.

분류: Conjectural. 5,302구역에서 저장된 단계 분율 중 0이 발견되어 양의 내부 Hessian 공식을 적용하려면 로그 분율·절단 경계 처리가 필요하다. 비이상 항의 곡률·물리 EOS·자체 GR 진화·scalar·관측 폐쇄는 남는다. 근거와 재현은 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 절단 조건과 실제 Coulomb 곡률

분류: Proven. 새 마스크의 원자 누락 확률 상계와 재사용 마스크의 실패 반례를 검산했다. 현재 로그 가중치 비율의 인증이 필요하다는 조건은 유지한다. 원소·분자 핵 수를 보존한 Fisher 계량의 저차원 곡률 투영식을 검산했다.

분류: Counterexample candidate. 5,735구역의 실제 Coulomb 값·미분·성분 곡률을 계산하고 독립 저장 EOS 21출력, 직렬 대조, 별도 Gram·고윳값 계산으로 감사했다. 이상 이온·전자+Coulomb의 최소 곡률 여유는 약 0.999892였다.

분류: Proven. 정확한 유리수로 해석한 저장 계수·재구성 양의 종에 한정한 성분 행렬 하한은 모든 구역에서 0.999624425469816267853739080232 이상이다. 네이티브 연속 오차나 전체 물리 EOS 하한을 뜻하지 않는다.

분류: Conjectural. 현재 0 마스크·교환·압력 이온화·분배/들뜸 항과 전체 물리 오차, 자체 GR 진화·scalar·관측 연결은 계속 수행한다. 결과와 경계는 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 보존한다.


## Request 33 실제 현재 이온화 마스크 인증

분류: Counterexample candidate. 별도 계측본에서 5,735구역의 최종 로그 가중치·마스크를 추출했고, 기존 21출력·원자·이온·분자 분율을 모두 비트 단위로 재현했다. 최종 호출은 모두 마스크 재사용이며, 수치 0인 307,605개 단계가 현재 마스크와 정확히 일치했다.

분류: Proven. 현재 저장 로그 가중치로 정의한 원자 분포의 누락 확률은 전 구역에서 약 1.673e-124 이하이다. 현재 최대 로그 차이를 정확한 유리수로 비교해 이전 마스크의 발생 이력과 무관하게 인증했다. 고정 마스크에서의 조건부 기울기·Hessian 전달식도 검산했다. 하드 절단 경계에서는 양의 작은 점프가 발생할 수 있어 원래 계산기의 전 구간 C2 성질을 주장할 수 없다.

분류: Conjectural. 현재 마스크 접근의 한계는 해소했으나, 로그 가중치 자체의 평가 오차와 결합 평형·비이상 항·물리 EOS의 인증은 남는다. 자체 GR 진화·scalar·관측 목표도 유지한다. 상세 결과와 정확한 적용 범위는 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 비선형 교환 성분 곡률

분류: Counterexample candidate. 마지막 내부 교환 입력과 EOS 반환 인수의 차이를 계측했다. 원래 반환 인수 대조와 괄호 순서만 바꾼 대조의 실패를 보존했고, 실제 내부 인수의 독립 재호출은 5,735구역에서 네 교환 출력을 비트 단위로 재현했다. 반환 인수에서의 Coulomb·전자·교환 결합 성분 곡률 여유는 최소 약 0.999892175였다.

분류: Proven. 엄격히 볼록한 grand-canonical 압력의 내부 Legendre 변환은 양의 전자 압축률을 갖는다. 별도로 저장된 이상 전자+교환 계수의 양의 부호를 정확한 유리수로 확인하여, 같은 고정 종과 Coulomb 행렬의 이전 하한 약 0.9996244254를 유지했다. 선언한 유한 계수 문제에 대한 결과다.

분류: Conjectural. 압력 이온화·분배/들뜸·결합 근과 전체 물리/연속 EOS, 자체 GR 진화·scalar·관측 및 투고본은 남는다. 가정·실패·입력 차이·독립 감사는 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 실제 MDH 압력 이온화 곡률

분류: Counterexample candidate. 원래 EOS 출력을 비트 단위로 유지한 별도 계측본에서 5,735구역의 MDH 모멘트·반경·열역학 입력과 출력을 추출했다. 실제 종별 모멘트와 압력·엔트로피·공급된 전체 도함수를 재구성했고 기존 전자·Coulomb·교환에 MDH를 더한 최소 성분 곡률 여유는 약 0.9990920164였다.

분류: Proven. MDH 다항식의 동차성·Hessian을 검산하고 독립 단항식 유리수 미분과 보존 제약 Gram을 사용했다. 선언한 고정 계수·종·반경·모멘트의 하한은 0.998908632203022086641605037458 이상이다. 원래 구현의 연속 도함수나 물리 EOS의 오차 보증은 아니다.

분류: Conjectural. 밀도 의존 들뜸/분배 함수, 누락 종·전체 평형·연속 및 물리 오차, 자체 GR·scalar·관측 목표는 남는다. 세부 가정과 독립 감사는 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 보존한다.


## Request 33 실제 들뜸 항과 자유에너지 항별 검산

분류: Counterexample candidate. 첫 들뜸 미분 대조는 39구역에서 실패했다. 실제 exp(600) 스케일을 너무 일찍 제거하면서 작은 기울기/Hessian이 소실된 것이 원인이었다. 원래 실패를 고정한 뒤 원본 스케일을 그대로 저장하는 별도 계측본에서 5,735구역 모두 같은 기준을 통과했다. 원래 EOS·MDH 출력과 행렬은 비트 단위로 유지됐다.

분류: Proven. 밀도 의존 분배 함수의 종별 Hessian과 보존 제약 투영식을 검산했다. 독립 종 쌍 미분·Schur 투영·정확한 유리수 계수 감사에서 선언한 고정 행렬 하한은 0.998716614014648890147575231282 이상이다. 감사기의 작은 양의 척도를 부동소수로 바꿔 0으로 만든 실패도 보존하고, 모든 양의 유리수에 정의되는 이진 척도로 해결했다.

분류: Proven. 지정 자유에너지 분해와 고정 온도·부피·양의 종·미분 가능한 분기를 가정하면 이상 종, 전자/교환, Coulomb, MDH, 들뜸 항이 비영 종 Hessian을 구성한다. 선형 결합 에너지·통계 가중치·Planck–Larkin 항의 고정 온도 Hessian은 0이지만 온도 미분·평형 잔차는 생략할 수 없다. 실제 소스의 PL 합산에서는 정확히 한 번만 남는다.

분류: Conjectural. 네이티브 들뜸 도함수 정확도·절단/연속 영역·전체 평형 잔차와 물리 EOS, 자체 GR·scalar·관측·투고본 인증은 남는다. 상세 결과와 실패 이력은 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 실제 원자 평형의 독립 잔차 감사

분류: Counterexample candidate. 미정 비활성 PL 슬롯을 읽은 최초 실패를 보존하고 활성 원소의 정의된 슬롯만 계측했다. 5,735구역의 실제 최종 이온화 상수와 자유에너지 기울기를 비교했고 133,117개 양의 원자 종 간 비교가 통과했다. 독립 재호출에서 원래 EOS 21출력·원자 분율과 이전 들뜸 계측 결과도 비트 단위로 유지했다.

분류: Proven. 저장된 원시 계수의 유리수 기울기와 바깥쪽 반올림 로그 구간으로 계산한 원자 반응 방향의 무차원 평형 잔차는 전 구역에서 2.002165257944909382995405630549291458193E-12 이하이다. 이는 분자를 고정한 양의 저장 원자 집합의 유한 계수 문제다. 비교 방향이 없는 2,310구역의 0은 전체 평형 인증을 뜻하지 않는다.

분류: Conjectural. 분자 반응 방향·누락 종·네이티브 원시 값 및 연속 결합 근의 오차, 물리 EOS·자체 GR·scalar·관측·투고본은 남는다. 전체 목표를 위의 유한 잔차로 대체하지 않는다. 상세 정의와 실패 이력은 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 실제 분자 평형과 양의 종 집합 전체 잔차

분류: Counterexample candidate. H2/H2+의 실제 분배 함수·평형 상수·절단 전 로그 밀도를 추출하고 5,735구역을 감사했다. 2,352개 양의 분자 생성 방향 대조가 통과했고, 원래 EOS 21출력·원자/분자 분율·이전 원자 계측은 전수 비트 일치했다.

분류: Proven. 핵 수 (1,1,2,2)를 사용한 수소 종의 모든 양의 지지 집합과 기준 선택에서 보존 방향의 완전성을 기호 검산했다. 4,283개 실제 방향의 독립 구간 잔차 상계는 9.577988737788866157797082220439514904151E-13이다. 이전 원자 항과 합친 선언된 양의 종 집합의 Fisher 쌍대 잔차 노름 제곱 상계는 8.729429563898607759394223039664804087338E-30이다. 고정 계수 문제에 한정한다.

분류: Counterexample candidate. 분자 모형이 꺼진 구역은 4,559개다. 실행된 분자 분기에서는 반환 0인 분자가 없어 언더플로 누락 상계가 생성되지 않았다. 꺼진 분자 모형의 오차가 0이라는 뜻이 아니다.

분류: Conjectural. 분자 모형 전환, 누락 종·분배 함수·원시 도함수와 연속 근·물리 EOS, 자체 GR 진화·scalar·관측·투고본은 남는다. 자세한 경계는 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 있다.


## Request 33 새 GR 구조·원천 연결 완료와 분자/복사 경계

분류: Counterexample candidate. 수정 EOS로 같은 바리온·핵종 재고의 GR 구조 계산을 완료했다. 구조 접합 잔차는 3.12639e-13, 바리온 질량은 3.92510144657137e+32 g이다. 독립 재호출로 5,735구역의 EOS 21출력이 모두 비트 일치했고 엔트로피 역산의 최대 점수는 0.910353로 기준 1 이하이다. 새 반응·전자·불투명도 입력 대조를 통과했다. 이는 모형 교체에 따른 구조 재계산이며 시간 진화가 아니다.

분류: Counterexample candidate. 독립 수송 계산의 상대 차이는 1.04142e-15, 닫힌 경계의 합산 에너지 잔차는 3.83023e-20였다. 두 공간 차분이 함께 음의 국소 부력을 표시한 질량 비율은 0.014474, 가장 짧은 해당 고유시간 척도는 22.1966 s다. 전역 안정성·대류·유체·계량 진화는 아직 풀리지 않았다.

분류: Proven. 이상 혼합의 새 종 x가 자유에너지에 kT*x*(ln x-1)로 들어오고, 양의 제공 종을 쓰는 생성 방향이 가능하며 다른 항의 방향 미분이 유한하면, [F(x)-F(0)]/x는 x→0+에서 -∞로 간다. 따라서 x=0을 강제한 평형은 전체 모형의 국소 최소가 아니다. 두 제한/비제한 최소가 전환 온도 근방에서 연속이면 갑작스러운 종 제거가 자유에너지의 양의 점프를 만든다. 이는 명시한 가정 아래의 정리이며 모든 네이티브 근의 가정을 자동 인증하지 않는다.

분류: Counterexample candidate. 실제 EOS의 마지막 분자 차단은 T>=1e6 K에서 실행된다. 별도 차단 해제 모형에서 원본 출력 전수 재현 후 비교한 압력의 최대 상대 변화는 5.66519e-09였다. 경계에 접근하는 유한 수열은 비에너지·엔트로피·자유에너지의 남는 차이를 보였다. 기존 Taylor 분배 함수의 고온 연장을 물리 보정이나 전구간 미분 인증으로 채택하지 않는다.

분류: Counterexample candidate. 복사 Rosseland 불투명도를 전도와 분리했다. 유효/복사 불투명도 비율의 최솟값은 0.00349988, 복사 수송 길이/적색편이 온도 척도의 최댓값은 1.67419이다. 저장 경계 밖 광학 깊이는 미지다. 광자 확산만으로 표면 대기와 관측 신호를 완성했다고 볼 수 없다.

분류: Proven. 양의 두 주파수 군에서 동일 가중치 조화 평균이 같은 불투명도 (1,1)과 (2/3,2)는 산술 평균이 각각 1과 4/3이다. 따라서 Rosseland 평균만으로 스펙트럼 불투명도나 일반적인 흡수 방출 계수를 유일하게 복원할 수 없다. 실제 주파수 자료 또는 추가 폐쇄 가정이 필요하다는 명시적 비유일성 대조다.

분류: Conjectural. 전체 물리 EOS·연속 원시 도함수·결합 근의 오차, 실제 비선형 열수송과 유체·계량·대류·복사 대기의 결합, 비영 scalar 구동·전하, 전체 비선형 관측 추론과 최종 투고본은 남는다. 전체 목표를 유한 감사로 종료하지 않는다. 자세한 수치·재현 명령·경계는 [GR 연결 노트](../notes/REQUEST33_DIRECT_EOS_GR_KO.md) 및 [EOS 인증 노트](../notes/REQUEST33_FERMI_UNIFORM_KO.md)에 보존한다.


## Request 33 첫 비선형 열수송의 전역 보존 실패

분류: Counterexample candidate. 새 GR 상태의 닫힌 고정 계량·밀도·조성 열수송을 실제 EOS와 불투명도를 매 반복 갱신하여 실행했다. 첫 단계의 국소 열역학 에너지 잔차는 1.30432124104e-14로 사전 기준 2e-12보다 작았으나, 전역 에너지 결손/절대 교환 에너지는 **0.401342351162**로 기준 1e-8에 실패했다. 최대 logT 변화는 0.0103515545799, 이산 엔트로피 증가는 7.20461227563e+24 erg/K였다. 엔트로피 증가와 국소 Newton 수렴으로 전역 보존 실패를 구제하지 않는다. 2/4단계 경로는 아직 실행되지 않았다.

분류: Proven. 이 실패 이후 정의한 제조 산술 대조에서 초기 에너지 (2^60,4)에 정확한 증분 (1,-1)을 더하면 증분 합은 0이나, binary64 끝점을 저장한 뒤 초기값을 빼면 (0,-1)이 되어 합이 -1이다. 큰 에너지 오프셋의 차감이 작은 보존 증분을 지울 수 있음을 보여 준다. 이는 실제 실패의 원인을 모두 확정한 진단이나 사전 등록된 양성 대조가 아니다.

분류: Counterexample candidate. 원래 실패의 계획·코드·로그·잔차는 SHA로 동결했다. 원래 실행은 실패 직전의 전체 끝점 배열을 저장하지 않았으므로, 같은 계산의 Python 프레임을 읽기 전용으로 관찰하는 별도 재실행에서 이를 수집 중이다. 물리식·Newton 절차·허용 기준을 바꾸지 않고 원래 실패의 재현과 구역별 온도/에너지 저장 정밀도 기여를 확인한다. 진단이 끝나기 전 정밀도를 원인으로 단정하거나 수정 성공으로 보고하지 않는다.

분류: Conjectural. 에너지 증분을 보존하는 열역학 좌표와 전체 EOS 연결의 오차, 새 비선형 시간 경로, 유체·계량·대류·복사 대기, 물리 EOS 및 실제 scalar·관측 폐쇄를 계속 해결해야 한다. 이번 실패 때문에 전체 목표를 축소하거나 완료로 표시하지 않는다. 근거는 [GR 연결 노트](../notes/REQUEST33_DIRECT_EOS_GR_KO.md)와 `outputs/direct-eos-gr33/gr-nonlinear-thermal/failure-manifest.json`이다.


## Request 33 열수송 실패 재현과 작은 증분의 보존 좌표

분류: Counterexample candidate. 물리식·알고리즘·기준을 바꾸지 않은 읽기 전용 끝점 재실행에서 전역 에너지 결손 0.401342351162를 원래 단계 기록과 정확히 동일하게 재현했다. 필요한 온도 증분이 저장 logT 간격의 절반보다 작은 구역과 실제 저장 온도가 변하지 않은 구역이 각각 4,181개였다. 이 구역들이 전체 절대 에너지 잔차의 78.9747%를 차지했다. 전체 초기 내부에너지 끝점 한 ULP 예산을 합하면 교환 에너지의 15.9139배다. 모든 잔차를 하나의 원인으로 단정하지 않고 나머지 연결 오차도 계속 검증한다.

분류: Proven. 밀도·핵종을 고정하고 경로에서 du/dlnT=cv*T가 성립하면, 에너지 증분은 큰 절대값 차감 대신 delta*integral_0^1 (cv*T)(lnT_old+delta*x) dx로 계산할 수 있다. 2점 Gaussian 식의 3차 다항식 정확성을 검산했다. 4점 Gaussian 식은 용량의 8차 미분 상계 M8이 있으면 정확한 수학적 나머지가 |delta|^9*M8/1778112000 이하이다. 실제 EOS의 M8, 원시 평가 및 절점 반올림 오차를 아직 제공한 것은 아니다.

분류: Counterexample candidate. 별도 증분 좌표 구현은 작은 delta를 절대 logT와 합쳐 저장하지 않고, 실제 EOS 용량의 4점 적분·2점 독립 비교·원래 에너지 허용 기준·1/2/4단계 시간 비교를 사용하도록 고정했다. 경계 유속도 작은 두 온도 증분의 차이를 별도로 유지한다. 영 증분은 기존 경계 값/기울기를 재현했고, logT=16에 더하면 사라지는 1e-18 증분을 분리한 식은 0이 아닌 유속을 보존했다. 준비 대조를 전체 새 시간 경로의 통과로 세지 않는다.

분류: Conjectural. 적분 좌표와 실제 EOS 자유에너지 사이의 연속 도함수·분기/절단 오차, 새 비선형 열수송 결과, 유체·계량·조성 반응·대류·복사 대기, 비영 scalar 및 관측·투고본 폐쇄는 계속 남는다. 원래 절대 에너지 끝점 시험의 실패는 철회하지 않는다. 재현 코드와 원시 끝점은 [GR 연결 노트](../notes/REQUEST33_DIRECT_EOS_GR_KO.md)의 산출물에 있다.


## Request 33 열 증분 1단계 감사와 초기 GR 열유속 제약

분류: Counterexample candidate. 증분 좌표의 닫힌 고정 밀도·계량·조성 열수송 1단계가 원래 기준을 통과했다. 전역 에너지 잔차/교환량은 5.38972e-17, 국소 잔차는 1.23966e-16, 2점/4점 적분 차이/교환량은 4.96097e-10이며 이산 엔트로피는 증가했다. 2/4단계 시간 경로는 별도 완료 판정을 기다린다. 원래 절대 에너지 끝점 시험의 0.401342 보존 실패는 유지한다.

분류: Proven. 저장된 질량·적색편이 가중치와 확장 정밀도 증분을 정확한 유리수로 변환해 감사했다. 1단계 전역 에너지 잔차/교환량의 상계는 5.389487862268683e-17, 국소 잔차 상계는 1.239650802338870e-16이다. 저장값의 최종 온도·에너지 증분을 재현했고 사후 1% 에너지 오염 대조는 거부됐다. 네이티브 EOS 원시 평가, 연속 적분 및 물리 EOS 오차의 인증은 아니다.

분류: Proven. 부호 규약 K_ij=-(1/(2N))*partial_(ct) gamma_ij, 구면 polar areal 좌표, 초기 물질 속도 0에서, 고유 열유속 q가 있으면 운동량 제약은 K^r_r=A=4*pi*r*a*q를 요구한다. 각방향 외재곡률이 0이면 추가 Hamiltonian 항은 항등적으로 0이다. dln(a)/dt=-c*N*A, dm/dt=-4*pi*r^2*N*q*c/a이며, 초기 Hamiltonian 전파 항등식을 기호 검산했다. q는 이 식에서 기하학 단위이다.

분류: Proven. 처음 작성한 물질 에너지율 식에는 가속도 항이 빠졌다. 초기 v=0은 v_dot=0을 뜻하지 않는다. Q=고유 열유속/c, v=물질 속도/c이면 E_normal_dot=epsilon_dot+2*Q*v_dot이므로 물질 비에너지율에 -(2*Q/rho_B)*v_dot을 더해야 한다. 초기 정역학 상태의 운동량식은 Q_dot+(epsilon+P)*v_dot=2*c*N*A*Q이다. 원래 자료를 보존하고 velocity-derivative-correction.json과 독립 Lorentz 텐서 검산을 추가했다. 기존 밀도율·초기 제약은 유효하나 원래 저장 물질 온도율은 v_dot=0에만 조건부로 유효하다.

분류: Counterexample candidate. 실제 5,735구역의 비영 열유속에 초기 외재곡률을 연결한 점별 운동량 잔차는 3.17226e-16이다. 최대 |dlnrho/dt|는 5.06399e-24 /s, 최대 |dm/dt|는 8.96794e-15 cm/s다. 기존 공간 Hamiltonian 이산 오차나 후속 유체·계량 진화를 인증한 것은 아니다.

분류: Proven. 국소 관성계의 균질 선형 baryon-frame 열유속 모형 w*V_dot+f_dot=0, tau*f_dot+f=-(K*T/c^2)*V_dot에서 비영 극은 -1/(tau-tau_min), tau_min=K*T/(w*c^2)이다. tau=0이면 성장하며, 이 모드의 감쇠에는 tau>tau_min이 필요하다. tau=2*tau_min과 3*tau_min은 동일한 평형 EOS/전도율에 서로 다른 감쇠 극을 준다. 정적 미시입력으로 과도 수송 계수를 유일하게 정할 수 없다. 이는 제한한 모형의 필요조건이며 모든 유체 이론·공간 파수·비선형 인과성의 정리는 아니다.

분류: Counterexample candidate. 실제 저장 상태에서 균질 필요 시간의 범위는 1.95690e-25–1.76525e-12 s다. 별도 고정 유체 Cattaneo 열방정식의 아광속 조건 K/(C_volume*c^2)는 1.48037e-19–2.05650e-4 s다. 두 값은 서로 다른 제한계의 진단이며 실제 열유속 완화시간 또는 거시적 chi 완화시간의 측정값이 아니다.

분류: Conjectural. 결합 유체의 유한 파수·인과성, 실제 열유속 폐쇄 계수와 복사 대기, 전체 EOS 연속 원시/근/도함수 오차, 유체·계량·반응의 자체 진화, 비영 scalar 및 전체 비선형 관측·최종 투고본은 계속 남는다.

재현: verification/gr_heat_initial_constraints.py verify, verification/gr_heat_fluid_frame.py verify, verification/gr_heat_relaxation_boundary.py verify, verification/verify_gr_caloric_increment.py verify. 자료는 outputs/direct-eos-gr33/의 같은 이름 하위 폴더에 보존한다.


## Request 33 결합 열유체의 국소 전파속도와 시간 분할 미달

분류: Proven. 국소 평형 Cattaneo 유체의 바리온·에너지·운동량·열유속 4변수 특성행렬과 모든 비영 실수 파수에 대한 Hurwitz 조건을 기호 검산했다. w=epsilon+P, b=rho*cv*T/w, r=P_lnrho/w, d=P_lnT/w, e=(P-rho*u_lnrho)/w, h=K*T/(w*c^2), g=r+d*e/b로 둔다. b,h,r>0, tau>h, g>0, 1-e-d+b*g>0, g*(1-e-d+b*g)-r>0이면 선언한 균질 평형 선형계의 모든 비영 파수 모드는 감쇠한다. 정확한 열역학 e=d에서 마지막 인자는 (b*r-d*(1-d))^2/b다. 퇴화하여 0인 경우를 엄격 감쇠 정리로 포함하지 않는다.

분류: Proven. 같은 계의 전파속도 제곱 z=(speed/c)^2는 A*z^2-C*z+r=0, A=b*(lambda-1), C=lambda*(b*r+d*e)+1-e-d, lambda=tau/h를 만족한다. A,C,r>0, 판별식>=0, C<=2*A, A-C+r>=0이면 두 z 모두 (0,1]에 있다. 정적 입력으로 정해지는 조건을 만족하는 서로 다른 tau를 만들 수 있으므로 안정성과 전파속도 검증도 물리 완화시간을 유일하게 결정하지 않는다.

분류: Proven. 실제 저장 도함수에서 e=d를 강제하지 않고 5,735구역을 정확한 유리수로 검사했다. 사전에 정의한 두 완화시간 후보가 모든 구역에서 위 엄격 감쇠 및 아광속 조건을 통과했다. 저장 계수의 산술 명제이며 네이티브 도함수/물리 불확실성까지 둘러싼 구간은 아니다. 후보 시간의 표시 범위는 각각 2.96074e-19–4.11300e-4 s와 4.44111e-19–6.16950e-4 s, 최대 전파속도 표시는 각각 약 0.707107c와 0.577351c다. 이는 실제 국소 열유속 시간이나 거시적 chi 시간의 측정이 아니다.

분류: Proven. 초기 비영 열유속 Q까지 포함하면 에너지 주부에 2*(Q/w)*v_dot, 운동량 주부에 2*c*(Q/w)*v_x가 추가된다. 두 후보×5,735구역의 11,470개 특성 사차식 각각에 대해 (-1,1) 안의 서로 겹치지 않는 네 유리수 구간에서 부호 변화를 확인했다. 중간값 정리와 차수에 따라 저장된 각 주부 행렬의 네 근은 모두 실수이며 광속보다 작다. 근삿값은 구간 제안에만 쓰고 정확한 끝점 다항식 값으로 인증했다. 비평형 배경의 하위차 항을 포함한 모든 파수 안정성이나 비선형 GR 안정성의 인증은 아니다.

분류: Proven. 명시한 covariant Cattaneo 법칙을 초기 v=0에서 축약하면 tau*Q_dot+w*h*v_dot=N*(Q_F-Q)다. 운동량식 Q_dot+w*v_dot=R, R=2*c*N*A*Q와 함께 풀었다. 기존 이산 열유속을 초기 Fourier 표적으로 쓰는 Q_F=Q 조건에서 Q_dot=-h*R/(tau-h), v_dot=tau*R/[w*(tau-h)]이며 두 식의 정확 저장값 잔차는 0이다. v_dot의 최대 크기 표시는 2.33091e-35 /s, 빠졌던 온도율 보정은 최대 9.65405e-37 /s다. 초기 미분이므로 아직 유한 시간의 유체·계량 경로가 아니다. Q_F=Q는 독립 연속 온도 기울기나 대기 검증을 대신하지 않는다.

분류: Counterexample candidate. 2단계 열 증분 경로의 에너지 보존과 양의 이산 엔트로피, 2/4점 적분 대조는 통과했다. 전 경로의 저장값 정확 산술 에너지 잔차/교환량 상계는 7.159308822931545e-17이다. 그러나 1단계 대비 최대 logT 차이는 **0.0010524272195028663**으로 사전 시간 분할 기준 **0.0001**에 미달했다. 에너지 보존 통과를 시간 정확도 통과로 해석하지 않는다. 원래 결과·기준을 유지하며 4단계 경로를 계산한다.

분류: Conjectural. 실제 과도 수송 계수의 보정, 유한 시간의 자체 유체·계량·반응·대류·대기 진화, 물리 EOS 및 연속 평가/도함수/근/시간 오차, 비영 scalar 구동·전하와 완전한 관측·최종 투고본은 계속 남는다. 국소 모형의 정확 산술 인증을 이 목표의 대체물로 삼지 않는다.

재현: verification/gr_heat_characteristics.py verify, verification/gr_heat_initial_tangent.py verify, verification/verify_gr_caloric_increment.py verify. 각각의 plan, symbolic, result, SHA manifest와 유리수 근 구간을 outputs/direct-eos-gr33/에 저장했다. GR 투영의 배경은 [Gourgoulhon의 원문](https://arxiv.org/abs/gr-qc/0703035), 열유속 불안정성의 선행 연구는 [Hiscock–Lindblom 원문](https://ccom.ucsd.edu/~lindblom/Publications/24_PhysRevD.31.725.pdf)이며, 여기서 사용한 제한계의 식과 조건은 별도로 유도·검산했다.


## Request 33 이차 엔트로피 모형의 전체 열유속 항

분류: Imported from prior work. [Maartens 원문](https://arxiv.org/pdf/astro-ph/9609119)의 식 (2.17), (2.20), (2.22), (2.24)를 대조했다. 열유속만 남긴 이차 엔트로피 전류와 그 전체 수송식을 사용하며 점성·열 교차 결합과 회전을 제외한다. 별도 물리 보정 없이 이를 실제 항성의 완전한 미시 이론이라고 주장하지 않는다.

분류: Proven. c=1 단위에서 S^mu=s*n*u^mu+q^mu/T-[tau/(2*K*T^2)]*q^2*u^mu를 가정하고, 보존식·미분 가능한 Gibbs 관계를 쓰면, 전체 열유속식에 -(K*T^2/2)*q^mu*div[(tau/(K*T^2))*u]를 포함했을 때 div S=q^2/(K*T^2)>=0이다. 이 항을 생략하면 추가 부호 미정 항 -(q^2/2)*div[(tau/(K*T^2))*u]가 남는다. 선언된 국소 미분값에서 음의 엔트로피 생성률을 만드는 대수적 대조를 검산했다. 기존 항성 자료에서 엔트로피 위반을 관측했다는 뜻은 아니다.

분류: Counterexample candidate. 초기 구역의 전파속도 조건으로 만든 최대 기준 시간에 각각 2와 3을 곱한, 공간·시간에 일정한 두 proper 시간 모형을 정의했다. 약 4.11300e-4 s와 6.16950e-4 s이며 이전의 공간 의존 후보와 구분한다. 물리 충돌시간 또는 거시적 chi 시간으로 보정된 값이 아니다. 열전도율은 기존 복사/전도 표에서 매 상태 평가되는 동일 K를 쓴다.

분류: Proven. 이 두 모형에서 전체 항이 바꾸는 시간 주부를 다시 유도하고, 실제 비영 열유속을 포함한 11,470개 저장 특성 사차식의 네 근을 정확한 유리수 부호 구간으로 모두 (-1,1) 안에 인증했다. Q_dot·v_dot·T_dot을 동시에 푼 저장 계수 방정식의 잔차는 정확히 0이다. 전체 항을 넣은 초기 v_dot 최대 표시는 2.57797e-7 /s, 물질 온도율의 가속도 보정 최대 표시는 1.06773e-8 /s다. 일정한 tau 아래 K와 T의 미분을 포함했으며 이전의 단순 Cattaneo 초기율로 대체하지 않는다.

분류: Counterexample candidate. 전체 항의 크기를 Q 항과 비교한 저장 상태의 최댓값은 **0.919130**이다. 따라서 이번 후보에서 그 항을 작은 보정으로 생략할 근거는 없다. 기존 단순 완화모형의 인증은 그 선언된 모형에만 남겨 두고, 엔트로피 전류 가정이 추가된 모형의 결과를 별도로 보존한다.

분류: Conjectural. 이차 엔트로피 근사 자체의 물리 유효범위·과도 계수 보정·실제 EOS의 Gibbs/도함수 오차, 비평형 전 파수 안정성과 비선형 해의 존재·시간 오차, 유한 시간의 자체 GR 유체·계량·반응·대류·대기 및 비영 scalar·관측 연결은 계속 남는다. 양의 엔트로피 항등식과 초기 특성 구간을 실제 전체 진화의 완료로 표시하지 않는다.

재현: verification/gr_heat_entropy_closure.py verify. 원문 PDF·계획·기호식·정확한 일정 시간·근 구간·초기 변화율과 SHA manifest를 outputs/direct-eos-gr33/gr-heat-entropy-closure/에 보존한다. 기존 시간 분할 미달도 유지한다.


## Request 33 비선형 구면 GR 질량수지와 닫힌 표면 영 대조

분류: Proven. 임의의 시간 의존 구면 polar-areal 4차원 계량에서 Ricci 텐서를 직접 계산했다. c=G=1, m=r*(1-a^-2)/2이면 G^t_t=-2*m_r/r^2, G^r_t=2*m_t/r^2, G^r_r=-2*m/r^3+2*N_r/(r*N*a^2)다. 따라서 m_r=4*pi*r^2*E, m_t=-4*pi*r^2*N*J/a, N_r/N=a^2*(m/r^2+4*pi*r*S)를 얻는다. 여기서 E,J,S는 normal 관측자의 에너지·운동량·방사 방향 응력이며, 정적 계량이나 영 물질 속도를 가정하지 않는다.

분류: Proven. 물질 표면 R_dot=N*v/a에서 Lorentz 변환된 텐서는 v*E-J=-P*v-Q를 정확히 만족한다. 따라서 d m(t,R(t))/dt=-4*pi*R^2*N/a*(P*v+Q)다. 표면 P=Q=0이고 다른 에너지·응력 통로가 없다면 비선형 내부 운동·열 재분배·내부 전환이 있어도 총 중력 질량은 일정하다. 연결된 구면 진공 외부의 m_r=m_t=0 및 lapse 적분은 시간 좌표 재정의 후 일정 질량 Schwarzschild 계량을 준다. 이는 알려진 구면 진공 경계의 직접 검산과 이 열유속 모형에 대한 적용이며 새로운 GR 정리라는 주장은 아니다.

분류: Proven. 따라서 닫힌 순수 GR 구면 부문에서 내부 열완화 상태만 추가해 외부 질량 단극 변조나 비영 scalar 전하를 얻는 경로는 실패한다. 이 경계를 벗어나려면 경계 에너지/일, 외부 결합, 비구면·다중극, 추가 장 또는 응력 중 어떤 가정을 바꾸는지 명시해야 한다. 임의 외부장에 대한 보편 SEP 정리로 확대하지 않는다.

분류: Conjectural. 저장된 외곽 EOS 경계는 아직 P=0인 물질 표면과 진공 외부에 접합됐다고 인증되지 않았다. 계산상의 두 열 경계를 닫는 것으로 실제 표면·복사 대기 조건을 충족했다고 보지 않는다. 이 정리를 실제 표면 인증이나 완전한 관측 예측으로 사용할 수는 없다.

분류: Counterexample candidate. 1/2단계 시간 분할 미달을 보존하고 동일한 직접 EOS 알고리즘·시간 길이·열 연산자·Newton 기준·에너지 기준·온도 차이 기준 1e-4로 별도 8/16/32단계 경로를 시작했다. 새 물리식이나 허용 기준 변경으로 실패를 구제하지 않는다. 4단계 결과는 완성 후 첫 8단계 감사에 원래 SHA를 확인하여 가져오며 새 계산 결과로 세지 않는다. 시작·준비 검산은 세 경로의 완료 또는 시간 정확도 인증이 아니다.

분류: Conjectural. 연속 EOS/원시/도함수/근/시간 오차, 과도 수송의 물리 보정, 유한 시간의 자체 GR 유체·계량·반응·대류·대기 경로, 비영 scalar 구동·전하와 실제 관측·최종 투고본은 계속 열려 있다. 순수 GR의 닫힌 질량 영 대조와 추가 결합이 필요한 관측 후보를 혼동하지 않는다.

재현: verification/gr_spherical_mass_balance.py verify. 추가 시간 경로: verification/gr_caloric_refinement.py run, 완료 경로의 저장값 감사: verification/gr_caloric_refinement.py audit N. 같은 run 명령을 실행 중인 작업에 중복 실행하지 않는다. 코드·계획·기호식은 outputs/direct-eos-gr33/gr-spherical-mass-balance/와 gr-caloric-refinement/에 연결되어 있다.


## Request 33 시간 분할 4단계와 직접 EOS 구역 적분

분류: Counterexample candidate. 기존 1/2/4단계 계산이 끝났다. 4단계 경로도 에너지·엔트로피·유한 2/4점 적분 대조는 통과했으나 2단계와의 최대 logT 차이는 **0.000657217536851674**로 시간 기준 0.0001에 미달했다. 정확한 저장값 유리수 연산에서 전 경로 에너지 잔차/교환량 상계는 6.413669804637097e-17이었다. 별도 8/16/32단계 계산은 같은 알고리즘·기준으로 계속한다.

분류: Counterexample candidate. 기존 1/2단계 시간 기준에 실패한 곳은 가장 바깥의 세 구역이다. 바리온 질량 비율은 3.11728e-14이지만 최대 오차 기준의 실패는 유지한다. 첫 구역의 Rosseland 수송 길이/온도 척도 비는 1.67419이며 저장 경계에서 첫 중간점까지의 광학 깊이는 0.0250585다. 이를 무한대에서 측정한 깊이나 완성된 대기 모형으로 해석하지 않는다.

분류: Counterexample candidate. 저장 점 압력의 3점/5점 차분과 TOV 압력 기울기를 비교한 최대 상대 차이는 각각 0.169957와 0.045477이었다. 여기서 생기는 수치적 가속도 불균형은 열유속에 의한 초기 가속도보다 컸다. 점 밀도를 동일 셀의 평균으로 놓으면 바리온 질량이 최대 2.43576% 달라진다. 따라서 현재 점 표본을 그대로 보존형 유체 적분기의 셀 평균으로 사용하거나 차분 불균형을 실제 가속도로 읽는 것은 정당화되지 않는다.

분류: Proven. 고정 엔트로피·조성의 제1법칙으로 d²[rho*(C_X*c²+u)]/drho²=(1/rho)*(dP/drho)_s를 얻는다. 이 값이 밀도 구간에서 양수이면 Jensen 부등식에 따라 비균일 등엔트로피 셀의 평균 에너지는 동일 평균 밀도·엔트로피의 균일 EOS 에너지보다 크다. 바리온 질량·부피·에너지·원래 엔트로피를 모두 보존하는 단일 균일 상태는 일반적으로 존재하지 않는다. 밀도 1과 3의 두 같은 부피 gamma=2 폴리트로프에서 이 차이를 정확히 검산했다. 실제 EOS의 전 구간 볼록성이나 평균화로 생긴 엔트로피를 물리적 열생성으로 인증하는 정리는 아니다.

분류: Counterexample candidate. 셀 내부 구조를 유지하기 위해 0,1,2,2972 구역에서 각 RHS마다 직접 EOS 엔트로피 역산을 수행하는 압력 좌표 GR 적분을 실행했다. 저장 lapse와 등엔트로피 엔탈피 관계로 경계 압력을 복원하고 반지름·중력질량·바리온 질량·고유 부피·내부에너지·압력 적분을 함께 구했다. 총 747회 엔트로피 역산, 3,018회 EOS 평가의 최대 원래 잔차 점수는 0.508167로 1 이하다. 외곽 세 구역의 적분 바리온 질량은 기존 재고와 상대 약 1.05e-9 이내, 내부 대조 구역은 약 1.31e-12 이내로 일치했다. 두 DOP853 허용 오차 대조는 연속 공간 오차의 구간 인증이 아니다.

분류: Counterexample candidate. 첫 구역의 직접 적분 평균 밀도는 1.85766e-9 g/cm³이며, 같은 평균 밀도·원래 엔트로피의 균일 EOS 상태는 적분 평균보다 비에너지가 약 3.95175e10 erg/g, 압력이 약 1.01113% 작았다. 첫 구역의 직접 적분 중력질량 두께는 3.05004e-10 cm인데 저장 면 질량의 차감값은 2.25555e-10 cm다. 총질량의 극소 차분에서 생기는 산술 오차와 실제 적분 오차를 다음 대조에서 분리한다. 이 차이를 곧바로 별 전체 질량의 같은 비율 오차로 확대하지 않는다.

분류: Conjectural. 전 구역의 보존형 재구성·평형을 보존하는 공간 적분, 중력질량의 작은 증분 처리·연속 오차, 실제 대기·수송 계수·물리 EOS·자체 GR 시간 경로와 scalar·관측 연결은 남는다. 중간점 구역 해와 실제 보존형 유체 상태를 구분하며 원래 시간 미달도 유지한다.

재현: verification/gr_spatial_preflight.py verify, verification/gr_direct_cell_integrals.py verify, verification/gr_cell_average_identity.py verify, verification/verify_gr_caloric_increment.py verify. 원시 경계 복원·각 직접 적분·평균화 대조·실패 구역과 SHA는 outputs/direct-eos-gr33/의 해당 하위 폴더에 보존했다.


## Request 33 평형을 보존하는 유체 보존식과 중단된 계산 복구

분류: Proven. 임의의 매끄러운 구면 polar-areal 계량에서 혼합 stress tensor의 공변 발산을 직접 계산했다(c=G=1). 정상 관측자 계수 E,J,S와 물질의 횡방향 압력 P에 대해 바리온 식 partial_t(a*r²*D)+partial_r(N*r²*D*v)=0, 운동량 식 partial_t(a*r²*J)+partial_r(N*r²*S)=-r²*E*N_r+2*N*r*P-r²*J*a_t를 얻었다. D=rho_B*W다. 에너지 투영에 Einstein 제약을 함께 적용하면 partial_t(r²*E)+partial_r(r²*N*J/a)=0이므로 작은 구역 중력질량을 공유 경계 유속으로 진화시킬 수 있다.

분류: Proven. 정수압 기준해 (N0,E0,P0)가 P0_r=-(E0+P0)*N0_r/N0를 만족하면 그 운동량 원천의 적분은 r²*N0*P0의 경계 차이와 같다. 이 기준 성분을 같은 공유 경계 값으로 평가하고 나머지 변화량만 구적하면, 정확히 표현된 영 열유속 정수압 기준에서 수치 운동량 갱신은 대수적으로 0이다. 이는 기준해의 재구성 오차가 0이라는 주장이 아니다. 실제 비영 J=Q와 계량 시간 항을 기준해에서 빼서 없애지 않으며 물리적 열 가속도를 유지한다.

분류: Proven. 물질 경계 속도가 R_dot=N*v/a일 때 움직이는 구역의 바리온 재고는 일정하다. 같은 경계의 중력질량 유속은 F_M=4*pi*r²*N/a*(P*v+Q), 운동량 유속은 F_Pi=4*pi*N*r²*(P+v*Q)다. 이차 열 엔트로피 모형의 구역 엔트로피는 Sigma=4*pi*integral[a*r²*W*(rho_B*s-tau*Q²/(2*K*T²)+v*Q/T)dr]이고, 경계 유속은 F_S=4*pi*N*r²*Q/(W*T), 부피 생성률은 4*pi*integral[N*a*r²*Q²/(K*T²)dr]다. 이들은 선언된 연속 보존식·Gibbs 관계·전체 열수송 모형 아래의 항등식이며 유한 단계의 엔트로피 안정성 보증은 아니다.

분류: Proven. 이웃 구역이 같은 경계 값을 공유하면 중력질량 유속은 전체 합에서 내부 경계별로 상쇄된다. 열 경계를 닫아도 움직이는 유한 압력 경계의 P*v 일 항은 남는다. 기존 닫힌 고정 구조 열 연산자의 보존을 움직이는 실제 별의 질량 일정성으로 옮기지 않는다. 작은 구역 질량 자체를 저장하는 표현과 셀 내부 열역학 모멘트가 필요하다.

분류: Counterexample candidate. 외곽 세 구역의 기존 RK4 덧셈을 그대로 재현하고 각 제안 질량 증분·실제 저장 증분을 보존하는 대조를 실행 중이다. 별도 작은 증분 표현 및 독립 압력 좌표 EOS 구역 적분과의 비교를 사전 등록했다. 통과 후 같은 바리온·조성·엔트로피와 원래 중심/표면 매개변수로 5,735구역 전체의 4/8 공간 세분 경로를 재계산하도록 준비했다. 현재는 준비·실행 단계이며 질량 오차 원인 확정이나 전체 재계산 통과를 주장하지 않는다.

분류: Proven. 실행 환경 재시작 뒤 원래 과학 계산 PID들이 사라졌고 새 WSL 가동 시간이 확인됐다. 남아 있던 8분할의 첫 네 승인 단계에서 누적 온도 증분 연쇄와 저장 에너지 잔차를 비트 단위로 재현했다. 원래 알고리즘 본문에 승인된 단계 초기화·반복 시작·기존 완료 경로 보존의 세 명시적 치환만 적용하고, 이를 되돌리면 원본 소스 및 AST가 완전히 같음을 검사했다. EOS·Newton·불투명도·구적·보존·시간 허용 기준은 바꾸지 않았다.

분류: Counterexample candidate. 네 단계에서 복구한 8/16/32 시간 계산과 아직 결과 저장 전이던 질량 덧셈 대조를 재실행했다. 이전 1/2 및 2/4 시간 미달은 그대로 유지한다. 재시작으로 잃은 미승인 계산을 완료 단계로 세지 않는다.

분류: Conjectural. 실제 전 구역 보존형 재구성·유체/열 primitive 역산·자체 계량/물질 경계 시간 경로·대기와 물리 수송 계수, 전체 EOS 및 연속 오차 인증·비영 scalar 구동/전하·완전한 관측 추론과 최종 투고본은 계속 남는다. 방정식 검산과 계산 복구를 이 전체 목표의 완료로 대체하지 않는다.

재현: verification/gr_balanced_conservation.py verify, verification/gr_material_conservation.py verify. 복구 근거는 outputs/direct-eos-gr33/gr-caloric-refinement/recovery-1/의 원본/복구 함수·승인 단계 스냅샷·원시 증분·SHA에 있다. verification/gr_caloric_resume.py run과 질량 대조는 실행 중 중복 시작하지 않는다. 전체 별 후보의 진입 검사와 기준은 gr-increment-structure/plan.json에 저장했다.


## Request 33 열유속을 포함한 비선형 원시 변수 역산과 질량 덧셈 손실

분류: Proven. 고정 조성과 comoving 열유속 Q에서 D=rho_B*W, J-Q=v*(E+P), epsilon=E-v*(J+Q)를 기호 검산했다. 따라서 D,E,J,Q가 주어지면 EOS와 함께 밀도·온도·속도를 연결할 수 있다. Q가 별도 미지수인데 세 보존량만으로 이를 결정할 수 있다는 주장은 아니다.

분류: Proven. w=epsilon+P, b=rho_B*cv*T/w, r=P_lnrho/w, d=P_lnT/w, e=(P-rho_B*u_lnrho)/w, j=Q/w라 놓으면 (lnrho,lnT,v)에서 (D,E,J)로 가는 Jacobian 행렬식은 rho_B*w²*W⁵*[b-v²*(b*r+d*e)+2*j*v*(b-d)]다. v=0에서 rho_B²*cv*T*w>0이며, 미분 가능한 EOS 아래 국소 역함수 조건을 준다. 실제 저장 계수 5,735개 각각에서 |v|<=1/2의 동결 다항식 하한을 정확한 유리수로 검사했고, b로 나눈 최소 하한 표시는 0.999996218590158이었다. EOS 계수가 변하는 실제 상태 근방이나 전역 역산의 유일성 인증은 아니다.

분류: Proven. 정지질량 크기의 두 값을 빼는 대신 E-D*C=D*C*(W-1)+rho_B*u*W²+P*v²*W²+2*v*Q*W², C=C_X*c²를 사용하고 W-1=v²*W²/(W+1)로 계산한다. 명시적 1e-25 열에너지 대조는 보존됐으며 큰 정지에너지에 더한 뒤 빼는 연산은 이를 소실했다. 이미 반올림되어 사라진 입력 정보를 되살리는 방법은 아니다.

분류: Counterexample candidate. 실제 EOS의 아홉 구역에 밀도·온도 변화와 v=-1e-4,0,+1e-4를 적용한 27개 점 상태를 만들고 보존량에서 원시 변수를 역산했다. 원래 hybrid 반복은 25/27 통과했고 내부 2972 구역의 두 속도 사례에서 정체했다. 동일 식·초기값·문턱에 단위 변수 척도만 명시한 별도 대조는 23/27로 전체 문제를 해결하지 못했다. 초기 Jacobian의 독립 두 해상도 차분은 모두 문턱 1e-5를 통과했으며 최대 행 상대 차이는 6.41532e-8, 초기 행렬 조건수 최댓값은 1.04241이었다. 작은 초기 조건수만으로 유한 비선형 반복의 수렴을 보장하지 않는다.

분류: Counterexample candidate. 같은 27개 보존 목표를 비트 단위로 유지하고, 최대 척도 잔차가 감소하는 단계만 받는 결합 Newton 계산을 별도로 수행했다. 27/27 모두 원래 밀도·온도·속도 및 잔차 기준을 통과했다. 최대 logT 오차는 9.41469e-13, 최대 척도 잔차는 9.41514e-13이었고 저장 밀도·속도 대조 차이는 0이었다. 최대 26회의 승인 반복이 필요했다. 완전 Newton 단계가 잔차를 증가시키는 사례와 감쇠 선택을 모두 기록했다. 이것은 원래 hybrid 내부 정체의 모든 원인을 특정했다는 주장이 아니며, 초기 성공/실패와 단위 척도 실패를 보존한다.

분류: Proven. 별도로 기존 RK4 질량 덧셈을 세 외곽 구역에서 재현했고 원래 step 및 저장 경계와 비트 일치를 확인했다. 첫 구역의 614회 덧셈 중 579회에서 질량 저장값이 바뀌지 않았다. 제안된 감소량 합은 약 3.05004065e-10 cm, 실제 저장 감소는 약 2.26528147e-10 cm이며, 제안량 대비 누락률은 정확히 104701738873659778876441/406933193777317072552985, 약 25.7295%다. 이는 저장 숫자의 덧셈 단계만의 정확 산술 판정으로 초기 seed 변환과 최종 단위/곱셈 반올림은 별도다. 다음 두 구역은 각각 약 1.85673% 과대 감소, 0.265284% 과소 감소였다. 이 상대 오차를 별 전체 질량 오차로 확대하지 않는다.

분류: Counterexample candidate. 작은 질량 증분을 별도로 누적하는 4/8 세분 대조는 실행 중이며, 원래 전체 별 상태를 바꾸지 않았다. 이 대조의 완료·독립 EOS 적분 기준 통과 후 전체 5,735구역 GR 공간 재계산을 진행한다. 복구한 열 시간 경로도 기존 기준으로 계속한다.

분류: Conjectural. 점 상태 역산은 비균일 셀의 모든 보존 모멘트를 만족하는 재구성이 아니다. 다음은 셀 내부 기준 구조를 유지한 보존형 역산과 경계 유속, 자체 유체·열·계량·대기 시간 적분의 결합이다. 물리 EOS·연속 평가/도함수/근/시간 오차·실제 비영 scalar 구동/전하·관측 추론과 최종 투고본은 미완료다. 27개 점 상태의 회복을 이 전체 목표로 대체하지 않는다.

재현: verification/gr_heat_primitive_inverse.py verify, verification/gr_heat_primitive_scaling.py verify, verification/gr_heat_primitive_newton.py verify, verification/gr_numerical_audit.py verify. 원래 실패·새 반복 궤적·동일 보존량·독립 차분·정확 질량 덧셈 감사와 SHA를 outputs/direct-eos-gr33/의 대응 하위 폴더에 보존했다.


## Request 33 비균일 셀 보존량 재구성

분류: Proven. 지정한 비균일 기준 셀 재구성에서 양의 구적 가중치와 rho,cvT,enthalpy 아래 기준점의 보존 모멘트 Jacobian 행렬식 B0*Cv_coordinate*enthalpy_proper는 양수다. 네 셀의 저장 인자를 정확 유리수로 검사했다. 단일 균일 EOS 치환의 Jensen 경계는 유지한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 열 8분할 감사와 전체 공간 재계산 착수

분류: Counterexample candidate. 열 8분할은 정확 저장값 에너지 감사 상계 1.13337e-11로 보존 기준을 통과했으나 4/8 시간 차이 0.000378246은 원래 0.0001 기준에 미달했다. 동일한 16/32 계산을 계속하며 질량 증분 보정의 두 공간 대조 통과 후 전체 5,735구역 정수압 재계산을 시작했다. 분류: Conjectural. 실제 GR 시간 진화·물리 EOS·관측 완료는 아직 입증되지 않았다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전체 물질 기준과 계량 정지에너지 변화

분류: Proven. 물질 측도의 정지에너지 가중치 kappa=average_B(1/a)의 계량 변화 항과 소실을 피하는 유한 증분식을 검산했다. 분류: Counterexample candidate. 전체 5735구역의 4세분 정수압 연결과 중심 포함 여섯 구역의 비균일 구적 사전 대조가 통과하여 전 구역 구적을 시작했다. 분류: Conjectural. 자체 GR 진화·물리 EOS·연속 오차·관측 폐쇄는 계속 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전체 공간 대조 완료와 공유 열 경계

분류: Counterexample candidate. 전체 5735구역의 4/8 정수압 연결과 면 좌표 차이 6.95940e-10<1e-8이 통과했고, 원래 물질 재고와 저장 구역 질량의 정확 산술을 독립 감사했다. 분류: Proven. 공유 경계의 바리온 선형 광도 보간은 규칙적인 중심에서 Q=O(r)를 준다. 분류: Counterexample candidate. 전체 공유 면 및 여섯 비균일 셀의 초기 열/곡률 대조를 통과했다. 분류: Conjectural. 전체 물리 EOS·자체 GR 시간 경로·관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 현재 EOS 정수압 기준의 선형 scalar 연결

분류: Counterexample candidate. 현재 EOS 정수압 기준을 기존 beta=-4 선형 scalar 식에 연결하고 선형/Riccati 및 두 ODE 허용 오차 대조가 1e-8 기준을 통과했다. 저장 응답 alpha_A/phi_infinity는 약 -4.00035656이다. 실패한 EOS 비트 일치 진입 검사와 두 호출의 유한 에너지 차이를 보존했다. 분류: Conjectural. 정적 선형 전하 계수는 실제 궤도 구동·내부 chi_A 완화·유한 backreaction 또는 관측 추론의 완료가 아니다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 밀도·온도 결합 증분과 실제 EOS 분기 실패

분류: Proven. 밀도·온도 결합 EOS 증분과 분기 점프가 필요한 경계를 유도했다. 분류: Counterexample candidate. 실제 native EOS는 54개 일반 대조를 통과하지만 1e6 K 경계 경로는 구적 일치에도 끝점 검사에 실패했다. 기존 분자 유지 비교는 같은 55개 검사를 통과하나 물리 EOS 인증이 아니다. 분류: Proven. 저장 열 경로 86025구간은 이 경계를 건너지 않지만 정수압 EOS 표 8개 셀은 경계를 포함한다. 분류: Conjectural. 물리 분배 함수·연속 오차·자체 GR 및 관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 분자 자료 연결과 고온 연장식의 배제 경계

분류: Imported from prior work. 공개 H2 원표 20000행과 준위 348개를 확보했다. 분류: Counterexample candidate. native Q와 원표는 최대 약 10.62% 차이가 있고 원표 Cv 정의도 별도 확인이 필요하다. 분류: Proven. 선언한 고온 Taylor 식은 1e6 K 주변 전체 구간에서 양의 고정 스펙트럼 분배 함수 조건에 위배된다. 별도 348준위 합산의 0–10차 미분 구간을 얻었다. 분류: Conjectural. 준위 완전성·H2+·플라스마 점유·전체 EOS 오차, 자체 GR 및 관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 H2+ 원자료 확보와 실제 포함 준위 감사

분류: Imported from prior work. 별도 CDS 경로로 H2+ 에너지/전이 원표를 확보했다. 분류: Counterexample candidate. 284168개 전이 행에는 고유 준위 337개만 들어 있어 원 논문의 423개보다 86개 적다. 분류: Proven. 명시한 337개 부분 스펙트럼의 10차까지 온도 미분 구간과 양의 내부 열용량을 얻었다. 분류: Conjectural. 누락 준위·물리 오차·플라스마 점유와 전체 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 H2+ 누락 86준위 원자료 확보와 423준위 합산

분류: Counterexample candidate. MOL-D에서 H2+ 에너지 423개를 실제 확보해 기존 누락 86키를 해소했다. 전체 키는 독립 논문의 각 v별 Nmax와 일치한다. 같은 모형의 1e6 K Q에서 새 86개의 기여는 약 40.61%다. 분류: Proven. 선언한 423준위 모형의 10차 미분 구간과 양의 내부 열용량을 얻었다. 분류: Conjectural. 물리 에너지/점유 오차·전체 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 423준위 분배 함수의 native EOS 연결

분류: Counterexample candidate. 실제 H2+ 423준위 분배 함수와 두 온도 미분을 공통 native EOS 루틴에 연결한 별도 라이브러리를 만들었다. 16개 온도의 직접 재현, 다섯 상태 재고/전하 및 기존 55개 결합 증분 검사를 통과했다. 원래 실패와 GR/시간 상태는 보존한다. 분류: Conjectural. H2 고온 함수·물리 점유/에너지/전체 EOS 오차, 자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 H₂ 302준위 원자료 계산 및 GR 역산 표현 한계

분류: Counterexample candidate. 원본 H2SPECTRE 7.4를 빌드해 H₂ 302준위를 계산했다. 초기 294개 통과와 8개 실패를 보존하고, 실패 8개 및 양성 대조 3개의 범위/간격 대조를 통과했다. 분류: Proven. 선언한 302준위 모형의 10차 미분 구간과 양의 내부 열용량을 얻었다. GR 역산 실패점에는 원래 예산을 만족하는 binary64 점이 없는 국소 입력 구간이 있음을 조건부로 보였다. 분류: Conjectural. native H₂/에너지 기준 연결, 물리 EOS·정밀 역산·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 H₂/H₂⁺ 분배 함수와 화학 에너지 순환의 native 연결

분류: Counterexample candidate. H₂ 302준위와 H₂⁺ 423준위의 공통 native EOS 및 해리/이온화 기준을 연결했다. 16개 온도 직접 대조·다섯 상태 물질/전하·같은 55개 결합 증분이 통과했다. 분류: Proven. 선언된 해리 에너지에서 두 분자 이온화 에너지를 유도하는 화학 순환 항등식을 검증했다. 분류: Conjectural. H₂⁺ 원자료 혼합·물리 점유/전체 EOS 오차, 정밀 역산·새 GR 구조·자체 진화·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 GR EOS 정밀 연산 후보의 엔트로피 역산

분류: Counterexample candidate. 원래 GR EOS의 128비트 연산 후보를 만들었다. 동일한 실패 목표와 네 알려진 근을 기존 예산으로 풀어 통과했다. 온도를 binary64로 되돌려도 새 평가기는 통과하므로 내부 연산 반올림의 영향도 드러났다. 기존 평가기의 실패는 보존한다. 분류: Conjectural. 실패 구역 재연결·평가/물리 EOS/연속 오차·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 새 분자 EOS의 전체 항성 기준과 16단계 시간 대조

분류: Counterexample candidate. 새 두 분자 EOS로 같은 물질·밀도·온도의 5,735개 구역을 재평가해 기준 엔트로피를 확보하고 유한 재고 검사를 통과했다. 새 GR 기하는 아직 맞추지 않았다. 분류: Proven. 16단계 경로의 저장 에너지 대수 검사를 통과했다. 분류: Counterexample candidate. 8/16 시간 대조는 2.04165e-4로 기존 1e-4 기준에 미달했다. 분류: Conjectural. 새 GR·32단계/엄밀 시간 오차·물리 EOS·관측 폐쇄를 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 정밀 EOS를 연결한 실패 구역의 전체 재검증

분류: Counterexample candidate. 실패 블록의 128개 구역과 3,072개 노드 역산을 기존 문턱으로 재검증해 통과했다. 한 실패 노드에만 같은 물리식의 정밀 평가기를 사용했으며 원 실패는 보존한다. 원 실행의 실제 완료 3,840개 구역과 미시작 1,767개 구역을 구분하고 후자를 별도 재개했다. 분류: Conjectural. 전체 집계·연속/물리 EOS 오차·새 GR 구조·자체 진화·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 고정 선형 scalar 해의 전구간 정칙성

분류: Proven. 고정 선형 scalar 모형의 5,736개 다항 구간을 정확한 유리수로 검사했다. Volterra 노름 κ≈1.50093e-4와 양의 외부 정규화 하한 0.99983295를 얻어 중심부터 외부까지 영점 부재, 지정 경계 정규화의 유일성 및 비자명 영 경계 정적 해의 부재를 증명했다. 이전 계수는 조건부 엄밀 구간 안이다. 분류: Conjectural. 물리 EOS/GR 계수 오차·유한 진폭·실제 동역학/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 새 분자 모형 전용 정밀 EOS와 GR 재구성 시작

분류: Counterexample candidate. 같은 두 분자 스펙트럼 소스의 전용 정밀 EOS에서 다섯 상태·세 온도 오프셋의 15개 알려진 근을 기존 문턱으로 회복했다. 새 GR용 17점 표를 재계산 중이며, 최초 NumPy 논리값 저장 오류를 보존한 별도 실행을 사용한다. 분류: Conjectural. 전체 새 GR·연속/물리 EOS·자체 진화·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 scalar 수치 계수의 전구간 잔차 인증과 독립 수정

분류: Proven. 고정 모형의 중심·전 구간·무한 외부를 포함한 잔차 인증으로 alpha_A/phi∞를 [-4.000356564041216,-4.000356563994881]에 묶었다. 분류: Counterexample candidate. 이전 전역 IVP 값은 이 구간 밖이며, 모든 다항 구간을 통과한 두 독립 적분은 -4.000356564018176으로 구간 안이다. 이전 유한 비교 판정은 보존한다. 분류: Conjectural. 물리 EOS/GR 계수 오차·유한 진폭·실제 구동/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 유한 scalar 진폭의 비선형 진공 계량 연결

분류: Counterexample candidate. 비영 scalar 에너지를 포함하는 비선형 진공 외부 매칭을 구현하고, 무한대를 포함하는 좌표의 14개 독립 적분을 해석적 Just 해와 대조해 통과했다. 분류: Proven. 분리한 질량/lapse 증분의 좌표식 및 Einstein/Jordan 점근 계수의 차이를 기호 검산했다. 분류: Conjectural. 실제 비선형 항성 내부·대기·자체 진화·구동/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 독립 복사 진화를 위한 물질 EOS 분리

분류: Proven. 광자 자유에너지와 모든 필요한 열역학 미분을 일관되게 제거하는 식을 기호 검산했다. 분류: Counterexample candidate. 동일한 새 분자 EOS의 native 복사 제외 옵션과 15개 대조를 통과했고, 저장된 5735개 상태의 물질 압력·열용량·등온 압축률 양성을 확인했다. 분류: Conjectural. 실제 흡수·산란·복사 수송, 대기와 자체 GR 진화는 계속 연결해야 한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 복사·물질·시간 의존 GR 계량의 보존식 연결

분류: Proven. 공변 텐서 계산으로 시간 의존 구면 GR의 에너지·운동량·바리온 식과 질량 제약 전파를 검산했다. 움직이는 복사의 4-force, 물질 표면 질량 수지 및 LTE 물질/광자 재결합을 연결하고 Tolman/null 복사 해석 대조를 통과했다. 분류: Conjectural. 실제 opacity·각분포·대기·중성미자/scalar 및 자체 항성 진화는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 혼합 정밀 EOS의 전 구역 native 검증 완료

분류: Counterexample candidate. 원래 성공 3840구역, 실패 블록 재계산 128구역과 미실행 1767구역을 합쳐 5735구역의 정확 coverage와 137640개 native 내부 노드 대조를 통과했다. 원래 실패는 보존했다. 분류: Conjectural. 연속 오차·물리 EOS 인증, 새 분자 GR와 실제 시간 진화·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 Planck 흡수 자료와 항성 상태의 연결

분류: Counterexample candidate. 실제 OP Planck 흡수 단면적 29563개 상태를 독립 Fortran/Python reader로 대조하고 새 분자 EOS의 3043개 항성 구역에 연결했다. 고밀도 2692구역(바리온 질량 94.3696%)과 Li/Be/B/F 자료는 미지원임을 측정했다. 분류: Conjectural. 전체 흡수/산란·각분포·대기 및 자체 GR 진화는 계속 연결해야 한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 고밀도 항성의 흡수·산란 스펙트럼 연결

분류: Counterexample candidate. 다섯 실제 구역의 12개 TOPS 흡수/산란 스펙트럼과 입력·178800개 합산 행을 검증했다. cutoff-OFF 평균은 직접 적분과 일치했다. cutoff-ON은 원시 배열을 바꾸지 않고 평균을 변경하므로 규약 확인이 필요하다. 외곽 서버 실패는 보존했다. 분류: Conjectural. 전체 구역·비-LTE 대기·산란 kernel·물리 EOS와 자체 진화/관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 외곽 평균 흡수와 고밀도 plasma 차단 규약 검증

분류: Counterexample candidate. 같은 외곽 입력의 gray 자료를 확보하고, 중심 아홉 그룹의 cutoff-OFF 평균을 독립 적분으로 재현했다. ON의 여섯 저에너지 그룹은 1e10 차단값이며 높은 세 그룹은 동일했다. 비상대론적 plasma 문턱 위 조건부 적분은 ON 평균에 근접하지만 물리 인증은 아니다. 분류: Conjectural. 외곽 스펙트럼/비-LTE 대기·매질 광자 EOS·자체 진화/관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 외곽 단색 적분 실패와 원천 그룹 적분의 분리

분류: Counterexample candidate. 원래 외곽 두 온도를 나눠 단색 자료를 확보했지만 Planck 구적은 미달했다. 원 실패를 보존하고 실제 원천 16/32그룹을 재합산해 평균과 최대 5.53e-5 이내 일관성을 확인했다. 정확한 차단 표지와 큰 opacity를 구별하도록 별도 reader를 수정했다. 분류: Conjectural. 선 모양·비-LTE 대기·전체 수송과 자체 GR/관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 매질 광자의 자유에너지·굴절·전 구역 진단 연결

분류: Proven. 선언한 횡파 준입자 모형의 자유에너지·상태 의존 gap 열역학·구면 GR 굴절 Hamiltonian을 검산했다. 분류: Counterexample candidate. 5735개 실제 상태의 광자 진단과 일곱 독립 적분을 통과했으며 중심에서 전자 축퇴/상대론 보정 필요성을 확인했다. 분류: Conjectural. 물리 plasma EOS·분극/수송·대기·자체 GR 진화와 관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 축퇴 전자 상태와 상대론 plasma moment 연결

분류: Counterexample candidate. 5735개 native 전자 상태의 Fermi 밀도 재구성과 독립 아홉 moment 대조를 통과했다. 중심 plasma frequency는 기존 비상대론식보다 약 3.91% 작다. 분류: Proven. 영온 moment·특성 속도 횡파 근의 유일성과 급수 꼬리를 검산했다. 분류: Conjectural. 일반 분산·공통 물리 EOS·자체 GR 진화/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 일반 선도 차수 횡파·종파 분산과 양의 급수 경계

분류: Proven. 원 분극 적분의 양의 moment 급수·절단 꼬리·횡파 및 시간꼴 종파 근의 유일성을 얻었다. 분류: Counterexample candidate. 5735개 상태의 분산 및 원 로그 적분 54개 독립 대조를 통과했다. 분류: Conjectural. 구간 수치 보증·공통 물리 EOS·자체 수송/GR 진화와 관측을 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 plasma moment 값·혼합 미분과 분산 kernel 연속 오차 보증

분류: Proven. 실제 전자 상태를 포함하는 상자에서 새 정확 Gauss 구적법의 moment 값·2차 혼합 미분 157성분을 인증했다. 분산 kernel 상계는 횡파 7.484e-16, 종파 3.867e-14 이하며 독립 유리수 전달을 통과했다. 분류: Conjectural. 기존 native 전구간·전체 물리 EOS·자체 진화/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 Coulomb 자유에너지 중복 및 유한 종파 경계의 직접 매칭

분류: Proven. 현재 EOS의 DH 항과 고전 ring 자유에너지가 같음을 검산하고 중복 계상·유한 종파 끝점 미분 경계를 밝혔다. 분류: Counterexample candidate. 실제 별의 바리온 질량 중 DH 영역은 약 6.09%이며 나머지는 보간/수정 OCP 영역이다. 분류: Conjectural. 공통 강결합 자유에너지·자체 수송/GR 진화·관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 전 구역 분산근 63085개 구간·반올림 인증

분류: Proven. 저장된 5735개 전자 상태에서 63085개 선도 차수 분산근을 지향 구간으로 검산했다. 독립 유리수 예산까지 포함한 h 오차 상계는 최대 3.931e-14다. 원 소수 예산 불일치는 보존하고 별도 반올림 예산으로 감쌌다. 분류: Conjectural. 실제 매개변수/상태 미분·공통 물리 EOS·자체 GR 진화/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 재고의 강결합 자유에너지 매칭과 공개 EOS 미분 수정

분류: Proven. 공개 양자 자유에너지의 압축률을 다시 유도하고 전자–이온 차폐식의 누락된 몫 미분 항을 수정했다. 분류: Counterexample candidate. 거의 완전 이온화된 3206구역(바리온 질량 97.5604%)의 실제 비이상 자유에너지를 매칭하고 같은 기준의 미분 재검산을 통과했다. 원 실패를 보존하며 PIMC 보정 표 아홉 값의 최대 차이는 약 0.088%였다. 분류: Conjectural. 물리/연속 EOS 오차·부분 이온화 평형 교체·자체 GR 진화/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 양자 혼합 선도계수와 열역학 오차의 엄밀 경계

분류: Proven. 고정 균일 전자 배경의 양의 점이온 혼합 WK 선도계수와 중성 조성 미분을 유도했다. 선언 양자 적합식의 전 R>0 절단 상계 및 실제 3206구역 22442개 열역학 평가의 반올림을 감싼다. 최대 평가 오차 2.578e-19로 기준을 통과했다. 분류: Conjectural. 실제 다체계의 고차 물리 오차·차폐/부분 이온화·전체 EOS 미분·자체 GR 및 관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 점이온 Hamiltonian의 비섭동 양자 자유에너지 경계

분류: Proven. 고정 균일 전자 배경의 양의 점이온 Boltzmann Hamiltonian에서 bulk 양자 자유에너지 보정이 0과 WK 선도항 사이임을 경로 적분과 양의 열 kernel로 증명했다. 유한 상자의 감김 항을 남겼으며 급수 수렴은 가정하지 않는다. 실제 3206구역의 상계는 이온당 kBT의 1.96915e-4 이하이다. 자유에너지 값 부등식의 미분은 허용되지 않는다. 분류: Conjectural. 전자 응답·고전 상관·부분 이온화 및 물리 미분·자체 GR/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전자 차폐 Hamiltonian의 자기항과 상태 미분 연결

분류: Proven. 원 전자 차폐 Hamiltonian을 쌍 항과 밀도/온도 의존 자기항으로 정확히 분해했다. 같은 장파장 κ와 양의 정적 유전함수도 다른 자기 에너지와 음의 실공간 꼬리를 허용하므로 기존 점이온 경계를 자동 확장할 수 없다. 상태 의존 Hamiltonian의 내부에너지·압력·열용량 미분을 명시했다. 분류: Counterexample candidate. 실제 3206구역의 이상 전자 κ₀ 및 독립 적분 대조를 통과했다. 분류: Conjectural. 유한 파수 전자 응답·상관/부분 이온화·물리 EOS·자체 GR/관측은 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 상태의 유한 온도·유한 파수 전자 응답

분류: Proven. 영온 Jancovici 매질 응답의 유한 온도 변환과 양의 정적 커널을 도출했다. 같은 저장 화학 퍼텐셜에서 생략한 열적 양전자의 상대 응답을 모든 파수에 대해 1.32e−158 미만으로 상계했다. 분류: Counterexample candidate. 3206구역·41678개 응답과 독립 고정밀 대조가 사전 기준을 통과했다. 두 구현 실패는 보존했다. 분류: Conjectural. 전 파수 자기 에너지/상태 미분·상관/전체 물리 EOS·자체 GR/관측 폐쇄를 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전자 편극 자기항의 전 파수 적분과 엄밀 꼬리 경계

분류: Proven. 이상 전자 매질 응답의 전 파수 모멘트와 자기 에너지의 전역·꼬리 상계를 도출했다. 3206상태의 생략 꼬리는 정규화 적분 기준 2.80e−11 미만으로 인증했다. 분류: Counterexample candidate. 실제 상태의 자기항 내부 적분과 두 외부 구적 대조를 통과했다. 이미 포함된 EOS에 중복 가산하지 않았다. 분류: Conjectural. 내부 구적/상태 미분·전체 물리 EOS·자체 GR/관측 폐쇄를 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 물질 재고의 전자 중성 근과 역미분 구간 인증

분류: Proven. 완전 이온화 선언 모형에서 실제 3206구역 물질 재고의 이상 전자 중성 근과 밀도/온도 1·2차 역미분을 MPFR 및 기존 연속 Fermi 경계로 인증했다. 출력 끝점 예산의 원 불일치는 보존하고 구간을 바꾸지 않은 정확한 감사로 보완했다. 분류: Conjectural. 자기항 열역학 연결·전체 상관/부분 이온화 EOS·자체 GR/관측 폐쇄를 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 참 중성 근 이동과 편극 열역학 계수 오차의 전 파수 상계

분류: Proven. 편극 자기항의 참 중성 근 이동과 역미분 계수 오차를 전 파수에서 상계했다. 분류: Counterexample candidate. 넓은 근 박스의 CV/PDT 상계는 사전 예산에 미달하며 원 결과를 보존했다. 인증된 근을 포함하는 64배 좁은 박스로 같은 미분을 재평가한다. 분류: Conjectural. 내부 구적/전체 EOS·자체 GR/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 동일 편극 자기항의 화학 퍼텐셜과 조성 미분

분류: Proven. 같은 자기 자유에너지의 화학 퍼텐셜·혼합 온도미분·대칭 rank≤2 조성 Hessian과 전하 보존 방향의 정확한 선형성을 도출했다. 분류: Counterexample candidate. 실제3206×26종 연결식과 명목 반응 계수 대조가 통과했다. 분류: Conjectural. 전체 상관/부분 이온화 EOS에의 일관된 결합·내부 오차·자체 GR/관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 좁은 중성 근 구간으로 열역학 역미분 오차 예산 충족

분류: Proven. 같은 물질 입력과 구적 규칙에서 η 박스만 64배 좁혀 전3206구역을 재인증했다. 자기항 일곱 열역학 출력의 참 근 이동/역미분 계수 오차 상계가 모두 기존 2e−7 예산 안이다. 원 넓은 박스 실패는 보존한다. 분류: Conjectural. 내부 구적/전체 EOS·자체 GR/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전자 응답 로그 특이점의 정확한 구적 규칙과 오차 계수

분류: Proven. 전자 응답의 log|1−u| 항에 사용할 보통/로그 가중 4점 Gauss 규칙을 정확한 유리수 대수와 양의 가중치·근 격리로 인증했다. 8차 미분 나머지 계수와 독립 지수 함수 대조를 확인했다. 분류: Conjectural. 실제 상태별 미분 경계와 내부/외부 적분 전체 오차·전체 EOS·GR/관측의 결합은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 32단계 열진화의 저장 보존식 감사와 시간 세분화 미달

분류: Proven. 고정 밀도/계량/조성 32단계 caloric 경로의 저장 수치를 정확한 유리수로 감사했다. 에너지·엔트로피·endpoint 재생 조건을 통과했다. 분류: Counterexample candidate. 16/32 온도 차이는 원 1e−4 기준에 미달하여 같은 알고리즘의 64단계 경로를 시작했다. 분류: Conjectural. 연속 시간/전체 EOS·자체 GR/관측 보증은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 편극 자기항의 중성 열역학 미분과 엄밀 꼬리

분류: Proven. 전자 중성 조건을 따라 편극 자기 자유에너지의 압력·내부에너지·엔트로피·열용량과 압력 미분을 도출하고 무한 파수 꼬리를 구간 상계했다. 분류: Counterexample candidate. 실제3206상태의 일곱 출력과 중성 근을 다시 푸는 유한차분 대조가 사전 기준을 통과했다. 분류: Conjectural. 내부 구적/연산 오차·전체 자유에너지 결합·자체 GR/관측을 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 중성 상태의 로그 응답 창과 고차 미분 구간 인증

분류: Proven. 실제3206상태×12파수에서 로그 특이점 창의 여섯 응답 편미분을 복소 Cauchy 경계와 정확한 가중 구적으로 감쌌다. 38472개 창 모두 출력 끝점을 포함한 1e−11 예산을 통과했다. 분류: Counterexample candidate. 독립70자리 대조72개가 구간에 포함됐다. 분류: Conjectural. 나머지 운동량 영역/외부 파수·전체 EOS·자체 GR/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 로그 창 밖 운동량 꼬리의 여섯 편미분 인증

분류: Proven. 실제3206상태×12파수에서 로그 창 밖 무한 운동량 꼬리와 여섯 편미분의 상계를 인증했다. 최대6.941e−20으로 사전1e−11 예산을 통과했고 올림 binary64 절단점도 정확 검산했다. 분류: Counterexample candidate. 독립70자리 대조72값이 상계 안에 있었다. 분류: Conjectural. 유한 운동량·외부 연속 파수·완전 EOS·자체 GR/관측 접합은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 복소 파수 응답과 유전함수 분모의 비영 영역

분류: Proven. η≤24, β>0의 모든 양의 실수 파수 주변에서 복소 응답 해석 영역과 여섯 상계를 유도하고 실제3206상태·38472개 분모 비영 원판을 인증했다. 원판은 보수적이며 전체 외부 적분 완료가 아니다. 16점 Gauss 규칙도 정확 인증했다. 분류: Counterexample candidate. 54개 독립 복소 대조 통과, 유한 운동량 전체 계산 진행. 분류: Conjectural. 전구간 접합·완전 EOS·자체 GR·관측 폐쇄는 계속한다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 연속 상태 상자에서 양의 복소 응답과 공통 비영 영역

분류: Proven. β∈[3e−6,0.006], η∈[−17,24] 연속 상자와 모든 q0>0에서 |Q−q0|≤q0/256의 응답 실수부가 양수임을 증명했다. 양의 커널 위상과 같은 파수 감소율의 내부/꼬리 경계를 결합했고 실제3206근 구간의 포함을 확인했다. 모든 B>0의 유전함수 분모가 비영이다. 분류: Conjectural. 실제 외부 적분·완전 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 자기항 외부 적분의 영파수 끝점 인증

분류: Proven. Q=0 끝점의 알려진 상수 기여와 세제곱 나머지를 분리해 실제3206중성 상태의 일곱 자기항 열역학량을 인증했다. 최대 오차2.585e−13 미만으로 사전2e−7 기준을 통과했다. 참 근 Hessian·조성·상수·출력 반올림을 포함한다. 분류: Conjectural. 두 파수 끝점 사이 실제 적분·완전 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 중성 상태의 전체 운동량 응답값과 다섯 편미분 인증

분류: Proven. 실제3206상태×12파수에서 유한 운동량·로그 창·무한 꼬리를 합쳐 p∈[0,∞)의 여섯 응답 성분(값과 다섯 편미분)을 인증했다. 유한 구간1e−11·전체3e−11의 출력 오차 예산을 모두 통과했다. 분류: Counterexample candidate. 독립70자리 대조72값 포함, 사전12개 재현 일치. 분류: Conjectural. 실제 외부 파수 적분·완전 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 자기항의 공동 복소 해석성과 직접 미분 상계

분류: Proven. 기존 상태 상자에서 Q 상대1/256·η 반지름0.01·lnT 반지름0.001의 공동 해석/양의 응답 영역을 인증했다. |G|≤1.000031과 파수 공통의 값/다섯 편미분 Cauchy 상계를 직접 얻어 큰 중간 응답 상계로 인한 손실을 줄였다. 여섯 포화 단항식 검산 통과. 분류: Conjectural. 실제 외부 적분·물리량 오차 합산·완전 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 고파수 응답에 재사용할 열 모멘트와 다섯 편미분 인증

분류: Proven. 실제3206중성 상태에서 17개 열 모멘트의 값과 다섯 편미분, 총327012성분을 인증했다. 최대 출력 오차1.335e−15 미만으로 사전1e−12 기준 통과. 분류: Counterexample candidate. 독립70자리 대조54값 전부 포함. 분류: Conjectural. 고파수 비선형 자기항·중간 파수 적분·완전 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 고파수부터 무한대까지 비선형 자기항과 일곱 열역학량 인증

분류: Proven. 실제3206상태에서 Q≥4P의 비선형 자기항 값과 다섯 편미분을 무한대까지 적분하고 일곱 중성 물리량으로 전파했다. 최대 출력 오차2.138e−18 미만으로 사전2e−7 기준 통과. 분류: Counterexample candidate. 세 독립 중첩 적분 대조 통과. 분류: Conjectural. 중간 파수 접합·완전 EOS·자체 GR·실제 관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 로그 특이점을 제거한 양의 응답과 일정 폭의 해석 영역

분류: Proven. 원래 로그 커널을 특이점이 제거된 양의 적분식으로 정확히 바꾸고 |ImQ|≤sqrt(β)/16의 해석성·양의 실수부를 인증했다. 분류: Counterexample candidate. 독립 원시함수/적분 구조로 계산한72성분 모두 기존 구간에 포함. 분류: Conjectural. 공동 매개변수·분모 영역과 실제 중간 파수 적분·완전 EOS·GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 영파수 부근까지 공동 비영 영역과 실제 외부 구적 계획

분류: Proven. 넓어진 Q/η/lnT 공동 영역에서 실제3206상태의 양의 응답 하계와 screening 크기의 분모 비영 영역을 인증했다. 공통|G|≤2와 여섯 미분 상계로 세 대표 상태의 정확한 외부 패널·구적 나머지2e−8을 확보했다. 분류: Conjectural. 실제 절점 계산·끝점 접합·완전 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 재사용 응답 보간·로그 곱 적분·연분수와 실제 절점 검증

분류: Proven. 32점 H 보간과 로그 곱 적분·엄밀한 연분수, 실제3206상태의 전체 H 꼬리 상계1.645e−20 미만을 확보했다. 분류: Counterexample candidate. 독립54모멘트·실제48절점/288성분 대조 통과. 초기 오차 정체 실패는 보존했다. 분류: Conjectural. 실제 H 계수·중간 적분·완전 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 재사용 H 구간 계수·독립 합성곱과 기준 적분 실패 진단

분류: Proven. 세 기준 중심의 실제 H 계수와36개 전체 운동량 응답을 인증했다(중간값 오차5.537e−17 미만). 분류: Counterexample candidate. 독립 합성곱216성분 최대 차이5.537e−17 미만. 기존 scalar 기준 적분의4개 정밀도 실패와 무시된 수렴 상태는 보존한다. 분류: Conjectural. 참 근 물리량의 중간 파수 적분·전체 EOS·자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 비영 열유속의 전 파수 감쇠와 실제 파수 절점 접합

분류: Proven. 기존5735구역·두 상수 열완화시간의 고정 계수 rank-one 모형11470개가 모든 유한 비영 파수에서 감쇠함을 정확 인증했다. 재사용 H의 실제80비트 파수48절점과 참 중성 물리량의 오차 접합 정의도 검증했다. 분류: Conjectural. 전체 중간 적분·변수 계수/실제 수송·완전 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).
