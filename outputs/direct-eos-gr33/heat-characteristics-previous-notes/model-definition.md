# Model definition — unified revision

The current scientific statement is [the unified manuscript](../paper/manuscript.md), Sections 3–4. Earlier runtime proposals remain in the dated REQUEST10 notes; they are not instructions to reopen runtime work.

Status: Counterexample candidate. The dimensionless prescribed-drive benchmark is

```math
\tau_\chi\dot\chi+\chi=\alpha F(t),\qquad
q(t)=c_YF(t)+c_\chi\chi(t),\qquad \beta=\alpha c_\chi,\quad\tau_\chi>0.
```

Status: Proven. With dimensionless F and chi, alpha, c_Y, c_chi and beta are dimensionless; tau_chi has units of time. For real F0 and `F=F0 cos(omega t)`, the settled response is

```math
\chi_{\rm ss}=\frac{\alpha F_0(\cos\omega t+\omega\tau_\chi\sin\omega t)}{1+\omega^2\tau_\chi^2},\qquad
G(i\omega)=c_Y+\frac{\beta}{1+i\omega\tau_\chi}.
```

Status: Proven. The full solution includes `chi_h exp(-(t-t0)/tau_chi)`. All periodic carrier results assume that term is absent or decayed. The stored templates do not estimate its independent amplitude.

Status: Counterexample candidate. The original mass readout `m_A/m0=1+q` and the timing pair potential `V_pj=-G m_p m_j(1+Delta0+q)/r_pj`, with fixed inertial masses, are separate realizations. The latter specifies the pairwise force modification for a prescribed external drive.

Status: Proven. The prescribed pair potential gives equal-and-opposite pair forces and allows energy exchange through its explicit time dependence. A varied position-dependent drive or inertial mass requires additional gradient/backreaction terms. The ODE alone does not derive these terms or a matching between the two realizations.

Status: Counterexample candidate. [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md) supplies a conditional force-level realization with an independent scalar charge Q_p, potential V(Q_p), inertia I and Rayleigh damping Gamma. Reciprocal pair forces follow from -(m_A m_B+Q_A Q_B)/r_AB; total orbital-plus-state energy loss is -Gamma Qdot_p^2. This is a leading Newtonian reduction, not a complete timing theory.

Status: Proven. On a stable branch kappa=V''(Q0)>0, equal fixed white-dwarf charge/mass ratios a_w give tau=Gamma/kappa and deltaDelta=B deltaU/(1+tau d/dt), B=a_w^2/(kappa m_p). For F=deltaU/Ustar, beta=B Ustar. A small inertial-error bound is needed for the one-pole approximation. Unequal companion ratios obstruct a common pair modulation.

Status: Conjectural. Actual EOS-to-body matching must still determine the coefficients and validate inertia, nonlinear response, feedback, radiation and transients for the real system. The conditional action does not supply those numerical inputs.

Status: Proven. A real local derivative comparator is `P_N(d/dt)F`, with coefficients shared across carriers. A finite-dimensional nuisance model still needs a rank test; finite dimension alone never guarantees pole observability.

Status: Imported from prior work. Explicit dynamical internal modes already appear in compact-body EFT; see Chakrabarti et al. (2013), Steinhoff et al. (2016), and Khalil et al. (2022), cited in the manuscript. The present claim is the specified comparator boundary, not priority for internal states.

Status: Counterexample candidate. Request 11.4 specifies a reciprocal overdamped comparator Gamma qdot+Kq=bF, readout b^Tq+c0F, with symmetric positive-definite matrices and all rates >=Lambda. It yields positive relaxation weights with tau_j<=1/Lambda. This physical comparator premise is not established EOS matching; nonconjugate readout, negative residues or slow/oscillatory modes need a different class.

Status: Proven. Request 11.5 shows that the six periodic response columns plus a static response cannot determine an initial-state exponential response: D*product(D^2+omega_k^2) annihilates the known inputs but not the exponential. A dedicated validated response is needed; an exponential coupling is not automatically an exponential TOA residual.


## Request 12 follow-through

Status: Proven. A positive reciprocal relaxation measure admits a two-frequency moment-variance equality identifying one observable relaxation time, with exact calibrated response; it does not count hidden physical states. Status: Imported from prior work. Request 12 constructs the corrected unequal-amplitude leading drive and dedicated exponential timing responses at 2, 52 and 500 days. Status: Conjectural. EOS coefficients and the fast-rate gap are not matched to J0337.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Counterexample candidate. Request 13 now fixes and numerically solves a concrete SLy/massless-DEF model (beta=-4, zero background scalar, gravitational mass 1.4378144085 solar masses). Its regular stellar scalar response and outgoing pole are calculated; it is not a choice of freely adjustable damping. The response is linear about the zero-scalar branch, where fluid and metric perturbations decouple. The published beta=-5 prototype is a separate positive control. See [stellar matching](../notes/REQUEST13_STELLAR_DERIVATION.md).

## Request 14 validated flow

Status: Proven. A separate interval implementation encloses the frozen internal GR 4-body 1PN IVP and its initial-state variational flow at recorded epochs, with four independent fractional-mass tangent columns in a local augmented run. Its input rejects nonzero dynamic-SEP and non-GR coefficients. This is a conditional numerical certificate for the GR comparator, not a completed matched scalar-tensor signal model. The physical timing-parameter-to-IVP map remains outside the certificate. See [certificate boundary](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).


## Request 15 후속 검증

분류: Proven. 활성 초기화의 질량 중심 위치 항과 회전 공변성을 깨는 속도 결합을 수정했다. 기존 물리 초기 상태를 재현하는 궤도·스핀 좌표를 별도로 계산했으며, 이는 새로운 힘이나 동적 관측량이 아니다. 질량 함수와 Kepler 근의 구간 미분은 인증했으나 초기 상태 전체로의 미분 연쇄는 미완료다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).


## Request 16 다체 관측식과 영 구동 경계

분류: Proven. 추가 천체의 GR Einstein 적분 항과 Shapiro 항을 실제·모의 관측의 공통 경로에 연결했고 천구 좌표 변경 시 보존된 전체 상태의 회전도 수정했다. 분류: Proven. 지정된 DEF 영 scalar 가지에서 모든 천체가 비스칼라화 상태이고 scalar 초기·입사 자료가 영이며 고전 초기값 문제가 유일하면, 궤도 운동에도 scalar=0인 GR 해가 유지된다. 외부 scalar에 대한 감수율이나 pole만으로 비영 구동을 얻지 못한다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).


## Request 17 비영 구동과 동반성 응답

분류: Counterexample candidate. DEF beta=−4, 배경 scalar=1e−5, SLy 중성자별과 차가운 mu_e=2 전자 축퇴 EOS의 백색왜성 두 개를 지정했다. 영 배경의 기준 별로부터 바리온 수를 고정한 비영 가지를 수치 계산했다. 선도 작은 배경 모형은 Lq=phi*1, L=diag(1/chi)−K, K_AB=1/r_AB이며 모든 전하가 서로 응답한다. 실제 J0337의 EOS·배경 값을 결정한 것은 아니다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).


## Request 18 열 백색왜성 구조와 응답 경계

분류: Counterexample candidate. 공개 MESA 열·외피 진화 자료에서 안쪽 WD의 광학 후보 이력 행과 그 전후의 내부 구조 두 개를 확보했다. 회전속도·각속도와 셀 질량·밀도로 자료 단위를 복원했다. 두 구조의 환산 Newtonian 질량은 약 0.198040 태양질량으로 타이밍 기준보다 0.2548% 높다. 영 배경·고정 밀도·평탄 시공간의 선형 scalar 모형을 별도로 정의했으며 완전한 열 GR 별로 취급하지 않는다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).


## Request 19 열 진화 재현과 질량 보정 대조

분류: Proven. 공개 MESA 7624 사진에 필요한 ca40 핵종과 실행 중 유지되던 확산·대기 경계 설정을 명시적으로 복원했다. 원래 구조의 3,095개 셀·29개 구조 및 조성 열, 다음 한 단계, 18000→19000의 1,000단계 주요 이력과 최종 내부 구조가 일치했다. 원래 배포 상수도 직접 확인했다.

분류: Counterexample candidate. 실제 외피 제거와 후속 열 진화로 Newtonian GM 목표에 맞춘 구조를 계산했다. 사전 질량·광학 탐색 기준을 함께 만족한 행은 11개다. 기록한 최적 모델 19057은 Teff=15828.041340 K, logg=5.74162539, GM/기준 GM_sun=0.19753638530700이다. 이는 GR 중력질량 및 형성 이력까지 검증한 matching이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).


## Request 20 열 구조 민감도와 질량 정의

분류: Counterexample candidate. 같은 초기 사진의 제거율 두 가지 및 공간·시간 분해능 변형 네 가지를 고정 나이 구간에서 검사했다. 실제 질량·광학 조건 통과 행은 fast 0개, slow 11개, mesh 14개, time 0개이며 유한 변형 판정은 미통과다. 독립 형성 모형이나 GR 질량 일치로 승격하지 않는다.

분류: Counterexample candidate. 원래 구간에서 실패한 두 실행을 대상으로 별도 등록한 다음 냉각 구간 검사에서는 fast 13개, time 25개의 실제 후보를 확보했다. 원래 실패 판정은 유지한다. 나이 이동과 양립하는 후보 회복이며, 고정 나이의 수렴 인증은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).


## Request 21 열 EOS와 GR 질량 매칭

분류: Counterexample candidate. FreeEOS 3.0.0 EOS1을 실제로 재평가하는 구대칭 TOV 모형을 구성했다. 중성 원자 에너지 기준, 고유 바리온 부피, 영압 외곽을 명시했다. 고정 T(P)·조성 조건에서는 GR 질량 일치만 통과하고 광학 조건은 실패했다. 별도 온도 배율 1.019968012525의 후보는 목표 질량과 조건부 Teff를 맞춘다. 불소를 제외해 재규격화한 조성, 정해 둔 T(P), 수학적 외곽 대기는 모형 가정이다.

세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).


## Request 22 GR 후보의 열수송과 물질 보존

분류: Counterexample candidate. Request21의 조정 후보는 원래 재규격화 프로필보다 총 바리온 질량이 0.069465%, 수소 재고가 5.42011% 작다. 같은 X(P)는 같은 물질을 보존하는 변환이 아니다. 별도 모형족의 조건부 GR 질량 해는 유지하지만, 원래 항성의 보존적 GR 변환이라는 해석은 통과하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).


## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성

분류: Counterexample candidate. Request23은 원래 5735개 구역의 바리온 질량, 불소 제외 재규격화 조성과 FreeEOS 기준 엔트로피를 물질 좌표에 고정하고 TOV를 다시 풀었다. 원래 내부 물질을 보존한 해는 지정 중력질량 목표보다 0.0698907% 높다. 별도 모형족에서 모든 구역 질량을 0.999301574445배로 조정하면 목표 GR 질량을 맞춘다. 이는 원래 총 물질을 보존한 변환과 다른 모형족이며 온도·엔트로피 조정 인자는 없다.

세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).


## Request 24 조성 변화의 반응 에너지 기준

분류: Proven. Request24는 조성 의존 에너지 기준을 연결했다. e_B=C_X c²+u_W를 유지하면서 u_Q=u_W+g(X), e0_Q=e0_W−g(X)로 바꾸면 총 에너지는 그대로다. g는 배포 원자량과 MESA 표준 Q를 만드는 질량초과의 차이이다. 따라서 핵 가열항도 dg/dτ만큼 함께 바뀌어야 한다. 22종 네트워크의 72항목 중 닫힌 반응식 64개에서 이 변환을 검산했으며 나머지 8개 보조율은 독립 반응식 인증 대상에서 제외했다. 이 기준 변환은 Request23의 GR 질량을 다시 적합하거나 바꾸지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).


## Request 25 새 GR 상태의 미세물리 평가

분류: Counterexample candidate. Request25는 기존 MESA 실행 파일의 모형 입력·초기 출력 경로로 새 GR 물질 상태의 핵 가열·반응 및 열적 중성미자·불투명도를 재평가했다. 출력의 밀도·온도·조성·구역 질량이 입력과 일치하는지 먼저 확인했다. FreeEOS가 정한 TOV 상태를 MESA 고유 EOS 보조량으로 평가한 혼합 모형이며 같은 EOS의 GR 진화 해는 아니다.

세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).


## Request 26 남은 폐쇄 조건의 실행 검산

분류: Counterexample candidate. Request26은 실제 네트워크를 63개 반응 그룹으로 조정하여 가열·중성미자·미분을 재구성하고, 소스 Q 정의로 22종 조성 변화 방향과 에너지 장부를 복원했다. 이는 native 전체 변화율의 독립 출력 검증과 구분한다. 초기 불소가 없어도 양의 생성 방향이 있어 19종 불소 생략 모형의 정확한 시간 불변성을 주장할 수 없다.

세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).


## Request 27 직접 반응 벡터와 명시적 PP 상태

분류: Counterexample candidate. Request27은 변경하지 않은 native 반응 함수에서 22종 변화율과 Jacobian을 직접 추출해, 기존 채널 복원 벡터를 독립 검증했다. 생성된 불소를 포함하여 같은 native EOS·반응망에 조성을 되먹이는 고정 밀도·온도 대조를 수행했다. 기존 축약 PP 망은 H=0 경계에서 바깥 방향 변화율을 보였고, 중간 핵종 4개를 명시한 26종 망은 같은 경계 대조를 통과했다.

세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).


## Request 28 반응 에너지 보존과 scalar 미분 구간

분류: Counterexample candidate. Request28은 26종 전체 조성·온도·누적 중성미자 에너지를 결합한 고정 부피 대조를 수행했다. 기본 EOS의 에너지 역산 실패와 온도를 분리한 HELM의 시간 정밀화 실패를 보존했다. 별도 고온 HELM 결합 대조는 같은 엄격 에너지·조성·온도 기준과 추가 중성미자 적분 기준을 통과했다. 물리적 공통 EOS나 GR 항성 진화의 인증은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).


## Request 29 공통 EOS와 GR 열 경로의 관측 연결

분류: Counterexample candidate. 불소를 복원한 26종 핵 재고를 보존하고, 이온 수·전하 수 보존 사상과 조성별 원자 결합에너지 기준 이동으로 FreeEOS 입력을 연결했다. 모든 초기 5,735구역의 유한 미분 대조는 통과했다. 미지원 원소의 부분 이온화·혼합 엔트로피와 물리 EOS 오차 인증은 남는다.

분류: Proven. 고정 조성의 합성 자유에너지 F_B=C F_FE(Cρ_B,T,ε)+ΔI는 원래 자유에너지의 열역학 항등식을 보존한다. 조성 변화에는 화학 에너지 항이 필요하다.

세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).


## Request 30 EOS 역산과 약반응 미분 수정

분류: Counterexample candidate. 같은 광도 입력과 고정 초기값 FreeEOS로 57,350개 저장 원천 끝점 방정식을 기존 허용량 안에서 재역산했다. 미지원 원소의 물리 EOS 인증과 전체 GR 재진화는 별도다. 실제 약반응 표의 선형 보간·열·중성미자 기여를 독립 재구성한 배정밀도 원천 함수를 정의했다.

근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).


## Request 31 조성 미분과 보존형 GR 재적분

분류: Counterexample candidate. 같은 EOS·반응 함수를 조성 25방향과 온도·밀도에서 재평가했다. 직접 미분한 초기 연산자는 명시적인 고정 적분 연산자로 사용하며 각 시간 단계의 비선형 원천과 EOS는 새로 계산한다. 새 1·2·4단계 고정 광도 GR 경로의 에너지 기준 일괄 판정은 실패다. 가장 미세한 4단계의 에너지 기준은 통과, 마지막 조성·온도 시간 대조는 실패다. 거친 경로의 실패는 미세 경로의 판정과 구분하여 보존한다. 별도 구조 피드백 잔차 보정 1단계의 에너지 점수는 9.7476097e-07로 통과했으며, 그 보정 방법의 시간 수렴은 아직 검증하지 않았다.

근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).


## Request 32 직접 원소 EOS 중간 검증

분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.

분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.

근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 32 구조 피드백 엔탈피 시간 대조

분류: Counterexample candidate. 고정 기준 압력의 엔탈피를 열 좌표로 사용하고, 매 단계의 예측 상태와 잔차 보정 상태에 각각 새 EOS 역산·GR 접합을 수행했다. 이전 실패 종점을 새 예측 상태로 대체하지 않았다. 1단계 에너지 9.9518558e-07 (통과); 2단계 에너지 5.0487311e-07 (통과); 4단계 에너지 4.2294414e-07 (통과). 마지막 조성 시간 점수 0.019911438, 로그 온도 차이 6.1622707e-10로 시간 대조는 통과다. 별도 에너지 단위 엔트로피 역산 GR 재투영의 점수는 8.2800499e-07로 통과다. 이 별도 재투영은 강화한 역산법의 전체 시간 재적분이 아니다.

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
