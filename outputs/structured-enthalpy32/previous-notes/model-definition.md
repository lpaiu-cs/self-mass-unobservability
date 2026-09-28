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
