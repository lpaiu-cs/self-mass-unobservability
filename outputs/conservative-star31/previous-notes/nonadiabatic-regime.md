# Nonadiabatic regime — exact scope

Status: Imported from prior work. Request 11.2 validates approximately 95% interval coverage within a prespecified full-nuisance, estimated-Fourier-covariance linear model (minimum 94.59% over tested cells). A diagonal covariance can under-cover severely even when the pole survives projection. See [the result and limits](../notes/REQUEST11_2_COVERAGE_RESULT.md).

Status: Proven. A passive first-order pole has no resonance peak. Its tau=Gamma/kappa need not be an inverse particle mass. Close to a scalar-charge instability, increasing susceptibility alone does not ensure negligible inertia, companion feedback or nonlinear response. The explicit physical validity inequalities are in [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md).

Status: Proven. For nonzero real beta and positive tau_chi, the settled response `c_Y+beta/(1+i omega tau_chi)` has a nonzero pole residue. No finite polynomial equals it on an open frequency interval.

Status: Proven. Finite sampling is different. A real degree-N polynomial with unrestricted shared coefficients matches K distinct positive carriers exactly if and only if `N>=2K-1`. Necessity follows from the 2K conjugate roots of `(1+tau z)(P-c_Y)-beta`; sufficiency follows from the explicit polynomial division in Theorem 3 of [the manuscript](../paper/manuscript.md).

Status: Proven. Three distinct positive carriers exclude real degree at most four; degree five suffices to match them. For unrestricted complex coefficients on positive-frequency samples alone, degree `K-1` suffices. Repeated or zero-amplitude carriers do not add independent information.

Status: Proven. Exact noninterpolation does not give a positive lower bound on the measurable residual. Nearby frequencies, small lag, and nuisance projection can make that residual small or zero. For whitened data, the required rank increment is `T_tilde^T (I-P_J) T_tilde>0`.

Status: Proven. An arbitrary independent complex projection at each carrier absorbs the response pointwise. A finite shared projection can also absorb it if its tangent directions span the signal. Rank five instead of six for a particular projection ansatz is not sufficient by itself to identify beta or tau_chi.

Status: Proven. Linear two-frequency forcing and linear readout produce only the input frequencies. Quadratic drive/readout terms can generate sum and difference sidebands, but local nonlinear static responses also generate sidebands. A sideband argument must specify and exclude that comparator.

Status: Counterexample candidate. The observable target is a shared pole relation against a bounded comparator with a justified drive and projection. The three-carrier timing benchmark is a conditional application, not a detection or a claim excluding every finite-order EFT.

Status: Imported from prior work. The stored J0337 finite-grid statistics report no detection; nuisance choices materially change conditional beta intervals. This revision reuses those artifacts without new runtime work.

Status: Proven. For known phases, the coefficient-space residual against a shared derivative comparator is invariant under block phase rotations. Multiplying it by the squared smallest singular value of the full-nuisance whitened carrier map gives a continuous all-phase information lower bound. Allowing the comparator to fit arbitrary independent quadratures is different and restores exact absorption.

Status: Imported from prior work. The item-4 map condition numbers are 98.19–124.82 under the specified covariance scenarios. Nonzero rank can coexist with a roughly 589-fold standard-error increase when fourth-order derivatives are allowed. The mathematical distinction alone therefore does not establish useful empirical precision.


## Request 12 follow-through

Status: Imported from prior work. Dedicated transient response columns pass registered amplitude-convergence gates and add at most 0.523 percent to beta standard errors in the nine tested cases. Status: Proven. Two arbitrarily close positive poles approach one pole quadratically, preventing uniform finite-precision state-count separation without additional separation/weight assumptions.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Counterexample candidate. Request 13 solves the outgoing scalar frequency response and one complex pole for a specified SLy star; changing the exterior contour and radius changes that pole by about 1.5e-8 relative. The low-frequency elastic scattering normalization and a constant-coefficient oscillator fitted to the pole are different reductions. No exact one-pole model over all frequencies or internal state count is inferred from this computation.

## Request 14 validated flow

Status: Proven. Direct interval variational integration removes finite-difference truncation from the frozen GR initial-state Jacobian calculation on its certified time domain. This does not certify dynamic-SEP response columns, physical state count or a complete nonadiabatic observable. The entire timing map still needs its own error propagation. See [certificate boundary](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).


## Request 15 후속 검증

분류: Proven. 다중 정밀도 전파, 초기화 수정·재매핑, 질량·Kepler·대수 관측식의 미분 검증을 수행했다. 이번 계산은 동적 SEP가 없는 GR 기준선에 조건부이며, 전체 기간 타이밍 인증과 비단열 검출 주장을 완성하지 않는다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).


## Request 16 다체 관측식과 영 구동 경계

분류: Proven. 지정된 DEF 영 가지의 물질·궤도 변화는 선형 scalar 방정식의 독립 외력을 만들지 않는다. 시간 의존 계수도 영 초기자료의 영 해를 보존한다. 안정성이나 다른 scalarized 가지의 부재를 증명한 것은 아니다. 이번 결과는 추가 GR 지연 수정 및 연속 보간 오차·영 구동 경계의 정리 진전이며 새로운 비단열 관측량의 확립은 아니다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).


## Request 17 비영 구동과 동반성 응답

분류: Proven. 비영 배경으로 궤도에 따라 변하는 전하는 생기지만, 이번 정적 해는 V_eff=−phi²*1^T L^-1 1/2로 정확히 제거된다. 관성을 생략한 선도 복사 모형의 빠른 집단 완화도 일 단위 상태를 확립하지 않는다. 분류: Conjectural. 비영 가지의 결합 모드·관성·고차 복사와 실제 구동·광자 전파를 일관되게 연결해야 비단열 관측 후보로 승격할 수 있다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).


## Request 18 열 백색왜성 구조와 응답 경계

분류: Proven. 고정 밀도의 약한 중력 열 WD scalar 모형에서 궤도 주파수 응답은 정적 값의 2.1e−10 이내이므로 큰 비단열 응답을 공급하지 않는다. 상반평면의 독립 scalar 성장 해도 축약 조건으로 배제된다. 하반평면의 모든 pole이나 유체·회전·중성자별 모드를 배제하는 정리는 아니다. 분류: Conjectural. 비영 배경의 유체·metric 결합과 실제 상호 구동은 추가 matching이 필요하다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).


## Request 19 열 진화 재현과 질량 보정 대조

분류: Proven. 광학 후보 구조를 보간 없이 복원해도 영 배경·고정 밀도·평탄 시공간의 독립 scalar 모형은 궤도 주파수에서 정적 응답과 극히 가깝다. 실제 진화와 내부 구조의 수치 재현을 확인했지만 유체·metric 결합이나 중성자별 모드의 비단열 응답을 배제하지 않았다.

세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).


## Request 20 열 구조 민감도와 질량 정의

분류: Proven. 네 변형 실행의 선택 구조에 기존 영 배경·고정 밀도·평탄 시공간 scalar 응답과 궤도 주파수 상계를 다시 적용했다. 별의 열 진화 시점 민감도와 독립 scalar 모형의 정적 근사는 다른 결과다. 유체·회전·metric 결합의 궤도 시간척도 상태를 배제하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).


## Request 21 열 EOS와 GR 질량 매칭

분류: Conjectural. 이번 계산은 열 EOS를 가진 GR 평형의 질량·부피 경계를 진전시켰다. 이 구조의 유체·metric·scalar 결합 동역학과 궤도 주파수의 관측 전달함수는 아직 계산하지 않았다. 정적 GR 질량 일치를 새로운 위상 지연이나 동적 관측량의 증명으로 세지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).


## Request 22 GR 후보의 열수송과 물질 보존

분류: Counterexample candidate. 원래 스냅숏의 순 핵반응 열원은 표면 광도의 약 284배이고 큰 음의 중력열 항을 가진다. 끝 상태 출력만의 에너지 잔차도 보존했다. 주변 모델에서 추정한 반지름 변화 시간 약 11617년은 진화 기울기이며 동적 응답의 pole이 아니다.

분류: Conjectural. 다음 경계는 바리온 좌표의 물질·조성 보존, 새 상태의 실제 반응·수송 계수, GR 열·조성 진화와 안정성·유체·metric·scalar 전달함수를 순서대로 연결하는 것이다.

세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).


## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성

분류: Counterexample candidate. 직접 FreeEOS 엔트로피 역산과 독립 DOP853 적분은 두 모형의 연결 잔차를 1e-8 이내로 재확인했다. 이 검사는 반응·열수송을 끈 지정 물질의 정수압 재구성이다.

분류: Conjectural. 다음 경계는 조성이 변할 때의 정지질량·핵반응·내부에너지 기준을 중복 없이 연결하고, 새 상태의 반응·손실·전도·복사·대류와 GR 열·조성 진화를 계산하는 것이다. 안정성 및 유체·metric·scalar 응답은 그 다음이며 정적 재구성을 궤도 완화나 관측 추론으로 확대하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).


## Request 24 조성 변화의 반응 에너지 기준

분류: Proven. Request22의 준정적 GR 광도 식을 반응하는 물질로 확장할 때 dL∞/dB=e^(2ν)[q_nuc−q_thermalν−T ds_B/dτ−Σ μ_i^th dY_i/dτ]로 화학 항을 명시해야 한다. q_nuc는 같은 에너지 기준에서 반응 중성미자를 이미 뺀 가열률이다. 바리온과 함께 움직이는 물질, 종별 확산 유속 없음, 투명한 중성미자 가정이며 확산 시 화학 에너지 유속도 필요하다.

분류: Conjectural. 다음 단계는 새 상태에서의 실제 반응률·중성미자·복사/전도/대류와 GR 열·조성 진화이다. 현재 시험은 반응량을 지정했으며 속도나 시간을 계산하지 않았다.

세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).


## Request 25 새 GR 상태의 미세물리 평가

분류: Counterexample candidate. 기존 핵 가열 미분의 직접 차분 대조는 사전 문턱을 실패했다. 작은 간격의 실제 평가 함수 차분으로 896개 구역의 국소 미분 자료를 별도 구성했다. 이 영역은 절대 핵 가열의 약 99.999808%를 포함하며 특정 상태 SHA에 결속되어 다른 상태로 사용할 수 없다. 유한 차분 수렴은 전역 미분 오차 보장이나 미래 GR 동역학의 검증이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).


## Request 26 남은 폐쇄 조건의 실행 검산

분류: Counterexample candidate. PP 분기비 미분 누락을 실제 채널 출력으로 검산하여 지정 상태의 기존 온도 미분 오차를 설명했다. 보정 후 차분 잔차는 최대 약 2.74e-8이다. 비선형 복사 경계를 가진 GR 동결 계수 열수송 대조는 에너지 수지를 통과했지만 국소 온도 변화가 약 14.86%여서 물리적 시간 진화 인증으로 세지 않는다.

분류: Proven. 열 변수를 제거한 scalar 응답에는 교차 결합의 곱이 들어간다. 긴 열 시간만으로 관측 가능한 느린 scalar pole을 주장할 수 없다.

세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).


## Request 27 직접 반응 벡터와 명시적 PP 상태

분류: Counterexample candidate. 26종 망의 명시적 중간 핵종에서 일 단위의 핵 가열 응답을 계산했다. 지정 두 주파수에서 4상태→1상태 축약의 상대 오차는 약 1.49e-6과 7.32e-9이며 빠른 모드 오차식 안에 든다. 전 주파수의 동등성이나 유한 carrier에서 자유로운 정적 고차 계수로 흡수되지 않는다는 주장은 하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).


## Request 28 반응 에너지 보존과 scalar 미분 구간

분류: Counterexample candidate. 외부 열욕을 제거한 별도 HELM 고정 부피 계산에서 26종 반응과 온도를 함께 진행했다. 2→4분할 조성 오차/허용량은 0.02855, 온도 ln 차이는 2.35e-11, 중성미자 손실 상대 차이는 1.25e-5다. 이는 지정 국소 근사의 시간 대조이며 실제 항성의 열·유체 모드나 scalar 힘의 비선형 전달함수는 아니다.

세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).


## Request 29 공통 EOS와 GR 열 경로의 관측 연결

분류: Counterexample candidate. 5,735개 물질 구역에서 1.6294일의 26종 반응·열·동일 바리온 TOV 경로를 계산했다. 고정 압력 4분할은 에너지와 온도 대조를 통과했으나 조성 시간 오차/허용량 3.37977로 실패했다. 별도 외삽은 Be7 오차를 줄였지만 H1·He4의 수 ulp 차이로 엄격 조성 기준을 실패했다.

분류: Conjectural. 고정 광도·준정적 계량 경로는 자체 수송·유체·계량 진화가 아니다. 조건부 scalar 읽기와 실제 주파수 구동의 phase lag 또는 pole 식별은 구별한다.

세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).


## Request 30 EOS 역산과 약반응 미분 수정

분류: Counterexample candidate. 같은 배정밀도 약반응 함수와 EOS 보조 입력 대조에서 896개 핵 가열 구역의 벡터·열 미분 기준 1e-3을 모두 통과했다.

분류: Proven. 지정된 네 약반응의 고정 밀도·자유전자·조성 모형에 한해 h=5e-5의 미분 계산 및 차분 절단 오차를 구간 보증했다. 정규화 상계 최대는 4.59802e-11 미만이다. 전체 EOS·반응·관측의 연속 미분 보증은 포함하지 않는다.

근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).
