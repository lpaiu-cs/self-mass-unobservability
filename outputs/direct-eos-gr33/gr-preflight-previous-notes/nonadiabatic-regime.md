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


## Request 31 조성 미분과 보존형 GR 재적분

분류: Counterexample candidate. 누적 조성 증분, 정지질량 변화, EOS 엔탈피, 직접 적분한 총 중성미자 손실과 광도 발산을 연결한 별도 원천 적분을 실행했다. 최대 Be7 시간 차이는 구조 재조정의 반응 피드백으로 약 99.995% 재현됐다. 기준 압력 엔탈피 좌표와 전체 함수의 잔차를 사용한 별도 보정도 실행했다. 새 1·2·4단계 고정 광도 GR 경로의 에너지 기준 일괄 판정은 실패다. 가장 미세한 4단계의 에너지 기준은 통과, 마지막 조성·온도 시간 대조는 실패다. 거친 경로의 실패는 미세 경로의 판정과 구분하여 보존한다. 별도 구조 피드백 잔차 보정 1단계의 에너지 점수는 9.7476097e-07로 통과했으며, 그 보정 방법의 시간 수렴은 아직 검증하지 않았다. 물리적 수송과 외부 구동의 응답 pole을 이 결과에서 추정하지 않는다.

근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).


## Request 32 직접 원소 EOS 중간 검증

분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.

분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.

근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 32 구조 피드백 엔탈피 시간 대조

분류: Counterexample candidate. 구조 피드백을 명시적으로 다시 평가한 엔탈피 보정법으로 1·2·4단계를 적분했다. 1단계 에너지 9.9518558e-07 (통과); 2단계 에너지 5.0487311e-07 (통과); 4단계 에너지 4.2294414e-07 (통과). 마지막 조성 시간 점수 0.019911438, 로그 온도 차이 6.1622707e-10로 시간 대조는 통과다. 별도 에너지 단위 엔트로피 역산 GR 재투영의 점수는 8.2800499e-07로 통과다. 이 별도 재투영은 강화한 역산법의 전체 시간 재적분이 아니다. 고정된 광도 면 입력과 초기 행렬을 사용했으며, 이 결과는 물리 수송이나 관측 완화 pole의 검출이 아니다.

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
