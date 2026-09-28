# Failure ledger — dynamic chi

Status: Proven. Request 11.1 isolates an omitted-nuisance bias mechanism: with truncated estimator l and omitted orthonormal residual basis D, an unbounded D-component can cancel any fitted beta signal whenever D.T l is nonzero. At fixed nominal Gaussian width, a finite multiplier cannot ensure uniform coverage. Full nuisance fitting removes this specific mean-bias mechanism in the fixed linear model. See [the registered audit result](../notes/REQUEST11_1_NUISANCE_AUDIT_RESULT.md).

Status: Imported from prior work. Request 11.2's full-space diagonal-noise interval reaches only 58.25% pointwise coverage under its specified extra-Fourier stress. Fitting the specified covariance by REML gives minimum coverage 94.59% across the follow-through conditions. This is conditional linear-model evidence; it neither validates all noise nor demonstrates failure of the historical K=10 grid envelope. See [coverage results](../notes/REQUEST11_2_COVERAGE_RESULT.md).

The unified revision supersedes the early runtime proposals formerly summarized here. Dated REQUEST10 notes and raw artifacts preserve the historical gate record; no frozen verdict is rewritten.

| Status | Condition | Exact failing step | Minimal missing ingredient |
| --- | --- | --- | --- |
| Proven | Reciprocal fast-spectrum assumptions absent | Positivity or the rate-gap quadrature inequality need not hold. | Independently justify conjugate readout, positive dissipation and the gap. |
| Proven | Comparator fits both quadratures freely at every carrier | All six periodic response columns are spanned. | Independently known drive phases/amplitudes or constrained comparator. |
| Proven | Periodic columns reused as an initial-state response | A causal differential annihilator is invisible on them but acts on a transient. | A dedicated transient column or a validated restricted forward operator. |
| Imported from prior work | Coarse three-phase grid treated as a supremum | It misses known narrow long-lag extrema from the historical grid. | Preserve seed extrema, refine and retain the analytic all-phase upper envelope. |
| Imported from prior work | Uniform single-cycle assurance across all gaps | The 223.37-day gap costs only 0.153 in chi-square under the a=1 stress. | A justified noise model and nonlinear pulse-number reconnection; no claim of an actual slip. |
| Imported from prior work | Numerical QR/SVD agreement treated as derivative accuracy | It certifies the supplied matrix, not its physical finite-difference error. | Projected error vectors and an error budget at the weak singular scale. |
| Proven | Inertia not small over the measured band | The charge response has denominator kappa-I omega^2+i Gamma omega, not a single pole. | Control epsilon_I at every carrier and treat homogeneous modes. |
| Proven | Unequal white-dwarf charge/mass ratios | Pair modulation differs by (a_i-a_o)deltaQ_p/m_p. | Pair-specific response columns or justified equality. |
| Proven | Responsive companions near small kappa | Feedback shifts stiffness by -sum C_j/r_pj^2. | Coupled-state treatment or a small-feedback bound. |
| Proven | Historical auxiliary potential drive | Its phase closure is zero instead of the derived 3.11837 radians. | Correctly prescribed amplitudes/phases and a new matched inference; rescaling old beta fails. |
| Proven | Tau treated as a Compton period or resonance | A monotone relaxation pole supplies neither identification. | Matching of inertia, damping and restoring force in a specified theory. |
| Proven | Zero lag, settled forcing | chi=alpha F gives a static coefficient shift. | Finite relaxation with nonzero readout. |
| Proven | Zero frequency after settling | No periodic quadrature remains. | Time-varying drive; an unmodeled transient is a separate signal. |
| Proven | beta=0 in the settled solution | No driven pole enters the readout. | Nonzero alpha c_chi; initial-state amplitudes must be treated separately. |
| Proven | Low-frequency band, finite error tolerance | The Taylor derivative residual is bounded by abs(beta) rho^(N+1). | Precision below that bound or a different sampled band. |
| Proven | One carrier, free F and dot F | Both quadratures are fit exactly. | Shared restrictions or more carriers. |
| Proven | K positive carriers, real degree N>=2K-1 | A shared polynomial interpolates the pole exactly. | Lower order/prior restriction, more carriers or an open frequency band. |
| Proven | K positive samples, complex degree N>=K-1 | Complex interpolation is exact. | A smaller comparator class. |
| Proven | Three distinct carriers, unrestricted real degree five | The explicit P5 in Appendix C fits all six conjugate samples. | A justified derivative-order ceiling. |
| Proven | Zero/repeated drive carrier | A nominal frequency provides no new independent response datum. | Additional nonzero distinct drive. |
| Proven | Independent complex projection per carrier | Lambda_k=O_k/(G_k F_k) absorbs every response. | Calibrated or constrained projection. |
| Proven | T lies in the specified nuisance span | The whitened residual information is zero. | A demonstrable rank increment; finite dimension alone is insufficient. |
| Proven | Singular, blind or pole-cancelling readout | Deprojection fails or removes the pole. | A nonzero nonsingular response channel. |
| Proven | Linear forcing and readout | Superposition creates no sum/difference carriers. | Nonlinearity plus exclusion of nonlinear static competitors. |
| Proven | Arbitrary initial chi ignored | The exponential transient is absent from periodic templates. | Settled initial condition or a fitted transient amplitude. |
| Counterexample candidate | Identifying mass readout with SEP pair coupling | ODE algebra does not derive the physical force law. | Action/force matching including gradients, inertia and backreaction. |
| Imported from prior work | K_dyn=10 treated as a calibrated uncertainty | K is an assumed width multiplier in stored Gaussian intervals. | Independently justified likelihood/noise and interval validation. |
| Imported from prior work | Full nuisance directions omitted from headline | Stored intervals widen by factors about 2.10–17.35. | Report both constructions; justify any truncation or prior. |
| Proven | beta interval relabeled as peak Delta | Carrier filtering and the independent instantaneous term are lost. | Drive normalization, amplitude definition and joint covariance. |
| Proven | Finite-grid maximum called arbitrary-phase marginalization | The sampled origins cover only the registered grid/domain. | Specify the domain or perform a separately designed phase analysis. |

Status: Proven. The analytic outcome is theorem progress: the exact real-carrier interpolation boundary, low-frequency remainder and nuisance-rank condition are explicit. The linear MVP's sideband attempt fails at superposition.

Status: Counterexample candidate. The remaining physical outcome is loophole progress: the shared pole is a candidate only against a restricted comparator and an observable projection that preserves it. Dynamic EFT is established prior art, and A4 is not the unique possible assumption boundary.

Status: Counterexample candidate. Request 11.3 now supplies a force-level EFT realization with explicit validity conditions. Its historical physical-drive interpretation fails a phase gate. This is a recorded boundary result, not a successful empirical scalar-tensor exclusion.

Status: Conjectural. Numerical EOS-to-body matching and a calibrated astrophysical timing likelihood remain outside the completed conditional analysis. They require additional physical inputs and cannot be supplied by prose or a unit-drive rescaling.


## Request 12 follow-through

Status: Imported from prior work. Request 12 half-step timing derivatives change individually by at most 0.6652 percent but rotate weak nuisance directions almost orthogonally; no physical derivative-error certificate follows. A weak 223-day pulse-count compensation demands e>1 at the proposed full-displacement fractions; local admissible probes do not rescue it. Status: Conjectural. Numerical EOS matching, omitted timing-force control and complete nonlinear pulse/noise inference remain incomplete.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Imported from prior work. Request 13 replaces the prior absence of numerical EOS work with a specified SLy stellar structure and outgoing scalar mode calculation. An initially 0.878% independent mass discrepancy was traced to sparse-table enthalpy construction; resampling the same pressure-energy continuum reduces it to 2.73e-7 relative. Timing bounds-check order is corrected and an algebraically equivalent residual evaluation is tested to reduce cancellation. Tightening integration or interpolation settings alone does not provide a rigorous derivative enclosure. Constrained nonlinear timing/noise fits now run inside physical eccentricity domains, but local optimizer output is not a global pulse reconnection or calibrated physical-signal inference. Outstanding completion gates are maintained in [Request 13 plan](../notes/REQUEST13_REMEDIATION_PLAN.md).

## Request 14 validated flow

Status: Imported from prior work. The new CAPD calculation supplies short-time continuum state and initial-state Jacobian enclosures, plus local independent mass derivatives. Direct Cartesian C1 and Hermite-Obreshkov representations exceed the declared Jacobian-width ceiling near six days; exact unit scaling extends this to about 32 days. Generic Jacobi reconstruction initially fails the remainder-inclusion test because exact common-position cancellation is represented numerically. Eliminating that position dependence algebraically restores validated stepping, with the width ceiling reached at about 45.78 days. None completes the approximately 2990-day forward span. See [raw bounds and remaining failures](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).

Status: Conjectural. The missing links to D2 are a useful whole-span variational enclosure, certified physical parameter initialization including masses, and all timing readout/inverse-time terms. No finite-step agreement, short-time success, or midpoint restart can replace these links. Broad interval bounds diagnose this enclosure method; they do not prove physical orbital instability.


## Request 15 후속 검증

분류: Proven. 추가 실패 경계를 기록한다. (1) 수정 전 초기화에는 질량 중심 항의 차원 오류와 회전 비공변 속도 항이 있다. 수정·재매핑 후에도 기존 수치 결과를 새 물리 제약으로 재해석하지 않는다. (2) 4체 운동에 3체 Einstein·Shapiro 지연 API가 연결되어 있으므로 전체 4체 관측 모형이라는 승격은 허용되지 않는다. (3) 질량·Kepler·호출 지점 관측식의 부분 미분 인증은 전체 초기화·누적 적분·보간·역변환의 결합 오차를 대신하지 못한다. 다중 정밀도 실행의 실패·종료 원인은 새 보고서에 개별 보존한다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).


## Request 16 다체 관측식과 영 구동 경계

분류: Proven. Request 15의 4체 운동/3체 지연 불일치는 격리 실행본의 공통 Einstein·Shapiro 합산으로 수정했다. 전체 4체 상대론 정확도나 scalar 관측 완성을 뜻하지 않는다. 지정된 영 scalar 가지에서는 외부 구동을 가정한 항성 산란 응답에서 실제 동반성 구동으로 넘어가는 단계가 성립하지 않는다. 최소 추가 조건은 비영 scalar 초기·입사 자료 또는 배경·scalarized 천체와 그에 맞는 전체 힘/관측 도출이다. 분류: Conjectural. 전 기간 IVP·전체 초기화·시선·누적 지연·보간 산술·역시간 연쇄 및 전역 pulse/noise 추론은 여전히 미완료다. 첫 1 cm 회전 검사 실패와 모의 지연 감사의 중복 단위 변환 오류는 보고서에 보존했다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).


## Request 17 비영 구동과 동반성 응답

분류: Proven. Request 16의 영 구동 경계는 비영 배경을 명시한 지정 후보에서 벗어났으며, 별 세 개의 유한 배경 평형과 선도 상호 구동을 계산했다. 그러나 작은 전체 되먹임이 작은 변동 신호의 상대오차를 보장하지 않아 고정 동반성 근사가 실패했다. 정적 전하 제거와 rank-one 복사의 조건부 빠른 완화는 새로운 느린 관측량을 공급하지 않는다. 차가운 안쪽 백색왜성 모형은 반지름 약 0.0212 태양반지름으로, Kaplan 등의 광학 반지름 0.091±0.005 태양반지름을 재현하지 못한다. 실제 항성 matching에는 열·외피 구조가 필요하다. 분류: Conjectural. 항성·동역학의 전 기간 구간 인증, 전체 매개변수 초기화·광자 전파, 새 쌍별 신호의 비선형 likelihood·pulse/noise 검증은 미완료다. 압력 좌표 실패·수치법 변경·축약 모형의 가정을 보고서에 보존한다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).


## Request 18 열 백색왜성 구조와 응답 경계

분류: Proven. 공개 열 진화 이력에는 광학 후보가 있지만 해당 실행의 저장된 내부 구조 중 사전 질량·광학 기준을 함께 만족하는 것은 0개였다. 최적 이력 행의 전후 구조는 온도 16409 K와 14577 K로 각각 기준을 벗어난다. 단위 보정 뒤 질량도 목표보다 0.2548% 높다. 임의 보간이나 질량 재규격화를 실제 EOS matching으로 승격하지 않았다. 조건부 고정 구각 모형의 산술 구간·주파수 상계는 확보했으나 실제 열 별의 구간 인증과 물리 timing·전체 비선형 추론은 미완료다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).


## Request 19 열 진화 재현과 질량 보정 대조

분류: Proven. 초기 재시작 실패는 기본 cno_extras의 21핵종과 사진의 22핵종 불일치였다. ca40을 추가한 뒤에도 대기·확산 제어의 재시작 누락으로 다음 단계 반지름이 0.343% 달랐다. 이 설정을 원래 사용자 코드에 맞게 복원한 후 상태·한 단계·1,000단계 대조가 통과했다.

분류: Counterexample candidate. 사전 질량·광학 탐색 기준을 함께 만족한 행은 11개다. 기록한 최적 모델 19057은 Teff=15828.041340 K, logg=5.74162539, GM/기준 GM_sun=0.19753638530700이다. 한 외피 제거율과 한 초기 이력만 검사했다. GR 질량 정의, 이력·제거율 및 격자 의존성, 실제 별의 전구간 오차 보장, 비영 배경의 전체 force/readout과 비선형 관측 추론은 별도 미완료다.

세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).


## Request 20 열 구조 민감도와 질량 정의

분류: Counterexample candidate. 고정 구간의 네 변형 검사 판정은 미통과이며 통과 행은 fast 0개, slow 11개, mesh 14개, time 0개다. 검사 밖 나이 구간이나 다른 형성 이력에 대한 결론으로 확장하지 않는다.

분류: Proven. 원래 MESA mass_correction은 배포 원자량과 바리온 조성의 비율이며 내부에너지·결합에너지까지 포함한 ADM 질량이 아니다. 정확한 소스와 조성으로 이를 대조했다. 절대 에너지 기준·고유 부피·구조 반작용이 누락된 단계가 GR 질량 완성의 경계다.

분류: Conjectural. 열 EOS의 절대 기준과 GR 구조를 먼저 완성한 후, 전체 force/readout 및 엄밀한 미분·관측 추론을 연결해야 한다.

분류: Counterexample candidate. 원래 구간에서 실패한 두 실행을 대상으로 별도 등록한 다음 냉각 구간 검사에서는 fast 13개, time 25개의 실제 후보를 확보했다. 원래 실패 판정은 유지한다. 나이 이동과 양립하는 후보 회복이며, 고정 나이의 수렴 인증은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).


## Request 21 열 EOS와 GR 질량 매칭

분류: Counterexample candidate. 고정 온도 경로의 실제 EOS/GR 재구성은 약 16566 K로 원래 광학 컷을 실패했다. 별도 온도 배율 조정으로 얻은 후보는 그 실패를 대체하지 않는다. 같은 EOS의 Newtonian–GR 비교에서 온도 차이는 약 3.24 K이므로 원래 MESA 후보와의 전체 차이를 GR 효과로 해석하지 않는다.

분류: Conjectural. 남은 경계는 원래 22개 동위원소 전체와 진화·열수송·회전을 일관되게 연결하는 실제 항성 모형, 안정성, 유체·metric·scalar 동역학, 전구간 미분 오차와 전체 관측 비선형 추론이다. 이번 수치 질량 허용오차는 EOS 정확도나 관측 질량 오차 인증이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).


## Request 22 GR 후보의 열수송과 물질 보존

분류: Counterexample candidate. 실패 단계는 질량 적분 이후에 원래 X(P)·광도·열 저장 항을 그대로 옮겨 동일한 실제 항성으로 해석하는 부분이다. 총 질량을 맞춰도 수소 재고와 열수송 수지가 남는다. 미정 엔트로피 변화율이나 온도 배율을 조정하는 것만으로 실제 진화 해결로 세지 않는다.

분류: Proven. 직접 압력 좌표의 독립 Gauss 재고 검산은 중심 끝점 때문에 실패했다. 실패 수치를 보존하고 sqrt(ln Pc-ln P) 좌표로 정칙화한 독립 적분이 수소 약 2.03e-8, 총 바리온 약 1.22e-12의 상대 차이로 통과했다.

분류: Conjectural. 최소 추가 조건은 물질 좌표와 각 반응·손실·수송 계수를 닫은 준정적 GR 진화다. Request21 질량 매칭과 고정 온도 광학 실패는 그대로 보존하며, 전체 항성·미분 인증·관측 추론 완료를 선언하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).


## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성

분류: Counterexample candidate. Request22에서 실패한 압력 고정 조성의 보존 해석은 원래 후보의 실패로 보존한다. 새 바리온 좌표 해는 원래 내부 재고와 선택한 EOS 엔트로피를 보존하지만 목표 중력질량과 일치하지 않는다. 별도로 등록한 균일 바리온 질량 조정 모형은 원래 총 재고를 약 0.0698426% 줄이는 다른 후보이며, 이를 원래 물질의 보존 변환이라고 부르지 않는다.

분류: Conjectural. 보존한 엔트로피는 원래 P,T에서 정의한 FreeEOS 값으로 원래 MESA 전체 EOS의 절대 엔트로피 인증이 아니다. 수학적 외곽 대기의 추가 바리온, 구역별 상수 템플릿, 불소 생략과 열·조성 진화 미완료를 계속 명시한다.

세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).


## Request 24 조성 변화의 반응 에너지 기준

분류: Proven. 조성 변화에 고정 조성 T ds 항등식을 그대로 적용하거나, 총 에너지에 정지질량 감소를 포함하면서 핵 가열을 다시 더하는 연결은 실패한다. 전자는 화학적 조성 항, 후자는 일관된 내부/총 에너지 식 선택이 필요하다. 배포 원자량과 표준 Q 표의 작은 차이도 조성 의존 기준 보정으로 명시했다.

분류: Counterexample candidate. 19종 FreeEOS 조성으로 지정 CNO 반응량의 국소 검증을 통과했다. 22종 전체 EOS, 실제 weak Q/중성미자 재평가와 72항목 네트워크 시간 적분은 미완료이다. 모든 prior 실패, 원래 물질 보존 질량 불일치 및 별도 바리온 배율 후보 구분을 유지한다.

세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).


## Request 25 새 GR 상태의 미세물리 평가

분류: Counterexample candidate. 기존 온도 미분과 실제 평가 함수 미분의 최대 차이는 작은 간격에서도 약 0.780312%로 남았다. 스크리닝 제거만으로 사라지지 않았고 weak_rate_factor=0도 반응 중성미자를 모두 제거하지 않아 순수 강반응 분리 대조로 인정하지 않는다. 근본 원인은 아직 특정하지 않았으며 기존 native 미분 통과로 승격하지 않는다.

분류: Conjectural. 같은 입력 상태에서도 MESA/FreeEOS 압력 비는 전체 구역에서 약 0.99172–1.02410이다. 공통 EOS 보조량·전체 조성 변화율·미분 오차 인증과 GR 열수송/시간 적분은 미완료다. 정적 질량 매칭과 실제 반응 원천 평가를 완전한 진화·동적 관측량으로 세지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).


## Request 26 남은 폐쇄 조건의 실행 검산

분류: Counterexample candidate. 단독 반응 합산과 개별 PP 비율 조정의 실패를 보존하고, 공유 PP 분기를 묶은 전체 네트워크 대조를 별도로 통과했다. EOS의 작은 간격 차분은 수렴하지 않았고 HELM 전체 영역 대조도 실패했다. 이를 연속 EOS의 제1법칙 위반으로 확정하지 않으며 전역 미분 보증으로도 사용하지 않는다. 두 미분법의 국소 부력 진단에서 음의 영역을 확인했지만 고유모드 인증은 아니다.

분류: Conjectural. 공통 EOS·독립 native 조성 벡터·열/유체/metric 진화·비영 배경의 교차 결합과 완전한 관측 전방 모형이 남는다. 모든 현재 레버를 실행 결과 또는 정확한 누락 조건으로 원장에 기록했으며, 접근 가능한 검사의 소진을 연구 전체의 완료로 바꾸지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).


## Request 27 직접 반응 벡터와 명시적 PP 상태

분류: Counterexample candidate. 직접 반응 벡터 검증의 미해결 항목은 지정 초기 상태에서 해소했다. 엄격 상대 잔차, 반올림만 포함한 잔차, 후진 Euler의 조성별 시간 정밀화 실패를 각각 보존했다. EOS 차분은 큰 간격에서도 실패했고 REAL 입력 경계에 걸친 출력 점프를 실제 측정했다. 이는 물리적 연속 EOS의 불연속성이나 제1법칙 위반의 증명이 아니다.

분류: Conjectural. 생성 불소와 PP 중간 핵종을 포함한 공통 EOS의 정칙성·오차 제어, 보존된 항성 초기화·열/유체/metric 진화 및 비영 scalar 배경의 구동/전하 연결이 남는다. 국소 핵 가열 상태를 완전한 관측 추론으로 승격하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).


## Request 28 반응 에너지 보존과 scalar 미분 구간

분류: Counterexample candidate. 기본 EOS의 온도 역산은 엄격 에너지 기준에 실패했고 4분할 첫 단계에서 별도 진단 한계도 넘어 중단했다. HELM 분리 적분도 조성 시간 기준에 실패했다. 이 실패를 보존하고 조성·온도·손실의 결합 적분을 별도로 통과했다.

분류: Proven. 고정 GR 보간 모형의 지정 scalar 미분은 구간 연산으로 보증했으므로 이 좁은 항목은 해소되었다.

분류: Conjectural. 물리적 공통 EOS와 에너지 영점 교정, 보존된 26종 GR 초기화·동적 진화, 실제 scalar 구동 및 질량 정규화 전하·관측 추론은 남는다.

세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).


## Request 29 공통 EOS와 GR 열 경로의 관측 연결

분류: Counterexample candidate. 공통 EOS의 초기 유한 미분 검사는 통과했지만 반응 함수의 반환 온도 미분은 EOS 보조 입력을 고정해도 벡터 척도 약 0.107의 차이를 보였다. 입력 좌표 차이는 0이므로 단순 입력 정규화로 설명하지 않는다. 두 eta 미분 인자 대조는 변화가 없어 순서를 판별하지 못했다.

분류: Counterexample candidate. 고정 부피 4분할의 첫 원천 단계는 에너지 역산 +19 erg/g에서 중단했다. 완료된 1·2분할은 에너지·시간 기준을 실패했다. 별도 고정 압력 4분할과 외삽도 엄격 조성 시간 기준을 실패했다. 저장 단계의 독립 EOS 역산 잔차는 전역 질량 에너지 점수와 별도로 판정한다. 원래 실패를 변경하지 않는다.

분류: Proven. 등방 GR 단극 질량 손실만으로는 횡방향 힘을 만들지 않는다.

분류: Conjectural. 최소 추가 연결은 미지원 원소의 물리 EOS·연속 오차 보증, 정확한 반응 미분과 보존적 시간 오차 제어, 자체 수송·유체·계량·대기, 비영 scalar 구동 및 전체 관측 전방 모형이다.

세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).


## Request 30 EOS 역산과 약반응 미분 수정

분류: Counterexample candidate. Request29의 큰 독립 역산 잔차에는 exp(nu)**2와 exp(2nu)의 광도 차감 차이가 주로 기여했다. 동일 강제 입력에서는 최대 잔차/허용량이 1.71875 이하로 줄었고, 별도 고정 초기값 역산은 최대 0.5625로 통과했다. 기존 실패 배열은 보존한다.

분류: Counterexample candidate. 약반응 혼합 가중치만 수정한 벡터 기준은 실패했다. 실제 경로는 약반응 온도·밀도 미분을 모두 0으로 덮으며 Q 및 Qnu의 미분도 전달하지 않는다. 단정밀도 선형 보간을 실제 표로 재현하고, 별도 배정밀도 기여와 완전한 부분 미분으로 초기 유한 기준을 통과했다. 실패한 설정 끝점 대조와 중간 미분 재구성을 보존했다.

분류: Conjectural. 전체 조성 Jacobian, 물리 EOS·연속 오차, 새 GR 시간 경로·자체 수송·유체·계량·대기, 비영 scalar 구동과 전체 관측 추론은 남는다.

근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).


## Request 31 조성 미분과 보존형 GR 재적분

분류: Counterexample candidate. He3의 원래 간격 실패와 온도 혼합 경계를 넘은 중성미자 미분 실패를 보존하고 별도 작은 간격을 대조했다. 가열 구역 마스크 밖의 H2 실패도 별도 기록한다. 음의 광도에서 1 erg/s로 분모가 잘리는 native 출력은 그대로 대류 비율로 해석할 수 없다. 중심 질량 간격을 생략한 수송 출력 재현 실패와 실제 경계를 반영한 대조를 모두 보존한다.

분류: Conjectural. 물리 EOS·전구간 미분, 자체 열수송·대류·대기·완전 GR 진화, 비영 scalar 구동·전하와 전체 관측 추론의 완료는 아직 주장하지 않는다.

근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).


## Request 32 직접 원소 EOS 중간 검증

분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.

분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.

근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 32 구조 피드백 엔탈피 시간 대조

분류: Counterexample candidate. 1단계 에너지 9.9518558e-07 (통과); 2단계 에너지 5.0487311e-07 (통과); 4단계 에너지 4.2294414e-07 (통과). 마지막 조성 시간 점수 0.019911438, 로그 온도 차이 6.1622707e-10로 시간 대조는 통과다. 별도 에너지 단위 엔트로피 역산 GR 재투영의 점수는 8.2800499e-07로 통과다. 이 별도 재투영은 강화한 역산법의 전체 시간 재적분이 아니다. 이전 단계의 시간·에너지 실패는 그대로 보존한다. 내부 확산의 고정 계수 대조에서 명시적 양성 시간 간격은 약 4.64e-5초다. 별도 암시적 밴드 풀이 대조가 통과해도 실제 표면 열손실·대류·대기·비선형 EOS와 GR 결합의 통과를 뜻하지 않는다.

분류: Proven. 수·전하 보존 원소 치환의 Coulomb 전하 제곱합 불일치를 명시했다.

분류: Conjectural. 물리 EOS 및 연속 오차, 전체 GR 진화와 실제 구동·관측 폐쇄는 아직 완료되지 않았다.

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
