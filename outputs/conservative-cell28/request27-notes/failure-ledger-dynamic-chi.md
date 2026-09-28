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
