# Observable targets — unified paper

## 단계 35 — native 결합의 실행 범위 (2026-09-17)

분류: Counterexample candidate. native 분자 EOS의 대표 5상태에서 상호적 scalar 구동·무구동 source 적분과 약 2차 시간 수렴을 확인했다. 무구동 drift와 부호 반전 대칭 오차는 저장 정밀도에서 0이었다. 고정 부피의 가역 source와 인위적인 scalar spring을 쓴 대조이므로 이 위상이나 진폭을 항성 전달함수로 읽지 않는다.

분류: Conjectural. 물리적인 scalar 외부 구동, 전 구역 물질·열수송·계량의 동시 진화, 정적 비교를 제거한 전하 응답은 남는다. 새 EOS의 전 셀 보존 초기화는 완료되어 아래의 해당 미완료 기록을 갱신하지만 관측 신호의 완성은 아니다. [단계 35](../notes/REQUEST35_NATIVE_COUPLING_KO.md).

## 단계 34 — 저장 GR 경로의 읽기 판정 (2026-09-17)

분류: Proven. 선언된 정확한 보간 계수 모형에서 70/140/280 경로의 초기 대비 끝점 `alpha/phi_infinity` 변화 절댓값은 각각 5.626e-13/9.134e-13/7.927e-13 이하이다. 유리수 다항식 및 구간 정규화로 얻은 상계이며, 사전 1e-9 읽기 척도에 미달한다. 계수 생성 오차·중간 시각·다른 구동·관측 민감도는 이 상계 밖이다. 닫힌 구면 GR의 내부 열 이동만으로 총 질량 변조를 만들 수 없다는 경계도 유지한다.

분류: Conjectural. 다음 관측 후보는 실제 비영 구동을 가한 경로와 같은 배경의 무구동 경로를 비교하고, 사전 지정한 순간·미분·nuisance 방향을 제거한 잔여 응답이다. 자유 열 이력의 정적 읽기는 전달함수나 검출이 아니다. [단계 34](../notes/REQUEST34_DRIVEN_GR_MILESTONE_KO.md).

## Current GR result — 2026-09-17

Status: Counterexample candidate. The original 70/140/280 paths completed the full coordinate interval 0.42117120910640804 s and passed the frozen five-field endpoint refinement gate. Measured orders for ln rho_B, ln T, v/c, Qtotal/w0 and Qrad/w0 are 1.974655, 2.035766, 2.003874, 2.062330 and 2.062330. Native residuals, exact BDF state/flux replay, conservation, sampled characteristic cones and the original 24-record iteration limit passed. This closes the full-duration endpoint refinement bottleneck as loophole progress. This is a numerical consistency result, not detection of relaxation, scalar charge or an astronomical observable.

Status: Conjectural. Rigorous continuous-time error bounds, physical EOS and conservative initialization of the new EOS, physical exterior/reaction/forcing and nonlinear observational closure remain open. The approved comparison is complete; no further long computation is launched automatically.

[Final result and preserved failures](../notes/GR_EXPERIMENT_REDESIGN_KO.md). Earlier GR progress entries below retain their historical verdicts; their pending full-duration refinement status is superseded by this dated result.


Current nuisance audit: [Request 11.1 result](../notes/REQUEST11_1_NUISANCE_AUDIT_RESULT.md). The full stored 90-direction space is the primary baseline for subsequent coverage work; the rank-71 result remains a sensitivity comparison. No physical prior justifying hard removal of the other directions has been established.

Status: Imported from prior work. [Request 11.2](../notes/REQUEST11_2_COVERAGE_RESULT.md) tests full nuisance and estimated Fourier covariance, obtaining minimum K=1 U coverage 94.59% within its specified family. This does not establish a universal astrophysical interval.

Status: Proven. [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md) derives a leading coplanar scalar-potential drive with unequal carrier amplitudes and phase closure 3.11837 radians for the stored parameters. The historical auxiliary physical-drive family has closure zero; no common time shift fixes it. Its beta_phys numbers are withdrawn as constraints on this realization. The unit-drive beta benchmark remains defined, but does not cover this physical drive family.

Status: Counterexample candidate. The primary target is a common transfer relation across known nonzero drives, `G(i omega)=c_Y+beta/(1+i omega tau_chi)`, tested after a specified nuisance projection. The instantaneous c_Y direction is co-fitted.

| Status | Proposed signature | Required comparison or limitation |
| --- | --- | --- |
| Proven | One-frequency quadrature | Exactly reproducible by freely fitted F and dot F coefficients. |
| Proven | Multi-frequency pole relation | K carriers exclude only real polynomial degree N below 2K-1, with known deprojection. |
| Proven | Sidebands | Absent in the linear MVP; nonlinear static comparators must also be excluded. |
| Proven | Nuisance-projected residual | Identifiability requires positive residual information, not merely a finite parameter count. |
| Counterexample candidate | Prescribed pairwise SEP readout | Requires physical matching before identification with a field-dependent mass or a universal SEP parameter. |

Status: Proven. The beta coefficient is tied to the normalized drive. Its relaxation-only carrier amplitude is `|beta d_k|/sqrt(1+omega_k^2 tau_chi^2)`. A peak is bounded by the sum of these amplitudes. A bound on total Delta also needs Delta0, c_Y and their covariance.

Status: Imported from prior work. The stored finite-origin-grid Gaussian construction gives, at tau_chi=2 days and assumed width inflation K_dyn=10, `U_beta=1.6795275e-9` for truncated nuisance and `3.5336060e-9` for the full stored nuisance construction. At the five tabulated lags the full/truncated ratio grows to about 17.35. No coverage calibration or drive-independent exclusion follows.

Status: Proven. A maximum over a registered finite origin grid is neither Bayesian phase marginalization nor a supremum over every possible relative phase. The outer orbital period is not a common period of the incommensurate carriers.

Status: Proven. A tidal `I2` correction scales as r^-7 in acceleration for a point source; a constant Nordtvedt correction scales as r^-2. Their coefficients cannot be equated without a theory-specific matching relation.

Status: Proven. The reciprocal fast-response comparator obeys S(omega_l)<=R_Lambda*S(omega_h), where S=-Im(H)/omega and R_Lambda=[1+(omega_h/Lambda)^2]/[1+(omega_l/Lambda)^2]. A positive slower pole violates it. The projected witness and a continuous known-phase information bound are derived in [Request 11.4](../notes/REQUEST11_4_COMPARATOR_RESULT.md). Free independent quadratures instead absorb all six periodic columns.

Status: Imported from prior work. After full nuisance and the declared covariance, allowing real derivatives through order four retains as little as 2.87821e-6 of the instantaneous-only information in the tested cases. Expanded independent-phase local refinement raises the unit-drive U by at most 3.72% in the tested cells; an analytic all-phase upper envelope is also reported, without treating local optimization as a certified supremum. See [Request 11.5](../notes/REQUEST11_5_PHASE_STATE_RESULT.md).


## Request 12 follow-through

Status: Proven. Inverting a joint confidence region in six carrier coefficients supports continuous phase/lag inference within its declared mean model. Status: Imported from prior work. Independent validation gives 95.17--95.43 percent inclusion in four specified covariance conditions. An omnibus carrier excess is not a relaxation detection; each tested physical lag section includes beta=0.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Imported from prior work. The specified SLy candidate has computed inner/outer carrier phase lags about 6.19e-9 and 3.08e-11 radians in the stated scattering convention. These are model response calculations, not measured timing lags. Request 13 also performs full live 28-parameter constrained timing/noise fits for alternative pulse assignments. Their stationary-point, coverage, matched-force and rigorous numerical-error gates are distinct and remain open; no relaxation detection is promoted.

## Request 14 validated flow

Status: Proven. The interval flow supplies explicit initial-state Jacobian error bounds and an outward-rounded width bound for the linear geometric-delay readout at recorded epochs. These do not yet bound a full timing residual, the 28 timing-parameter derivatives or the normalized nuisance projector. No empirical rank, exclusion or lag claim is promoted. See [validated variational certificate](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).


## Request 15 후속 검증

분류: Proven. 실제 호출 입력에서 기하 지연 및 Shapiro·수차 지연의 구간 편미분을 검증했다. 방출시각과 Einstein 변환의 연쇄법칙 및 잔차 기반 오차식을 도출했다. 분류: Conjectural. 전체 기간의 운동·Einstein 누적 적분·스플라인·Tempo2 의존성을 연결하고 네 번째 천체의 생략 지연을 포함하거나 상계로 정당화해야 한다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).


## Request 16 다체 관측식과 영 구동 경계

분류: Proven. 고정 궤도·매개변수·펄스 번호의 12,474개 TOA에서 추가 GR 지연의 잔차 효과는 RMS 1.3273 ns, 최대 2.6057 ns였다. 이 값은 새 적합이나 검출 통계가 아니다. 128개 모의 관측의 Einstein 성분 차이는 최대 4.4865 ns였다. 분류: Conjectural. 실제 scalar 신호에는 비영 배경·동반성 전하·초기 또는 입사 구동과 그에 맞는 상호 힘·광자 전파의 도출이 더 필요하다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).


## Request 17 비영 구동과 동반성 응답

분류: Proven. 지정된 선도 모형에서 동반성 전하 응답을 포함하면 안쪽 쌍의 내궤도 변동은 고정 동반성 근사보다 약 36.6배, 바깥쪽 쌍의 외궤도 변동은 약 18.6배 커진다. 절대 진폭은 각각 약 4.88e−17, 7.07e−17의 무차원 힘 결합이다. 이는 TOA 잔차나 검출값이 아니다. 공통 두 쌍 결합 템플릿과 새로운 세 쌍 모형을 동일시하지 않는다. scalar 부호 반전은 물질 관측량을 보존하며 영 배경에서 통상적인 phi 선형 Fisher 근사가 퇴화한다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).


## Request 18 열 백색왜성 구조와 응답 경계

분류: Proven. 지정된 두 열 WD 구각 모형에서 l=0 scalar 산란 응답의 정적 값은 약 1169.86 m다. 안쪽 궤도 주파수에서 정적 값과의 상대 차이는 2.1e−10 미만이라는 조건부 해석 상계를 얻었다. 이것은 해당 독립 scalar 모형의 산란 응답이며 삼중계 힘·타이밍 잔차 또는 관측 검출값이 아니다.

세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).


## Request 19 열 진화 재현과 질량 보정 대조

분류: Proven. 원래 저장되지 않았던 광학 후보 모델 18969의 실제 내부 구조를 재실행으로 확보했다. Teff=15786.146896 K, logg=5.82750698이다. 원래 단위의 GM 환산 질량은 기준보다 0.254509% 높다. 질량을 보정한 별도 진화의 탐색 결과는 다음과 같다. 사전 질량·광학 탐색 기준을 함께 만족한 행은 11개다. 기록한 최적 모델 19057은 Teff=15828.041340 K, logg=5.74162539, GM/기준 GM_sun=0.19753638530700이다. 관측 우도 적합이나 독립 반지름 측정으로 취급하지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).


## Request 20 열 구조 민감도와 질량 정의

분류: Counterexample candidate. Request 20의 실제 진화 행에 동일한 질량·Teff·logg 선택 규칙을 적용했다: fast 0개, slow 11개, mesh 14개, time 0개. 모든 교차와 통과 행을 보존했다. 같은 물리적 나이의 표본 차이와 온도 교차를 구분했으며, 보간값을 실제 내부 구조나 관측 우도로 쓰지 않았다.

분류: Counterexample candidate. 원래 구간에서 실패한 두 실행을 대상으로 별도 등록한 다음 냉각 구간 검사에서는 fast 13개, time 25개의 실제 후보를 확보했다. 원래 실패 판정은 유지한다. 나이 이동과 양립하는 후보 회복이며, 고정 나이의 수렴 인증은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST20_THERMAL_ROBUSTNESS_KO.md).


## Request 21 열 EOS와 GR 질량 매칭

분류: Counterexample candidate. 별도 조정 후보의 GR 질량은 0.197536385339 GM_sun 단위, 반지름은 0.0993123229 배포 R_sun, 조건부 Teff는 15799.999991 K, logg는 5.73936077다. Teff와 반지름은 유지한 모형 광도를 조건으로 조정에 사용한 값이며 독립 관측 예측이 아니다. 두 적분 좌표와 직접 EOS 재평가가 질량 수치 허용오차를 통과했다.

세부 근거: [한글 보고서](../notes/REQUEST21_GR_MASS_MATCHING_KO.md).


## Request 22 GR 후보의 열수송과 물질 보존

분류: Counterexample candidate. 원래 프로필을 정확히 재현해 추가 열수송 자료를 복원했다. 조정 GR 구조에 원래 불투명도·열원을 동결 이식하면 확산 광도비 중앙값은 1.08982이고, 중력열을 포함한 적색편이 순 광도는 0.274013 L_sun으로 유지한 0.552221과 맞지 않는다. 이는 동결 계수 실험이며 새 상태의 실제 반응·수송 계산이나 독립 광학 예측이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST22_THERMAL_CLOSURE_KO.md).


## Request 23 바리온 좌표와 엔트로피 보존 GR 재구성

분류: Counterexample candidate. 별도 바리온 질량 조정 후보의 ADM 질량은 0.197536385307002 GM_sun 단위, 광학 반지름은 68819648.589 m, 유지한 광도를 조건으로 한 Teff는 15834.371 K, logg는 5.743136이다. 이 단계에서 광학량을 적합하지 않았으나 물질 템플릿과 광도가 외부 조건이므로 독립적인 실제 항성 예측이나 새로운 동적 관측량으로 세지 않는다.

세부 근거: [한글 보고서](../notes/REQUEST23_BARYON_ENTROPY_KO.md).


## Request 24 조성 변화의 반응 에너지 기준

분류: Counterexample candidate. 지정된 CNO 순반응량의 직접 EOS 시험에서 기준 보정 누락은 열수지 6.21335 ppm, 화학적 조성 항을 뺀 엔트로피 적분은 0.481113% 오차를 냈다. 보정한 유한 에너지 역산과 독립 미분 적분은 일치한다. 이들은 에너지 장부 검증이며 실제 반응 시간, 궤도 완화시간이나 새로운 관측량은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).


## Request 25 새 GR 상태의 미세물리 평가

분류: Counterexample candidate. 같은 GR 바리온 질량과 적색편이 가중치에서 재평가 핵 가열은 219.464816 Lsun으로 기존 값을 옮긴 156.792765 Lsun보다 약 39.9713% 높다. 열적 중성미자를 뺀 원천은 유지한 광도의 약 397.42배다. 이는 저장·팽창·조성 진화를 풀어야 한다는 순간 원천 진단이며 관측 광도 예측이나 에너지 보존 위반 판정이 아니다.

세부 근거: [한글 보고서](../notes/REQUEST25_FRESH_MICROPHYSICS_KO.md).


## Request 26 남은 폐쇄 조건의 실행 검산

분류: Counterexample candidate. 물질 표본·수학적 대기를 연결한 GR 계량/trace 보간 배경에서 scalar 정적 감수율 약 1166.852694 m를 얻었다. 두 지정 일 단위 주파수의 곡률 제거 monopole 산란 위상/ω는 약 3.8922 마이크로초다. 다체 힘·광자·TOA 전달함수나 단일 완화 pole의 증명은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST26_REMAINING_CLOSURE_KO.md).


## Request 27 직접 반응 벡터와 명시적 PP 상태

분류: Counterexample candidate. 별도로 초기화한 국소 반응 부분계에서 외부 온도 구동→핵 가열 응답의 느린 단일 상태를 분리했다. 최대 가열 구역의 고유시간은 약 5.248988일이며, 고정 GR 배경의 좌표시간으로 약 5.249074일이다. 입출력 미분을 native 재평가로 대조했다. 이는 핵 가열 전달함수이며 scalar 전하·자유낙하 힘·표면 광도 전달함수는 도출하지 않았다.

세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).


## Request 28 반응 에너지 보존과 scalar 미분 구간

분류: Proven. 정적 감수율의 변분은 δχ=−∫[φ²δv+(φ′)²δp]dr이며, 고정 φ∞에서 α_A=−φ∞χ/M의 변분에는 질량 정규화 항도 들어간다. 보편 결합의 영 결합에너지 한계에서 Q=αM이면 내부 에너지 재분배로 δ(Q/M)이 생기지 않는다. 고정 GR 보간 모형의 지정 퍼텐셜 모양 미분은 [251.36919,251.70105] m로 구간 보증했다.

세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).


## Request 29 공통 EOS와 GR 열 경로의 관측 연결

분류: Proven. 등방 에너지 방출의 단극 운동량 변화가 u^μ에 평행하면 정규화 투영 후 횡방향 가속도는 0이다. 순수 GR의 이 단극 경계에서 핵 가열·질량 손실만으로 새 자유낙하 힘을 얻지 못한다.

분류: Counterexample candidate. β=−4의 지정 GR 보간 모형과 고정 광도 열 경로에서 질량 정규화 scalar 읽기를 계산했다. 4분할의 δ(α_A/φ∞)는 약 −5.02184e−11이다. 유한 Wronskian 식은 외부 꼬리·계량·질량 정규화를 포함한다. 실제 구동·비선형 역반응·관측 likelihood 및 이 미소 잔여항의 물리적 오차 인증은 아니다.

세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).


## Request 30 EOS 역산과 약반응 미분 수정

분류: Counterexample candidate. 실제 반응 계산기의 weaklib 경로가 온도·밀도 미분을 0으로 지정함을 끝점 대조로 확인했다. 보정된 원천 함수의 전체 초기 상태 유한 미분 대조가 통과했다.

분류: Conjectural. 이 원천 함수에서 새 GR 시간 경로·질량 정규화 scalar 읽기·실제 구동·관측 모형을 다시 계산해야 한다. 이전 조건부 scalar 수치를 새 미시물리 모형의 결과로 소급하지 않는다.

근거: [한글 보고서](../notes/REQUEST30_CLOSURE_PRECISION_KO.md).


## Request 31 조성 미분과 보존형 GR 재적분

분류: Counterexample candidate. 반환된 조성 편미분과 EOS·조성 보조량을 포함한 전체 미분을 구분했다. 실제 값과 대조하지 않은 Jacobian을 관측 오차 보증으로 사용하지 않는다.

분류: Conjectural. 새 상태의 자체 수송·대기·유체·계량, 실제 scalar 구동과 전체 비선형 관측 추론은 여전히 연결해야 한다.

근거: [한글 보고서](../notes/REQUEST31_CONSERVATIVE_STAR_KO.md).


## Request 32 직접 원소 EOS 중간 검증

분류: Counterexample candidate. Li·Be·B·F를 직접 포함한 24원소 EOS 후보를 별도로 만들었다. 5,735개 기존 원소 호환성 대조는 비트 단위로 일치하고 6개 핵종의 희박 고온 전자수 대조도 통과했다. 실제 전체 조성의 미분 검사는 한 구역에서 실패했으며, 기존 EOS에서도 같은 실패를 재현했다. 실제 전자 교환 인자의 Cody–Thacher η=1 근사 경계 통과를 확인했다. 작은 차분의 국소 통과를 연속 오차 보증으로 취급하지 않는다.

분류: Conjectural. 직접 Fermi 적분 대조, 동위원소 및 물리 EOS 오차, 자체 GR 수송·유체·계량 진화와 실제 구동·관측 연결은 계속 수행해야 한다.

근거: [단계 32 보고서](../notes/REQUEST32_STRUCTURED_ENTHALPY_KO.md).


## Request 32 구조 피드백 엔탈피 시간 대조

분류: Counterexample candidate. 같은 최종 시점의 보정 GR 상태와 누적 조성 증분을 독립 대조했다. 1단계 에너지 9.9518558e-07 (통과); 2단계 에너지 5.0487311e-07 (통과); 4단계 에너지 4.2294414e-07 (통과). 마지막 조성 시간 점수 0.019911438, 로그 온도 차이 6.1622707e-10로 시간 대조는 통과다. 별도 에너지 단위 엔트로피 역산 GR 재투영의 점수는 8.2800499e-07로 통과다. 이 별도 재투영은 강화한 역산법의 전체 시간 재적분이 아니다.

분류: Conjectural. 자체 수송과 방사 대기, 비영 scalar 구동·전하, 전체 비선형 관측 추론은 계속 연결해야 한다.

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

분류: Conjectural. 비균일 셀 모멘트 역산 통과는 실제 항성의 광도·scalar 전하·궤도 신호를 예측한 결과가 아니다. 전 구역의 공유 경계·계량/물질 이동·수송/대기를 연결한 뒤 관측식에 전달해야 한다.

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


## Request 33 첫 실제 전 파수 자기항·중성 열역학 인증과32개 H 표

분류: Proven. 첫 대표 상태의0≤Q<∞ 자기항과 일곱 중성 물리량을13632개 실제 중간 절점·두 끝점·참 근 이동·출력 반올림까지 합쳐 인증했다. 최대 오차 상계 약9.035e−10으로 원2e−7 기준 통과. 전체 확장용 H 표32개도 기준 통과. 분류: Conjectural. 나머지 상태·전체 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 두 번째 전 파수 인증·비선형 GR 보존식과 복사 표면 연결

분류: Proven. 두 번째 전 파수 자기항 상태도 원 오차 기준을 통과했고 첫 상태의 두 적분 경로는 일곱 최종 구간이 겹친다. 비선형 GR 보존·제약 전파와 조건부 표면 광도식을 검산했다. 저장된 양압·영열유속 경계는 무표면층 순수 Vaidya 접합 필요조건에 미달한다. 분류: Conjectural. 전체 상태·물리 EOS·대기 접합·유한 자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 세 대표 상태의 전 파수 자기항 인증 완료

분류: Proven. 세 대표 상태의 전 파수 자기항을51904개 실제 중간 절점과 끝점/참 중성 이동까지 합쳐 인증했다. 일곱 물리량씩21개 출력의 최대 오차 상계는4.932e−8 미만으로 원2e−7 기준을 통과했다. 분류: Conjectural. 전체3206상태·물리 EOS·대기·유한 자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 네 번째 전 파수 사전 인증과 전체3206상태 실행

분류: Proven. 새 사전 상태까지 서로 다른 네 상태의 전 파수 자기항28개 출력을 인증했다. 원2e−7 기준을 유지했다. 분류: Counterexample candidate. 동일 소스로 전체3206상태 실행을 시작했다. 분류: Conjectural. 전체 완료·물리 EOS·대기 접합·유한 자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 인증된 차수 절단과 모든 방사방향 부스트의 고정 계수 감쇠

분류: Proven. Legendre 고차 생략의 전 실수축 오차 상계를 얻고 원 응답 문턱을 유지했다. 양의 잔여계수·서로 다른 아광속 특성근을 가진11470개 고정 열수송 모형은 모든 상수 방사방향 Lorentz 부스트와 유한 비영 파수에서 감쇠한다. 분류: Counterexample candidate. 짝지은 실제 파수 대조의 중앙 속도 비율1.60–1.98을 확인하고 전 파수 재검증을 시작했다. 분류: Conjectural. 전체 물리 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전체 표 첫 추가 묶음과 여덟 전 파수 상태 인증

분류: Proven. 전체 표의 첫 추가 네 상태를 원 오차 문턱과 실제 원시/H 해시로 검증했다. 서로 다른 여덟 전 파수 상태, 총120064개 실제 중간 절점과56개 물리량 출력이 인증되었다. 최대 오차 상계4.932e−8 미만이다. 분류: Conjectural. 전체3206상태·물리 EOS·자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전역 원시 상태 복원의 충분조건·실제 두 해 반례와 운동학적 괄호

분류: Proven. 연결된 양의 온도 가지 위의 원시 상태 유일성 충분조건과 같은 보존량의 두 해를 갖는 명시적 열역학 반례를 증명했다. 저장5735계수는 모든 아광속 속도의 고정 Jacobian 양성 검사를 통과했다. 분류: Counterexample candidate. 실제 EOS27개 목표/복원값은 보존량 기반 운동학적 괄호 안이었다. 분류: Conjectural. 실제 EOS 연속 가지·비균일 셀·열변수·대기·유한 자체 GR·관측은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 적응 차수 전 파수 인증과 두 번째 운동량 경로 대조

분류: Proven. 적응 차수의 세 전 파수 상태에서51904개 파수점·21개 최종 물리량이 기존 오차 문턱을 통과했다. 분류: Counterexample candidate. 두 번째 상태의19024개 정확 파수·57072개 중심 보정 응답 대조도 통과했다. 두 경로의 공유 구적·끝점 의존성은 유지한다. 분류: Conjectural. 전체 상태·물리 EOS·자체 GR·대기·관측 연결은 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 계량 결합 EOS 보존 복원·반지름 좌표 진단과 세 번째 운동량 대조

분류: Proven. 비균일 셀 계량 재구성과 전체 연쇄 도함수를 검산했고 세 번째 기준 상태의 모든 파수/일곱 물리량 경로 대조를 통과했다. 분류: Counterexample candidate. 원12개 역산 실패를 보존하고 누적 상자만 제거한 동일36개 목표가 원 오차 기준을 통과했다. 전체5735셀의 반지름/부피 일관성을 유한 진단했다. 분류: Conjectural. 연속 EOS/공간 오차·진화 Q·자체 GR·대기·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 EOS 도함수 좌표 정정·독립 열 구동과 전체 H 표 인증

분류: Proven. 압력/밀도 도함수 변환과 일반 열 구동 초기식을 검산하고3206개 저장 중심의 H 표 오차를 인증했다. 분류: Counterexample candidate.137640개 변환 비열은 양수이고9개 native 좌표 대조 및12개 수정 척도 잔차가 통과했다. 원 잘못 명명한 비열 해석을 정정하되 산출물은 보존했다. 전체5735셀의 수정 열 구동 평가를 완료했고 저장 Q와의 불일치는 남아 Q_F=Q는 별도 가정이다. 분류: Conjectural. 물리 EOS·연속 미분 오차·실제 구동의 자체 GR 진화·대기·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 독립 열 구동의 초기 변화율 보정과 경계 적합성

분류: Proven. 독립 Q_F−Q가 만드는 초기 변화율 차이를 완전한 시간 행렬에서 풀었고 단열 물질 경계의 필요조건 hR=NQ_F를 얻었다. 초기 정수압이면 열 구동도0이어야 하며 τ만 조정해 해결할 수 없다. 원래 면 반지름에 대응하는 보간 매개변수의 근 구간도 유리수로 감쌌다. 분류: Counterexample candidate. 기존 여섯 셀의288개 보정이 원 열 구동 필드를 그대로 보존하면서 정확한 차분 잔차0을 만족했다. 원래 외곽 반지름의 비영 기울기는 셀 밖 다항식 연장을 포함하며 중심 잔차는 절점 수에 따라 부호가 바뀐다. 물리 경계 인증으로 세지 않는다. 분류: Conjectural. 물리 경계·연속 EOS/미분 오차·유한 자체 GR·대기·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 계량 결합 역산 실패의 목표 기반 예측 진단

분류: Counterexample candidate. 원 전체 역산에서24개 반복 제한 실패를 보존했다. 대표 두 실패의 초기/마지막 Jacobian 유한 대조는 원 기준을 통과했다. 같은 시작점·보존 목표에서 바리온/운동량 예측 뒤 원 전체 Newton 함수를 실행한 별도 후보가 두 실패와 영 목표 대조를 원 복원·잔차 기준으로 풀었다. 실패 뒤 선택한 세 사례이며72개 확대 대조와 전체 역산 완료는 별도다. 분류: Conjectural. 물리 EOS·연속 오차·유한 자체 GR·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 전체 구역의 비선형 열·유체·계량 시간 경로

분류: Counterexample candidate. loophole progress. 정지질량 상쇄를 물질 내부에너지 항등식으로 먼저 처리하여 실제5735구역에서 열·유체·조성 이류·질량·계량을 함께 갱신하는 유한 시간 경로를 실행했다. 동일 종료 시각의4/8/16단계 비교가 지정한 유한 수렴 기준을 통과했다. 원 시도 오류와 이전24개 역산 실패를 보존했다. 분류: Conjectural. 공간/연속 오차, 물리 대기·외부, 장시간 핵반응 결합·물리 EOS·구동/관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 결합 진화의 초기 보존과 두 공간 차분 수정

분류: Counterexample candidate. 초기 질량 편차를 약1.04e-12로 줄이고 열 면 유속과 바리온 원시 변수 차분의 두 불일치를 실제 적분기에서 수정했다. 수정5735구역 4/8/16 전체 경로와 원 유한 시간 수렴 기준을 통과했다. 원 오류 경로·국소 에너지 결함은 보존했다. 분류: Conjectural. 공간 보존·장시간·물리 경계·EOS·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 암시적 결합 진화의 반복 수렴과 인과 영역 이탈

분류: Counterexample candidate. 실제 암시적 결합 단계의 반복 정체는 동일 문턱으로 해소했다. 그러나 같은1.542ms에서 두 세밀한 경로가 최외곽1.074c/1.135c를 보여 초기 아광속 조건의 전체 진화 연장이 실패했다. 원 경로를 보존하고 사전 선언된 두 번째 일정 열시간으로 같은 초기 자료·좌표 격자의 결합 적분을 시작했다. 분류: Conjectural. 전 구간 인과성·경계·물리 EOS/반응/관측 연결은 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 외곽 열 유입 경로 차단과 실제 결합 재적분

분류: Counterexample candidate. 원 마지막 구역값의 면 복사가 초기 경계에5.46e33 erg/s의 열 유입을 만들었다. 공통 반사 면 유속의 별도 계산 경계로 원 EOS·완화시간·시간 절점의 실제 결합 재적분을 수행했다. 같은 두 번째 해상도의1.542ms 상태가 원1.074c에서0.715c로 바뀌었고16개 내부 단계도 모두 아광속이었다. 분류: Conjectural. 반사벽은 실제 대기가 아니다. 전체 시간/공간 대조·물리 EOS/반응/관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 복사·전도 시간 분리와 실제 결합 적분

분류: Counterexample candidate. 반사 경계 장시간 경로도 0.014986초에서 인과성 실패. 두 열유속을 분리해 복사 수송시간의 상태 미분을 포함한 실제 5735구역 적분을 시작했고 첫 단계·native 재생을 통과했다. 분류: Conjectural. 전체 경로·시간 수렴·물리 EOS·비평형 복사·대기·관측 폐쇄는 미완료.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 결합 진화의 공변 운동량 수정과 전체 경로 재개

분류: Proven. 실제 공변 운동량 행에서 누락된 a_t*S를 수정해 계수 2를 복원했다. 분류: Counterexample candidate. 기존 RHS의 오류 재현·수정 대조와 실제 5735구역 새 경로의 첫 단계들을 통과했다. 기존 잘못된 경로는 보존하고 수정된 전체 1/2/4 시간 세분화에 자원을 집중한다. 분류: Conjectural. 전체 종료·시간 수렴, 공간 에너지 보존·물리 EOS·외부/관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 총 에너지와 물질 보존을 실제 GR 결합 진화에 연결

분류: Proven. 바리온·총 에너지·운동량·핵종의 공유 면 수지와 BDF2 항등식을 연결했다. 분류: Counterexample candidate. native EOS를 넣은 새 보존형 5735구역 진화의 실제 두 단계를 통과했고 두 번째 국소 에너지 결함은 초기 열용량 대비 1.2242e-11이다. 에너지 드리프트가 큰 원시 변수 경로를 보존·일시정지하고 보존형 전체 경로에 자원을 집중한다. 분류: Conjectural. 전체 시간 수렴·연속 오차·물리 EOS·외부·반응·관측 폐쇄는 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 보존형 결합 진화의 실제 선 탐색 병목 해소

분류: Counterexample candidate. 보존형 원 경로의 다섯 번째 선 탐색 실패를 보존했다. 같은 실패 상태에서 매 반복 접선을 갱신하는 새 경로가 실제 31개 잔차 점수 0.0037066과 에너지·바리온·26종 보존 및 특성속도 문턱을 통과해 다음 단계로 진행했다. 물리식·시간 간격·수락 문턱은 그대로다. 분류: Conjectural. 전체 시간 수렴·연속 오차·물리 EOS·외부·반응·관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 실제 EOS 조성 응답을 포함한 보존형 GR 비선형 진화

분류: Counterexample candidate. 접선 갱신과 donor 경계 수정만으로는 여섯 번째 실제 단계를 풀지 못한 원 실패를 보존했다. native EOS의 수송 조성 응답을 반복 행렬에 연결한 새 경로가 같은 시간·물리식·문턱에서 여섯 번째와 일곱 번째 단계를 수락했다. 실제 31개 잔차 점수는 각각 0.263285와 0.120397이며 보존·특성속도 기준을 통과했다. 두 방향 근사는 반복 행렬에만 적용되고 실제 26종과 EOS는 유지된다. 분류: Conjectural. 전체 시간 수렴·연속/미분 오차·물리 EOS·복사/외부/반응·관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 보존형 native GR 결합 진화의 첫 전체 시간 경로

분류: Counterexample candidate. 보존형 native GR 결합 진화의 첫 전체 39단계·0.421171초 경로와 저장 상태 재검증을 완료했다. 최대 잔차 점수 0.64653, 최대 국소 에너지 수지 결함/초기 열용량 1.77381e-09, 최대 특성속도 0.47291c로 기존 문턱을 통과했다. 26종은 조성 이류이며 핵반응은 아직 포함하지 않는다. 분류: Counterexample candidate. 이후 같은 초기 상태의 78/156단계와 모든 저장 상태의 native 잔차·보존 수지·특성속도 검증도 완료했다. 동결 1/2/4 시간 대조에서 밀도·온도·속도의 관측 차수는 각각 -0.61169, -0.47868, -0.70413이고 두 열유속은 2.05478이다. 전체 시간 수렴 판정은 미달이며 원래 결과를 보존한다. 분류: Conjectural. 시간 수렴 실패의 해소, 물리 EOS·연속/미분 오차·실제 외부/반응/구동/관측 폐쇄는 미완료다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 새 분자 EOS의 보존 초기 상태 생성 연결

분류: Counterexample candidate. 새 분자 EOS의 GR 해에서 셀별 바리온·에너지를 적분하고 그 값을 보존하는 native 초기 상태와 열유속을 생성하는 코드를 구현해 전체 5735구역 계산을 시작했다. 다섯 실제 EOS 온도 역산과 기존 기호 검사를 통과했다. 분자 GR-4/8의 접합 잔차는 각각 2.6291e-13, 3.1264e-13 미만이며, 반지름·질량·logP의 유한 격자 대조 최댓값 7.04258e-10이 동결 문턱 1e-8을 통과했다. 분류: Conjectural. 전체 보존 초기화 및 새 EOS의 시간 진화는 아직 미완료이며, 물리/연속 EOS·외부/반응/구동/관측 폐쇄도 남는다.

세부 근거: [Request33 한글 보고서](../notes/REQUEST33_DIRECT_EOS_GR_KO.md).


## Request 33 GR 시간 수렴 실패의 원인 추적 (2026-09-14)

분류: Counterexample candidate. 중심부 유체 대조가 실제 GR의 세분도 간 밀도·온도·속도 차이를 최대 0.1071% 이내로 재현했다. 같은 비균일 격자의 일정 매질 음향 대조에서 압력 기울기와 속도 발산의 비호환 보간이 약 10.2313 s⁻¹의 인공 성장률을 만들었고, 가중 수반 관계를 만족시킨 대조에서는 실수부가 4.96e-14 s⁻¹ 미만이었다. 분류: Proven. 해당 선형 대조에서 같은 면 가중치가 남기는 부호 부정의 이차 에너지 항과 호환 가중치의 상쇄를 증명했다. 분류: Conjectural. 전체 GR의 일관된 공간 이산화 수정, 내부 조성 경계 및 전체 시간 수렴 재검증은 남는다. 원래 39/78/156단계의 실패 판정은 유지한다.

세부 근거: [GR 시간 수렴 원인 추적](../notes/GR_TIME_CONVERGENCE_CAUSE_KO.md).


## Request 33 GR 실험 재설계와 실제 결합 예비실험 (2026-09-14)

분류: Counterexample candidate. 수정된 GR 공간 구현은 사전 지정한 두 비영 섭동의 1/2/4 시간 대조를 통과했다. 39개 공통 시각의 차수 최솟값은 1.678, 가장 미세한 시간 해의 상대 오차 최댓값은 0.01395%다. 이는 같은 유한 공간 연산자의 시간 대조이며 새로운 chi 관측량이나 실제 관측 추론 결과는 아니다. 분류: Conjectural. 전체 native 예비실험·전체 기간 및 관측 연결은 남는다.

세부 근거: [GR 실험 재설계](../notes/GR_EXPERIMENT_REDESIGN_KO.md).


## Request 33 완료 경로 재사용과 정확한 EOS 캐시

분류: Counterexample candidate. 전체 예비실험의 두 열유속 시간 차수가 1.3502로 미달하여 관측 승격은 없다. 정확한 EOS 캐시와 동일 방정식의 완료 경로 재사용은 계산 비용을 줄이는 구현 변경이며 새 관측량이나 EOS 인증으로 세지 않는다. 분류: Conjectural. 추가 시간 대조와 전체 GR 시간 수렴·물리 EOS·외부·관측 폐쇄는 남는다.

세부 근거: [GR 실험 재설계 및 가속 검토](../notes/GR_EXPERIMENT_REDESIGN_KO.md).


## Request 33 추가 64단계 시간 판정

분류: Counterexample candidate. 추가 16/32/64 경로를 원 실행 계획의 시각으로 재검증했다. 재사용 32단계의 세 시각에서 생긴 1 ULP 재분할 차이를 출처 검증으로 해소했으며 적분 자료·수락 문턱은 변경하지 않았다. 전 경로의 native·보존·특성속도 검사는 통과했지만 온도 차수 1.10825 < 1.5로 최종 판정은 실패다. 두 열유속 차수는 1.91196이고 모든 끝점 최대 차이는 감소했다. 온도 최대 차이는 같은 외곽 셀 5733에 있다. 원 8/16/32 및 원 39/78/156 실패도 보존한다.

분류: Conjectural. 외곽 온도 시간 수렴과 전체 기간 GR 수렴·물리 EOS·외부·관측 폐쇄는 남는다.

세부 근거: [GR 실험 재설계 및 최종 시간 판정](../notes/GR_EXPERIMENT_REDESIGN_KO.md).


## Request 33 예비실험 32/64/128 수렴 통과

분류: Counterexample candidate. 128단계 완료 뒤 32/64/128의 엄격한 native 재생·보존·특성속도 및 시간 문턱을 모두 통과했다. 다섯 끝점 차이가 모두 감소했고 밀도/온도/속도/총 열유속/복사 열유속 차수는 1.97680/1.73409/2.05741/2.03979/2.03979다. 온도와 두 열유속의 기준은 그대로 1.5다. 검증한 기간은 0.0016452000355719064초의 전체 격자 예비실험이며 원 세 실패는 보존한다. 실제 예비실험의 온도 시간 수렴 병목을 넘어선 loophole progress다.

분류: Conjectural. 전체 0.42117120910640804초의 수정 GR 시간 수렴과 연속 오차 상계·물리 EOS·외부·반응·관측 폐쇄는 남는다.

세부 근거: [GR 예비실험 최종 통과](../notes/GR_EXPERIMENT_REDESIGN_KO.md).


## Request 33 전체 기간 GR 연결

분류: Conjectural. 전체 0.421초 native GR 시간 수렴을 다음 문턱으로 삼아 실행했다. 끝점 다섯 변수와 전 구간 공통 시점의 중심/외곽 차이를 확인할 계획이며, 외부 구동 또는 관측 검출 결과로 해석하지 않는다.

세부 근거: [전체 기간 본실험](../notes/GR_EXPERIMENT_REDESIGN_KO.md).


## Request 33 전체 기간 이류 분기 실패

분류: Counterexample candidate. 전체 기간 42단계의 실제 실패를 재현해 영속도 근처 donor 전환과 기존 분기 접선의 불일치를 분리했다. 조성 응답·계량 미분 대조와 두 donor 교체 수정의 실패를 보존한다. 고정 donor 대조의 통과는 실제 해의 통과가 아니다.

분류: Conjectural. 원 방정식·문턱을 유지하는 실패 구간 분할 적분을 검증 중이며 실제 복구와 전체 시간 수렴은 미완료다. [원인과 근거](../notes/GR_EXPERIMENT_REDESIGN_KO.md).

## Request 33 두 절반 복구의 실패 경계

Status: Counterexample candidate. 동일 native 방정식으로 나눈 첫 절반은 0.04760797602936204148초에서 잔차 0.6866074325480188 및 보존·특성속도 문턱을 통과했다. 두 번째 절반은 잔차 17.560264824407337로 반복 감소에 실패했다. 이 수락 상태는 보존하지만 실패 구간 전체의 복구로 세지 않는다. 동시 donor 및 조성 분기 미분 후보도 실패했고 생산 진화에 적용하지 않는다. 정확한 실패 단계와 근거는 notes/GR_EXPERIMENT_REDESIGN_KO.md에 있다.

Status: Conjectural. 실제 분기 전환에 일관된 반복 연산자의 수렴, 전체 기간 시간 수렴, 물리 EOS·외부·관측 폐쇄는 남는다. 이번 결과는 실패 메커니즘과 부분 수락 상태를 좁힌 loophole progress이며 물리적 no-go나 완전한 진화 인증이 아니다.

## Request 33 native 분기 복구와 원 전체 격자 재개

Status: Counterexample candidate. donor·조성 미분의 분기를 맞추고 원 native 조성 기준값을 유지한 반복 수정이 같은 42단계·같은 시간 간격에서 잔차 18.768935482897138을 0.002270952319166062로 낮췄다. 실제 31개 방정식·보존·특성속도·정확 BDF 유속 재생을 통과했다. 실패한 고정 donor 대조나 두 절반 경로를 해로 대체하지 않는다. 이 복구와 기존 42/64/128 수락 상태로 원 70/140/280 전체 기간 계산을 재개했다. 세부 기준점 차이·실패·체크포인트는 notes/GR_EXPERIMENT_REDESIGN_KO.md에 있다.

Status: Conjectural. 전체 시간 수렴·연속 오차·물리 EOS·외부·관측 폐쇄는 남는다. 현재 결과는 실제 실패 단계 복구와 후속 진화 연결이라는 loophole progress다.


## Request 36 전 구역 scalar 동시 결합

분류: Counterexample candidate. 새 분자 EOS의 5,735셀·26핵종·두 열수송 채널을 scalar와 계량의 같은 비선형 진화에 연결하여 3.213281319e-6초의 구동·무구동 첫 단계 대조를 완료했다. 최종 방법은 에너지 보정이 있는 scalar midpoint 갱신식이며, native 정규화 잔차 3.65633e-4, scalar 상대 에너지 수지 1.32535e-19, 기존 물질 보존 및 특성속도 기준을 통과했다. 두 저장 상태의 native 인자와 재구성 배열도 정확히 재생했다. 최초 가짜 열원에 의한 실패와 물리 교환만 적용한 scalar 에너지 실패는 보존한다.

분류: Conjectural. 보정 없는 원 midpoint 식의 통과나 연속 방정식 인증이 아니다. 운동량 0 회전점에서 수치 보정의 정칙성, 시간·공간 합치성, 해상 가능한 외부 파동과 실제 동반천체 구동, 관측 연결은 남는다. 추가 장기 native 계산 전에 scalar–계량 보존 이산화의 이 병목부터 해결한다. [수정·실패·예산·최종 근거](../notes/REQUEST36_SPHERICAL_COUPLING_KO.md).

분류: Counterexample candidate. 외곽 phi=1e-6*t/h ramp와 0의 비교에서 최대 delta ln T=0.0842616을 얻었지만, 이는 미해상 경계 파동을 포함한 유한 이산 대조다. 실제 항성 온도 응답·전달함수·관측 검출로 해석하지 않는다.


## Request 37 회전점 정칙 scalar–계량 적분

분류: Proven. 두 시간 상태의 이산 질량 제약을 정확히 차분하고 방사 방향 수반을 취하여 운동량·상태 변화량·에너지 결함으로 나누지 않는 scalar discrete gradient를 구성했다. 양의 정칙 가지에서 외곽 질량 보존, 운동량 0의 정칙성, 고정 유한 격자에서 매끄러운 해의 시간 2차 합치성이 성립한다. 공간식은 u=r phi의 중점 운동 에너지와 선형 u의 정확한 dual-cell gradient 에너지를 쓰는 별도 구적으로 정의한다. 매끄러운 내부의 radial mass/lapse 전개와 flat 경계 고유함수 항등식을 검산했다.

분류: Counterexample candidate. 진공 자기중력의 수정 파동 대조 8개가 동결한 기준을 통과했다. 비선형 시간 차수는 phi/Pi=2.00002/1.99882, 공간 차수는 2.00094/1.99527, 최대 끝점 상대 질량 차이는 1.10203e-18이었다. 정확히 Pi=0인 초기 상태와 이후 부호 전환을 포함한다. 저장 끝점의 질량·정확해 오차·차수 재계산도 일치했다. 첫 공간식의 수렴 실패와 수정 실행기의 기본 인자 누락 실패를 보존했다. 이는 theorem progress 및 진공 부문의 loophole progress다.

분류: Conjectural. 새 kernel의 native 열유체 적용, 비평형 열 상태의 정준 폐쇄, 균일한 비선형 PDE 오차 보장, 실제 외부 구동과 관측 연결은 미완료다. Phase36의 보정 native 첫 단계와 이번 무물질 파동 시험을 합쳐 완전한 새 결합 진화라고 하지 않는다. [정의·증명·원 실패·수렴·계산 예산](../notes/REQUEST37_REGULAR_SCALAR_KO.md).

분류: Conjectural. 닫힌 진공 상자의 질량 보존과 약 2차 파동 수렴은 동반천체에 의한 비정적 관측량이 아니다. 회전점에서 정칙한 수치 전달을 확보했으며, 실제 열유속·scalar 교환을 일관되게 결합한 뒤 물리 경계 구동과 정적 비교 제거를 진행해야 한다.


## Request 38 동적 비평형 열 상태와 scalar 교환

분류: Proven. 평형 엔트로피와 열유속을 고정하는 확장을 필수로 요구하지 않는다. 바리온·운동량·계량 및 scalar 일·두 완전한 열수송식을 함께 풀면 비평형 엔트로피 생산 항등식과 v=0,Q≠0의 정칙한 국소 상태 변화가 성립한다. determinant 조건과 역운동량 없는 유한 직접 교환을 기호 검산했다. 같은 로그 계량 일을 정칙한 계량 secant로 쓰면 Phase37의 방사 방향 수반에 물질의 E+R 항과 trace source를 포함할 수 있다.

분류: Counterexample candidate. 실제 native EOS·불투명도·26종 조성을 쓰는 대표 셀의 두 열유속·scalar 동시 국소 적분을 수행했다. 8/16/32단계 및 무열유속·가역 대조가 원 기준을 통과했다. 고정 척도 최대 노름의 시간 차수는 2.01053, 최대 에너지 수지 / 초기 열용량 척도는 5.660e-15였다. 다섯 저장 경로를 fresh native EOS로 재구성했다. 본래 초기 5,735셀에서도 국소 미분 행렬의 정칙 조건을 확인했다. 이는 국소 결합이라는 loophole progress이며 전체 구면 계산은 아니다.

분류: Conjectural. 새로운 수반·공유 면 에너지 유속을 native 구면 운동량·열수송·조성 이류와 동일 잔차에 조립하는 작업이 남는다. 실험의 초기 열유속·외부 계량 변형·scalar 진동자와 시간 척도는 시험용이다. 물리적 열수송 계수·전체 비선형 PDE 오차·실제 구동·관측 폐쇄를 주장하지 않는다. [방정식·조건·native 시험·실패 복구·예산](../notes/REQUEST38_HEAT_COUPLING_KO.md).

분류: Conjectural. 이 native 국소 대조의 약 4.71e-14초는 초기 복사 완화시간으로 정한 시험 단위다. 항성 모드·궤도 구동·비정적 관측량을 얻었다고 해석하지 않는다. 실제 공간·경계의 결합과 정적 비교 제거가 선행되어야 한다.

## Request 39 정칙 scalar와 전 구역 native 열유체 조립

분류: Proven. 동일한 방사 질량 수반에 물질의 로그 계량 일, reciprocal trace source와 공유 면 수송을 함께 넣으면 총질량 증가가 scalar 경계 일·물질 경계 유속·수반 가중 물질 잔차의 합이 된다. 양의 정칙 가지에서 유한 항등식을 기호 검산했으며 역운동량 보정을 쓰지 않는다. 전체 Einstein 제약의 연속 수렴 정리는 아니다.

분류: Counterexample candidate. 작은 16셀 대조 뒤 실제 5,735구역의 native EOS·불투명도·두 열유속·26종 조성 이류·scalar·계량을 같은 첫 단계 잔차에서 풀었다. 구동·무구동 두 경로가 227.69초에 통과했다. 구동 native 정규화 잔차 0.0577209, scalar 잔차 4.186e-24, 질량–경계–잔차 항등식 상대 오차 6.028e-20였다. 무구동 끝점은 새 잔차를 이미 만족해 재사용했고 구동에는 native 호출 13,434회와 세 Newton 보정이 필요했다. 이는 첫 단계 조립이라는 loophole progress다.

분류: Counterexample candidate. 직접 trace 물질 항의 최댓값은 에너지 잔차 허용치의 0.0219279배이며, 이 항을 저장 상태에서 제거해도 native 판정이 바뀌지 않는다. c*h는 외곽 셀 폭의 0.0261배다. 따라서 총 에너지 수지 통과를 직접 물질 피드백의 분리 검출이나 해상된 항성 파동으로 해석하지 않는다.

분류: Conjectural. 다음 병목은 물질 내부와 겹치는 해상된 scalar 구동을 작은 보존 구면 모형에서 만들고 무교환 대조와 분리하는 것이다. 비영 이전 scalar 상태의 다단계 식·시간 수렴, 물질 포함 공간 수렴, 실제 외부 구동·관측 폐쇄는 남는다. [식·기준·실행·한계·예산](../notes/REQUEST39_REGULAR_SPHERICAL_KO.md).

## Request 40 온도 대조의 통과와 전하 읽기의 분리

분류: Counterexample candidate. 초기 바리온 가중 RMS Jordan logT의 구동−무구동 끝점 응답은 3.37240e-4, 구동−beta=0 응답은 2.11865e-4다. 48/96/192단계의 측정 시간 차수는 1.63671/1.64396, 응답/경험적 시간 오차 지표는 26.185/25.852다. beta=0은 물질 frame 결합·source·힘을 함께 없애고 scalar 응력의 중력을 유지한 대조다. 온도 차이를 scalar 전하나 자유낙하 관측량으로 바꾸어 부르지 않는다.

분류: Proven. sin^8 펄스 종료 시 경계 구동의 7차까지 도함수가 0이다. 따라서 구동 jet가 0일 때 증분이 0인 경계장·4차 미분의 고정 다항식만으로는 비영 종료 후 증분을 표현하지 못한다. 내부 전파장·음향 상태·임의 차수 nuisance를 포함하는 모형에 대한 no-go가 아니다.

분류: Conjectural. 선언 쌍성의 목표는 외부 매칭을 거친 delta(alpha_A/phi0)의 정적 비교 제거 잔여 1e-9와 결합 오차 1e-10이다. 이는 관측 민감도가 아니다. 다음은 정규화 전하 민감도 또는 상계를 구하는 작업이며 실제 관측 추론은 남는다. [전체 정의와 판정](../notes/REQUEST40_RESOLVED_RESPONSE_KO.md).

## Request 41 외부 전하와 고정 재고 열 미분

분류: Counterexample candidate. 새 분자 EOS의 24셀 보존 재고와 고정 반경에서 phi_infinity=0.001의 구속 상태를 만들고 비선형 Just 외부장에 매칭했다. alpha/phi_infinity=-4.000350085060503이며, 바리온·scalar·매칭 잔차를 통과했다. 물질의 정수압·열평형이나 실제 궤도 진화는 아니다. native 중앙 차분은 선도 열 미분과 상대 1.09e-6로 일치했다. 계획한 부호 방향이 균일 방향과 같아 6개 중복 상태가 있었으므로 독립 native 미분 방향은 하나로 센다.

분류: Proven. 정확한 저장 계수 모형에서 고정 셀 바리온·조성·반경의 선도 열 미분 delta(alpha/phi_infinity)=sum h_i delta ln T_J,i를 구간 연산으로 계산했다. sum abs(h_i)<1.144e-5다. 질량 정규화의 미분까지 포함한 조건부 선형 상계이며, 실제 EOS 미분 오차와 유한 진폭의 나머지항을 포함하지 않는다.

분류: Conjectural. 선언 쌍성의 1e-9 목표에 이 고정 재고 열 경로만으로 도달하려면 정규화 물질 증폭 약 4.64e8 이상이 필요하다. 실제 증폭·반경 이동·바리온 재분배·조성 변화는 미계산이므로 전체 모형의 no-go가 아니다. 다음은 이 자유도를 포함한 동적 전하 미분과 비영 배경의 정상성 판정이다. [식·상태·조건부 상계·중복과 미완료](../notes/REQUEST41_NORMALIZED_CHARGE_KO.md).

## Request 42 물질 이동·반경 변화와 준정적 전하 읽기

분류: Proven. B=a rho W V의 보존 미분에 부피·계량·면 반경 변화를 포함했다. 움직이는 열유속 물질의 trace는 -E+R+2P=-epsilon+3P이며, 정지 배경에서 delta E와 delta R에는 2Q delta v가 들어간다. 정지 상태용 -E+3P를 움직이는 물질에 쓰지 않는다. 기호 검산을 통과했다.

분류: Counterexample candidate. 셀 바리온·온도·면 반경·속도의 전하 미분을 구현하고 내부 반경 이동·전체 팽창·총 바리온 보존 재분배를 native EOS 차분과 대조했다. 실제 초기 열유속의 속도 미분은 검출 바닥 아래이며 별도 비영 열유속 운동학적 양성 대조만 분해됐다. 저장된 9개 결합 진화 끝점의 같은 재고 준정적 전하 투영에서 구동-무구동은 -1.17210e-8, 구동-물질 결합 제거는 -7.43760e-9이고 경험적 시간 차수 약 1.52/1.50으로 사전 기준을 통과했다. 온도만의 +3.58e-10은 전체와 부호가 반대다. 바리온 재분배가 지배하고 선형 바리온·온도 합의 잔여는 첫 대조에서 약 1.78%다. 새 시간 적분은 없었다.

분류: Conjectural. 이번 연결은 움직이는 물질 상태의 준정적 전하 함수다. 실제 빠른 펄스의 복사 전하, 비영 배경의 유체 증폭, 자유 표면 진화, 실제 궤도 주파수·관측 연결은 미완료다. 반경 미분의 존재를 자유 표면 진화로 세거나, Eulerian 바리온 재분배와 같은 이동의 Lagrangian 반경을 중복 합산하지 않는다. 원 벽·정수압 기준 보정과 과거 결과를 보존한다. [보존 미분·native 대조·수치 연결·남은 의존성](../notes/REQUEST42_FLUID_CHARGE_KO.md).

## Request 43 비영 배경의 동적 외부장과 이동 표면

분류: Proven. 정적 scalar 진공에서 비영 주파수의 Einstein 제약은 delta m=r^2 b Phi delta phi를 준다. 이를 소거한 u=r delta phi의 파동 퍼텐셜은 V=2N^2(m/r^3-Phi^2)다. 같은 incoming 조건에서 움직이는 표면의 나가는 진폭 차이는 (Delta phi_contrast-xi_contrast Phi)/h(R)이며, 표면 이동을 전하로 오인하는 좌표 항을 제거한다. 기호 검산을 통과했다.

분류: Counterexample candidate. 비영 Just 배경의 복소 표면 경계와 무한대 나가는 진폭 전달을 구현했다. 실제 표면의 첫 세 궤도 조화 주파수 및 강한 진공·평탄 공간 대조의 18개 외부 해가 정적 극한·유속·시작점/적분 허용치 비교를 통과했다. 물질 시간 적분이나 새 EOS 호출은 없다. 이것은 동적 외부 연산자의 진전이며 결합 항성의 전하 완성이 아니다.

분류: Conjectural. 실제 비영 배경의 내부 유체·열·조성·계량·scalar, 자유 표면/대기, 동반천체 incoming 구동과 정적 비교 제거를 같은 해로 연결해야 한다. 물리적 배경과 결합 오차가 미완료여서 전체 목표를 닫지 않는다. [유도·구현·검증](../notes/REQUEST43_RADIATIVE_EXTERIOR_KO.md), [전체 완료 판정](dynamic-charge-completion.md).


## Phase44 — 비영 배경의 비정상 내부 결합

분류: Counterexample candidate. 같은 비영 배경의 odd=(plus−minus)/2와 even=(plus+minus)/2−undriven을 구분했다. 온도·속도·내부 scalar 응답 모두 사전 시간 기준에 미달했다. 가장 촘촘한 온도 odd RMS 1.32906e-7은 마지막 시간 차이 2.89846e-7보다 작다. 이 값을 수렴된 신호·검출·상한으로 채택하지 않는다. 무구동 배경의 log 온도 RMS 변화 0.002153은 별도 기록한다.

분류: Imported from prior work. 수치·구현·실패·예산·재현 범위는 [단계 44 보고서](../notes/REQUEST44_NONSTATIONARY_INTERIOR_KO.md)에 기록했다.


## Phase45 — 지수 전파와 물질 일

분류: Counterexample candidate. 같은 126단계 대조에서 큰 직접 파동을 포함한 내부 scalar는 시간 기준을 통과했다. 마지막 차이/응답은 2.57e-9다. 온도·속도는 각각 응답/마지막 차이 1.94433/4.09092로 원 기준 5에 미달했다. 전체 응답 및 미세한 물질 매개 scalar 전하를 수락하지 않는다.

분류: Imported from prior work. 세부 근거와 원 실패는 [단계 45 보고서](../notes/REQUEST45_EXPONENTIAL_COUPLING_KO.md)에 기록했다.


## Phase46 — 비영 내부 결합의 시간 수렴 통과

분류: Counterexample candidate. 같은 초기 상태·진폭·24셀·종료점·원 문턱을 유지한 126단계 비교에서 온도·속도·scalar 읽기가 모두 통과했다. 가장 촘촘한 odd RMS는 1.3135347e-7, 6.1365241e-9, 1.9227544e-6; 시간 차수는 1.72446, 1.92589, 2.07680; 응답/마지막 시간 차이는 6.78333, 13.72977, 3.77986e8이다. 이들은 고정 벽 내부 읽기이며 큰 직접 파동을 포함한다. 미세한 물질 매개 전하·무한대 복사 전하·관측 검출을 뜻하지 않는다.

분류: Imported from prior work. 구현·원 문턱·재생·비용과 제한은 [단계 46 보고서](../notes/REQUEST46_CENTERED_MATTER_KO.md)에 기록했다.


## Phase47 — 같은 재고의 비영 기계적 배경

분류: Counterexample candidate. 새 비영 기계적 배경의 형식적 외부 계수 alpha_A/phi_infinity는 -4.000356686614, Einstein 광구 반경은 68,818,637.74977m다. 같은 재고 GR 배경보다 약 699.70m 작다. 이는 반경·물질·scalar가 함께 조정된 정적 초기 상태의 값이며 관측 질량·실제 구동 응답 또는 동적 전하가 아니다. 유한 압력 절단의 Just 정규화와 실제 물질-진공 접합을 구분한다.

분류: Imported from prior work. 방정식·native 대조·실패 복구·비용·재현 범위는 [단계 47 보고서](../notes/REQUEST47_HYDROSTATIC_BACKGROUND_KO.md)에 기록했다.


## Phase48 — native 외층과 조건부 자유 표면

분류: Counterexample candidate. 조건부 기체 자유 표면은 Einstein 반경 69,193,585.01253m에 있으며, 생략 물질원의 반경 반폭 상계는 2.954e-6m이다. 이는 동결한 배경 및 명시한 EOS 극한에 조건부인 정적 값이다. 전체 EOS·배경 오차 또는 관측 반경 오차로 읽지 않는다.

분류: Conjectural. 작은 외층 정적 질량은 유체 증폭·이동 경계·복사 전하 오차의 상계가 아니다. 실제 내부 응답과 동적 외부를 같은 해로 연결한 뒤 직접 장과 물질 매개 응답을 분리해야 한다. 실제 benchmark의 신호 판정은 아직 없다.

분류: Imported from prior work. 결과·조건부 정리·원 실패·예산·재현 범위는 [단계 48 보고서](../notes/REQUEST48_MATERIAL_SURFACE_KO.md)에 기록했다.


## Phase49 — 자유 표면의 단열 내부·외부 결합

분류: Counterexample candidate. 단계 48의 표면 위치·구간은 EOS 질량 환산 계수 오류 때문에 철회한다. 수정 반경은 69,193,584.98351425m, 물질원만의 조건부 반폭은 2.954e-6m다. 전체 물리 오차가 아니다.

분류: Counterexample candidate. 단열 결합 해의 물질 운동 outgoing scalar 계수는 고정 배경 ADM 질량 정규화로 약 -6.251e-9/단위 입사 scalar다. 공통 최대 구동 척도로 환산한 delta alpha/phi0는 약 1.178e-15다. 같은 행렬의 운동량만 고정한 비교와 차이 방정식으로 분리했으며 자유 진화 항성 간 비교는 아니다. 전하 공간 차이 2.039–2.043%는 사전 2%에 미달한다. 실제 Fourier 진폭·전체 질량 정규화·nuisance 제거·관측 또는 no-go 상계로 해석하지 않는다.

분류: Imported from prior work. 정정·방정식·원 실패·계산값·범위와 재현 근거는 [단계 49 보고서](../notes/REQUEST49_FREE_SURFACE_RESPONSE_KO.md)에 기록했다.


## Phase50 — 반응·화학 에너지와 자유 표면 준정적 연결

분류: Counterexample candidate. 단계 50의 초기 반응 방향을 자유 표면 준정적 구조에 연결했다. 반경 변화율의 격자 상대 차이 2.964e-5와 질량–중성미자 수지 상대 차이 1.453e-4는 원 기준을 통과했다. 실제 시간 진화의 표면 속도나 관측 응답으로 승격하지 않는다.

분류: Counterexample candidate. scalar 전하와 ADM 질량 분모를 함께 미분한 전하 변화율 후보는 -1.199e-19/s지만 격자 차이 22.63%로 2% 기준에 미달했다. 이 값은 검출·상계·no-go에 사용하지 않는다. 단계 49의 별도 기계적 전하 기준 미달도 유지한다.

분류: Imported from prior work. 식·원 실패·자원·재현 범위는 [단계 50 보고서](../notes/REQUEST50_REACTIVE_FREE_SURFACE_KO.md)를 따른다.


## Phase51 — 인과적 중성미자 수송과 동일 원천 연결

분류: Counterexample candidate. 단계 51의 고정 배경/초기 일정 방출률에서는 한 광행시간 0.230805초에 누적 중성미자 에너지의 약 7.74%만 물질 표면을 통과했다. 부피 방출량을 그 시각의 표면 유출 또는 무한대 ADM 질량 변화로 동일시하지 않는다. 이는 중성미자 수송 결과이며 동적 scalar 전하나 관측 신호가 아니다.

분류: Conjectural. 복사 에너지·운동량·압력을 총 응력과 물질 반대 원천에 연결한 같은 시간 해에서 전하를 읽어야 한다. 단계 49/50의 전하 공간 기준 미달과 실제 구동·정적 비교·nuisance·오차 요구사항은 유지한다.

분류: Imported from prior work. 가정·해석해·실행·제한은 [단계 51 보고서](../notes/REQUEST51_CAUSAL_NEUTRINOS_KO.md)에 기록했다.


## Phase52 — 반응·물질·scalar·계량·중성미자 시간 결합

분류: Counterexample candidate. 원 native 초기 반응/중성미자의 결합 부분 응답에서 물질 속도 및 Eulerian scalar의 질량 가중 RMS 시간 기준을 통과했다. 한 광행시간 끝의 값은 6.842e-8 m/s 및 5.986e-26이다. 두 RMS는 수치 읽기이며 무한대 동적 전하나 관측 신호가 아니다.

분류: Conjectural. 초기 핵반응 부분 응답의 작은 값으로 전체 열 매개 효과를 제한하지 않는다. 전하 공간 기준 미달, 실제 구동·물리 열 경계·정적 비교·nuisance와 오차 요구사항은 유지한다.

분류: Imported from prior work. 식·대조·예산·범위는 [단계 52 보고서](../notes/REQUEST52_REACTIVE_CAUCHY_KO.md)에 기록했다.


## Phase53 — 물리 열수송 구성식의 식별 경계

분류: Proven. 같은 정적 열 계수의 양의 소산 모형들이 서로 다른 시간 전달을 가질 수 있다. 따라서 정적 EOS·수송 일치와 시간 지연 검출만으로 해당 지연의 실제 물리적 기원을 식별할 수 없다.

분류: Conjectural. 이번 구성식 비유일성은 실제 항성의 신호 부재 또는 no-go가 아니다. 실제 구동·동적 전하를 예측하려면 수송 모형의 물리 입력과 오차를 먼저 특정해야 한다.

분류: Imported from prior work. 정확한 반례·물리 입력 계약·출처·접근 실패는 [단계 53 보고서](../notes/REQUEST53_THERMAL_CLOSURE_KO.md)에 기록했다.


## Phase54 — 미시적 전도 연결과 저에너지 적용 경계

분류: Counterexample candidate. 20% 정상 전도도 호환 통과는 관측 지연의 확인이 아니다. 저장한 상관 후보의 tau_eff와 그 값으로 정규화한 동적 대조는 저운동량 발산 때문에 수락하지 않는다.

분류: Proven. 선언한 탄성 연장의 저주파 s^(2/3) 경계층은 유한 완화시간 전개를 막는다. 이를 실제 항성 신호나 정적 비교를 벗어난 관측량으로 승격하지 않는다.

분류: Imported from prior work. 원식·실패·복구·수치·증명·예산은 [단계 54 보고서](../notes/REQUEST54_ELECTRON_COLLISION_RESPONSE_KO.md)에 기록했다.


## Phase55 — 에너지를 교환하는 전자 충돌 연산자

분류: Counterexample candidate. 대표 세 상태의 유한 기저 첫 시간 모멘트는 1.858–2.010e-17 s다. 선언 쌍성의 세 번째 조화에서 국소 충돌 메모리만의 상대 변화 상계는 2.24–2.42e-20이다. 이 값은 공간 열확산·전역 유체·광자·scalar 전하나 관측 신호의 상계가 아니다.

분류: Conjectural. 같은 EOS 및 열유속식에 상태별 충돌 계수를 연결하고 적용 오차를 판정해야 실제 동적 전하 읽기로 나아갈 수 있다.

분류: Imported from prior work. 식·모형·원 실패·복구·수치·예산은 [단계 55 보고서](../notes/REQUEST55_ELECTRON_ENERGY_EXCHANGE_KO.md)에 기록했다.


## Phase56 — 내부 전도와 GR 시간 경로의 연결 및 수렴 미달

분류: Counterexample candidate. 최신 전도 성분의 scalar 질량 RMS는 약 5.606945e-30이며 지정 시간·계수·외부 경계 대조를 통과했다. 속도 성분의 시간 차수는 미달이다. 이 국소/유한 영역 내부 읽기를 무한대 동적 전하, 실제 궤도 구동의 신호나 관측 no-go로 바꾸지 않는다.

분류: Conjectural. 큰 반응 성분으로 정규화한 합성 응답의 통과는 작은 전도 기여의 정확도를 보증하지 않는다. 물리적 내부-외층 열유속 연결과 전체 오차·정적 비교·관측은 남는다.

분류: Imported from prior work. 식·범위·실패·수치·예산·다음 결정은 [단계 56 보고서](../notes/REQUEST56_CORE_CONDUCTION_COUPLING_KO.md)에 기록했다.


## Phase57 — 축퇴 경계를 넘는 전자 전도 연산자

분류: Counterexample candidate. 단계 57의 국소 전도 응답은 강한 축퇴 선택선 바깥에서도 유한 양의 pole을 갖고 지정 구적·기저 대조를 통과했다. 분류: Conjectural. 이 짧은 국소 시간 모멘트를 별 전체의 열 기억이나 관측 가능한 이탈 신호의 no-go로 해석하지 않는다. 실제 구동·동일 재고 정적 비교·nuisance·전하와 관측의 판정 조건은 유지한다.

분류: Imported from prior work. 모형·원 실패·수정·수치·예산·남은 연결은 [단계 57 보고서](../notes/REQUEST57_WARM_CONDUCTION_KO.md)에 기록했다.


## Phase58 — 조성 전이층을 해상한 전도 응답의 반경 연결

분류: Counterexample candidate. 단계 58에서 보간된 전도의 정상 계수·시간 모멘트와 omega*tau=0.01,1,100의 복소 응답을 분리 검증했다. 최종 독립 두 상태의 최대 차이는 0.03010%, 0.06164%, 정규화 응답 0.03951%다. 분류: Conjectural. 이를 관측 가능 신호나 전 주파수 오차 인증으로 바꾸지 않는다. 실제 구동·정적 비교·전하·관측 연결은 미완료다.

분류: Imported from prior work. 원 실패·수정·독립/개발 표본·예산·채택 파일은 [단계 58 보고서](../notes/REQUEST58_RADIAL_CONDUCTION_KO.md)에 기록했다.


## Phase59 — 같은 EOS의 광자 분리와 원자료 수락 경계

분류: Proven. 양의 불투명도 분포의 정확한 산술/조화 모멘트와 양의 하한이 주어지면 단일 조화 pole과 복소 손실 응답의 차이를 조건부로 제한할 수 있다. 분류: Counterexample candidate. 실제 미세 native 그룹은 같은 구간의 분할 가법성 대조에서 Rosseland 평균이 0.740286% 달라 원 기준 미달이다. 따라서 얻은 형식적 오차식을 실제 항성 응답·검출·no-go로 승격하지 않는다.

분류: Imported from prior work. 식·수치·원 실패·수정·예산·채택 및 거부 입력은 [단계 59 보고서](../notes/REQUEST59_PHOTON_INPUTS_KO.md)에 기록했다.


## Phase60 — 원 광자 격자 복원과 물질 흡수 교환

분류: Counterexample candidate. 원 격자 복원으로 기존 두 온도와 새 독립 두 온도의 정적 평균을 동시에 회복했고, 같은 물질 EOS의 주파수별 흡수 교환을 시간 진화했다. 초기 1 K 섭동은 성분 검증용이며 실제 항성 구동이나 검출 결과가 아니다. 분류: Conjectural. 실제 구동·정적 비교 모형·전하 공간·관측으로의 연결은 미완료다.

분류: Imported from prior work. 원인·독립 자료·식·수치·한계·재현은 [단계 60 보고서](../notes/REQUEST60_NATIVE_PHOTON_MEASURE_KO.md)에 기록했다.


## Phase61 — 광자 공간 수송과 Compton·물질 교환

분류: Counterexample candidate. 공간 열 응답과 에너지·유속·압력·물질 힘을 같은 시간 상태에서 읽었다. Compton 제거는 지정 모드의 초기 노름 기준 1.5371% 차이를 만들었다. 분류: Conjectural. 이 성분 검증 입력은 실제 궤도 구동·동적 전하 관측량·정적 비교나 nuisance 제거를 대체하지 않는다. 실제 불투명도와 미분 산란 입력의 물리 오차가 다음 우선 병목이다.

분류: Imported from prior work. 식·원 실패·유한시간 대조·민감도·예산은 [단계 61 보고서](../notes/REQUEST61_PHOTON_SPATIAL_COUPLING_KO.md)에 둔다.


## Phase62 — 유한 주파수 이동 광자 커널의 결합

분류: Counterexample candidate. 동일 자유 전자 조건에서 유한 이동 커널과 확산식의 응답 차이는 초기 노름의 1.51638%, 유한 이동 커널과 탄성식만의 차이는 0.03342%다. 단계 61의 큰 확산식 스펙트럼 변화를 검증된 물리 신호로 채택하지 않는다. 분류: Conjectural. 현재 성분 입력은 실제 궤도 구동·동적 전하 관측량·정적 비교 및 nuisance 제거를 대체하지 않는다.

분류: Imported from prior work. 식·근사 판정·수치 대조·보존·예산·남은 경계는 [단계 62 보고서](../notes/REQUEST62_FINITE_JUMP_PHOTONS_KO.md)에 둔다.


## Phase63 — EOS 점유수와 일치하는 집단 광자 결합

분류: Counterexample candidate. 같은 셀 적분의 독립 전자 경로와 집단 산란 경로의 지정 응답 차이는 초기 노름의 8.7835063e-06다. 원 이산 속도 구적 차이가 효과 크기와 비슷한 경로는 미수락으로 보존했고 연속 유전 스펙트럼의 결합 대조로 수정했다. 분류: Conjectural. 이 성분 효과는 실제 궤도 구동·동적 전하·정적 비교 및 관측 신호의 완료를 뜻하지 않는다.

분류: Imported from prior work. 식·원 실패·표현 수정·수치 판정·예산·남은 경계는 [단계 63 보고서](../notes/REQUEST63_COLLECTIVE_PHOTONS_KO.md)에 둔다.


## Phase64 — 실제 온도 흡수와 좁은 선 적분

분류: Counterexample candidate. 직접 온도 원자 후보와 기존 TOPS 후보의 끝점 응답 차이는 초기 에너지 노름의 0.0348483782, 물질 온도 차이는 +0.00137241374K다. 이는 원자 모형·점유수·선 표현을 함께 바꾼 후보 비교이며 TOPS 온도 보간 오차, 실제 항성 신호 또는 자유낙하 전하 검출이 아니다. 분류: Conjectural. 실제 구동·정적 비교·관측 연결은 그대로 남는다.

분류: Imported from prior work. 원 실패·선 적분 식·수치 대조·물리 경계·예산·재현은 [단계 64 보고서](../notes/REQUEST64_ACTUAL_PHOTON_OPACITY_KO.md)에 둔다.


## Phase65 — EOS 준위 공급과 광학 상세평형

분류: Counterexample candidate. 90개 H·He II 전이의 조건부 상세평형 잔차는 최대 3.99e-14다. 이는 선언한 전이 법칙의 검사이며 불투명도 전체나 관측 응답의 오차가 아니다. 새 광자 결합·항성 진화를 실행하지 않았고 기존 관측 미완료 판정도 유지한다.

분류: Imported from prior work. 구현·원 실패·수치·조건부 경계·다음 연결은 [단계 65 보고서](../notes/REQUEST65_EOS_OPTICAL_POPULATIONS_KO.md)에 둔다.


## Phase66 — 공통 원자 자유에너지와 광학 상태

분류: Counterexample candidate. 공통 준위 후보에서 저장 외곽 상태의 압력 변화는 -7.47e-7이지만 C I/N I/N III/O III의 재고 변화는 각각 +21.78%/+48.19%/-14.15%/-18.54%다. 총량의 근접성을 근거로 이전 흡수 스펙트럼을 그대로 채택할 수 없다. 분류: Conjectural. 이것은 원자 모형 변화이며 자유낙하 신호나 관측 오차가 아니다. 실제 단면적·구동·정적 비교를 연결해야 한다.

분류: Imported from prior work. 구현·수치·원 실패·모형 경계·예산은 [단계 66 보고서](../notes/REQUEST66_SHARED_ATOMIC_FREE_ENERGY_KO.md)에 둔다.


## Phase67 — 공통 원자 단면적과 역방출

분류: Counterexample candidate. 실제 전이 진폭을 사용하는 공통 원자 흡수·역방출 공급기를 구현하고 저장 상태의 점별 상세평형을 통과했다. 분류: Conjectural. 이는 전체 스펙트럼 오차, 새 광자 결합 응답 또는 관측 신호 검출이 아니다. 소멸/누락 성분과 선폭·공명 적분을 닫은 뒤 같은 물질·광자 시간식에서 응답을 판단한다. 실제 구동·정적 비교·관측 매칭 기준은 유지한다.

분류: Imported from prior work. 식·결과·원 실패·예산·남은 경계는 [단계 67 보고서](../notes/REQUEST67_COMMON_ATOMIC_RATES_KO.md)에 둔다.


## Phase68 — 새 공통 입력의 실제 결합 진화

분류: Counterexample candidate. 실제 결합 끝점에서 원자 입력 제거 대조와의 차이는 초기 정규화 노름으로0.0692391이며, 시간 차이는1.474e-6이다. 입력이 실제 시간식에 작용했음을 확인한다. 분류: Conjectural. 누락 채널과 EOS가 다른 단계64 온도 응답과의 차이를 신호나 오차 상계로 해석하지 않는다. 궤도 구동·정적 비교에서 살아남는 전하·관측 판정은 미완료다.

분류: Imported from prior work. 수치·예산·단위 정정·원 실패·완료 경계는 [단계 68 보고서](../notes/REQUEST68_COMMON_INPUT_COUPLED_KO.md)에 둔다.


## Phase69 — H/He II와 자유-자유 흡수의 실제 결합

분류: Counterexample candidate. 같은 EOS·초기 조건에서 H/He II와 자유-자유 채널을 추가하자 끝점 온도 교란이0.97902174 K에서0.81368576 K로 바뀌었다. 추가 채널의 결합 효과이며 개별 채널 분리나 관측 전하 신호는 아니다. 분류: Conjectural. 실제 구동·정적 비교·관측으로 가는 완료 경계는 유지한다.

분류: Imported from prior work. 실제 연결·결과·모형 경계·다음 시간 폐쇄는 [단계 69 보고서](../notes/REQUEST69_HHE_COUPLED_KO.md)에 둔다.


## Phase70 — 유한 이온 점유의 실제 결합

분류: Counterexample candidate. 유한 이온 점유의 실제 국소 응답과 순간 LTE의 상태 차이는 초기 노름 대비0.04684768이며 수치 대조보다 크다. 이 차이는 지정 ms 실험의 모형 효과이며 관측 전하·궤도 pole·정적 비교 제거를 입증하지 않는다.

분류: Imported from prior work. 열역학·실제 경로·원 실패와 비용은 [단계70 보고서](../notes/REQUEST70_POPULATION_COUPLED_KO.md)에 둔다.


## Phase71 — 실제 GR 반경 입력 재연결과 실패 지속

분류: Counterexample candidate. 반경 전도 입력을 적용한 실제 GR 후보의 속도 및 scalar 시간 차수가 모두 미달했다. 끝점 후보값을 수락된 전하·관측 신호로 사용하지 않는다. 큰 전체 원천으로 기존 실패 부근의 국소 시간 대조를 가리지 않았다.

분류: Imported from prior work. 실제 경로·원 기준·비용과 남은 원인은 [단계71 보고서](../notes/REQUEST71_GR_RADIAL_RECONNECT_KO.md)에 둔다.


## Phase72 — 실제 GR 과도응답 수정 후보의 기각

분류: Counterexample candidate. 고차 시간법의 전체 속도/scalar 차수는 -0.355131/1.293933으로 원1.5 기준에 미달한다. 국소 읽기 개선이나 작은 절대 크기로 전체 성분의 실패를 구제하지 않는다. 관성 후보의 scalar 계수 대조6.0424%도 원2% 기준 미달이며 물리 신호로 채택하지 않는다.

분류: Imported from prior work. 실제 경로·원 판정·비용·거부 이유는 [단계72 보고서](../notes/REQUEST72_GR_TRANSIENT_REPAIR_KO.md)에 둔다.


## Phase73 — 원 GR 직접 전파의 국소 실패

분류: Counterexample candidate. 직접 전파의 전체 속도와 scalar는 지정 주파수 절단 대조를 통과했으나 기존 절단면 속도의37.4005% 차이와 음의 차수는 미달이다. 전체 RMS로 작은 국소 성분의 실패를 숨기지 않는다.

분류: Imported from prior work. 실제 적용·비용·원 판정은 [단계73 보고서](../notes/REQUEST73_GR_DIRECT_PROPAGATION_KO.md)에 둔다.


## Phase74 — 제약 보존 전파의 재판정

분류: Counterexample candidate. 정확한 열 원천 맵과 원 위치 제약의 동적 좌표 제한을 구현했으나 축약계 전파는 기각했다. Laguerre 시간 전개를 원 GR 식에 실제 적용한 경로도 기존 절단면 속도 등에서 미수렴이다. 원 실패·문턱을 유지하고 추가 경로를 실행하지 않았다. 분류: Conjectural. 다음은 원 연속식의 에너지 형태와 압력·열 원천·공간 이산화의 합치성이다. 현재 공간식의 인과적 책임이나 전체 물리·전하·관측 완료를 주장하지 않는다.

분류: Imported from prior work. 실제 경로·기각 이유·비용은 [단계74 보고서](../notes/REQUEST74_GR_MATRIX_PROPAGATION_KO.md)에 둔다.


## Phase75 — 실제 GR 전파 통과와 공간 속도 실패

분류: Counterexample candidate. 같은 반경 전도 입력의 실제 네 읽기 모두 고정 격자 전파·계수·외곽 대조를 통과했다. 기존 절단면 속도의 마지막 전파 차이는0.00003337%다. 그러나 속도 공간 차이5.78–29.46%가 남아 검증된 물리 관측량이나 원 실패 완전 해결로 승격하지 않는다. 광자·전하·구동·정적 비교·관측 연결 요구를 유지한다.

분류: Imported from prior work. 실제 경로·수치·원 실패·예산은 [단계75 보고서](../notes/REQUEST75_GR_CANONICAL_EVOLUTION_KO.md)에 둔다.


## 단계76 — 실제 공간 수정 후보의 미수락

분류: Counterexample candidate. 단계76에서도 동일한 65개 시각의 전체 속도/scalar 및 기존·새 경계 속도 읽기를 유지했다. 직접 적분의 scalar와 새 경계 성분 기준 통과를 다른 방법의 속도 통과와 합성해 수락 해로 만들지 않는다. 전체 속도와 기존 경계 속도의 시간 차수는 미달했다.

분류: Conjectural. 네 성분 동시 전파 수락과 이미 선언한 공간 대조가 남는다. 국소 성분의 개선은 관측 신호나 정적 모형 밖 식별성의 증거가 아니다.

분류: Imported from prior work. 상세 식·원 실패·실행 예산·판정은 [단계76 보고서](../notes/REQUEST76_GR_SPATIAL_REPAIR_KO.md)와 연결된 저장 근거에 둔다. 이번 분류는 loophole progress이며 연구 완성은 아니다.


## 단계77 — 직접 속도 차이 감소와 전체 수렴 미수락

분류: Counterexample candidate. 같은65개 시각의 네 읽기를 유지했다. 4차 직접 시간법의 마지막 상대 차이는 전체 속도0.00636476%, scalar1.46644e-7%, 기존 경계 속도0.0100882%, 새 경계 속도2.60418e-7%였다. 기존 경계만 차수 기준도 통과했다. 방법마다 통과한 성분을 합치지 않으며 작은 단계 차이를 참해 오차로 부르지 않는다.

분류: Conjectural. 실제 관측량 연결에 앞서 네 성분 공통 오차 제한과 고정 공간 대조가 남는다.

분류: Imported from prior work. 식·원 실패·실측 예산과 판정은 [단계77 보고서](../notes/REQUEST77_GR_COUPLED_REPAIR_KO.md)에 보존한다. 이번 분류는 loophole progress이며 연구 완성이 아니다.

## Phase78 — 네 성분의 산술·전파 병목 분리

분류: Counterexample candidate. 네 native RMS 정의를 유지했다. 약한 읽기의 기존 단계 차이는 같은 단계 산술 변경에 민감했고 두 속도의 차이는 그 변경보다 훨씬 컸다. 작은 상대 차이를 참해 오차 상계로 부르지 않는다. 네 성분의 공통 수렴은 미통과이며 공간·계수·외곽·구적 대조는 시작하지 않았다.

분류: Conjectural. 다음 수정 대상은 같은 원천의 전체 행렬 해·투영 직접 해·투영 고유분해 해의 전달 오차로 분리한다. 단계 수·기저 수를 먼저 확대하지 않는다.

분류: Imported from prior work. 실제 코드·예산·실패 판정·저장 재생은 [단계78 보고서](../notes/REQUEST78_FOUR_COMPONENT_CONVERGENCE_KO.md)에 둔다.

## Phase79 — 공통 시간 전파 통과와 새 경계 공간 실패

분류: Counterexample candidate. 같은 네 native RMS의 시간 전파 기준을 함께 통과했다. 작은 차이를 참해 오차 상계 또는 관측 가능한 신호로 부르지 않는다. 가장 거친512계수의 약64% 차이에서 얻은 큰 비율을 엄밀한 방법 차수로 해석하지 않는다. 전체 목표와 original_failure_resolved는 미완료로 유지한다.

분류: Conjectural. 다음 수정 대상은 열 입력 종료 면 부근의 공간 파동 표현과 원천 점프 처리다. 통과한 전파를 유지하며 자동 격자 확대 없이 계획·비용·중단 기준을 재평가한다.

분류: Imported from prior work. 실제 경로·원 실패·비용·독립 재생은 [단계79 보고서](../notes/REQUEST79_GR_JOINT_PROPAGATION_KO.md)에 둔다.

## Phase80 — 네 성분의 공통 수렴 기준 통과

분류: Counterexample candidate. 전체 속도·scalar·기존 절단면·새 지지 경계의 네 native RMS가 한 경로에서 시간·공간 및 원 조건부 대조를 함께 통과했다. 이 수치 통과를 새 관측 신호나 실제 물리 오차 인증으로 확대하지 않는다.

분류: Imported from prior work. 원 기준·실제 결과·예산·완료 경계는 [단계80 보고서](../notes/REQUEST80_GR_COMMON_CONVERGENCE_KO.md)에 둔다.


## Phase81 — 온도 되먹임의 실제 결합

분류: Counterexample candidate. 같은65시각의 네 native 응답에 대한 온도 되먹임 보정은 최대3.99064e-11이며 원 응답 대비1e−4 판정 문턱 아래다. 이것은 GR 상태 읽기의 제한된 영향 판정이고 새로운 관측량·위상 지연의 비흡수성·전하 검출을 뜻하지 않는다.

분류: Imported from prior work. 식·원 계획·비용·독립 검사·범위는 [단계81 보고서](../notes/REQUEST81_GR_TEMPERATURE_FEEDBACK_KO.md)에 둔다.


## Phase82 — 광자·물질의 실제 반경 GR 결합

분류: Counterexample candidate.18.5344μs 실제 결합에서 물질 온도 변화는 약±0.25352K, GR 역응답을 끈 대조와의 끝점 차이는 최대0.0007521K다. 이는 국소 상태 응답이며 자유 계수의 정적 모형을 벗어나는 관측량이나 동적 전하 검출이 아니다. 전파된 scalar 가중 RMS는 약6e−52로 상대율만으로 의미를 부여하지 않는다.

분류: Imported from prior work. 식·원 실패·보정·수치 판정·범위는 [단계82 보고서](../notes/REQUEST82_PHOTON_MATTER_RADIAL_GR_KO.md)에 둔다.


## Phase83 — 실제 궤도 조화와 기계적 전하

분류: Counterexample candidate. 첫 세 조화의 기계적 전하 기여에서 동일한 정적 전달값을 뺀 진폭 합은 delta(alpha)/phi0=1.22546314e-19다. 자유로운 실수4차 미분 비교 후 여섯 전하 계수 노름은3.77561361e-26다. 실제 timing·잡음·nuisance는 미적용이다. 분류: Proven. 계수 제한 없는 실수5차 비교는 이 세 복소 성분을 보간하므로 그 비교 아래에서는 세 성분만으로 식별할 수 없다.

분류: Imported from prior work. 식·수치 판정·비교 경계는 [단계83 보고서](../notes/REQUEST83_ORBITAL_CHARGE_FEM_KO.md)에 둔다.


## Phase84 — 전체 단열 외부 응답과 질량 정규화

분류: Counterexample candidate. 전체 단열 별의 정적 표면 비교는 전파만으로 생기는 큰 차이를 포함한다. 같은 계량의 시험파를 비교에 넣으면 첫 세 실제 궤도 조화의 복사 전하 차이 진폭 합은2.84403183e-15다. 추가4차 미분 잔여2.80174e-24는 유한 대조 변화5.62897e-24보다 작아 비흡수성이 미확정이다. 전체 관측 상계나 no-go로 해석하지 않는다.

분류: Imported from prior work. 식·대조·판정은 [단계84 보고서](../notes/REQUEST84_FULL_ORBITAL_EXTERIOR_KO.md)에 둔다.


## Phase85 — 실제 궤도 구동의 전도·GR 결합

분류: Counterexample candidate. 실제 첫 세 궤도 조화에서 계산된 내부 전도가 추가하는 정규화 복사 전하 진폭 합은3.02319751e-29다. 단계84의 잔여에 합쳐도4차 미분 잔여는 대조 변화보다 작다. 이 전도 성분이 기존 수치 비흡수성 미확정을 뒤집지는 못하며 전체 물리 신호 상계는 아니다.

분류: Imported from prior work. 식·원식 검사·수치 판정·범위는 [단계85 보고서](../notes/REQUEST85_ORBITAL_CONDUCTIVE_GR_KO.md)에 둔다.


## Phase86 — 전 반경 내부 복사 수송과 입력 합치성

분류: Counterexample candidate. 회색 복사와 미시적 전도를 합친 세 조화의 정규화 복사 전하 기여 합은2.88409708e-26이다. 단계84 전체 단열 잔여와 합친4차 미분 비교 잔여는2.80080232e-24로, 기존 수치 대조 변화5.62896703e-24보다 작다. 현재 전체4차 비흡수성 미확정 판정은 유지한다. 이 전하 계수 노름은 관측 잡음 노름이나 실제 전체 신호 상계가 아니다.

분류: Imported from prior work. 식·원 실패·실측 비용·수치와 범위는 [단계86 보고서](../notes/REQUEST86_RADIAL_RADIATIVE_TRANSPORT_KO.md)에 둔다.


## 단계87 — 실제 열 배경 진화와 외곽 수지

분류: Counterexample candidate. 단계87에서 국소 온도 변화는 0.370053ms에 1%지만 해당 셀의 질량 비중은 1.03909e-14다. 이 값만으로 전체 전하 실패나 상계를 주장하지 않는다. 빛의 전파 시간을 통과하지 않은 경계 값을 outgoing 관측량으로 읽지 않았으며 실제 전하에 대한 배경 민감도가 남는다.

상세: [단계87 보고](../notes/REQUEST87_THERMAL_BACKGROUND_DRIFT_KO.md).


## 단계88 — 배경의 수송 변화를 실제 궤도 전하에 전달

분류: Counterexample candidate. 실제 배경의 유한 수송 갱신이 바꾼 세 조화 외부 전하 진폭 합은6.38986e-32로 계산됐다. 변화 자체의 공간 대조9.36~9.40%는 원2% 기준에 미달하므로 정밀 확정값이나 물리 상계로 채택하지 않는다. 갱신 후4차 잔여2.80080e-24는 기존 대조 변화5.62897e-24보다 작아 전체 비흡수성 미확정은 유지된다.

상세: [단계88 보고](../notes/REQUEST88_BACKGROUND_TRANSPORT_CHARGE_KO.md).


## 단계89 — native EOS 온도식과 실제 결합 정정

분류: Counterexample candidate. 정정한 실제 세 궤도 조화의 열적 외부 전하 진폭 합은2.92076668e-26다. 전체4차 비교 잔여2.80078666e-24는 기존 수치 대조 변화5.62896703e-24보다 작아 전체 비흡수성이 미확정이다. 분류: Proven. 자유 실수 미분 계수의 세 복소 주파수5차 보간 경계는 유지된다. 분류: Conjectural. 계산된 열 성분을 전체 관측 상계나 보편 no-go로 바꾸지 않는다.

분류: Counterexample candidate. 정정 공지: 단계81과85–88의 cp 기반 총 밀도 온도 및 열 되먹임 해석은 단계89로 대체한다. 원 소스·수치 판정·실패는 보존한다. 단계87의0.370ms 온도 교차와 단계88의11.4% 수송 변화 및6.39e-32 전하 변화는 현재 물리 추정으로 인용하지 않는다. 단계88 유한 수송 갱신 자체를 수정 온도로 재계산한 것은 아니다.

상세: [단계89 보고](../notes/REQUEST89_NATIVE_TEMPERATURE_CLOSURE_KO.md).


## 단계90 — native EOS 외향 복사 외피

분류: Counterexample candidate. 새 외피 광도2.1733737400e33 erg/s와 원 내부 확산 광도는0.584212768% 다르다. 원 궤도1주기의 방출 에너지/기존 ADM 에너지9.64475e-17은 조건부 에너지 진단이며, 전하·SEP·관측 상계로 해석하지 않는다. 수정된 전체 배경과 관측 nuisance를 연결하기 전 전체 비흡수성은 미확정이다.

상세: [단계90 보고](../notes/REQUEST90_NATIVE_RADIATIVE_ENVELOPE_KO.md).

## 단계91 — 외부 강제항과 최종 전하의 구별

분류: Counterexample candidate. 실제 외피 크기의 짧은 시험 방출로 생긴 추가 외부 scalar·계량 응답을 계산했다. 2R·3R에서 읽은 (r/R)δφ는 최대 약3.12e-35·3.63e-35이나, 표면 δφ를 고정한 부분 문제의 값이다. 내부 응답·질량 정규화·장시간 구동이 빠져 있으므로 정규화 전하·전체 신호·관측 상계로 해석하지 않는다.

분류: Conjectural. 최종 판정에는 새 EOS 배경과 표면 연속성, 그 배경의 물리 구동 및 정적 비교 모형이 모두 필요하다. 시험 펄스를 천체물리적 구동의 발견으로 세지 않는다.

상세: [단계91 보고](../notes/REQUEST91_NATIVE_PHOTON_EXTERIOR_KO.md).

## 단계92 — 재구성된 배경의 전하 읽기 경계

분류: Counterexample candidate. 동일 물질 재고로 재접합한 배경의 정적 alpha는약-0.0040003566866002다. 큰 기존 정적 값과의 작은 직접 차이를 새 동적 전하로 채택하지 않는다. 동일 매개변수를 사용한 정밀도 대조에서 표면 전하가 동일한 것도 전하 수렴의 증거가 아니다.

분류: Conjectural. 최종 차이는 새 배경의 열 재조정·물질·scalar·계량 결합 응답에서 직접 계산해야 한다. 현재 이 기여의 소멸이나 관측 상계는 미확정이다.

상세: [단계92 보고](../notes/REQUEST92_NATIVE_WHOLE_STAR_MATCH_KO.md).


## 단계93 — 새 배경의 실제 결합 진화와 남은 GR 공간 실패

분류: Counterexample candidate. 새 배경의3.434ms 경로에서 온도 대조는 통과했으나 속도/scalar 공간 최대 차이는192.49%/15.618%로 원3% 기준을 넘었다. 표면 scalar 생산값3.72518e-37을 최종 전하로 채택하지 않는다. 분류: Conjectural. 광자 방출의 작은 에너지 비율이나 수렴한 온도만으로 전하 오차를 제한할 수 없으며, 물리 경계·해당 GR 응답과 인과적 외부 읽기가 먼저 필요하다.

상세: [단계93 보고](../notes/REQUEST93_NATIVE_COUPLED_READJUSTMENT_KO.md).


## 단계94 — 측정점 정렬 시험과 이동 표면 보존식

분류: Counterexample candidate. 단계94의 scalar 시간·공간 차이는 각각0.03793%,0.8294%로 지정3% 기준을 통과했다. 속도 공간 차이는211.49%로 실패했고 별도 질량 적분법 대조도 미실행이다. 따라서 이 scalar 부분 통과를 전체 수렴, 최종 전하 또는 정적 비교 모형을 벗어난 신호로 해석하지 않는다.

상세: [단계94 보고](../notes/REQUEST94_NATIVE_WAVE_COLLOCATION_KO.md).


## 단계95 — 보존 연속 원천과 독립 고차 대조

분류: Counterexample candidate. 단계95의 일관 질량 p4/p6 대조는 scalar 공간 차이0.02224%, 온도0.04198%를 보였지만 속도6.9121%로 원3% 기준을 넘었다. 부분 scalar 통과는 최종 정규화 전하나 그 오차 상계가 아니다. 최초 GLL p2/p4 통과를 전체 수락으로 쓰지 않는다.

분류: Conjectural. 남은 유체 오차의 전하 영향과 이동 표면 접합을 판정해야 하며, 원천 보간에 대한 물리 불확실성도 별도로 남는다.

상세: [단계95 보고](../notes/REQUEST95_NATIVE_CONSERVATIVE_SOURCE_KO.md).


## 단계96 — 표면 상쇄 수정과 실제 이동 광자 성분의 결합

분류: Counterexample candidate. 기준 유속 상쇄를 수정한 실제 표면 이력의 이동 광자 추가 GR 응답을 계산했다. 이 구간 광자 압력 일은 약1.10447erg다. 분류: Conjectural. 작은 일/방출 비율이나 표면 scalar를 최종 정규화 전하 상계·정적 모형 비흡수·관측 검출로 바꾸지 않는다. 같은 차수의 전체 경계 항과 원 내부 공간 오차의 전하 영향이 아직 필요하다.

상세: [단계96 보고](../notes/REQUEST96_NATIVE_BALANCED_MOVING_RAYS_KO.md).


## 단계97 — 계량과 전 반경 수송의 동시 되먹임

분류: Counterexample candidate. 같은 해의 계량–열유속 되먹임 연결 후 끝점 delta L/L0는7.76075e-29에서2.74300e-29로 변했다. 내부 성분의 전후 최대값 정규화 차이6e-16 이하는 독립 신호 검출 또는 정확한0의 증거로 쓰지 않는다. 분류: Conjectural. 유지 압력만으로 외부 scalar 응력이 결정되지 않으므로, 이 상태에서 최종 정규화 전하 상계를 확정할 수 없다. 실제 대기 또는 그 영향의 통제가 필요하다.

상세: [단계97 보고](../notes/REQUEST97_NATIVE_METRIC_TRANSPORT_KO.md).


## 단계98 — native 비선형 진공 팽창과 보존 GR 원천

분류: Counterexample candidate. 보존식으로 얻은 자유 팽창의 정규화 scalar 원천 모멘트는2.18076e-32,73/145상태 차이는0.17413%다. 이것은 표면 alpha*A*N과 면적을 고정한 국소 원천 적분이며 최종 무한대 전하·관측 신호·정적 비교 비흡수성으로 읽지 않는다. 정지질량 구적 잔여를 그대로 넣으면약-2.02e-27의 가짜 원천이 생겨 보존량을 먼저 소거해야 한다.

상세: [단계98 보고](../notes/REQUEST98_NATIVE_VACUUM_RELEASE_KO.md).


## 단계99 — 구면 보존 기체 유동과 탄성 광자 교환

분류: Counterexample candidate. 구면 기체의 trace 적분 변화는-1.89043e24erg,896/1792셀 차이는1.61232%다. 이는 최종 정규화 전하가 아니다. 에너지 잔여를 해당 trace 변화로 나누면1.35650%이며 원 보존 판정의 분모와 별도 공개했다. 점별 희박 전면 유속의 수렴이나 엄밀한 전하 오차 상계를 주장하지 않는다.

상세: [단계99 보고](../notes/REQUEST99_NATIVE_RADIAL_RELEASE_KO.md).


## 단계100 — 보존 대기 유출의 직접 지연 전하

분류: Counterexample candidate. 원3.434431ms 이력으로 판독한3.430428ms 끝점의 직접 정규화 전하 변화는−2.5896961e−30이다. 두 저장 격자 파형 차이는0.101344%, 시간 보간 차이는0.026487%다. 최종 전체 전하·관측 신호가 아니며, 한 경계 해제 파형을 정적 EFT 비흡수성이나 물리적 완화 검출로 해석하지 않는다.

상세: [단계100 보고](../notes/REQUEST100_NATIVE_RELEASE_CHARGE_KO.md).


## 단계101 — 보존 대기 에너지와 직접 전하 재계산

분류: Counterexample candidate. 에너지 보존 유출 이력을 같은 직접 Green 계산에 적용한 끝점 정규화 전하는−2.5897362e−30이다. 두 격자 파형 차이는0.108729%, 이전 엔트로피 진화 파형과의 최대 차이는0.004026%다. 원 에너지 잔여를 제거해도 이 직접 성분은 남지만 전체 전하·정적 EFT 비흡수성·관측 검출을 확정하지 않는다.

상세: [단계101 보고](../notes/REQUEST101_NATIVE_CONSERVATIVE_ENERGY_KO.md).


## 단계102 — 안쪽 경계 수정과 실제 계량 원천 적용

분류: Counterexample candidate. 경계를 수정한 직접 지연 전하 끝점은−2.58976186e−30, 추가한 물질 계량 응력 성분은−1.08782202e−35, 층 내부 강제 질량 성분은−1.09146725e−37이다. 계산한 성분 합은−2.58977284e−30이지만 전체 GR 전하나 그 오차 상계가 아니다. 산란 광자의 지연 원천·scalar 되먹임·실제 내부·꼬리·화학을 해결하기 전에는 신호 생존이나 정적 EFT 탈출을 판정하지 않는다.

상세: [단계102 보고](../notes/REQUEST102_NATIVE_METRIC_RELEASE_KO.md).


## 단계103 — 보존 부피 되먹임과 퍼텐셜 응답 상계

분류: Counterexample candidate. 선언한1792셀 원천·배경의 전 차수 scalar 퍼텐셜을 포함한 끝점 정규화 계수는[−2.5925617092e−30,−2.5869839773e−30]으로 제한된다. 지정 퍼텐셜 항만으로 해당 성분이 소거되지는 않는다. 이는 격자·EOS·물리적 꼬리·물질 변위·광자·Bondi 질량·전체 관측 오차를 합한 구간이 아니며 최종 물리 전하의 생존 판정은 아니다.

상세: [단계103 보고](../notes/REQUEST103_NATIVE_CONSERVED_WAVE_KO.md).


## 단계104 — 실제 대기의 유한 화학과 순간평형 경계

분류: Counterexample candidate. 단계103 구간은 지정 LTE 유출 원천의 성분 결과로 보존한다. 그 원천의 재결합 요구가 도입한 반응 사건 예산을 초과하므로 이 구간을 최종 물리적 전하의 생존 판정으로 채택하지 않는다. 누락 중성화 반응을 모두 제한한 결과는 아니다. 새 화학·유출을 연결하기 전 동일 외부장 계산의 정밀도를 늘려도 이 결론 공백은 닫히지 않는다.

상세: [단계104 보고](../notes/REQUEST104_NATIVE_ION_CLOSURE_KO.md).


## 단계105 — 고정 화학 EOS 저온 실패 해결

분류: Counterexample candidate. 고정 화학 입력의 저온 평가 실패는 해소됐지만 실제 유출 및 전하에 아직 적용되지 않았다. ln(rho/rho0)=−6에서 같은 초기 엔트로피의 고정 화학 압력은 LTE의0.05167이므로 기존 LTE 전하 성분을 최종 물리 전하로 채택하지 않는 경계는 유지한다. 이 대조는 전하의 부호나 생존 판정이 아니다.

상세: [단계105 보고](../notes/REQUEST105_NATIVE_COLD_POPULATION_KO.md).


## 단계106 — 고정 화학 입력의 실제 유출·직접 전하 적용

분류: Counterexample candidate. 같은896셀 직접 전하는 LTE−2.5925218e−30, 고정 화학−2.5931668e−30이다. 파형 변화0.02680%는 고정 화학 격자 차이0.41636%보다 작아 화학 모형의 작은 차이를 분해한 검출이 아니다. 두 지정 모형에서 직접 성분은 같은 부호로 남으며 물질 재분배가 지배한다. 최종 물리 전하·정적 EFT 탈출·관측은 미완료다.

상세: [단계106 보고](../notes/REQUEST106_NATIVE_INVENTORY_FLOW_KO.md).


## 단계107 — 유한 수소·광자 교환의 실제 유출 연결

분류: Counterexample candidate. 실제 유한 수소 진화의 직접 전하 끝점은 -2.5932036022e-30이다. 고정 화학 대비 파형 변화0.0061509%는448/896셀 차이0.420584%보다 작다. 원 trace 공간 판정7.7739% 실패를 보존하므로 이 파형은 진단 읽기이며 수락된 최종 물리 전하나 관측 검출이 아니다.

상세: [단계107 보고](../notes/REQUEST107_NATIVE_HYDROGEN_EXCHANGE_KO.md).


## 단계108 — 실제 배경 온도의 반응 내부·대기 접합

분류: Counterexample candidate. 실제 배경 온도의 반응 접합은 원 trace 공간 차이1.9117%로2% 기준을 통과했다. 이 성공을 단계107 진단 전하에 자동 이식하지 않는다. 같은 저장 계량의3.434ms 스칼라 과거 광원뿔은 표면 아래약1029.61km까지 이르며, 이번 전체 물질 구간은400m 깊이다. 더 깊은 체적 열·화학 원천을 정하거나 제한해야 최종 전하를 판정할 수 있다.

상세: [단계108 보고](../notes/REQUEST108_NATIVE_REACTIVE_INTERFACE_KO.md).


## 단계109 — 깊은 내부의 비선형 수소·열·광자 결합

분류: Counterexample candidate. 깊은 내부의 실제 비선형 물질·광자 시간 결합은 완주했으나16셀 적분 trace는 최종 전하가 아니다. 지연·스칼라 공간 가중·질량 정규화가 없고, 움직이는 대기와 같은 광자 유속으로 연결하지 않았다. 거친 P1 공간 재구성에서 양의 모멘트 경계를 넘는 지표가 있어 물리적 광자장 해석을 수락하지 않았다. 이전 국소 전하와 직접 비율 비교하지 않는다.

상세: [단계109 보고](../notes/REQUEST109_NATIVE_CAUSAL_PHOTONS_KO.md).


## 단계110 — 양의 각도별 광자 진화와 미량 분자 제약 수정

분류: Counterexample candidate. 단계110의 표면 방출 스펙트럼은 시간1.5294%·각도1.2122% 차이로 원2% 기준을 통과했다. 그러나 전체 내부 trace 시간 차이4.8889%는 실패했고 전체 해를 수락하지 않았다. surface-port-candidate.npz는 후보 입력일 뿐 움직이는 대기와의 공동 진화나 최종 전하가 아니다. 단계109의 다른 P1 시간 통과를 새 각도 모형에 이식하지 않는다.

상세: [단계110 보고](../notes/REQUEST110_NATIVE_ANGULAR_PHOTONS_KO.md).


## 단계111 — 비분할 물질·광자 시간 결합 통과

분류: Counterexample candidate. 단계111의 전체 이력 trace 시간 차이는0.0372244%, 표면 스펙트럼 시간 차이는0.00627956%로 원2% 기준을 통과했다. 각도 대조도 통과했다. 끝점 체적 trace −5.5704101095e32erg는 비정상 초기장의 고정 배경 과도응답이며, 지연·GR 되먹임·질량 정규화가 포함된 최종 전하나 궤도 완화 관측량이 아니다. 표면 입력 파일은 공간/주파수·실제 대기 연결이 남은 후보로만 저장했다.

상세: [단계111 보고](../notes/REQUEST111_NATIVE_UNSPLIT_PHOTONS_KO.md).


## 단계112 — 실제 양방향 광자와 움직이는 대기

분류: Counterexample candidate. 단계112의 두 전체 시간 경로에서 내부 trace·대기 trace·표면 스펙트럼 대조를 통과했고 대기에서 내부로 돌아오는 실제 복사와 기체 냉각이 생겼다. 이 과도 체적 trace는 최종 정규화 전하나 궤도 완화 검출이 아니다. fine 경로의 전체 공간 판정은 저온 EOS 중단으로 미완료다.

상세: [단계112 보고](../notes/REQUEST112_NATIVE_TWO_WAY_ATMOSPHERE_KO.md).


## 단계113 — 저온 EOS 병목을 고친 실제 결합 완주

분류: Counterexample candidate. 수정 입력으로 원128단계 결합 진화를 완주했다. 기준값 오류 수정 뒤 전체 구간의 명목 공간 차이는 대기 trace0.462530%, 외부 질량0.679781%로 원2% 안이다. 별도로 등록한 엄격한 재시작 중첩·독립 읽기 검사는 실패다. 이 결과를 최종 scalar 전하나 관측 완화 검출로 승격하지 않는다. 다음 원천 읽기는 보존 변수와 회복 오차를 사용한다.

상세: [단계113 보고](../notes/REQUEST113_NATIVE_COLD_COUPLING_KO.md).


## 단계114 — 실제 결합의 지연 전하와 광자 질량 분모

분류: Counterexample candidate. 실제 결합 원천의 지연 직접 성분2.2368e-27과 출사 광자의 질량 정규화4.1530e-27의 합은 선언 모형에서6.3897e-27로 양수다. 초기 출사 점유를 고정한 질량 성분만 뺀 값도5.7626e-27이다. 완전한 정적 비교 항성, nuisance 흡수, 전체 GR 전하 또는 궤도 관측 검출이 아니다. 내부 기계·계량·반경/주파수·각도 오차를 닫아야 한다.

상세: [단계114 보고](../notes/REQUEST114_NATIVE_COUPLED_CHARGE_KO.md).


## 단계115 — 내부 운동과 직접 전하 상쇄

분류: Counterexample candidate. 내부 질량 재배치를 지연 원천에 넣자 기존 직접 성분+2.2368e-27에−2.2569e-27이 더해져 약99percent 상쇄되고 직접 합은−2.0151e-29다. 기존 출사 광자 질량 정규화까지 더한 합은+4.1328e-27이다. 작은 잔여의 부호·크기는 내부 반경 오차가 닫히기 전 물리적 검출로 해석하지 않는다. 이전의 양수 성분 해석은 이 새 조건부 결과와 함께 읽는다.

상세: [단계115 보고](../notes/REQUEST115_NATIVE_INTERIOR_MOTION_KO.md).


## 단계116 — 내부 물질과 광자의 실제 양방향 진화

분류: Counterexample candidate. 실제 양방향64/128 경로에서 직접 전하 끝점은+1.14392e-27, 광자 질량 정규화까지 더한 합은+5.29688e-27다. 시간 차이는 각각0.943percent와0.207percent다. 단계115의 약99percent 상쇄는 새 모형에서 유지되지 않는다. 이동 산란 끝점 규칙도 동시에 고쳤으므로 차이 전부를 단일 되먹임 항의 효과로 식별하지 않는다. 최종 물리 전하·정적 EFT 이탈 검출은 아니다.

상세: [단계116 보고](../notes/REQUEST116_NATIVE_INTERIOR_FEEDBACK_KO.md).


## 단계117 — 실제 물질·광자 원천의 GR 기여

분류: Counterexample candidate. 영역 내 계량 응력/질량 기여는 각각−1.4640e-32,−1.3986e-34다. 임의 비음수 외향 외부 광자 분포와 scalar 퍼텐셜 모든 반복을 허용한 추가 GR 상계 뒤에도, 고정된 실제 원천의 관측 직접 변화 하한은+1.14170e-27다. 이 조건부 하한은 원천 공간/주파수/내부 각도·초기 GR 및 전체 재진화 오차를 포함하지 않는다.

상세: [단계117 보고](../notes/REQUEST117_NATIVE_FEEDBACK_GR_KO.md).


## 단계118 — 실제 공유 물질 경계와 인공 음향 확산 수정

분류: Counterexample candidate. 수정된 공유 물질 진화의 직접 전하 변화는1.01758836e-27, 실제 원천의 추가 GR·외부 광자/퍼텐셜 상계 후 조건부 직접 하한은1.01540309e-27다. 첫 upwind 시간 수렴 음수 결과는 압력 점프 질량 확산이 거의 전부를 설명하여 물리 결과로 채택하지 않았다. 새 양수도 원천 공간·초기 Einstein 제약·동적 GR·정적 비교 모형이 닫히지 않은 조건부 성분 결과다.

상세: [단계118 보고](../notes/REQUEST118_NATIVE_MATERIAL_JOIN_KO.md).


## 단계119 — 실제 native 경계층과 국소 공간 대조

분류: Counterexample candidate. 같은 관측 시계에서 국소16→19셀 차이는 직접 전하0.003351%, 직접+광자 질량 합0.133996%다. 새 실제 GR 원천·외부 광자/퍼텐셜 상계 후 조건부 직접 하한은1.01537118e-27로 양수다. 전체 공간 오차막대나 최종 항성 전하로 승격하지 않는다.

상세: [단계119 보고](../notes/REQUEST119_NATIVE_BOUNDARY_LAYER_KO.md).


## 단계120 — 실제 셀 재고와 초기 제약 연결

분류: Counterexample candidate. 초기 배경의 정적 질량·K 보정은 동적 관측 전하가 아니다. 새 초기 입력의 셀별 재고·native 대조는 통과했지만 새 실제 결합 궤적과 전하 판정은 없다. 단계119의 조건부 양수 하한을 새 배경에 승계하지 않는다. 작은 계량 상대 보정만으로 훨씬 작은 동적 신호의 영향이 작다고 추정하지 않는다.

상세: [단계120 보고](../notes/REQUEST120_INITIAL_CONSTRAINTS_KO.md).


## 단계121 — 수정 초기값의 실제 결합 진화

분류: Counterexample candidate. 새 배경의 직접 전하 변화는+1.02100041e-27, 지정 광자 질량 항과 합친 명목값은+5.16706005e-27다. 원 readout passed=True. 기존 양수 GR 하한은 승계하지 않는다. 최종 물리 전하·정적 EFT 이탈·관측 검출은 미확정이다.

상세: [단계121 보고](../notes/REQUEST121_PROJECTED_COUPLED_EVOLUTION_KO.md).


## 단계122 — 비등방 GR·스칼라 변화량 진화

분류: Counterexample candidate. 특성 경계를 분할한 새 직접 끝점은+1.02100007e-27, 표현한 비등방 GR 항 포함 끝점은+1.02098476e-27다. 독립 Jordan 반경 적분과1.810e-14로 일치한다. 원 호환 실패는 보존하며 최종 전하·정적 EFT 이탈·관측 검출은 미확정이다.

상세: [단계122 보고](../notes/REQUEST122_ANISOTROPIC_DYNAMIC_GR_KO.md).


## 단계123 — 동적 lapse와 실제 보상 광자 수송

분류: Counterexample candidate. 실제 출사 광자의 질량 감소와 반경 압력을 함께 lapse에 연결하고 같은 내부·대기 광자 수송의 기하 응답을 완주했다. 수송의 reference-frequency 에너지 변화는 최종 전하가 아니다. 추가 충돌·유체·외부 scalar 폐쇄 뒤에 다시 전하를 읽어야 한다.

상세: [단계123 보고](../notes/REQUEST123_DYNAMIC_LAPSE_TRANSPORT_KO.md).


## 단계124 — 광자·물질 에너지/H의 동시 응답

분류: Counterexample candidate. 단순 기하 광자 응답을 실제 물질 에너지/H 응답에 연결하고 지정 시간/배경 대조를 통과했다. 다음 GR용 추가 에너지·반경 압력·trace를 저장했지만 아직 적용하지 않았다. 작은 성분 응답을 전체 물리 전하의 상계나 관측 검출로 사용하지 않는다.

상세: [단계124 보고](../notes/REQUEST124_MONOLITHIC_COLLISION_RESPONSE_KO.md).


## 단계125 — 실제 보존 물질 응답

분류: Counterexample candidate. 추가 물질 이동의 실제 시간/배경 비교를 통과했다. 작은 질량 이동도 정지질량 에너지 source를 만들므로 전달 열만으로 생략하지 않는다. 새 원천의 GR 재적용과 최종 물리 전하는 아직 판정하지 않는다.

상세: [단계125 보고](../notes/REQUEST125_CONSERVED_MATERIAL_RESPONSE_KO.md).


## 단계126 — 실제 물질·광자·GR 귀환

분류: Counterexample candidate. 새 물질·광자 원천의 추가 compact 전하 성분은1.645727905e-44다. 이는 이전 표현 성분에 대한 추가 응답이며 최종 물리 전하·관측 검출·양의 하한 승계가 아니다.

상세: [단계126 보고](../notes/REQUEST126_MATTER_PHOTON_GR_RETURN_KO.md).


## 단계127 — 외부 광자와 무한원 scalar 전하

분류: Counterexample candidate. 외부 scalar 두 성분의 합은-5.300536877e-33이다. 지정 내부 원천에서 외부 각도 자유도와 전 차수 퍼텐셜 등을 포함한 조건부 scalar 하한은+1.019915831e-27로 남는다. 원천·EOS·결합 오차와 정적 비교 비흡수성을 인증하지 않았으므로 실제 관측 신호 검출이 아니다.

상세: [단계127 보고](../notes/REQUEST127_GLOBAL_SCALAR_CLOSURE_KO.md).


## 단계128 — 수락 단계의 보존 이력과 GR 전하

분류: Counterexample candidate. 실제 단계 원천을 적용한 내부 scalar 자유 성분은1.021135582147e-27이다. 새 외부·퍼텐셜·질량 오차와 조건부 삭제 물질 직접 원천 상계를 포함한 하한은1.020019540300e-27이다. 전체 EOS·시간 연속·결합 오차막대나 관측 검출은 아니며, 이전17시점의 명목 광자 질량 합을 새 결과로 승계하지 않는다.

상세: [단계128 보고](../notes/REQUEST128_STAGE_ENERGY_GR_SOURCE_KO.md).


## 단계129 — 보존 native 압력의 실제 전하 적용

분류: Counterexample candidate. native 보존 압력을 실제 trace/반경 응력에 넣은 내부 자유 scalar 끝점은+1.021133892155e-27이며 기존 대비-0.000165501percent다. 이 국소 구성 관계 안정성은 최종 물리 전하·관측 식별성이나 정적 비교 비흡수성의 증명이 아니다.

상세: [단계129 보고](../notes/REQUEST129_CONSERVATIVE_NATIVE_EOS_KO.md).


## 단계130 — native 열·반응 수정의 실제 결합 진화

분류: Counterexample candidate. 수정된 열·반응 입력의 실제64/128 경로를 GR 전하에 적용한 fine 자유 성분은+1.021135080724e-27이며 이전 상대 변화는-4.910446e-7이다. 지정 초기 연산자/원천과 외부·손실 조건에서 새 scalar 구간은[1.020019038924e-27,1.022251122525e-27]이다. 전체 시간·EOS·수송 오차나 관측 식별성까지 포함한 최종 물리 구간이 아니다.

상세: [단계130 보고](../notes/REQUEST130_NATIVE_THERMOCHEMISTRY_EVOLUTION_KO.md).


## 단계131 — 실제 대기 끝점의 native 보존 역산

분류: Counterexample candidate. 새 실제 fine 끝점223셀의 native D/S/K/H 역산과 동일 광자 충돌 판독을 통과했다. 압력 차이최대1.436e-8, 대기 끝점trace 변화1.272e-8, 순H교환 차이0.00766percent다. 원574상태 전체 예산 실패를 보존하고223셀로 줄여 같은180s 안에서 완료했다. 전 중간 이력·EOS 미분·새 retarded 전하는 인증하지 않는다.

분류: Conjectural. 같은 대기 보간을 이유로 궤적을 다시 돌리지 않는다. 다음은 수정 원천의 GR 귀환과 실제 물질·광자 결합/삭제 물질 효과이며, 전체 대기 이력과 작은 포트 공간 오차는 열린 조건이다. [단계131 보고서](../notes/REQUEST131_NATIVE_ATMOSPHERE_INVERSE_KO.md).


## 단계132 — 수정 원천 GR의 실제 광자·열·수소 귀환

분류: Counterexample candidate. 수정 EOS 궤적의 공간 GR/lapse를 실제531셀 광자 반경·각도·주파수 수송과 이동 충돌·물질 에너지/H 동시 방정식에 적용해 원 세 경로를 완주했다. 시간 차이최대0.112567%, 배경 경로 차이0.162445%, 에너지/H 잔차최대4.42e-12로 원 기준을 통과했다. 전체 각도 분포·물질 에너지/H·충돌 전달량·반경 포트를 저장했다. 입력 확인만이 아니라 실제 응답 적분의 결과다.

분류: Counterexample candidate. 비트 단위 재생은 실패로 남기며 최초1e-12 재생 기준은 통과했다. 원17시점 GR의 전체 이력 대조0.332625% 실패도 보존하고, 등록된33시점 대조로0.0427522%에 도달했다. 직렬 응답 예산 실패 후 같은 세 경로를3CPU·합6GiB 상한으로 명시적으로 재설계해207.61s/650s 안에 마쳤다. 물리 격자·시간 간격·기간·수락 기준은 늘리거나 완화하지 않았다.

분류: Conjectural. 새 에너지/H/운동량 전달을 추가 자유 유체 운동과 GR로 다시 반환하여 결합 잔차를 판정해야 한다. 전체 EOS 이력·삭제 물질 수송·작은 포트 공간 오차·비선형 GR·최종 전하·정적 비흡수성 및 관측 연결은 미완료다. 기존 양의 조건부 구간을 전체 귀환의 새 인증으로 승격하지 않는다. [단계132 보고서](../notes/REQUEST132_UPDATED_GR_RETURN_KO.md).


## 단계133 — 수정 전달량의 자유 유체·GR 귀환

분류: Counterexample candidate. 새 광자 충돌 전달량을 실제 공유 물질 유속에 적용하고 세 경로의 바리온·운동량·에너지/H를 전체 구간 진화했다. 직접 압력/trace와 광자 응력을 compact GR에 적용해 추가 전하1.64551e-44, 기존 성분 대비1.61145e-17을 얻었다. 물질/응력의 시간·배경 대조, 보존, pressure/directional 및 독립 GR 기준을 통과했다. 이전 누락된 내부 반경 운동 응력2T도 계량 일에 넣었다.

분류: Counterexample candidate. 앞선 광자 풀이와 자유 유체 반환 이력의 에너지/H 불일치는1.02237/0.439882다. 작은 compact 보정만으로 결합이 닫혔다거나 전체 누락 효과가 작다고 결론 내리지 않는다. 최종 전하·결합 고정점·외부/삭제 물질·전체 EOS 이력·비선형·관측은 미완료다.

분류: Conjectural. 다음은 새 물질 운동/재고/기계 수송을 이동 광자 충돌에 되돌려 동시 진화하고 잔차를 다시 판정하는 일이다. [단계133 보고서](../notes/REQUEST133_UPDATED_MATERIAL_GR_RETURN_KO.md).


## 단계134 — 물질·광자 재결합과 새 GR 반환

분류: Counterexample candidate. 실제 물질→광자/에너지/H→물질→compact GR 반환 후 공통17시점의 절대 에너지/H 불일치는 이전의0.0331454%/0.254817%로 줄었다. 추가 compact 전하는 기존 단계130 성분의1.61144e-17배지만 전체 누락 응답의 상계가 아니다. 이 불일치 감소나 양의 조건부 전하를 정적 모형 비흡수성·관측 검출로 승격하지 않는다. [단계134 보고서](../notes/REQUEST134_UPDATED_JOINT_FEEDBACK_KO.md)


## 단계135 — 새 GR의 실제 귀환과 물질 단계 경계 수정

분류: Counterexample candidate. 새 GR 보정의 실제 수송·물질·GR 반환 성분은1.35048e-64, 단계130 성분의1.32253e-37배다. 이 유한 구동의 작은 비율은 모든 누락 효과의 상계나 균일 수축성·관측 제약이 아니다. 이번 작은 보정의 에너지/H 불일치0.978809/0.020705와 기존 큰 이력의0.000334498/0.00160665는 서로 다른 정규화다. 기존 큰 반환의 구간 경계 수정이 다음 우선순위다. [단계135 보고서](../notes/REQUEST135_COMPENSATED_GR_RETURN_KO.md).


## 단계136 — 기존 큰 물질 반환의 수정과 GR 적용

분류: Counterexample candidate. SSP 종료 원천 구간과 영 가중치 probe 수정을 기존 큰 물질 응답에 적용하고 실제 압력·trace·compact GR까지 연결했다. 물질 압력 시간 차이는0.084372%에서0.001907%로 개선됐으나 광자와의 에너지/H 불일치는0.089750%/0.143735%로 남았다. 수정된 추가 전하1.64549994e-44와 변화-1.69971e-50는 유한 성분 결과이며 전체 오차 상계가 아니다. 분류: Conjectural. 다음은 수정된 큰 물질 운동을 실제 광자 동시 풀이에 반환하는 일이다. 고정점·전체 EOS/미분·외부/삭제 물질·비선형·관측·최종 전하 조건은 유지한다. 원 예산 중단과 재평가는 [단계136 보고서](../notes/REQUEST136_CORRECTED_LARGE_RETURN_KO.md)에 보존했다.


## 단계137 — 수정된 운동의 실제 광자·물질·GR 반환

분류: Counterexample candidate. 실제 수정 운동→광자→물질→GR 반환 후 에너지/H 불일치가약117배/153배 감소했고, 절대 잔차는 시간 간격 대조 차이의0.389652/0.0105163배다. 이번 추가 compact 성분1.64550663e-44와 이전 대비 변화6.68498e-50는 유한 응답이며 관측 검출·정적 비흡수성·전체 오차 상계가 아니다.

분류: Conjectural. 다음은 저장된 잔차를 같은 보존 압력/trace와 retarded GR 전하에 연결하고 증폭을 포함한 오차 한계를 구성하는 일이다. 작은 waveform 반복을 자동 추가하지 않는다. 상세 수치·원 예측 실패·예산·조건은 [단계137 보고서](../notes/REQUEST137_CORRECTED_MOTION_FEEDBACK_KO.md)에 보존했다.


## 단계138 — 실제 잔차의 전하 연결과 조건부 상계

분류: Proven. 선언한 저장 E/H 원천 집합의 compact 전하 상계는1.446615046085e-52이며, 최종 E/H 원천 오차가A배 집합 안에 있으면 영향은A배 상계 이하다. 분류: Counterexample candidate. 이 크기는 기존 조건부 하한의1.418224e-25배지만, 전체 결합 증폭A와 다른 원천 오차는 인증하지 않았다. 실제 잔차 끝점은 시간 경로에 따라 부호가 달라져 음의 물리 전하나 관측 신호로 부르지 않는다.

분류: Conjectural. 같은 작은 E/H waveform의 반복보다 주 전하의 EOS·시간·경계 원천 오차와 실제 결합 증폭을 제한하는 일이 다음 우선순위다. 전체 물리·관측 완료 요구사항을 유지한다. 상세 전제·증명·실패·예산은 [단계138 보고서](../notes/REQUEST138_RESIDUAL_CHARGE_ENVELOPE_KO.md)에 둔다.


## 단계139 — 수정 EOS 실제 저장 이력의 전하 적용

분류: Counterexample candidate. 실제 수정 궤적17시점의 심부323개·대기3,017개 보존 상태를 native 역산하고 초기 offset을 유지한 압력·반경 응력·trace 교정을 같은 GR에 적용했다. 주 자유 compact 끝점은1.021135080724e-27에서1.021134154346e-27로 변했고 보정은 기존 신호의9.072049e-7배였다. 역산·9/17시점 원천·4/8차 구적·합산 원천의 직접 적용과 독립 raw 검산을 통과했다. 새 물리 궤적은 계산하지 않았다.

분류: Conjectural. 이 결과는17개 저장 시점의 제한된 원천 교체다. 전체129상태·연속 EOS/미분·초기 GR 계수·결합 증폭·외부/삭제 물질·비선형·정적 비흡수성/관측 및 최종 전하를 인증하지 않는다. 다음은 저장 native 상태와 실제 광자장의 충돌/반응률 이력 연결이다. 원 실행예측 부적격·정확 호출 캐시·원 예산 내 검산 재배분 및 수치는 [단계139 보고서](../notes/REQUEST139_CONSERVED_EOS_HISTORY_KO.md)에 보존한다.


## 단계140 — native 충돌의 실제 유한시간 응답과 물질 반환 경계

분류: Counterexample candidate. 저장17시점의3,340개 native 보존 상태를 실제 충돌 입력으로 연결하고,531셀·3.434431ms의 광자/열/H 결합64/128 응답을 완주했다. 최대 시간 대조 차이는0.1737491%, 최대 에너지 잔차는1.3961e-10이다. 부호 있는 산란 차분, 심부 H 차분 괄호 정정과 에너지/운동량을 보존하는 명시적 광자 수 반올림 보정을 사용했다. 원 실패를 보존하며, j 미복원이 속도 오류라는 진단은 실제 Pi 기반 호출을 확인하여 철회했다.

분류: Counterexample candidate. 자유 물질의 첫 두 단계에서 심부11번 셀 H 상대 변화1.465821457e-6이 원 선형 기준1e-6을 넘었다. 실제 진폭 비선형 대체식도 시험했지만 유량 차분 해상도 지표0.002425358383이0.002 기준을 넘었다. 물질 전체 생산과 새 GR 전하는 실행하지 않았다. 원 양의 전하가 이번 충돌 교정까지 포함해 유지되는지는 미판정이다.

분류: Conjectural. 다음은 같은 실제 입력·진폭을 유지한 물질 유량/primitive 압력의 보상 차분이다. 이를 통과한 뒤 같은64/128 물질 경로와GR 반환으로 간다. 전체 EOS/미분·결합 증폭·외부/floor·비선형·정적 비흡수성/관측 및 최종 완료 기준을 축소하지 않는다. 상세 수치·수정·실패·예산은 [단계140 보고서](../notes/REQUEST140_NATIVE_COLLISION_RESPONSE_KO.md)에 보존한다.


## 단계141 — 보상 유량의 유한 물질 반환과 실제 GR 전하

분류: Proven. 고정 면 유량B와 그 발산을 함께 제거하면 보존율-div(F)+G가 그대로다. 분류: Counterexample candidate. 이 항등식을 심부 배경 압력과 공유 면에 적용해 단계140의 유량 차분 병목을 해결했다. 실제 진폭의 물질64/128 전체 경로와 압력/응력·GR 대조를 통과했다. native 충돌/물질 교정의 전하 끝점은+1.380775607e-34이며, 기존 native 압력 교정까지 합산한 원천을 직접 적용한 끝점은+1.021134292423e-27이다. 지정 계산에서 양의 전하가 유지되며 원 선형 small-state 실패는 보존한다.

분류: Counterexample candidate. 이 결과는 저장 초기 기하의free compact 첫 변분과 보간 EOS의 유한 물질 반환이다. 유량/압력16eps 지표를 엄밀한 EOS 인증으로 세지 않는다. 분류: Conjectural. 다음은 같은 native 입력을 유지한 새 물질 운동의 실제 광자 반환 및 결합 오차 판정이다. 전체 EOS/미분·퍼텐셜/외부/floor·비선형·정적 비교/관측과 최종 완료 요구사항을 유지한다. 수치·실패·예산은 [단계141 보고서](../notes/REQUEST141_COMPENSATED_FINITE_RETURN_KO.md)에 보존한다.


## 단계142 — native 유한 물질의 실제 광자·물질·GR 귀환

분류: Counterexample candidate. 단계141의 실제 바리온·운동량·재고·비충돌 수송을 광자/열/H 동시 방정식에 되돌리고, 동일 native 충돌 입력을 한 번만 유지한 새 전달을 유한 물질·압력·GR까지 적용했다. 원64/128 경로를 완주했고 합산 자유 compact 끝점은+1.021134292423207e-27이다. fine 에너지/H 잔차는4.7530e-12/7.7409e-12로, 절대 잔차가 시간 대조 차이의1.6981e-8/4.4552e-9배다. 원천·보존·출사·시간·독립 GR와 합산 직접 적용을 통과했다.

분류: Counterexample candidate. 유한 물질을 반환했지만 광자 primitive/충돌 지도는 여전히 선형이다. 끝점의 대기236셀 H 상대 변화4.15355e-4를 포함한 EOS/충돌 나머지를 인증하지 않는다. 전하 성분의 전후 차이-2.33748e-41은 시간 대조 규모보다 작아 물리적 음의 피드백 검출로 해석하지 않는다. 분류: Conjectural. 같은 파형의 추가 반복 대신 실제 EOS/미분·결합 증폭·외부/floor 오차를 전하 불확실성에 연결한다. 전체 퍼텐셜·비선형·정적 비교/관측과 최종 완료 조건은 유지한다. [단계142 보고서](../notes/REQUEST142_NATIVE_FINITE_MOTION_FEEDBACK_KO.md).


## 단계143 — 유한 충돌 잔여의 부분 계산과 자원 중단

분류: Counterexample candidate. 수락한 접선 연산자는 유지하고 유한 보존 역산·충돌 차분만 확장 정밀도로 계산했다. 저장17시점 중9시점에서 잔여/선형 충돌항의 에너지 가중 L1 최대1.87742e-7, 절반 진폭 잔여비0.25606–0.27035를 얻었다. 분류: Proven. 잔여 정의와 두 주파수 광자 수 보정의 에너지/동일 각도 운동량 보존 항등식을 symbolic 검산했다.

분류: Counterexample candidate. 원 binary64 실패와 예산 부적격을 보존한다. I/O 중복 읽기를 줄였으나 남은8시점의 실측 예측이 재배분 잔여 예산을 넘어 중단했다. 새 광자·물질·GR 생산은 미실행이며, 단계142의 양의 끝점을 이번 교정을 포함한 값으로 인용하지 않는다. source_complete, finite_remainder_propagated_to_GR, uniform_EOS_derivative_bound, uniform_nonlinear_remainder_bound, final_charge_solved, full_goal_complete는false다. 분류: Conjectural. 다음은 초기화/저장 비용을 포함해 남은 원천의 실행 가능성을 재평가하고 실제 전파까지 연결하는 일이다. 전체 완료 조건을 축소하지 않는다. [단계143 보고서](../notes/REQUEST143_FINITE_COLLISION_REMAINDER_KO.md).


## 단계144 — 유한 충돌 교정의 실제 광자·물질·GR 적용

분류: Counterexample candidate. 저장9시점을 재사용해17시점의 유한 충돌 잔여를 완성하고, 같은531셀·3.434431ms의64/128 광자 증분을 적분했다. 저장된 원 응답과 합친 총 전달량으로 유한 물질을 새로 풀어 압력·응력·GR까지 적용했다. 합산 질량 정규화 free compact 첫 변분 끝점은+1.0211342924230259e-27로 양수다. 광자 시간 대조 최대0.0842258%, 물질 보존 최대9.89709e-16, GR 시간 대조5.35748e-4와 독립 직접 적용을 통과했다. 분류: Proven. 같은 고정 선형 연산자의 응답 합산 및 정의한 잔여의 복원 항등식을 symbolic 검산했다.

분류: Counterexample candidate. 단계142 대비 변화는 시간 대조 규모의0.26835%로, 작은 음의 효과가 검출됐다는 결론은 아니다. 원 binary64·작은 상태·예산 실패를 보존하고 미사용 원천 예산만 옮겨 성공한 접두부부터 재개했다. finite_collision_remainder_applied_to_photons_material_and_GR와residual_below_temporal_comparison_scale는true다. uniform_EOS_derivative_bound, uniform_nonlinear_remainder_bound, coupled_fixed_point_verified, exterior_floor_feedback_closed, full_nonlinear_GR, final_charge_solved, full_goal_complete는false다. 분류: Conjectural. 다음은 합산 원천의 EOS/미분·결합·외부/floor/퍼텐셜 오차를 주 전하의 정확도 요구에 연결하는 일이다. 추가 작은 파형 반복을 자동 실행하지 않으며 전체 비선형·정적 비교/관측의 완료 조건을 유지한다. [단계144 보고서](../notes/REQUEST144_FINITE_COLLISION_GR_RETURN_KO.md).


## 단계145 — 최신 출사 이력의 무한원 전하와 질량 정규화

분류: Counterexample candidate. 수정 EOS의 실제 단계 광도와 단계144의 전체 native 광자 증분을 무한원 외부 응력·대응 질량 감소에 적용했다. 도착 광자의 질량 분모까지 포함한 명목 변화는+5.238124070132e-27이며, 그80.5058%가 질량 정규화 항이다. 원0.2% 적분·2% 시간 기준과 단계 포트·독립 원시함수·고정밀 정규화 검산을 통과했다. 분류: Proven. δα=(s+α₀ε)/(1−ε)는 주어진 분자·분모의 대수 항등식이다.

분류: Counterexample candidate. 모든 셀의 새 원천 노름·부호 있는 방출·정확 유리수 포트 표현으로 퍼텐셜/지정 질량 제약 오차를 합친 조건부 상한은3.75271e-32다. 이는 전체 물리 오차 구간이 아니며, 외부 시간 대조 약1.3%·EOS/미분·연속 결합·삭제 물질·비선형 오차를 포함하지 않는다. 질량 분모 변화를 동적 새 효과로 세지 않는다. actual_current_emission_at_null_infinity와actual_emission_mass_normalization는true, uniform_EOS_derivative_bound, full_source_error_enclosed, coupled_fixed_point_verified, nonlinear_GR, final_charge_solved, full_goal_complete는false다. 분류: Conjectural. 다음은 정규화 항을 분리한 scalar 잔여에 전체 구성/결합 오차와 동일 재고 정적 비교·관측을 연결하는 일이다. 전체 완료 조건과 기존 계산 예산 규율을 유지한다. [단계145 보고서](../notes/REQUEST145_NATIVE_CHARGE_NULL_INFINITY_KO.md).


## 단계146 — 삭제 물질의 수동 복사 오차

분류: Counterexample candidate. 수동 삭제 물질의 복사 효과 상계는 고정 질량 scalar에서9.57728e-32, 질량 분모까지 포함하면1.03285e-29다. 각각 기존 scalar 잔여의0.00938%,1.01147%다. 이 항만 적용한 명목값은 양수이나 전체 오차 구간이나 새로운 동적 관측 신호가 아니다. [단계146 보고서](../notes/REQUEST146_DISCARDED_RADIATION_KO.md).


## 단계147 — 저밀도 native 물질의 실제 결합 적용

분류: Counterexample candidate. coarse 경로의 삭제 바리온은236640.54g에서593.42g으로 줄었고 마지막 상태에3502.10g의 저밀도 물질이 양의 압력·운동과 함께 남았다. 이는 실제 원천의 변화이며 전하 오차가 같은 비율로 감소했다는 결론은 아니다. 새 GR 전하 값은 없으며 이전 배경의 양의 명목 전하를 이번 결과까지 반영한 값으로 인용하지 않는다. [단계147 보고서](../notes/REQUEST147_RETAINED_NATIVE_MATERIAL_KO.md).


## 단계148 — 저밀도 물질의 실제 전하 연결

분류: Counterexample candidate. 새 원천의 compact 전하 변화는+1.0211271551555008e-27, 무한원 질량 정규화 변화는+5.238117415710266e-27이다. 약80.51%는 질량 분모 항이다. 동일 원 thermal 배경과 비교한 저밀도 물질 교정은-7.92557e-33이나 그 자체의4.21767% 시간 대조는원2% 기준에 미달한다. 전체 전하 시간 대조 통과와 미세 교정 분해능을 구별하며 새 관측 신호로 세지 않는다. [단계148 보고서](../notes/REQUEST148_RETAINED_FINE_GR_KO.md).


## 단계149 — 현재 native EOS 원천의 실제 전하 적용

분류: Counterexample candidate. 현재 native 원천을 적용한 compact 변화는 1.0211262287752335e-27, 무한원 정규화 변화는 5.2381164893299992e-27다. 전 이력의 EOS 압력 교정은 주 scalar 잔여의 최대 9.07213428e-07이며 원 물질성 기준을 통과했다. 관측 검출이나 전체 오차 구간으로 부르지 않는다. [단계149 보고서](../notes/REQUEST149_RETAINED_NATIVE_SOURCE_KO.md).


## 단계150 — 실제 native 광자·물질·GR 반환

분류: Counterexample candidate. 실제 native 힘·충돌을 적용한 compact 교정은 8.26273103e-34, 실제 패킷의 무한원 정규화 교정은 3.40613247e-34다. 최종 명목 값 5.2381168299432462e-27는 양수이며 모든 등록된 최종 비교를 통과했다. 관측 식별성이나 전체 오차 인증으로 확대하지 않는다. [단계150 보고서](../notes/REQUEST150_RETAINED_NATIVE_RETURN_KO.md).


## 단계151 — 실제 물질의 광자·GR 재귀환

분류: Counterexample candidate. 물질 재귀환의 무한원 전하 증분은 3.2937902e-40, 새 명목 전하는 5.2381168299435755e-27다. 작은 교정 자체의 시간·구적 비교를 별도로 확인했다. 비영 값은 정적 nuisance를 제거한 관측 식별성의 증거가 아니다. [단계151 보고서](../notes/REQUEST151_RETAINED_MOTION_RETURN_KO.md).


## 단계152 — 현재 계량의 광자·물질 재적용

분류: Counterexample candidate. 현재 계량 반환의 정규화 증분은 1.64550691e-44, 새 명목 끝점은 5.2381168299435755e-27다. 증분 자체의 시간·구적 기준을 확인했다. 비영 전하를 정적 nuisance에서 분리된 관측 검출로 해석하지 않는다. [단계152 보고서](../notes/REQUEST152_RETAINED_METRIC_RETURN_KO.md).


## 단계153 — native 음향 미분의 실제 결합 적용

분류: Counterexample candidate. 이번 HLL 음향 보정의 정규화 증분은-1.34141139e-40, 새 명목 끝점은5.2381168299434413e-27다. 지정 수치 응답의 결과이며 같은 재고의 정적 비교에서 분리된 관측 신호가 아니다. [단계153 보고서](../notes/REQUEST153_RETAINED_NATIVE_ACOUSTIC_KO.md).


## 단계154 — 전 반경 순간 정적 비교

분류: Counterexample candidate. 같은 심부 공급 에너지 총량의 두 정적 분포가 서로 다른 전하를 주었다. 이 감도 예를 인과적 진화나 엄밀한 오차 범위로 읽지 않는다. 현재 자발적 완화의 명목 전하와 이번 불완전 정적 성분을 빼서 관측 잔여를 만들지 않는다. [단계154 보고서](../notes/REQUEST154_WHOLE_STAR_STATIC_RESPONSE_KO.md).


## 단계155 — 외부 스칼라 입력의 실제 결합 응답

분류: Counterexample candidate. 외부 입력의 1차 차등 응답을 실제 계산했으나 compact 물질 매개 반환에 한정한다. 직접 입사파의 전체 무한원 산란·외부 광자·같은 재고 정적 비교와 관측 비흡수성은 미완료다. 기존 무구동 완화의 명목 양의 전하와 구별한다. [단계155 보고서](../notes/REQUEST155_EXTERNAL_SCALAR_COUPLED_RESPONSE_KO.md).


## 단계156 — 물질·광자 상호 반환

분류: Counterexample candidate. 물질→광자→물질 상호 반환을 수락한 같은 이력에서 compact GR 응답을 판독했다. 직접 외부 입사파/광자의 완전한 무한원 판독과 정적 비교·관측 비흡수성은 아직 필요하다. [단계156 보고서](../notes/REQUEST156_RECIPROCAL_INCIDENT_RESPONSE_KO.md).


## 단계157 — 실제 자체 GR 반환

분류: Counterexample candidate. 실제 자체 GR 반환 뒤 compact 전하 보정을 분리 판독했다. 전체 외부 입사파/광자의 무한원 기여와 동일 입력의 정적 비교·관측 식별성은 남아 있다. 분류: Counterexample candidate. 동일한 정밀 GR128 입력에 조건부인 판정이며, 서로 다른 상위 GR 입력의 초기 비교 실패는 보존한다. 상위 입력 불확실성은 닫히지 않았다. [단계157 보고서](../notes/REQUEST157_INCIDENT_SELF_GR_RETURN_KO.md).


## Phase158 — 실제 입사 응답의 조건부 무한원 판독

분류: Counterexample candidate. 실제 단계156/157 결합 이력과 signed SDIRK 각도별 방출을 초기 GR 계수·광선의 무한원 판독에 연결했다. 같은 fine retained 질량 정규화 아래 직접 산란, 물질 매개 응답, 자체 GR 보정을 분리했다. 물질 끝점은 -2.53643242961e-51로 유지됐고 외부 광자·질량 연결의 상대 변화는 2.7729987e-9다. 시간 간격 대조는 0.001674247%다. 새 유체·EOS 적분은 없다.

분류: Counterexample candidate. 직접 산란 끝점은 +2.82690789279e-35다. 직접장의4/8차 차이는 물질 신호의 약5.48e10배, 유한 higher-return 추정 척도는 약7.43e8배다. 두 척도는 물리 오차 상계가 아니다. 같은 입력의 직접 성분을 공통으로 제거하는 bookkeeping은 동일 재고 정적 EFT 또는 관측 nuisance 비교가 아니다. 신호가 관측 가능하다고 판정하지 않는다.

근거: [단계158 보고](../notes/REQUEST158_INCIDENT_NULL_INFINITY_KO.md), `outputs/direct-eos-gr33/native-incident-infinity/final-result.json`. 원 출력·입력 SHA·기호 및 성분별 검사는 단계158 manifest가 소유한다.


## Phase159 — 선택 외부 항의 크기와 관측 판정 경계

분류: Counterexample candidate. 단계158 저장 물질 전하 −2.53643242961482e−51를 고정 기준으로 둘 때, 직접 방출 복사 계량에 의한 입사 primary·유한 첫 Born 장의 반응, 이 원천의 추가 정적 potential 반복 및 광자 도달 질량 기여의 합은 1.8115806123377687e−53 이하로 제한됐다. 조건부 구간은 [−2.554548235738198e−51, −2.518316623491442e−51]다. 원 점 추정의 수렴 실패는 그대로다.

분류: Counterexample candidate. 이는 열거한 항만으로 저장된 음의 물질 전하를 지우지 못한다는 결과다. 기준 물질 전하의 오차, 실제 공급 장의 보간 오차, 전체 직접장의 정확도, 내부 연산자·혼합 제약·변한 광선의 scalar 판독은 이 구간에 들어 있지 않다. 동일 재고 정적 비교·미분 nuisance 제거·궤도 및 관측 식별성은 여전히 미완료다.

근거: [단계159 보고](../notes/REQUEST159_EXTERIOR_INCIDENT_COUPLING_KO.md), `outputs/direct-eos-gr33/native-exterior-incident/final-result.json`.


## Phase160 — 내부 배경 연산자 반환과 실제 경계 접합

분류: Counterexample candidate. 선택한 내부 primary 반환은 −2.8122596136721473e−55로 저장 물질 전하의 0.01109%이며 같은 부호다. 그 내부 부문의 유한 Born·추가 정적 potential 상계는 4.514759127688983e−54, 즉 저장 물질 전하의 0.1780% 이하다. 기준 전하와 primary 점 추정을 고정한 조건부 구간은 음수다.

분류: Counterexample candidate. 기준 전하·원천 표현의 오차, 전체 외부·상호 혼합 GR·새 파동의 물질 반환은 이 구간 밖이다. 같은 재고 정적 비교, 직접장 제거 정확도와 관측 비흡수를 증명하지 않았다. 단계159 원 점 추정 실패도 유지한다.

근거: [단계160 보고](../notes/REQUEST160_INTERIOR_OPERATOR_RETURN_KO.md), `outputs/direct-eos-gr33/native-interior-incident/final-result.json`.


## 단계161 — 상호 혼합 GR 원천 적용

분류: Counterexample candidate. 내부 상호 혼합 GR 원천을 기존 outgoing 성분 판독에 적용했다. 보정은 −8.784791844740518e-55로 기존 물질 전하의 0.0346344%이며 같은 부호다. 단계160까지의 성분과 합하면 −2.5375921347606613e-51이다. 이 값은 외부 mixed ADM/도달 질량·물질 반환·전체 물리 오차가 해결된 관측 신호가 아니다.

분류: Conjectural. 동일 재고 정적·미분 nuisance 비교와 실제 관측 식별성 판정은 여전히 필요하다.

실행·방정식·범위: [단계161 기록](../notes/REQUEST161_RECIPROCAL_MIXED_GR_KO.md).


## 단계162 — 혼합 GR 장의 실제 수송 반환

분류: Counterexample candidate. 실제 혼합 GR 수송 반환을 signed 방출·도달 질량과 함께 고정 연산자의 무한원 전하로 읽었다. 추가 반환은 -2.0088782264382117e-79, 수정된 혼합 GR 성분 대비 2.285603836506345e-25이며 별도 열로 보존했다. 기존 선택 성분의 음의 부호가 유지된다. 전체 외부장 분해 정확도와 동일 재고 정적 비흡수성·관측 식별성은 미완료다. [단계162 보고서](../notes/REQUEST162_MIXED_GR_TRANSPORT_KO.md).


## 단계163 — 질량항의 공통 전하 정규화

분류: Counterexample candidate. 질량 정규화 보정+3.4829714644393073e-57을 실제 선택 전하에 반영해 끝점-2.5375890992998707e-51을 얻었다. 변화는 기존 크기의 약1.37255ppm이며 음의 부호가 유지된다. 전체 물리 전하의 오차 인증이나 정적/관측 비흡수성은 아니다. [단계163 보고서](../notes/REQUEST163_MATCHED_MASS_READOUT_KO.md).


## 단계164 — 외부 광자 혼합 원천의 전하 상계

분류: Counterexample candidate. 새 외부 광자 혼합 부문의 전하 상계1.36903e-57을 기존 외부 상계와 합하면 선택 끝점-2.53758910e-51의0.7139523% 이하다. 조건부 구간[-2.55570627e-51,-2.51947192e-51]은 음수다. 전체 물리 오차 구간이나 정적/관측 비흡수성 판정은 아니다. [단계164 보고서](../notes/REQUEST164_EXTERIOR_PHOTON_MIXED_BOUND_KO.md).


## 단계165 — 생성 scalar의 외부 혼합 질량과 에너지 흐름

분류: Counterexample candidate. 생성 scalar의 새 선택 상계2.28147e-57과 이전 외부 상계의 합은 선택 전하의0.7140422% 이하이며 조건부 구간[-2.55570856e-51,-2.51946964e-51]은 음수다. 새 내부 질량 접합 보정과 전체 물리/관측 오차는 이 구간 밖이다. [단계165 보고서](../notes/REQUEST165_GENERATED_SCALAR_EXTERIOR_KO.md).


## 단계166 — 질량 장부와 이산 적색편이 일

분류: Counterexample candidate. 실제 입력에서 누락된 공간 일0.3853807705erg는 기존 paired 질량 잔여 절댓값의1.253845배다. 앞선 선택 전하 부호 구간은 이 수송 수정의 결합 응답을 포함하지 않는다. 실제 전하를 다시 읽기 전 갱신하지 않는다. [단계166 보고서](../notes/REQUEST166_MASS_ENERGY_RECONCILIATION_KO.md).


## 단계167 — 실제 적색편이 수정 결합 예비 계산

분류: Counterexample candidate. 면 적색편이 수정 원천을 동시 광자·열·수소 응답에 실제 적용했으나 원2% 시간 게이트가 실패했다. 새 GR 전하는 아직 없으며 기존 선택 전하 부호를 전체 수정 결과로 확대하지 않는다. [단계167 보고서](../notes/REQUEST167_CONSERVATIVE_REDSHIFT_COUPLED_PREFIX_KO.md).


## 단계168 — 동시 Radau 적분의 실제 예비 대조

분류: Counterexample candidate. 같은 수정 원천의 실제 두 단계 동시 Radau 계산에서 물질 압력 대조는0.885%로 통과했지만 물질 에너지·수소·충격량은3.8–4.3%로 실패했다. 새 GR 전하를 수락하지 않는다. [단계168 보고서](../notes/REQUEST168_RADAU_TRANSFER_PREFIX_KO.md).


## 단계169 — 충돌 원천 적분의 실제 예비 대조

분류: Counterexample candidate. 물질 에너지·수소·충격량의 시간 차이는 원2% 기준에 들어왔으나 압력 차이는6.2624%로 실패했다. 수정된 동일 결합 해의 최종 전하 결론은 판정 불가다. [단계169 보고서](../notes/REQUEST169_COLLISION_SOURCE_PREFIX_KO.md).


## 단계170 — 수송 원천 적분 실패와 강한 감쇠 극한

분류: Counterexample candidate. 수송·충돌·출구를 함께 적분한 실제 대조는 여섯 시간 기준을 모두 실패했다. 수정된 동일 결합 해에서 최종 전하 결론은 아직 판정 불가다. [단계170 보고서](../notes/REQUEST170_KNOWN_FORCING_STIFF_LIMIT_KO.md).


## 단계171 — 실제 단계 원천 직접 적용

분류: Counterexample candidate. 실제 직접 원천 Radau의 다섯 시간 성분은 통과했으나 물질 압력11.4298%가 실패했다. 수정된 동일 결합 해에서 최종 전하 결론은 판정 불가다. [단계171 보고서](../notes/REQUEST171_DIRECT_SOURCE_RADAU_KO.md).


## 단계172 — 같은 결합 해의 압력 도달 구간 해상

분류: Counterexample candidate. 같은 실제 결합 해의 압력 시간 차이가11.4298%에서1.6284%로 줄어 여섯 성분 모두 원2% 기준을 통과했다. 수락 범위는 예비 구간이며 최종 전하 결론은 미판정이다. [단계172 보고서](../notes/REQUEST172_PRESSURE_ONSET_RESOLUTION_KO.md).


## 단계173 — 수정 결합 해의 전체 기간 계속

분류: Counterexample candidate. 전체 기간 광자·열·수소 시간 대조와 출구 이력이 원 기준을 통과했다. 자유 물질·GR 반환이 아직 없으므로 최종 전하의 기존 결론 유지 여부는 미판정이다. [단계173 보고](../notes/REQUEST173_COUPLED_CONTINUATION_PLAN_KO.md).


## 단계174 — 동일 수정 해의 실제 자유 물질 반환

분류: Counterexample candidate. 수정된 광자와 전체 자유 물질 해를 연결했지만 상호 반환과 GR 판독은 미완료다. 최종 전하의 기존 부호·크기가 유지되는지는 여전히 미판정이다. [단계174 보고](../notes/REQUEST174_SAME_SOLUTION_MATERIAL_PLAN_KO.md).


## 단계175 — 실제 상호 반환 완료와 수소 입력 실패

분류: Counterexample candidate. 단계175의 실제 상호 반환 뒤 M_H 입력 기준이 실패했으므로 예정한176 GR·최종 전하 판독은 실행하지 않았다. 작은 광자/물질 E/H 차이나 단계 원천 잔차를 전체 상호 결합 수락으로 대체하지 않는다. [실행 및 원 판정](../notes/REQUEST175_ACTUAL_MATERIAL_RECIPROCITY_KO.md), [미실행 GR 판독](../notes/REQUEST176_SAME_SOLUTION_CHARGE_READOUT_KO.md).


## 단계177 — 단계 내 H 수송과 실제 반환 실패

분류: Counterexample candidate. 최종 전하 결론의 유지 여부는 미판정이다. 단계 방정식의 수치 수락을 최종 관측량 성과로 대체하지 않는다. 실제 자유 물질 반환의 수송 불일치 때문에 전체 기간·GR 판독은 실행하지 않았다. [실행과 판정](../notes/REQUEST177_NATIVE_NEUTRAL_COUPLING_KO.md).


## 단계178 — 실제 충돌 이력 반환과 남은 물질 실패

분류: Counterexample candidate. 최종 전하 유지 여부는 미판정이다. 실제 단계 충돌을 같은 예비 해에 적용한 수소 수송 개선은 확인했지만 전체 물질·광자 상호 일치의 수락으로 대체하지 않는다. [실제 실행과 남은 판정](../notes/REQUEST178_ACTUAL_STAGE_COLLISION_RETURN_KO.md).


## 단계179 — 동일 광자·물질 단계의 실제 판정

분류: Counterexample candidate. 최종 전하 결론은 미판정이다. 실제6개 하위 단계에서 물질과 광자의 동시 식이 통과한 것을 전하 결론의 유지나 새로운 관측량 검출로 대신하지 않는다. [실행과 판정](../notes/REQUEST179_JOINT_NATIVE_FLUID_RADAU_KO.md).


## 단계180 — 실제 동시 해의 원 시간 대조

분류: Counterexample candidate. 최종 전하는 미판정이다. 바리온 시간 실패를 유지하고 전체 기간을 확대하지 않았다. 보정 응답을 과거 선택 전하에 더하는 방식은 원 물질 연산자의 분기에서 중첩이 보장되지 않아 최종 판정 경로가 아니다. [실행과 판정](../notes/REQUEST180_JOINT_FLUID_TIME_KO.md).


## 단계181 — 전체 입사 구동의 동일 단계 연결

분류: Counterexample candidate. 최종 전하의 기존 결론은 미판정이다. 전체 선언 입사 구동의 같은 해가 짧은 단계식·보존 기준을 통과했다. 시간 수렴·자기 GR·동일 해의 최종 전하가 없어 부호나 정적 모형 탈출을 주장하지 않는다. 별도 보정 해 또는 진단값의 사후 합산도 하지 않는다. [실행과 판정](../notes/REQUEST181_FULL_INCIDENT_JOINT_KO.md).


## 단계182 — 전체 구동의 원 시간 대조

분류: Counterexample candidate. 최종 전하 결론은 미판정이다. 원 시간 차이는 B5.89185%,물질 에너지3.67449%,물질 압력3.81585%로2% 기준을 실패했다. 전체 기간·GR로 확대하지 않았으며 기존 부호나 정적 모형 탈출을 계승하지 않는다. [실행과 판정](../notes/REQUEST182_FULL_INCIDENT_TIME_KO.md).


## 단계183 — 실제 물질 전선의 시간 분해

분류: Counterexample candidate. 최종 전하의 결론은 미판정이다. 동일한 전체 입사 구동 결합 해에 물질 전선의 시간 수정을 적용해 원 시간 실패를 해소했다. B차이는5.89185%에서0.672533%로 줄었다. 이 결과는 지정 구간의 실제 병목 수리이며 최종 전하·정적 모형 탈출의 결론을 대신하지 않는다. [실행과 판정](../notes/REQUEST183_NATIVE_MATERIAL_FRONT_KO.md).


## 단계184 — 동일 해를 보존한 실제 연속 적분

분류: Counterexample candidate. 최종 전하 결론은 미판정이다. 동일 해를 실제로 이어가 원 시간 기준을 통과했지만 이는 전체 기간·자체 GR·전하 판독의 중간 성과다. 이전 전하의 부호나 오차 한계를 계승하지 않는다. [실행과 판정](../notes/REQUEST184_SAME_SOLUTION_CONTINUATION_KO.md).


## 단계185 — 동일 해의 전체 기간 계산 등록

분류: Conjectural. 최종 전하 결론은 미판정이다. 이번 전체 기간의 물질·광자 응답은 그 해의 생성 GR과 질량·전하 판독을 위한 필수 입력이다. 전체 적분 통과를 최종 관측량 판정으로 바꾸지 않는다. [고정 계획](../notes/REQUEST185_FULL_INCIDENT_HORIZON_KO.md).


## 단계186 — 동일 결합 해에서 생성 GR로 연결

분류: Counterexample candidate. 최종 전하는 미판정이다. 수락된0.4293ms의 동일 해에서 실제 compact GR 응답을 계산했지만 입사파 전체의 무한대 전하나 자기 결합 완료로 해석하지 않는다. [근거](../notes/REQUEST186_SAME_JOINT_SOLUTION_GR_KO.md).


## 단계187 — 같은 해의 생성 GR 반환 입력

분류: Counterexample candidate. 최종 전하는 여전히 미판정이다. 작은 반환 입력의 단계식 적용은 통과했지만 그 원천을 포함한 실제 결합 진화와 동일 해 전하 판독이 남는다. 입력 크기만으로 최종 전하 상계를 선언하지 않는다. [근거](../notes/REQUEST187_SAME_SOLUTION_GR_RETURN_KO.md).


## 단계188 — 같은 단계 해의 GR 반환 시간 적분

분류: Counterexample candidate. 최종 전하는 미판정이다. 첫T/64의 실제 GR 반환 증분과 자체 에너지·경계 이력이 생겼지만 전체 기간·자기 결합·동일 해의 최종 전하 판독은 남는다. [근거](../notes/REQUEST188_SAME_SOLUTION_GR_RETURN_EVOLUTION_KO.md).


## 단계189–193 — 마지막 구간의 잔차 병목

분류: Counterexample candidate. 최종 전하는 미판정이다. 선형 시스템의 검사나 반환 분기 진단을 실제 전체 결합 해의 전하 수락으로 대체하지 않는다. [근거](../notes/REQUEST189_193_FINAL_INTERVAL_REPAIR_KO.md).


## 단계194–195 — 실제 단계 실패와 넉넉한 속행 예산

분류: Counterexample candidate. 최종 전하는 미판정이다. 실제 단계의 벡터 잔차 실패는 작은 물리 모멘트 오차로 면제하지 않았다. 수락된 같은 해의 전체 이력이 준비된 뒤 전하를 읽는다. [근거](../notes/REQUEST194_195_NATIVE_RESIDUAL_CONTINUATION_KO.md).


## 단계195–197 — 실제 잔차 보정과 결합 속행

분류: Counterexample candidate. 최종 전하는 미판정이다. 실제 속행은 failed 상태이며, coarse 수락 단계는 113개다. 원 수락 기준을 넘긴 결과를 완료로 승격하지 않았다. 반환 GR과 같은 해의 전하 판독까지 통과해야 전하 결론을 내릴 수 있다. [근거](../notes/REQUEST195_197_ACTUAL_NATIVE_CONTINUATION_KO.md).


## 단계198–199 — native 정밀도와 실제 진화

분류: Counterexample candidate. 최종 전하는 미판정이다. 속행은 coarse에서 종료됐고 coarse 수락 단계는 114개다. 전체 기간 및 최종 전하를 완료로 판정하지 않았다. 같은 해의 자기 GR·무한대 전하와 지배 오차를 함께 판정해야 한다. [근거](../notes/REQUEST198_199_NATIVE_PRECISION_CONTINUATION_KO.md).


## 단계204–206 — 실제 단계 이력과 초기 GR 원천 오차

분류: Counterexample candidate. 최종 전하는 미판정이다. 초기 원천의B·광자E·Pr가 원시간 기준을 실패해 새 실제 GR 반환 진화의 수락으로 연결하지 않았다. [근거](../notes/REQUEST204_206_STAGE_HISTORY_GR_SOURCE_KO.md).


## 단계207–208 — GR 시간 오차 전달과 광자 이력 완결

분류: Counterexample candidate. 최종 전하는 미판정이다. 초기U시간차25.8719%가 남고 같은 짧은구간compact읽기를 무한대 전하로 승격하지 않았다. 원 방정식에 맞는 누락 광자 이력을 확보했다. [근거](../notes/REQUEST207_208_GR_TIME_ERROR_AND_PHOTON_HISTORY_KO.md).


## 단계200–202 — 고정밀 실제 유속과 결합 진화

분류: Counterexample candidate. 최종 전하는 미판정이다. 속행은 coarse에서 종료됐고 두 경로의 수락 단계는 115/215개다. 전체 기간 통과로 판정하지 않는다. 같은 전체 해의 GR 원천·장 판독은 아직 수락하지 않았다. 같은 해의 자기 GR·무한대 전하 및 지배 오차를 함께 판정해야 한다. [근거](../notes/REQUEST200_202_ACTUAL_PRECISE_NATIVE_FLUX_KO.md).


## 단계209–210 — 운동량 수정과 실제 속행 입력 동결

분류: Counterexample candidate. 최종 전하는 미판정이다. 115개 수락 상태에서 실제 속행을 시작했으며 여기의210근거는 실행 입력 동결이다. 다음 단계나 최종 전하의 수락 결과로 세지 않는다. [근거](../notes/REQUEST209_210_ACTUAL_MOMENTUM_CONTINUATION_KO.md).


## 단계211–212 — 실제 Radau·입사 펄스의 GR 적용

분류: Counterexample candidate. 최종 전하는 미판정이다. 시간 표현 수정 후 실제U대조40.8632%로 원2%를 실패했다. 짧은compact끝점 부호를 무한대 전하나 관측량으로 해석하지 않는다. [근거](../notes/REQUEST211_212_RADAU_PULSE_GR_KO.md).


## 단계213 — 같은 장기 해의 광자 이력 속행 입력

분류: Counterexample candidate. 최종 전하는 미판정이다. 이미 수락된공통15/16기간의 물질 이력을 재사용해 누락 광자107/207단계 복원을 시작했다. 앞4/8단계 재사용·시각/가중치/상태 접두부·재시작 검사는 통과했다. 이 기록은 실행 입력 동결이며 종료·GR·최종 전하 수락이 아니다. 원 실패와 전체 완료 요건을 유지한다. [근거](../notes/REQUEST213_LONG_SAME_SOLUTION_PHOTON_HISTORY_KO.md).


## 단계213–215 — 동일 좌표·연산자 복원과 실제 단계 경계

분류: Counterexample candidate. 최종 전하는 미판정이다. 실제116단계는 원 기준을 통과했지만117선형 풀이와 저장 제안의 실제 비선형식은 실패했다. 넓힌 잔차 보정만으로 수락하지 않았다. 저장652단계의 정확 보존좌표와 원 구간별 모델 수명을 복구했고 실제 광자 재개에 적용해 기존 실패 쌍을 통과했다. 배열 비교 오류를 고친 뒤 기존 체크포인트에서 속행한다. 이 기록은 전체 복원·자기GR·최종 전하 통과가 아니다. 원 실패·기준·전체 완료 조건을 유지한다. [근거](../notes/REQUEST213_215_COORDINATE_RECOVERY_AND_STAGE_LIMIT_KO.md).


분류: Counterexample candidate. 214재개는원T/16누적반경출구 기준4.741e-12>1e-12로 종료됐다. 광자끝점·물질수지·각도출구와 개별 결합식/native동일성 통과로 이를 대체하지 않는다. 마지막coarse11/fine15복원 체크포인트를 보존했고 최종 전하는 미판정이다. [종료 근거](../notes/REQUEST213_215_COORDINATE_RECOVERY_AND_STAGE_LIMIT_KO.md).


## 단계216–219 — 보정기를 실제 결합 진화에 적용

분류: Counterexample candidate. 최종 전하는 미판정이다. 같은117선형식의 인접 표현 좌표 보정이7.9803e-15로 통과한 뒤, 이를 실제 진화에 적용해117단계를 원 비선형5.2715e-13와 물리 모멘트 기준으로 수락했다.116저장 상태·이력과 선형 우변·분기·잔차 일치를 검사했고 원 실패들은 보존했다. 광자 이력의 누적 반경 출구 실패는 원 누적 순서와 정확히 같은 시각으로도 유지되므로 전체 복원을 수락하지 않는다. 본 기록은 첫 실제 진전과 계속 실행할 입력의 증거이며 전체 기간·자기GR·전하의 완료가 아니다. [근거](../notes/REQUEST216_219_ACTUAL_ROUNDING_CONTINUATION_KO.md).


## 단계220 — 원 광자 입력과 추가 반복의 한계

분류: Counterexample candidate. 최종 전하는 미판정이다. 원 체크포인트의 광자 입력으로 실패 구간8블록을 다시 풀어도 반경 출구4.5472e-12는 원1e-12기준을 실패했다. 저장된 마지막 선형쌍만11회 더 보정한 뒤에도 출구4.5456e-12로 거의 바뀌지 않았다. 원 전체 결합식·native동일성·에너지/물질수지 통과로 이 경계 실패를 대체하지 않는다. 원 실패와 강화목표 실패를 모두 보존하고 동일 반복을 자동 확대하지 않았다. 실제 결합218의 고정 실행은 별개로 유지한다. [근거](../notes/REQUEST220_PHOTON_INPUT_AND_REFINEMENT_KO.md).


## 단계221–223 — 원 결합 광자 단계의 회수와 적용

분류: Counterexample candidate. 최종 전하는 미판정이다. 짧은 원 결합 구간8단계를 재생하며 누락된 실제 광자 단계를 저장했고, 모든 원 물리 배열·내부 끝점을 정확히 재현했다. 같은 이력의 반경 출구1.61324e-13은 원1e-12기준을 통과했다. 수정된16단계 체크포인트를 후속 복원에 실제 적용하여 두 경로가 추가 단계를 수락했다. 원 조건부 실패와 출처가 독립적이지 않았던221끝점 대조 실패를 보존한다. 전체 이력·GR 시간 오차·최종 전하 수락은 아직 아니다. [근거](../notes/REQUEST221_223_ORIGINAL_PHOTON_CAPTURE_KO.md).


## 단계224 — 배경 구간을 연결한 동일 이력 GR 판독

분류: Counterexample candidate. 최종 전하는 미판정이다. 수락된T/8 물질·광자·출구 이력을 실제 배경 구간과 원 입사 파형으로 GR에 연결했다. 배경 에너지 상쇄로 실패한 원 다항식 기준을80자리 동일식 계산으로 해결하고 원 수치 대조를 통과했다. 원 실패와 이전T/64 시간 실패는 유지한다. 전체 공통 이력이 수락되면 같은 판독기를 적용하며, 전 기간·자기GR·비선형·최종 전하의 수락은 아직 아니다. [근거](../notes/REQUEST224_DENSE_GR_BRIDGE_KO.md).


## 단계225–226 — 같은 실제 단계에 밀집 GR 반환

분류: Counterexample candidate. 최종 전하는 미판정이다. 동일T/8 이력의 수락된 밀집 GR을 물질·광자 방정식에 실제 반환했다. 첫 배경 전환의 저장 native 재현 실패를 원 모델의 구간별 재생성으로 수정했고, 거친15단계를 원 기준으로 완료했다. 첫8단계와 재시작 이력의 정확 동일성을 확인했다. 이전 국소 시간 실패와 원 실패는 보존한다. 짝 시간·자기GR·전 기간·비선형·최종 전하의 전체 수락은 별도다. [근거](../notes/REQUEST225_226_DENSE_GR_RETURN_KO.md).


## 단계223–231 — 실제 단계 GR 및 마지막 결합 구간

분류: Counterexample candidate. 최종 전하 미판정.226의 실제 두 경로는 완료했으나 시간 최대5.03484%로 실패했다. 동일 GR의 단계 시각 및 미분을 수정하여227의 실제 두 경로를 다시 완료했다. 원 시간 판정=True, 최대0.416827614%. 이전 실패·기준을 보존하며 자기GR 고정점이나 최종 전하로 확대하지 않는다. 공통15/16광자 이력은111/215단계 복원을 완료했다. 전 기간118단계는 선형 보정과 보존량 변환의 고정밀화를 각각 실제 풀이에 적용했으나 비선형 기준에서 계속 실패했다. 변환 정밀도만으로 해결됐다는 결론을 기각하고117수락 상태를 유지한다. [근거](../notes/REQUEST223_232_ACTUAL_STAGE_GR_KO.md).


## 단계233–235 — 열 좌표 산술 수정의 실제118단계 적용

분류: Counterexample candidate. 최종 전하 미판정. 이전 비선형 실패를 정확히 재현한 뒤, 고정밀 native 내부의 열 좌표 반올림을 제거하고 실제 같은 결합 해에 적용했다. 막혀 있던118단계가 원 비선형 잔차1.10532e-15와 물리 모멘트 기준을 통과했다. 원 실패·기준을 보존한다. 남은 단계와 다른 산술의 미세 경로 교차 평가, 같은 해의 전 기간·자기GR·최종 전하 수락은 별도다. [근거](../notes/REQUEST233_235_THERMAL_NATIVE_KO.md).


## 단계236–238 — 긴 구간 실제 GR 반환 및 마지막 운동량 단계

분류: Counterexample candidate. 최종 전하 미판정.235의 마지막119단계 실패를 정확히 재현한 뒤 운동량 수정과 수지를 실제238풀이에 적용했다. 마지막 단계가 원 잔차9.16708e-15로 통과하여 거친 경로의 원 전체 기간3.434431ms를119단계로 완료했다. 수정 산술의prefix 수지·미세 경로 교차 평가·짝 시간 판정은 별도다. 공통15/16의 원천·GR 시간 통과는 원111/215단계 실제GR반환에 연결했고 독립 장 계산 세 개를 별도 코어에서 병행한다. 원 실패·기준과 자기GR·최종 전하 미해결 조건은 유지한다. [근거](../notes/REQUEST236_238_PARALLEL_GR_MOMENTUM_KO.md).


## 단계239 — 동일 산술의 원 미세 경로

분류: Counterexample candidate. 최종 전하 미판정. 이전 산술의232미세 경로는220단계 수락 후221비선형 잔차에서 실패했다. 원 실패와 수락 상태를 보존한다.239는215단계 모든 저장 배열을 정확히 재시작하고, 거친119단계를 통과한 같은 열 좌표·B/S 산술로 원16개 미세 단계와 짝 시간 비교를 진행한다.238의 거친 저장prefix 물질 수지는 통과했으나 전체 벡터 균일 인증이나 최종 전하 완료가 아니다. [근거](../notes/REQUEST239_COMMON_ARITHMETIC_FINE_KO.md).


## 단계240–243 — 같은 전체 기간의 누락 광자 연결

분류: Counterexample candidate. 최종 전하 미판정. 거친 원119단계의 광자·물질·출구 이력을 모두 연결했다. 기존111단계와 실제118/119광자 캡처를 재사용하고 누락6단계만 원 방정식·정확 native 이력·끝점·수지·출구 기준으로 판정한다.240의물질 복원 버전,241의시각 캐시,242의정규화B잔차 실패를 보존한다. 실제199복원·원 저장 시각·모든 보존량이 정확히 같은B좌표 복원을 적용하고 수락 상태·실패 광자 쌍을 재사용한다.238의 두 저장prefix 수지는 통과했으나 균일 벡터 인증은 아니다. 원 미세 경로/짝 시간 판정과 실제GR반환, 전체 전하 완료 조건은 유지한다. [근거](../notes/REQUEST240_243_COMPLETE_COARSE_PHOTONS_KO.md).


## 단계244 — 실제 캡처로 전체 기간 원천 연결

분류: Counterexample candidate. 최종 전하 미판정.239의 실제 미세 완료와 원 짝 수락 이후, 기존 거친 광자 이력·미세215단계와 실제32개 새 캡처를 재사용하는 소비자를 연결했다. 모든 배열이 같은 기존 거친 완료를 재현하고 시각 불일치를 거부하는 대조를 통과했다. 대기 상태는 실제 전체 조립이나 물리 수락이 아니다. 원 기준·이전 실패·236실제GR반환과 전체 완료 조건을 유지한다. [근거](../notes/REQUEST244_FULL_CAPTURED_SOURCE_KO.md).


## 단계245 — 실제 반환 기하와 같은 해의 원천 판독

분류: Counterexample candidate. 최종 전하 미판정. 짧은 실제227반환 해의15/29단계 끝점을 그 해의 반환 기하·물질·광자·출구와 함께 판독했고, 원천 시간 최대0.635%와 압력/수지가 원 기준을 통과했다. 진행 중인 긴236이 실제 수락된 뒤 동일 판독을 적용한다. 높은/낮은 성분을 유지하며 최초 구동의 한 모드 표현을 반환 기하로 대체하지 않는다. 끝점 검사는 연속 시간·최종 전하·물리 폐쇄의 완료가 아니다. 원 실패와 전체 완료 범위를 유지한다. [근거](../notes/REQUEST245_RETURNED_JOINT_SOURCE_KO.md).


## 단계246–247 — 같은 실제 반환 해의 compact 전하 판독

분류: Counterexample candidate. 최종 전하 미판정. 짧은 같은 반환 해의 compact 전하는 음의 부호를 유지했고 낮은/높은 성분 비는1.63043e-27이었다. 시간 대조1.918%는 원2%기준에 가까우며 무한대 전하나 물리 부호 보장은 아니다. [근거](../notes/REQUEST247_SAME_RETURN_CHARGE_KO.md).


## 단계248 — 같은 인과적 GR 이력의 전체 기간 연장

분류: Counterexample candidate. 최종 전하 미판정. 기존535시각 중532개를 재사용하고 마지막3개와 접합2개만 계산해 독립 세 전체 GR 장을 모든 저장 배열에서 정확히 재현했다. 실제 원천 계수 변화는 거절한다. 전체244원천 수락 뒤 원535시각과 같은 과거 원천을 그대로 유지하며 추가 시각만 계산하도록 연결했다. 원 수락 기준과 초기 실패를 보존하고 전체 기간 실제 반환·자기GR·물리 오차·최종 전하 범위는 축소하지 않는다. [근거](../notes/REQUEST248_CAUSAL_GR_EXTENSION_KO.md).


## 단계249 — 종료점 미분을 반영한 실제 결합 해의 연장

분류: Counterexample candidate. 최종 전하 미판정. 기존 종료점이 내부 시각이 되면 원천 다항식의 미분이 바뀌므로14/16저장 상태에서15/16을 겹쳐 다시 풀고16/16까지 연결한다. coarse103단계·44개 저장 배열의0단계 복원은 정확했고 late native anchor도 원 기준을 통과했다.236/248실제 수락 뒤 같은 전체 입력으로 물질·광자를 이어 풀도록 연결했다. 원 해·실패·모든 수락 기준 및 전체 최종 전하 범위를 유지한다. [근거](../notes/REQUEST249_FULL_RETURN_CONTINUATION_KO.md).


## 단계250 — 전체 원 경로와 긴 실제 GR 반환 쌍 수락

분류: Counterexample candidate. 최종 전하 미판정. 원 물질·광자의 전체119/231단계는 시간 대조 최대0.00267%, 실제 GR 반환 공통111/215단계는0.01733%로 원2%기준을 통과했다. 원 전체 경로는 이미 완성됐고 누락된 비교 파일만 연결해 동일 audit를 다시 실행했으며 물리 이력 SHA는 그대로다.188초기 국소 실패·균일 오차 및 전체 최종 전하 범위를 유지한다. 같은 반환 해의 전하 판독과 전체 기간의 실제 반환으로 이어간다. [근거](../notes/REQUEST250_COMPLETED_PHYSICAL_PAIRS_KO.md).


## 단계251 — 같은 반환 해의 긴 전하 판독과 전체기간 연결

분류: Counterexample candidate. 최종 물리 전하 미판정. 공통15/16실제 해의 compact전하 대조 수락=True, 조건부 부호 유지=True. 전체 원천의 종단 시간 표현 때문에494개 정확히 같은 과거 GR출력만 재사용하고81개를 재계산한다. 원 source-prefix실패를 보존했다. 전체249반환을 동일 해의 전하 판독에 연결했으며, 짧은 해의 모든 원천·GR배열과 전하를 정확히 재현했다. 전체 EOS·미분·공간·경계·비선형/자기GR·정적/관측·무한대 범위와 기존 실패는 유지한다. [근거](../notes/REQUEST251_FULL_RETURN_CHARGE_ROUTE_KO.md).


## 단계252 — 같은 반환 해의 외부·질량 판독

분류: Counterexample candidate. 최종 물리 전하 미판정. 긴 같은 반환 해의 실제 Radau 방출과 동일 원천 질량을 고정 외부 연산자로 읽어 음의 부호를 유지했다. 작은 낮은 성분의 질량 정규화 변화도 따로 적용하고 시간·구적을 성분별로 판정했다. 물리 에너지 변환·시간 의존 외부·배경 질량 재정규화와 전체 EOS/미분/공간/비선형/관측 범위는 남는다. 전체249계량은 완료되어 실제 마지막 두 구간의 결합 진화에 적용됐으며,251후속에 같은 외부 판독을 연결했다. [근거](../notes/REQUEST252_SAME_RETURN_EXTERIOR_KO.md).


## 단계254 — 실제 반환의 후반 선형 풀이 속행

분류: Counterexample candidate. 최종 물리 전하 미판정.249의 완료 계량을 실제 반환에 적용했으나113단계의 선형 풀이가 실패했다.253은12회 보정으로 첫 선형계를 통과했으나 다음 계에서 잔차가 정체됐다.254는 실패한 호출만 더 큰 Krylov공간으로 이어 풀며 원 물리 방정식·수락 기준을 유지한다.111/112저장 상태 재현은 정확히 통과했다. 전체 기간 쌍이 통과하면 동일 해의 전하 및 고정 외부·질량 판독을 자동 수행하도록 연결했다. 전체 물리 전하·자기GR와EOS/미분/공간/비선형/정적·관측 조건 및 원 실패는 남는다. [근거](../notes/REQUEST254_ACTUAL_RETURN_KRYLOV_KO.md).


## 단계256 — 전체 기간 실제 반환 완료

분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 원119/231실제 반환을 완료하고 원 시간 대조를 통과했다.254fine229의 선형 실패는 기존 오른쪽 전처리와 저장 해 재사용으로 실제 단계에서 해소했다. 완료coarse와fine215를 재사용하고228단계 재생의 정확한 일치를 확인했다. 같은 해의 전하와 고정 외부·질량 판독을 연결했다.255물리 경계 에너지 입력은 실제 시계로 준비했지만 전파 일·물리 외부/질량 접합 및 실제 계량 적용은 미완료이며,옛 진단을 전하에 가산하지 않았다. 자기GR·EOS/미분/공간/비선형/정적·관측 조건과 원 실패는 그대로 남는다. [근거](../notes/REQUEST256_FULL_RETURN_RESIDUAL_KO.md).


## 단계257 — 물질 성분 정확도 수정과 동일 해 전하

분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 전체 119/231단계 실제 해의 fine 220단계 물질 성분 오차를 더 엄격한 실제 단계 수락으로 수정했다. coarse 119단계와 fine 215단계를 재사용하고 후반 16단계만 다시 풀어, 원 dense 1e−12와 실제 시간·compact 전하 기준을 통과했다. 조건부 compact 음의 부호는 유지됐다. 그러나 고정 외부·질량 감사는 low 동질 질량의 시간 차이 2.005740886%가 원 2%를 넘어 탈락했다. 차이는 주로 기체 비정지 에너지 항에 있으며 합산 반올림으로 설명되지 않는다. 원 256 실패와 연결 키 오류, 새 감사 탈락을 보존한다. 물리 외부 에너지/전파 중 에너지 변화·질량 접합과 실제 계량 반환, 자기 GR 및 EOS/균일 오차/비선형/정적/관측 조건도 남는다. [근거](../notes/REQUEST257_MATERIAL_ACCURACY_CHARGE_KO.md).


## 단계258 — 질량 압력일의 구간 끝 미분 수정

분류: Counterexample candidate. 최종 물리 전하는 미판정이다. 실제 저장 수송률700개를 정확히 재현하여 시간 차이가 명시적 계량 압력일에 집중됨을 확인했다. 닫는 Radau단계에서 다음 원천 구간의 미분을 읽는 문제를 수정한 입력으로 원119/231실제 결합 반환을 시작했다. 원257질량2.005740886%실패와 모든 기준을 보존하며 새 해의 질량·전하 판정은 아직 없다. 다른 해의 진단값이나 부분적분 보정값을 전하에 가산하지 않는다. 같은 새 해의 전하·질량 감사까지 후속 실행을 연결했고 물리 외부·EOS/균일 오차·자기GR·비선형/정적/관측 범위도 유지한다. [근거](../notes/REQUEST258_MASS_WORK_ENDPOINT_KO.md).


## 단계259 — 같은 결합 해의 반복 비용 수정과 속행

분류: Counterexample candidate. 최종 전하의 결론은 미판정이다. 원258미분 수정과 수락 기준을 유지한 동일 저장 선형계에서 잦은 참 잔차 갱신이 비용을15.677배 줄였다. 그 방식을 실제 결합 진화에 적용했고, 저장8단계의 모든 배열·native률을 정확히 복원한 뒤 새 canonical 구간이 기존 진행량을 따라잡았음을 확인해 기존 프로세스를 교체했다. 비용 비교를 전체 물리 정확도나 최종 성과로 확대하지 않는다. 같은 새 해의 질량·전하 자동 판독과 전체 물리 폐쇄의 미완료 범위를 유지한다. [근거](../notes/REQUEST259_ACTUAL_KRYLOV_COST_KO.md).


## 단계260 — 수정된 동일 해의 조건부 전하 수락

분류: Counterexample candidate. 원119/231단계를 완료한 수정 결합 해의 질량 시간 차이는 2.00574089%에서 1.32954088%로 줄어 원2%기준을 통과했다. 같은 해의 compact 및 고정 외부·질량 판독에서 기존 음의 전하 부호가 유지됐다. 원 실패·수락 기준은 보존했다. 전체 물리 전하의 미판정 사유는 이제 이 수치 수락 실패가 아니라 시간 의존 물리 외부·배경 정규화·자기GR 및 원 EOS/관측 폐쇄다. 이번은 동일 해의 전하까지 도달한 loophole progress다. [근거](../notes/REQUEST260_CORRECTED_SAME_SOLUTION_CHARGE_KO.md).


## 단계261 — 같은 이력의 동적 외부 광자 원천

분류: Counterexample candidate. 기존 수락 해의 방출과 같은 EOS 진공·입사장을 사용해 광자의 에너지·반경·방향 변화를 원 전체 기간에서 실제 전파했고 네 원 구적 대조와 high 방출 시간 기준을 통과했다. 실제 패킷 응력을 원 반경의 GR 질량·lapse 제약 특수해에 연결했다. 단계260의 조건부 음의 전하 수락은 보존하되, 새 원천에 대한 완전한 경계·물질 재적용과 최종 전하는 미판정이다. 반환 계량 기하 응답, 상호 스칼라 일·배경 연산자·질량 접합을 누락한 채 전하에 보정값을 가산하지 않는다. [근거](../notes/REQUEST261_PHYSICAL_EXTERIOR_PHOTONS_KO.md).


## 단계262 — 반환 계량의 외부 반경 연결

분류: Counterexample candidate. 같은 수락 이력의 반환 질량·광자 압력·outgoing 스칼라 경계를 실제 외부 반경으로 연장했고 기존 출구 lapse와 독립 미분 대조를 통과했다. 이 계량을 배경 광자 전파에 연결했으며 대표 비용 검사 후에만 원 네 구적 경로를 실행한다. 최초 비용 중단과 중복 초기화 오류를 보존했다. 새 외부 광자 전 기간, reciprocal 스칼라·배경 연산자·접합·새 물질 결합 및 최종 물리 전하는 아직 미판정이다. 기존 조건부 음의 전하 수락은 유지한다. [근거](../notes/REQUEST262_RETURNED_RADIAL_METRIC_KO.md).


## 단계263 — 같은 특성식의 누적량 전파

분류: Counterexample candidate. 시간 미분과 광선 교차 압력 항을 정확한 변수변환으로 누적 방출량·장 값에 옮기고, 원 허용오차의 실제 대표 cohort 세 개를 완료했다. 서로 다른 변환의 궤적·응력 대조가 원0.2%기준을 통과했다. 복원된 계량 일 항등식은 독립 일 적분 검증으로 세지 않으며 원 시간 미분식의 별도 실제 광선 대조를 전체 생산의 전제조건으로 둔다. 실측 보수적 경로 예상6.09시간을 근거로 원4시간 비용 부적합을 보존하고8시간 예산으로 원 네 경로를 연결했다. 물질·EOS를 다시 계산하지 않았다. 전체 물리 전하는 미판정이다. [근거](../notes/REQUEST263_PRIMITIVE_PHOTON_CONTINUATION_KO.md).


## 단계264 — 독립 광선 대조 완료와 전 기간 광자 전파

분류: Counterexample candidate. 단계263의 원 시간 미분식 대조가 1시간 관문 예산에서 멈춘 실패와 원 기록을 보존했다. 등록 바인딩 확인 뒤 저장된 packet 0를 재사용하고 packet 31만 같은 코드·허용오차로 완료했으며, 두 광선 대조 최대 1.154e-05로 원0.2%기준을 통과했다. 이어 등록된 네 구적 경로가 원 16개 출력 시각을 모두 완료했고 응력·광자 질량 원천·lapse 경계 대조가 원0.2%기준을 통과했다. 격자·기간·경로 수·기준은 바꾸지 않았다. 반환 원천 64/128 시계 대조와 reciprocal 스칼라·배경 연산자·물질 질량 접합, 결합 해 적용이 남아 최종 전하는 미판정이며 단계260의 조건부 음의 전하 결론을 유지한다. [근거](../notes/REQUEST264_DIRECT_REFERENCE_RETRY_KO.md).


## 단계265 — 1회 반환 외부 광자 경계를 적용한 같은 해의 전하

분류: Counterexample candidate. 입사 계량의 배경 광자 기하 lapse를 575개 적용 시각에서 계산해 반환 계량 outer lapse에 넣고, 원 119/231 결합 쌍을 t=0부터 다시 진화했다. 16개 매듭은 단계261과 비트 일치했고 계량·결합·판독의 원 기준을 통과했다. 첫 재진화는 단계257 내부 물질 기준의 H 성분(1e−13)이 반올림 바닥(1.00e−13)에 걸려 멈췄고, 사용자 승인으로 H만 2e−13으로 바꿔 다시 진화했다. 그 실행은 coarse 마지막 단계에서 double 선형 풀이가 long-double 수락 연산자를 대표하지 못해 멈췄다. 두 double 풀이가 모두 실패할 때만 쓰는 long-double flexible GMRES를 더하고, 15/16 구간을 정확한 재시작 검사 뒤 재사용해 마지막 구간을 다시 진화했다. 같은 해의 compact·고정 외부·질량 정규화에 배경 광자의 발사 에너지까지 더한 미세 시계 합은 -2.483473e-51로 음수로 유지됐다(단계260 대비 상대 -3.86e-08). 다섯 실패 시도(H 바닥, NameError, Newton 12회, double 선형 풀이, 대체 풀이 로그의 JSON 오류)와 작업자 스케줄 실패(컨트롤러 종료 시 WSL 세션 작업자 전체 종료)는 보존했다. 외부 스칼라 연산자 변분·자기GR·EOS/공간·관측 폐쇄가 남아 전체 물리 전하는 미판정이다. [근거](../notes/REQUEST265_ONE_RETURN_PHOTON_BOUNDARY_KO.md).


## 단계266 — 수락된 primary 전하의 깊이 분해와 지배 오차

분류: Counterexample candidate. 단계265에서 유지된 조건부 음의 전하의 high 성분을 원천 깊이 대역으로 선형 분해했다. 대역 합은 전체와 같았고 저장값을 비트 단위로 재현했다. 셀 11(276–345km)이 62%, 셀 12가 29%를 차지하며, 가장 깊이 보이는 셀 9·8은 양의 기여 -11%다. 대역별 시간 격자 차이는 지배 셀에서 1e−4 이하다. 반응이 약 3개의 미세분 68.75km 셀에 몰려 있으므로, 현재 지배 오차는 구동된 primary 이력의 내부 반경 해상도로 판단한다(Conjectural). 음·양 기여의 비는 약 10:1이다. 세분 격자에서 primary 이력을 다시 진화하는 결정적 시험이 다음이며, 최종 전하는 미판정이다. [근거](../notes/REQUEST266_PRIMARY_DEPTH_DECOMPOSITION_KO.md).


## 단계267 — 내부 셀을 세분한 primary 재진화와 끝점 전하

분류: Counterexample candidate. 단계266이 지배 오차로 지목한 내부 반경 해상도를 고쳤다. 셀 8–15를 2배 세분한 격자(539셀)에서 구동 primary(64 시계)를 같은 최종 방정식으로 t=0부터 끝점까지 진화했다. 끝점 compact 전하는 -2.3721e-51로 음이고, 원 격자 -2.4833e-51보다 크기가 4.5% 작다. 사전 등록 규칙에 따라 조건부 음의 전하 결론은 유지된다. 2% 기준을 넘으므로 이 수준의 해상도 수렴은 미달이다. 마지막 2.5 macro 단계(실제 5단계)는 벡터 선형·비선형 잔차가 표현 바닥(증폭 약 5×10⁸)에 걸렸다. 그래서 사용자 승인에 따라 물리 모멘트·물질 성분 1e−13 기준으로 수락했다. 모든 원 기준을 지킨 t=61/64·T 비교는 −5.4%다. 변화는 셀 11(−5.0%)과 셀 12(−3.1%)가 주도한다. 재생성 범위의 교훈: 구동 primary는 단계149–157 보정을 쓰지 않았다(입력 2배 교란에도 첫 단계 비트 동일). 존재 검사가 산술 경로를 바꾸는 경우도 감사해야 한다. 최종 물리 전하(128 시계, 1회 GR 반환·외부 광자 항의 세분 재계산, 자기GR·EOS·비선형·관측 폐쇄)는 미판정이다. [근거](../notes/REQUEST267_REFINED_PRIMARY_KO.md).


## 단계268 — 내부 셀 4배 세분의 해상도 수렴 시험

분류: Counterexample candidate. 셀 8–15를 4배 세분한 격자(555셀)에서 구동 primary(64 시계)를 같은 최종 방정식으로 t=0부터 끝점까지 진화했다. 끝점 compact 전하는 -2.3443e-51로 음이다(2배 -2.3721e-51, 1배 -2.4833e-51). 사전 등록 규칙에 따라 조건부 음의 전하 결론은 4배 격자에서도 유지된다. 2×→4× 크기 변화는 t₆₁ -1.42%, T -1.17%이며, 두 시각 모두 2% 이하이므로, 2배 결과를 2% 수준의 해상도 수렴으로 판정한다. 관측 차수는 t61 2.01, t62 2.01, t63 2.00, T 2.00다. 최종 물리 전하(128 시계, 1회 GR 반환·외부 광자 항의 세분 재계산, 자기GR·EOS·비선형·관측 폐쇄)는 미판정이다. [근거](../notes/REQUEST268_QUADRUPLE_REFINED_PRIMARY_KO.md).


## 단계269 — 전하 지배 층의 EOS 물리 민감도

분류: Counterexample candidate. 계산 줄기의 FreeEOS(option 11: Planck–Larkin + MDH)와 수소 준위 라이브러리를 PL만 끈 변형으로 다시 빌드했다(동일성 빌드는 원본과 비트 단위로 일치). 그 라이브러리로 내부 셀 0–26의 EOS·광학 배열을 다시 만들어 2배 격자 primary를 끝점까지 진화했다. 전하 층의 중성 분율은 최대 3배, 일부 광자 계수는 수백 배 바뀌었다. 그런데도 끝점 compact 전하는 -2.3721e-51로 PL/MHD 해와 상대 +9.2e-11만 다르다. 사전 등록 규칙에 따라 조건부 음의 전하 결론은 유지된다. 전하 원천은 정지질량(metric stress) 성분이 지배하고, EOS·광학에 민감한 열·광자 성분은 그 1e−7–1e−8 규모다. 배경 밀도 구조의 EOS·불투명도 의존성은 시험하지 않았다. 최종 물리 전하는 미판정이다. [근거](../notes/REQUEST269_EOS_PHYSICS_SENSITIVITY_KO.md).


## 단계270–271 — 끝점 전하의 물리적 구성과 자유낙하 재현

분류: Counterexample candidate. 4배 해의 끝점 compact 전하를 원천 성분·부분별로 정확히 분해했다(닫힘 1e-16). 전하는 바리온 질량 섭동의 상태 응답이 0.999996를 차지하고, 입사장×배경의 직접 결합은 1e-08이며, 열·광자 성분은 3e−9 이하다. 바리온 섭동은 질량을 보존하는 재배치이며, 끝점 전하의 99.9%는 펄스가 이미 떠난 층의 동결 변위(기억)에서 온다. 선언 이론(A=exp(−2φ²))에서 유도한 무압력 자유낙하 모형은 지배 셀의 끝점 δM을 4배에서 6e−6 이내로 재현했다. 그 전하의 연속 극한 -2.3351e-51는 계산 줄기의 극한 -2.3350e-51와 상대 +5e-05로 일치한다. 전하는 배경 ρ₀·φ₀와 입사 펄스의 명시적 범함수이며, 최종 물리 전하는 미판정이다. [근거 270](../notes/REQUEST270_CHARGE_COMPOSITION_KO.md), [근거 271](../notes/REQUEST271_FREEFALL_REPRODUCTION_KO.md).


## 단계272–273 — 배경 구조 민감도와 정적 EFT 붕괴 경계

분류: Counterexample candidate. 자유낙하 전하는 면 밀도에 정확히 선형이다. 4배 격자의 면 핵(합의 닫힘 2e-16)은 거의 모두 전하와 같은 부호이며(반대 부호 몫 6.1e-04), 깊이 240–465 km에 모인다. 부호를 뒤집으려면 면 밀도의 최대노름 상대 변화가 0.9988여야 하므로, 밀도가 양수인 한 어떤 봉투 분포도 부호를 뒤집지 못한다. 크기는 이 깊이의 밀도에 비례하여, 분포가 10 km 어긋나면 약 19%(바깥쪽)/-17%(안쪽) 바뀐다. 전하 이력(128시각)에서 동적 변위 부분의 비중은 모든 시각에 1−1e−8이며, 전하는 펄스가 표면을 떠난 뒤에도 더 깊은 층의 변위로 약 600배 자란다. 배경 항성의 가장 낮은 반경 단열 모드 주기는 178 s, 전하 층의 음향 차단 주기는 28–47 s다. 긴 파장(λ≫R) 단극 힘은 펄스 힘의 3.3e-08배이고, J0337 내측 궤도에서 (ω/ω₀)²=1.6e-06이다. 따라서 이 자유낙하 기억은 궤도 시간척도에서 정적 계수로 붕괴한다(no-go 경계, 조석 채널·소산은 미포함). [근거 272](../notes/REQUEST272_BACKGROUND_STRUCTURE_KO.md), [근거 273](../notes/REQUEST273_STATIC_EFT_BOUNDARY_KO.md).


## 단계274–275 — 조석·소산 경계와 광구 밀도 불확실성

분류: Counterexample candidate. 정적 구대칭 배경에서 compact 전하는 선형 차수로 구동의 l=0 성분에만 반응한다(Proven 선택 규칙). 선언 배경은 중심 밀도가 평균의 488배라 조석 근점 상수가 k₂=3.29e-04이고, J0337 내측 궤도에서 정적 조석의 상대 힘은 9.1e-12다. J0337 조석 구동은 l=2 g모드 차수 약 1453에 해당하며, 광학적으로 두꺼운 봉투만으로 잰 감쇠 깊이가 8.9e+06이라 이산 공명 없는 진행파 영역이다. 자유낙하 기억이 궤도 시간척도에서 붕괴해 들어가는 단극 정적 구조 계수는 |κ_struct|≤8.8e-09(직접 계수 β의 2.2e-09, φ_∞=1e−3)이고, 이를 지연 상한으로 써도 J0337 SEP 진동은 2.7e-18로 Paper B 한계 1.7e−9보다 6.3e+08배 작다. 따라서 no-go 경계를 조석·소산 채널로 확장한다. 전하 층은 광구다(회색 광학깊이 240 km 0.005, 362 km 1.06, 466 km 7.9; 열 시간 1초 이하). Kaplan et al. 2014의 log g 5.82±0.05, T_eff 15,800±100 K로 회색 대기를 정역학 상사 변환하면 전하 크기는 선언 대기의 3.4–3.7배(1σ 1.6–7.7배)이고 부호는 2σ 범위에서 유지된다. 전하 크기의 지배 불확실성은 표면중력이 정하는 광구 밀도다. [근거 274](../notes/REQUEST274_TIDAL_DISSIPATION_KO.md), [근거 275](../notes/REQUEST275_PHOTOSPHERE_DENSITY_UNCERTAINTY_KO.md).


## 단계276 — 최종 전하 판정과 관측 폐쇄

분류: Counterexample candidate. 최종 전하 결론은 유지된다. 지배 수치 오차를 해결한 같은 결합 해의 끝점 compact 전하는 연속 극한 -2.3350e-51(φ_∞=1e−3, η=1e−30; q/(ηφ_∞²)=-2.335e-15)로 음수다. 부호는 EOS·광학, 봉투 밀도 분포, 관측 광구(2σ)에 견고하다. 크기는 관측 표면중력에서 -7.85e-51–-8.73e-51(1σ -3.8e-51–-1.8e-50)로 수정된다. 자유낙하 닫힌 식의 기호 잔차는 0이다. 분류: Conjectural. 관측 폐쇄: Cassini 2σ로 |α₀|≤3.54e-03이면, 내측 백색왜성 단극 전하의 선형·수동 지연이 J0337 지연 한계 1.7e−9에 닿으려면 |κ_lag α_p|≥1.56e+03여야 한다. 정적 감수율 전체(|β|=4)가 완화되어도 4.4e-12|α_p|, 자유낙하 기작은 7.5e-21|α_p|다. 따라서 이 감도의 지연 신호는 백색왜성이 아니라 중성자별 전하 쪽이어야 한다. 분류: Conjectural. 미션 분류는 theorem progress(no-go 경계)이며, 이 최소 상태에 대해 A4가 유지된다. 실패 원장: 정확한 붕괴 단계는 ω≪ω₀, ω≪ω_ac, g모드 진행파 영역에서 자유낙하 기억이 정적 구조 응답(|κ_struct|≤8.8e−9)으로 바뀌는 단계다. 최소 누락 가정은 궤도 주기 완화 시간과 |κ_lag|≳1.6e3/|α_p|의 결합을 가진 내부 상태이며, 약한장 백색왜성에서는 성립하지 않는다. 남은 조건: GR 반환 1회(ADM·무한대 정규화 없음), φ_∞² 스케일은 유도, 끝점 프로토콜 T=2D, 회색 대기, 느린 자전, 원고 통합은 Pandoc·TeX 부재로 보류. [근거 276](../notes/REQUEST276_FINAL_CHARGE_CLOSURE_KO.md), [영문 절 초안](white-dwarf-free-fall-charge-section.md).


## 단계277–278 — 자기 일관 GR·비선형의 폐쇄와 전하의 전체 시간 이력

분류: Counterexample candidate. 자기 일관 GR 고정점 잔차는 4.5e-18 이하(부등식 Proven), ADM 교차 에너지는 6.8e-11, 무한대 꼬리는 4.2e-06 이하(측정 4.4e-15), 비선형은 1.0e-27(α=βφ 정확, Proven)다. 질량 정규화를 적용한 연속 극한은 -2.33502e-51다. 판독 원천의 셀 가중치가 α(φ₀)/r이므로 상태 부분은 φ_∞²에 비례하고, 짝함수 성질로 q/η는 φ_∞의 짝함수다. 선형 단열 반경 모형(전체 별, 중심 반사 포함)은 끝점 전하를 -0.81%로 재현했다. 끝점 값은 정적 창 값의 639배인 지연 단극장이다. 들어오는 구간(0–0.35 s)에서 전하는 음수로 약 12자릿수 자라고(중심 도달 −7.5e−41), 나가는 통과에서 부호가 진동하며 출사 시각 0.463 s에 최대 3.5e-36(입사 진폭의 1.5e−11)다. 펄스가 떠난 뒤 반경 p모드(주요 주기 44–59 s)로 진동하며(rms 8.1e-45), 시간 평균은 0이다. 양의 모드 감쇠에서 영구 전하는 η의 1차에서 0이다(Proven). 따라서 끝점 T=2D에 묶인 조건을 닫는다. 음의 끝점 결론은 0.35 s까지의 모든 끝점으로 넓어지며, 영구 전하는 0이다. [근거 277](../notes/REQUEST277_GR_NONLINEAR_CLOSURE_KO.md), [근거 278](../notes/REQUEST278_TRANSIT_LONG_TIME_KO.md).


## 단계279 — 관측 광구의 직접 재구성과 비회색 LTE 보정

분류: Counterexample candidate. 선언 EOS·Rosseland 표·중력 밀도 cx·ρ로 구면 회색 Eddington 봉투를 적분했다. 깊이 원점은 선언 봉투와 같은 P_gas=1 dyn/cm² 자름점이다. 이 봉투는 선언 배경의 밀도를 전하 층에서 핵 가중 0.25%로 재현했다(2% 관문 통과, 앞선 네 번의 실패 87%·8.6%·3.6%·5.0%는 보존). 질량을 고정한 관측 대기에서 끝점 전하는 선언 모형의 3.25배(Kaplan 중심), 1σ 1.63–6.17배, 2σ 0.78–11.2배다. 이는 단계275의 상사 변환 극한과 맞으며, 두 불투명도 극한 가정을 대체한다. 같은 H·He 연속 불투명도의 LTE 비회색 복사평형(흐름 일정성 5.5e−5, 표면 T₀/T_eff=0.64)은 전하를 1.7–2.8배 키운다. 그래서 관측 대기의 전하는 6.9배(1σ 3.9–11.7배)다. 부호는 모든 경우에 유지된다. 비회색 보정은 LTE, 선·금속 생략(단순 불투명도의 Rosseland 평균은 표의 0.55–0.74배), Eddington 닫힘에 조건부다. 따라서 크기는 회색 값과 LTE 비회색 값을 함께 적는다. [근거 279](../notes/REQUEST279_ATMOSPHERE_RECONSTRUCTION_KO.md).


## 단계280 — 남은 조건의 폐쇄와 최종 전하 진술의 개정

분류: Counterexample candidate. 끝점 T=2D의 compact 전하는 음수다. 질량 정규화 연속 극한은 -2.33502e-51이고, 관측 광구에서는 회색 -7.58e-51, LTE 비회색 -1.61e-50다. 다만 이 값은 입사 펄스의 단극 산란 신호 가운데 가장 이른 광구 부분이다. 신호는 들어오는 구간에서 음수로 자라고, 나가는 통과에서 최대 3.5e-36로 진동하며, 펄스가 떠난 뒤 영평균 반경 모드 진동이 된다. 영구 전하는 η의 1차에서 0이다(양의 감쇠, Proven). 단계276에 남은 조건을 모두 닫았다: GR 고정점 4.5e-18, ADM 6.8e-11, 무한대 4.2e-06, 비선형 1.0e-27, φ_∞ 짝함수·φ_∞² 가중(단계277), 전 시간 이력(단계278), 관측 광구 재구성과 비회색 보정(단계279). 남는 가정은 양의 모드 감쇠, 내부 Born 산란 무시, 비회색 보정의 LTE 연속 불투명도, 직접 결합 부분의 φ_∞ 지수, 느린 자전이다. 원고 통합은 Pandoc이 없어 보류했다. no-go 경계와 A4 유지 판정은 바뀌지 않는다. [근거 280](../notes/REQUEST280_OPEN_CONDITIONS_CLOSED_KO.md), [개정 영문 절 초안](white-dwarf-free-fall-charge-section.md).


## 단계282 — 독립 심사와 통합 되돌림

분류: Imported from prior work. 통합 원고 §4.6(커밋 a322c6462)을 gpt-6-astra(Codex, 읽기 전용), opus5.5, fable5.1이 서로 모른 채 심사했다. 셋 모두 주요 수정을 권고했다(astra는 차단 1건). 수치는 기록과 일치했고, 핵심 결론을 뒤집는 결함은 없었다. 합의된 지적은 다음과 같다: J0337 비교(계수 구간 대 진폭, 두 쌍 채널, 척도 범위, Cassini 재척도), 주장 표지, 판독 시각 기준 이력, 조석 서술, 기호 충돌. 사용자 결정에 따라 통합을 커밋 b40864a22로 되돌렸고 논문 검증은 통과한다. 수정 뒤 같은 세 심사자로 재심한다. [근거 282](../notes/REQUEST282_INDEPENDENT_REVIEW_REVISION_KO.md).


## 단계284 — 재심 응답(4판)과 결론 정정

분류: Imported from prior work. 개정 3판의 재심에서 gpt-6-astra는 주요 수정을 권고했다(차단: §5 감도 결론). opus5.5와 fable5.1은 경미 수정 후 수락이었다. 4판(`docs/white-dwarf-free-fall-charge-section.md`)은 다음을 고쳤다. 감도 배제 결론은 척도 비교로 낮췄다. 열 완화 세기 ≤ 𝒮_struct와 단일 완화는 계산하지 않은 가정으로 명시했다. 순간 응답 β_sδφ_mod(지연 없음, §4.3 조건)를 분리했다. 부호 문장은 '최대노름 99.87% 미만의 변화는 부호를 바꾸지 못한다'로 고쳤다. 새 결합 계산은 없다.

분류: Conjectural. 결론 정정: 궤도 시간척도에서 이 자유낙하 상태가 관측량을 만들지 않는다는 이전 진술과 A4가 유지된다는 진술(단계273·274·276·280)은 증명되지 않았다. 성립하는 것은 두 가정 아래의 척도 비교다. 그 아래에서 구조 변조의 척도가 §5 저장 척도보다 8자릿수 이상 작다. 두 비공통 채널의 타이밍 감도는 계산하지 않았다. 분류는 조건부 theorem progress다. 붕괴를 피하는 최소 추가 조건은 궤도 주기와 비슷한 완화 시간을 가진 상태가 있고, 그 완화 세기가 𝒮_struct를 넘거나 단일 완화가 아닌 것이다. 이전 노트의 정오표(부호 문장, 단계278 이력 값, 반올림)는 [근거 284](../notes/REQUEST284_REVISION4_RESPONSE_KO.md)에 있다. 같은 세 심사자로 다시 재심한다.


## 단계285 — 4판 재심 응답(5판)과 분류 정정

분류: Imported from prior work. 개정 4판의 재심에서 gpt-6-astra는 차단을 해제했지만 주요 수정을 권고했다. opus5.5와 fable5.1은 경미 수정 후 수락이었다. 5판(`docs/white-dwarf-free-fall-charge-section.md`)은 다음을 고쳤다. 순간 응답 β_sδφ_mod는 Section 3의 지연 없는 계수로 적었다. 이 항은 §4.3의 고정 동반성 축약이 빠뜨리며, 펄서–내측 쌍에서 약 1.24e−9 a_p² 이하다. 초록과 요약에는 두 완화 가정을 모두 적었다. 판독은 compact 부분 𝒬_c와 질량 정규화 𝒬로 나눴다. 새 계산은 없다.

분류: Conjectural. 분류 정정: 단계284 항목의 '두 가정 아래 no-go'는 틀렸다. 두 가정은 지연의 크기(≤|𝒮_struct|/2)를 묶을 뿐 지연을 없애지 않는다. 반례는 H=s/2+(s/2)/(1+iωτ)로, 두 가정을 만족하면서 궤도 진동수의 직교 성분이 s/4다. 분류는 '선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress'이며, 관측 배제나 A4 유지는 성립하지 않는다. no-go가 무너지는 정확한 단계는 열 완화 세기다. 단열 정적 계산은 이를 묶지 못한다. 빠진 최소 계산은 두 가지다: 깊은 층의 비단열 열 응답(세기와 극점 구조), 두 비공통 채널의 타이밍 응답. 깊은 층의 열 완화 상태는 세기를 계산하지 않은 loophole 후보로 남는다. 단계272–284 항목의 부호·no-go·A4 진술은 [근거 285](../notes/REQUEST285_REVISION5_RESPONSE_KO.md)의 정오표로 대체된다. 같은 세 심사자가 반영 여부를 확인한다.


## 단계286 — 반영 확인 재심 통과와 원고 재통합

분류: Imported from prior work. 5판의 반영 확인 재심에서 fable5.1은 수락, gpt-6-astra와 opus5.5는 경미 수정 후 수락이었다. 경미 지적을 반영한 6판을 통합 원고 §4.6으로 다시 넣었다. 6판은 순간 응답의 전체 상한과 정적 SEP 양립 조건(Cassini 수준 a_o에서 |a_p|≲4e−3), 외측 백색왜성의 같은 순간 응답, 구조 변조의 합계(각 2.13e−18, 합 4.3e−18), 초록의 t→∞ 한정을 담는다. 초록·§6·데이터 가용성 문장을 함께 넣었고, 참고문헌 세 항목을 복원했다. main.tex는 Pandoc 3.11로 다시 만들었다(고치기 전 원고에서 바이트 재현을 먼저 확인). PDF와 제출 zip은 §4.6 이전 판이다(TeX 없음). 저장소 분류는 선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress다. 관측 배제도 A4 유지도 아니다. [근거 286](../notes/REQUEST286_MANUSCRIPT_REINTEGRATION_KO.md).


## 단계287 — 깊은 층 열 완화 세기의 판정

분류: Imported from prior work. 사용자 승인(2026-09-28, 부분 진행)으로 열 완화 세기를 싼 계산으로 판정했다. 모형은 층별 열 완화다: 각 층의 완화 시간은 위쪽 층의 열 시간이고, 완화 극한은 등온 Γ_T=P_gas/P로 괄호를 쳤다. 단계274 정적 풀이기를 1e−3–1e13 s의 절단 161개로 다시 풀었다. 연산자는 모두 양정치였고, 단열 dq/ε는 그대로 재현됐다. 궤도 진동수의 Debye 가중 지연은 |𝒮_struct|의 4.0e−9(내측 궤도)와 3.3e−7(외측 궤도)이고, 창 상한으로도 ≤8.9e−5다. 수락 기준(1%)을 충족한다. τ_th=1/ω인 층은 깊이 2,600–6,600 km, 위쪽 질량 몫 1e−9–1e−7로 가볍다. 비공통 채널 타이밍과 §4.3 강성 조건은 해석 추정만 남겼다(강성 비 약 1e−13|β_p|). 원고는 고치지 않았다.

분류: Conjectural. 이 상태의 궤도 지연 결합은 층별 완화 모형에서 단열 구조 척도의 ≲3.3e−7이고, 쌍 인자로는 ≲7e−25다. 분류는 이 상태에 대한 계산된 정량 경계(theorem progress)다. 관측 배제나 일반 A4는 아니다. 남은 최소 계산은 두 가지다: 완전한 비단열 반경 확산 응답, 비공통 채널의 타이밍 응답. [근거 287](../notes/REQUEST287_THERMAL_RELAXATION_STRENGTH_KO.md).


## 단계288 — §4.6에 층별 열 완화 추정 반영(단독 검토)

분류: Imported from prior work. 사용자 결정(2026-09-28)으로 문장 수준 변경은 외부 심사 없이 단독 검토로 반영했다. 통합 원고 §4.6의 "Deeper layers" 항목, 요약, 가정 목록, 데이터 가용성에 단계287의 층별 열 완화 추정을 넣었다. 추정 지연은 |𝒮_struct|의 4.0e−9(내측)와 3.3e−7(외측)이고, 계산 결과는 Imported, 모형 한계는 Conjectural로 표지했다. 두 가정 상한, 초록, §6 문장은 여전히 참이라 그대로 두었다. main.tex는 Pandoc 3.11로 다시 만들었고 논문 검증을 통과한다. [근거 288](../notes/REQUEST288_RELAXATION_ESTIMATE_SECTION46_KO.md).


## 단계289 — Physical Review D 제출 준비

분류: Imported from prior work. 사용자 지시(2026-09-28)로 PDF를 만들고 투고처를 조사해 PRD(Regular Article)를 1순위로 정했다. 근거는 범위와 APS의 2026년 6월 AI 정책이다. 이 정책은 실질적 AI 사용을 논문 안에 공개하는 조건으로 허용한다. 대안은 CQG다. Tectonic 0.17.0을 설치해 30쪽 PDF를 빌드했다. 빌드 중 Pandoc이 수식을 놓친 §4.6의 `K≈934` 한 곳을 고쳤다. 원고에는 'Use of AI tools' 절을 넣고 Status 줄을 지웠다. 저널용 소스 zip(자체 컴파일 확인), 1쪽 커버레터, 평문 초록, 체크리스트를 `output/submission/`에 두었다. 남은 일은 저자 몫이다: AI 절 확인, 소속·ORCID, 공개 저장소 push(로컬이 1,137커밋 앞섬), 제출. `paper/package_revision.py`는 manifest를 덮어쓰므로 쓰지 않았다. [근거 289](../notes/REQUEST289_PRD_SUBMISSION_PREP_KO.md).


## 단계290 — 제출 전 최종 검토

분류: Imported from prior work. 원고 전문을 정독하고 자동 점검(참고문헌 실재·교차 참조·조판·표기)과 쪽 렌더링을 했다. 제출을 막는 오류는 없었다. 고친 것은 다음과 같다: 참고문헌에 인쇄되던 내부 메모 삭제, DOI 2건 추가, AI 모델 버전(Claude Opus 5.5, GPT-6-Astra), 미국식 철자·수식·en dash 정리, 저장소 내부 표현 삭제, 데이터 가용성의 공개 스냅숏 문구. 전체 이력(51.6 GB)은 GitHub 한도를 넘어, 저자 결정에 따라 대형 배열을 뺀 공개 스냅숏으로 올린다(`PUBLIC_SNAPSHOT.md`). 남은 필수 항목은 소속·ORCID뿐이다.

분류: Conjectural. PRD 게재 가능성은 약 15–35%(중심 약 25%)로 본다. 초기 반려 약 25–40%, 심사로 가면 약 35–55%다. 주된 약점은 분량·문체, 핵심 정리의 제한된 새로움, 비검출·조건부 결과다. 초점을 좁히고 새로움의 위치를 명시하면 가능성이 오른다. [근거 290](../notes/REQUEST290_FINAL_REVIEW_KO.md).


## 단계291 — 원고 초점 축소

분류: Imported from prior work. 사용자 지시로 PRD 제출 원고의 초점을 좁혔다. 본문(11쪽, 새 제목 *Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715*)은 식별성 경계(유한 반송파 보간 경계, 속도 간극 부등식, 두 표본 모멘트-분산 식별과 가까운 극의 한계, 유한 주기 기록의 과도 응답 불충분성), 감쇠 스칼라 전하 실현과 위상 폐합, J0337 조건부 적용만 남긴다. 알려진 사실은 한 줄이나 인용으로 줄였다. 이전 전체 원고는 제목·초록만 바꿔 보충 자료(`paper/supplement.md`, 30쪽)로 옮겼다. 새 계산·새 수치는 없고(본문 수치가 모두 이전 원고에 있다), 최종 전하 결론은 유지된다. [근거 291](../notes/REQUEST291_FOCUSED_MANUSCRIPT_KO.md).


## 단계292 — 최종 독립 검토와 근점 규약 정정

분류: Imported from prior work. Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra의 독립 심사가 모두 주요 수정을 권했다. Opus가 J0337 물리 구동의 η·κ 규약 오류를 찾았다. 코드는 η=e sin ϖ를 쓰는데 분석은 η=e cos ϖ로 읽었다. 이를 바로잡고 영향받는 계산을 다시 돌렸다(각 수 초). 이전 출력은 `outputs/research-completion/withdrawn-periastron-convention/`에 보존했다.

실패 기록: 실패한 단계는 `symbolic/physical_matching.py`와 `verification/physical_drive_completion.py`의 근점 각 계산이다. 빠진 최소 가정은 런타임 코드의 매개변수 규약 확인이다. 이전 결론 '모든 물리 지연 단면이 β=0을 포함한다'를 철회한다.

분류: Imported from prior work. 수정 뒤 여섯 단면이 모두 비어 있다. 기록된 반송파 초과(옴니버스 16.3525, 명목 p≈0.012)를 물리 구동의 순간+완화 응답이 재현하지 못한다. relaxation 검출은 아니며 원인은 분리하지 못했다. 백색왜성의 최종 전하 결론은 유지된다. 문헌 위치, 두 표본 조건, 한정어, 정의 등 나머지 지적도 반영했다(본문 13쪽). [근거 292](../notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md).


## 단계293 — 확인 심사와 경미 수정

분류: Imported from prior work. 단계292 개정판을 세 심사자(Opus 5.5, Fable 5.1, GPT-6-Astra)가 확인 심사했다. 모두 경미 수정 후 수락을 권했고, 세 명 모두 근점 규약 정정을 독립적으로 확인했다. J0337 결론의 범위를 평가한 여섯 지연과 순간항+단일 완화 모형으로 한정했다. 2일 단면의 근소한 기각을 정량화했고(보정 순서통계량 12.27–12.82, 표준오차 약 0.13), 보관된 도함수 조건을 명시했다.

분류: Proven. 통계량은 볼록 이차식이고 β=0 직선이 모든 단면에 들어 있다. 제약 없는 β̂은 18일을 빼면 음수다. 따라서 등전하 실현의 β≥0 아래에서는 여섯 단면 모두 최솟값이 14.59 이상이다. 백색왜성의 최종 전하 결론은 유지된다. [근거 293](../notes/REQUEST293_CONFIRMATION_REVIEW_KO.md).

## 단계295 — 투고 전 조건 정정과 공개 재생

분류: Imported from prior work. 공개 재생에서도 기존 J0337 여섯 단면과 비음수 beta 판정은 유지된다. 열 완화의 4.0e-9 / 3.3e-7은 열 시간 1e13 s까지 계산한 층의 부분합이며, 전체 지연 상한으로 해석하지 않는다.

분류: Proven. 같은 Debye 껍질 모형에서 미계산 총 세기 T_tail의 잔여항은 T_tail/(omega*tau_cut) 이하이다. 분류: Conjectural. T_tail 자체의 보장은 없으며, 모형 전체의 열 지연 인증과 전체 물리 목표는 완료하지 않았다. 실패한 단계는 부분합에서 전체 지연으로의 승격이고, 빠진 최소 조건은 깊은 층의 총 완화 세기 또는 그 잔여항의 보장이다. 기존 기록과 기준을 보존하고 원고 범위를 정정했다. [근거 295](../notes/REQUEST295_SUBMISSION_CORRECTIONS_KO.md).
