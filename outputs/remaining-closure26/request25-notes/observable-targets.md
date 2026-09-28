# Observable targets — unified paper

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
