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
