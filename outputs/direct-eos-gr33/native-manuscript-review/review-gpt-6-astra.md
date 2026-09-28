**1. 종합 권고: 주요 수정**

끝점 전하의 부호와 주요 수치는 저장 산출물에 대체로 충실합니다. 그러나 **조건부 계산을 궤도 시간척도의 no-go와 중성자별 기원이라는 결론으로 연결하는 논증**에는 중요한 공백이 있으며, §5의 관측량 정의와도 충돌합니다.

대상 HEAD는 `a322c6462`와 일치했습니다. 다른 심사 의견은 열람하지 않았고, 파일 수정이나 2분 이상의 계산은 하지 않았습니다.

**2. 지적 사항 — 심각도순**

아래에서 `O/`는 `outputs/direct-eos-gr33/`입니다.

**[차단] ① §5의 계수 구간을 백색왜성 전하의 관측 감도로 직접 사용할 수 없습니다.**

- **위치:** [원고 §4.6, 404–411행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/paper/manuscript.md:404), §6 추가 문단.
- **문제:** §5의 \(U\)는 지정된 세 carrier 구동에 대한 **계수 \(\beta\)**의 구간입니다. §5.4 자체가 이것을 SEP 진동 진폭으로 읽으면 안 된다고 설명합니다. 그런데 §4.6은 \(|\delta\Delta|\)와 \(U\)를 직접 비교합니다.
- **Proven:** 힘의 변화도 다릅니다. 내측 백색왜성의 전하만 변하면 \((pi,po,io)\) 쌍의 스칼라 결합 변화는 각각 \((\alpha_p\delta\alpha_i,0,\alpha_o\delta\alpha_i)\)입니다. §4.2–4.3의 저장 template는 펄서–동반성 두 쌍에 공통 변조를 가합니다. 외측 별에 대한 차등가속도 계수 \(-\alpha_o\delta\alpha_i\)만 일치시켜서는 timing template가 같아지지 않습니다.
- **근거:** 원고 260–315행, 474행, 484–494행, 518–529행. `O/native-final-closure/scripts/phase276-closure.py`는 `limit=1.7e-9`를 직접 대입하며 이 대응을 계산하지 않습니다.
- **수정:** 현재 비교를 조건부 크기 비교로 낮추고, 검출 가능성이나 기원 판정으로 사용하지 마십시오. 그 결론을 유지하려면 백색왜성 전하 변화의 힘·위상·진폭을 실제 추론 모형에 연결해야 합니다.

같은 산술 비교를 유지하더라도, \(1.68\times10^{-9}\)가 **71방향·\(K=10\)** 시나리오임을 밝혀야 합니다. 90방향 값은 \(3.534\times10^{-9}\), \(K\simeq934\) 시나리오는 \(1.57\times10^{-7}\)입니다. 필요한 \(|\kappa_{\rm lag}\alpha_p|\)는 각각 약 \(1.54\times10^3,\ 3.25\times10^3,\ 1.44\times10^5\)가 됩니다. 이 산술은 관측량 대응 문제를 해결하지는 않습니다.

**[주요] ② 양의 감쇠·정적 감수율만으로 A4 유지와 모든 지연의 배제를 결론낼 수 없습니다.**

- **위치:** §4.6 386행, 398–413행; 초록 추가문; §6 684행.
- **Proven:** 양의 감쇠는 무한시간에서의 소멸 조건이지, 궤도 시간척도보다 빠른 이완의 증명이 아닙니다. \(I\ddot x+\Gamma\dot x+kx=F\)는 같은 양의 \(I,k\)에서도 큰 \(\Gamma\)에 대해 느린 시간 \(\Gamma/k\)를 가집니다. 원고 §4.5가 이미 설명하는 경계입니다.
- **Proven:** 수동성만으로 주파수별 감수율을 정적 감수율 이하로 묶을 수도 없습니다. 수동적인 감쇠 진동자는 공명에서 정적 응답을 초과합니다. 단일 완화 또는 양의 상호적 완화 스펙트럼의 상계를 적용하려면 그 구조를 입증해야 합니다.
- **근거:** `O/native-tidal-photosphere/phase274/tides.json`의 `structural.kappa_struct_max`는 **단열 정적 계산**입니다. `phase274-tides.py`에는 비단열 전달함수나 감쇠율 계산이 없습니다. REQUEST274는 단일 완화 상한을 가정으로 남깁니다.
- **수정:** “정적 감수율의 계산”, “선택한 완화족에서의 진폭 상한”, “실제 별의 이완시간”을 분리하십시오. \(\alpha_p\)의 허용 범위도 명시해야 합니다. 현재 근거로는 이 자유낙하 기작의 조건부 작은 응답을 제시할 수 있지만, 지연이 반드시 중성자별에 있어야 한다는 필요조건은 성립하지 않습니다.

또한 **Counterexample candidate:** 계산한 것은 공간적으로 분포한 반경 변위와 다수 모드입니다. §3의 단일 상태·단일 pole로 축약한 결과는 아닙니다. 펄스 폭 \(cD\simeq515\) km는 별 반지름 약 \(69{,}100\) km보다 훨씬 작으므로, 짧은 펄스 결과가 점입자 EFT의 적용 영역 밖이라는 REQUEST270·273의 설명을 본문에도 복원해야 합니다.

**[주요] ③ 조석 감쇠 깊이의 근사 적용 조건이 충족되지 않습니다.**

- **위치:** §4.6 401행.
- **문제:** \(8.9\times10^6\)이라는 적분값 자체는 재현 가능한 저장값입니다. 그러나 사용한 준단열 감쇠식으로 이산 공명 부재와 자전 45분 경계를 확정하기 어렵습니다.
- **Counterexample candidate:** 저장 범위의 가장 보수적인 조합으로도
  \[
  \frac{Kk_r^2}{\omega}
  \ge
  \frac{6K_{\min}N_{\min}^2}{\omega^3R^2}
  \simeq1.54\times10^3.
  \]
  따라서 열확산을 작은 교란으로 다루는 조건 \(\ll1\)과 크게 어긋납니다.
- **근거:** `O/native-tidal-photosphere/phase274/tides.json`의 `thermal.K_thick`, `thermal.N_thick`, `R`, `drives["2 n_in (non-rotating tide)"].omega`; `phase274-tides.py:53–62`. 해당 감쇠식의 준단열 전제는 [Ahuir·Mathis·Amard의 원 논문](https://doi.org/10.1051/0004-6361/202040174)에도 명시됩니다.
- **수정:** 이를 강한 열확산의 지표로 제시하고, 공명 부재·45분 경계는 미확정으로 남기거나 적절한 비단열 문제로 검증하십시오. 구대칭·비회전 배경에서의 선형 \(l=0\) 선택 규칙은 **Proven**으로 유지할 수 있습니다.

**[주요] ④ GR·곡률·비선형 보정의 크기 추정이 증명된 전하 오차 상계로 승격되어 있습니다.**

- **위치:** §4.6 388–394행.
- **문제:** \(2\pi G\rho_{\max}T^2\)는 뉴턴 반경 자기중력 되먹임의 변위 노름 추정입니다. 이를 전체 물질·광자·계량 결합 사상의 수축 상수로 사용하고, 다시 상대 전하 오차로 옮기는 단계가 없습니다.
- **근거:** `O/native-closure-transit/scripts/phase277-bounds.py:18–32`는 밀도를 보간하고 \(L\), \(L/(1-L)\)를 계산합니다. 완전한 결합 연산자나 전하 판독의 오차 증폭은 계산하지 않습니다. REQUEST277 표도 적용을 **Conjectural**로 구분합니다.
- **Proven:** 평탄한 외부에서 들어오는 성분을 제거한 구형 파동의 \(r\delta\varphi\)가 지연시각만의 함수라는 명제는 타당합니다. 그러나 \(M/r_{\rm out}\)이라는 차수 추정만으로 곡률 꼬리의 계수 1짜리 상대 상계가 증명되지는 않습니다.
- **수정:** 보정별로 정확한 항등식, 조건부 부등식, 수치 대입, 물리적 적용 추정을 분리하십시오. \(\alpha=\beta\varphi\)의 선형성도 특정 이차항의 비교를 정당화할 뿐, 상쇄를 거친 전체 전하의 상대 오차까지 자동으로 보증하지 않습니다.

**[주요] ⑤ 관측 대기 보정의 방식과 수렴 범위가 생략되어 있습니다.**

- **위치:** §4.6 396행.
- **Counterexample candidate:** 회색 재구성의 0.25%는 **전하 핵 가중 오차**입니다. 전하 층 최대 오차는 0.376%입니다. 비회색 \(5.5\times10^{-5}\)는 **\(\tau_R\le100\)**에서의 흐름 오차입니다.
- **근거:** `O/native-atmosphere-reconstruction/phase279/combined.json`의 `gray_validation.kernel_weighted`, `max_abs_rel_charge_layers`, `cases["Kaplan central"].nongray_flux_error_tau_le_100`. `nongray.log`에서는 모든 경우 1500회까지 실행되었고 원래 평가 범위의 흐름 오차는 약 \(1.2–2.0\times10^{-3}\)입니다.
- **Conjectural:** 최종 ×6.9는 표 불투명도의 구면 회색 밀도에, 단순 H·He 불투명도의 평행평판 비회색/회색 밀도비를 곱하고 **기존 전하 핵을 고정하여** 얻은 보정입니다. 같은 관측 대기에서 물질·광자·GR을 다시 결합한 해는 아닙니다. `phase279-combine.py:4–5,23–24`도 이 합성을 Conjectural로 선언합니다.
- **수정:** 합성 절차, 고정한 핵, 깊이 원점, 제한된 수렴 구간을 명시하십시오. 표 불투명도와 단순 불투명도의 Rosseland 평균 차이 0.55–0.74도 남겨야 합니다. 1σ·2σ 범위는 주로 \(\log g\)를 개별 변동한 결과이며, 공동 확률분포에서 산출한 신뢰구간처럼 표현하면 안 됩니다.

**[주요] ⑥ 주장 표지와 초록·논의의 확정성이 본문 근거보다 강합니다.**

- **위치:** 369–376행, 386–394행, 413행, 초록, §6 추가문.
- **수정할 표지:**

| 대상 | 적절한 구분 |
|---|---|
| 자유낙하 닫힌 식이 명시된 ODE를 만족 | **Proven** 유지 |
| 그 식과 결합 계산의 \(5.2\times10^{-5}\) 일치 | 별도 **Counterexample candidate** |
| 양의 강성을 가진 유한 선형계의 영평균, 양의 모드 감쇠하 소멸 | 조건을 명시하면 **Proven** |
| 실제 별에 그 조건이 적용됨 | **Conjectural** |
| \(q/\eta\)의 배경장에 대한 짝함수 성질 | \(\eta\)의 일차 응답에 한정하여 **Proven** |
| 동적 성분의 \(\varphi_\infty^2\) 법칙 | 약한 배경장 전개의 선도항임을 명시; 현재 물리적 적용은 **Conjectural** |
| 보정 전체가 작다는 포괄 문장 | **Proven** 제거 |
| 마지막 경계 문단 | 허용되지 않은 **Boundary statement** 대신 **Conjectural** |

초록과 논의에서는 **선언 모형에서 펄스가 유도한 일차 동적 전하 변화**가 사라진다고 써야 합니다. 배경 스칼라 전하 자체가 0이라는 뜻은 아닙니다. 전체 별 이력의 끝점 일치 −0.81%는 유용하지만, 이후 모든 시각의 결합 모형 정확도나 비단열 안정성을 검증한 것은 아닙니다.

**[경미] ⑦ 통과 시간과 정적 창의 부호가 부정확합니다.**

- **위치:** 380–383행.
- **Counterexample candidate:** 중심 도달은 **0.230497 s**, 출사는 **0.462712 s**입니다. 0.35 s는 이미 나가는 통과이며, 음의 부호가 유지되는 구간과 입사 구간을 혼동했습니다.
- 정적 창은 \(+3.62246\times10^{-54}\), 지연 전하는 음수입니다. “1/640”은 절댓값 비이며 부호는 반대입니다.
- **근거:** `O/native-closure-transit/phase278/transit.json`: `center_arrival_s`, `exit_s`, `validation.static_window_T`, `validation.q_T`.
- **수정:** 입사·출사·부호 유지 구간을 구분하고 절댓값 비라고 적으십시오. 끝점에서 중심 부근까지의 증가는 약 10.5자릿수이고, 약 12자릿수는 0.35 s까지의 증가에 해당합니다.

**[경미] ⑧ 단위·기호·대표값의 정의를 보완해야 합니다.**

- **위치:** 358–376행, 407–411행.
- **Proven:** \(U=r\delta\varphi\)라면 \(q=-U/M\)에서 \(M\)은 기하학적 길이 \(GM_{\rm phys}/c^2\)여야 합니다. 앞의 태양질량 단위 \(M\)과 구분해야 합니다.
- \(r_0\)와 \(R_d\), \(a,B,x,u,\delta M_e\)를 정의하십시오. \(g(s)\)는 \(0<s<1\) 밖에서 0임을 써야 의도한 compact \(C^3\) 펄스가 됩니다.
- \(\beta=-4\)는 §3·§5의 응답 계수 \(\beta\)와 다릅니다. \(\beta_{\rm ST}\) 등으로 분리하는 편이 안전합니다.
- **Counterexample candidate:** 대표 전하 \(-2.335\times10^{-51}\)는 직접 4배 출력이 아니라 Richardson 외삽에 질량 정규화를 적용한 값입니다. 정규화 계수는 이전 1배 실행에서 가져왔다는 사실도 밝혀야 합니다.
- Cassini 적용값 \(|\alpha_0|\le0.00353556\)는 \(\varphi_\infty\le8.83889\times10^{-4}\)에 해당합니다. 계산 기준점 \(10^{-3}\)과 관측 상한에 맞춘 재척도를 구분하십시오.

**[경미] ⑨ 생성 소스는 맞지만 새 계산의 재현 안내가 부족합니다.**

[main.tex의 해당 절](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/paper/main.tex:741)에는 네 목록, 세 equation label, 새 인용들이 보존되어 있습니다. 확인한 범위에서 Markdown→TeX 내용 누락은 없습니다.

다만 README의 실행 안내는 주로 기존 원고 검증과 Request12용입니다. 새 절에는 개별 manifest·실행 스크립트·외부 입력의 대응표가 필요합니다. `verify_unified_paper.py`의 통과는 새 §4.6 물리 논증의 검증을 뜻하지 않습니다. PDF·제출 ZIP·`main.bbl`이 이전 판이라는 설명은 적절히 남아 있습니다.

**3. 직접 확인한 수치 목록**

아래 계산값의 분류는 **Counterexample candidate**입니다. “일치”는 산출물과의 일치를 뜻하며 물리적 상계의 증명을 뜻하지 않습니다.

| 원고 수치 | 대조 결과 | 근거 파일·JSON 키 |
|---|---|---|
| 관측 질량 0.1975, \(T_{\rm eff}=15800\pm100\), \(\log g=5.82\pm0.05\), \(R=0.091\pm0.005\) | 일치 | **Imported from prior work:** [Kaplan 원문 초록](https://arxiv.org/abs/1402.0407) |
| 배경 \(M\simeq0.198M_\odot,\ R=6.91\times10^9\) cm, \(T_{\rm eff}=15898\) K, \(\log g=5.740\) | 일치 | `native-tidal-photosphere/phase274/tides.json`: `M_sun`, `R`, `thermal.Teff`; `combined.json`: `cases.declared.logg` |
| \(\beta=-4,\varphi_\infty=10^{-3},\eta=10^{-30}\) | 선언 입력과 일치 | `phase274-tides.py`, `phase277-bounds.py`; `bounds.json.eta` |
| \(D=1.717\) ms, \(T=2D\) | 일치: 1.71721556, 3.43443112 ms | `native-closure-transit/phase277/bounds.json`: `D`, `T` |
| 전하 \(-2.335\times10^{-51}\), 정규화 \(-2.335\times10^{-15}\) | 일치, 외삽·질량 정규화 값 | 같은 파일 `continuum.mass_normalized`, `per_eta_phi2` |
| 관측 차수 2.00, 2배→4배 −1.2% | 일치: 1.997469, −1.173956% | `native-quad-refined-primary/result.json`: `rows["64"]` |
| 직접 기하 몫 \(10^{-8}\), 바리온 몫 0.999996 | 일치 | `native-charge-mechanism/result.json`: `composition_4x_T` |
| 지나간 층의 몫 99.9% | REQUEST270 기록과 일치; 원천 배열 재계산은 미실시 | REQUEST270 25행 |
| PL-off 변화 최대 \(2.2\times10^{-10}\) | 일치: \(2.17583\times10^{-10}\) | `native-eos-sensitivity/result.json.max_abs_relative_change` |
| 자유낙하 외삽값 차이 \(5.2\times10^{-5}\) | 일치: \(5.23309\times10^{-5}\) | `native-charge-mechanism/result.json.freefall.limits_relative` |
| 부호 여유 0.9988, 중심 깊이 363 km | 일치: 0.9987908, 363.402 km | `native-structure-eft-boundary/result.json.background_structure` |
| 광학깊이 1.06 | 일치: 1.05739 | `native-final-closure/closure.json.accepted.photosphere_tau` |
| 정적 창 \(3.6\times10^{-54}\), 1/640 | 크기는 일치; 부호 설명 불일치 | `transit.json.validation`; 실제 비 −0.00156397 |
| 끝점 재현 −0.81%, 변위 \(1.2\times10^{-6}\) | 일치: −0.805834%, \(1.17815\times10^{-6}\) | `transit.json.validation` |
| 입사 구간 0.35 s | 불일치: 중심 도달 0.230497 s | `transit.json.center_arrival_s` |
| 중심 부근 \(-7.5\times10^{-41}\) | REQUEST278 표와 일치; 저장 근접 표본은 약 \(-8.1\times10^{-41}\) | REQUEST278; `transit.json.transit` |
| 최대 \(3.5\times10^{-36}\), 입사 대비 \(1.5\times10^{-11}\) | 일치: \(|q|_{\max}=3.48940\times10^{-36}\) | `transit.q_max_abs`; 같은 파동 진폭으로 환산할 때 일치 |
| 주기 44–59 s, 진폭 \(10^{-44}\) | 일치: 주요 주기 44.62–59.25 s; RMS \(8.10957\times10^{-45}\) | `transit.json.relaxation` |
| \(\omega_0^2=1.25\times10^{-3}\), 최저 주기 178 s | 일치: 약 \(1.24876\times10^{-3}\), 177.803 s | `relaxation.omega_lowest`, `period_lowest` |
| GR \(4.5\times10^{-18}\), 곡률 \(4.2\times10^{-6}\), 비선형 \(1.0\times10^{-27}\), ADM \(6.8\times10^{-11}\) | 저장 산술과 일치; 상계 해석은 지적 ④ | `bounds.json`: `fixed_point`, `infinity`, `nonlinear`, `adm` |
| 회색 재현 0.25% | 가중 오차로 일치; 최대 오차라는 해석은 불일치 | `combined.json.gray_validation` |
| 회색 ×3.2, 1σ 1.6–6.2, 2σ 0.78–11.2 | 일치: 3.24512, 1.62512–6.17249, 0.776015–11.19981 | `combined.json.cases`, `one_sigma`, `two_sigma` |
| 비회색 ×2.1→×6.9, 1σ 3.9–11.7 | 일치: 2.11862→6.87519, 3.92020–11.67670 | 같은 파일 |
| 비회색 흐름 \(5.5\times10^{-5}\), \(T_0/T_{\rm eff}=0.64\) | 제한 구간에서 일치: \(5.50517\times10^{-5}\), 0.643682 | `cases["Kaplan central"]` |
| 긴 파장 힘 비 \(3.3\times10^{-8}\), 주파수비 제곱 \(1.6\times10^{-6}\) | 일치 | `native-structure-eft-boundary/result.json.static_eft` |
| \(k_2=3.3\times10^{-4}\), 감쇠 깊이 \(8.9\times10^6\) | 일치: 0.000328918, 8871071.93 | `tides.json.love`, `drives` |
| 자전 약 45분 경계 | 노트의 근사 추정과 일치; 검증된 경계는 아님 | REQUEST274 |
| \(\kappa_{\rm struct}\le8.8\times10^{-9}\), \(|\beta|\) 대비 \(2.2\times10^{-9}\) | 일치 | `tides.json.structural` |
| 구동 \(3.08\times10^{-10}|\alpha_p|\) | 내측 이심률 성분에 한해 일치 | `tides.json.j0337.dphi_modulation_per_unit_charge` |
| Cassini \(3.5\times10^{-3}\) | 명시한 2σ 변환과 일치: 0.003535556 | `native-final-closure/closure.json.closure`; [Cassini 측정 원문](https://doi.org/10.1038/nature01997) |
| 요구값 \(1.5\times10^3\), 상한 \(4.4\times10^{-12}\), \(7.5\times10^{-21}\) | 반올림 산술 일치; 관측 해석은 지적 ①–② | 같은 파일 `closure`; 첫 값은 원고의 1.68e−9로 재계산 |
| 핵 compactness 약 \(10^{-4}\) | REQUEST278의 근사 기록과 일치; 전체 배경에서 독립 재계산하지 않음 | REQUEST278 한계 문단 |

**Proven:** 자유낙하 식은 \(F_1'=g,\ F_2'=F_1\)을 사용해 두 번 미분하면 명시된 강제 ODE를 만족합니다. 구형 얇은 껍질의 평탄 시공간 지연 평균도 제시한 \(c/(2r_e)\) 창을 줍니다. 실제 별에 대한 적용 근사와 이 두 수학적 확인은 별개입니다.

새 참고문헌 세 항목의 서지는 원 출처와 일치했습니다. [Ransom 논문](https://www.nature.com/articles/nature12917)도 권·쪽·연도·DOI를 확인했습니다. 날짜 2026-09-27과 새 manifest 연결도 일치합니다.

**4. 확인하지 못한 것과 그 이유**

- **결합 진화·전체 별 통과·대기 계산의 독립 재실행:** 계산 제한 때문에 하지 않았습니다. 저장 JSON·생산 스크립트·로그·해시를 확인했습니다.
- **전체 별 이력의 공간·시간 수렴과 실제 비단열 감쇠율:** 제시된 끝점 대조와 모드 전파 대조만으로는 확인되지 않습니다.
- **완전한 GR 고정점·곡률 꼬리의 엄밀한 오차 상계:** 대응 증명이나 연산자 평가가 제시 자료에 없습니다.
- **관측 DA 조성에서의 완전한 비회색 결합 해:** 현재 자료는 혼합 조성의 대기 계산과 고정 핵 보정입니다.
- **PDF 배치·실제 인용 해소:** 현재 PDF가 대상 변경 이전 판이므로 판단하지 않았습니다. `main.tex` 소스 보존만 확인했습니다.
- **원격 저장소의 공개 가용성:** 로컬 revision과 산출물 연결만 검증했습니다.

해시 검사는 대상 manifest 10개와 그 아래 등록 산출물 523개 모두 통과했습니다. 종료 시 저장소 변경은 없었습니다.

