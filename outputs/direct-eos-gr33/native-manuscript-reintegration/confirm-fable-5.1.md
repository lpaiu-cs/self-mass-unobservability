## 반영 확인 재심 보고서 — fable5.1 (개정 5판 `docs/white-dwarf-free-fall-charge-section.md`, 2026-09-27)

### 1. 종합 권고

**수락.** 4판 재심에서 낸 주요 1건·경미 5건은 모두 해소되었다. 5판에서 새로 들어간 서술(순간 응답의 Section 3 귀속, 1.24e−9 a_p²와 |a_p|≳0.5 문턱, §4.3 고정 동반성 축약과 강성 조건, 𝒬_c/𝒬 분리, Debye 지연 상한, 두 가정을 담은 초록·요약·§6 문장, manifest 파일명과 판정 필드 대체 문장)은 원고 §4.3 원문, 저장 JSON, 인용된 노트 행과 대조해 모두 정확했다. 정오표 7–15는 인용 행·manifest 필드·저장 표본에서 전부 재현된다. 차단·주요 지적은 없다. 남는 것은 선택적 문구 조정 3건(§4 참조)이며, 어느 것도 주장을 바꾸지 않으므로 통합 커밋에서 처리하거나 그대로 두어도 된다. 통합 전제는 이전과 같다: bib 3항목 복원, README·revision manifest의 PDF 이전 판 기록, Pandoc 재생성, `verify_unified_paper.py` 실행, 문구를 바꾸면 초안 해시 재결속.

### 2. 이전 지적(4판 재심)별 판정표

| 이전 지적 | 판정 | 이유 |
|---|---|---|
| 주요 1 순간 응답의 §4.3 교차 참조와 척도 누락 | **해소** | 되먹임 조건 귀속 문장을 지우고 "a zero-lag coefficient, not an internal state"로 적었다. §4.3 원문(275행 "hold the companion charges fixed"; 315행 "Neglect requires this small relative to κ")과 정확히 맞는 두 문장으로 나눴다: 고정 동반성 축약이 동상 변조를 빠뜨린다는 것과, 같은 감수율이 강성을 옮기지만 그 조건은 평가하지 않았다는 것. 척도 1.24e−9 a_p²와 문턱 |a_p|≳0.5는 재현된다(§4). 해석 문장은 Conjectural로 옮겼고, 136행이 "either part in these two non-common channels"로 순간 항을 포함한다. 내가 제안한 되먹임 비율 ~4e−13은 싣지 않았는데, 이 절에서 κ가 정해지지 않는다는 이유는 타당하다. |
| 경미 2 초록·요약의 가정 수 | **해소** | 초록 "a thermal relaxation strength no larger than that term and a single relaxation time"; 요약 "the two stated assumptions on the thermal relaxation"; 가정 목록 "with a single relaxation time". |
| 경미 3 정오표 누락 (REQUEST274 공명·자전, REQUEST276 39행, 원장 2563행) | **해소** | 정오표 10(274 48·54·58·60행), 11(276 37·39·43·47행: 일반 수동 상한 4.4e−12, 필요조건 1.56e3, "도달할 수 없다", 중성자별 귀속 철회), 15(원장 2563행). 인용 행을 모두 열어 원문이 있음을 확인했다. |
| 경미 4 Born 비교의 노름 표현 | **해소** | "the stored maximum ratio of the first-return Born field to the incident field" — 저장 필드 `maximum_first_return_over_incident`=1.738e−10(phase267-born 로그)과 일치. |
| 경미 5 116행 표지 | **해소** | 쌍 인자 대수를 Proven, J0337 척도 비교를 Conjectural로 분리. |
| 경미 6 통합 체크리스트 | **해소** | 현 초안 SHA-256 7112a679…와 REQUEST285 2615c700…이 `paper/revision-manifest.json`(12138·12229행)에 묶였고, `native_section_revision_5` 항목이 "pending confirmation review"와 `superseded_verdict_fields`를 기록한다. 1.24e−9는 "at most"로 유지. |
| (4판 미확인) REQUEST278 표 불일치 원인 | **해소** | 정오표 13: 옛 값 −2.7e−43·−7.5e−41·−2.4e−39가 저장 표본 0.0982·0.2282·0.3482 s의 −2.672e−43·−7.526e−41·−2.376e−39와 정확히 같음을 transit.json에서 확인. 원인(직전 표본 판독)이 특정되어 4판 정오표 2의 "확인 못함"이 닫혔다. |

### 3. 새 지적

**차단: 없음. 주요: 없음.**

**경미 (모두 선택 사항, 주장 불변)**

1. **전체 별 모형의 −0.81% 비교 대상** (69·71행). 5판이 "it evaluates 𝒬_c … the mass term is not applied to its history"라고 명시했는데, transit.json `validation`은 모형의 𝒬_c(−2.3162e−51)를 질량 정규화된 연속 극한 𝒬(−2.33502e−51)와 비교해 −0.8058%를 얻는다. 𝒬_c끼리(−2.33497e−51) 비교하면 −0.8036%로 반올림이 −0.80%가 된다. 문장 자체("reproduces the endpoint readout", 곧 𝒬(t_e))는 참이고 차이는 2.2e−5 수준이라 무해하다. 원하면 "to −0.81% (−0.80% against 𝒬_c)"로 적는다.
2. **데이터 가용성 문장의 일반화** (174행). "The later steps read runtime arrays of the endpoint run … cannot be re-run from the deposit alone"은 270–273·277–279에는 맞지만(스크립트가 `readout268-quad64-work/gr/*.npz`를 읽고, 이 배열들은 `large_arrays_sha256`에 해시만 있음), 276(phase276-closure.py는 저장 JSON만 읽음)과 REQUEST284 산술은 저장물만으로 재실행된다. REQUEST285 표 자체가 이를 "없음(저장 JSON)"으로 바르게 적고 있으므로 "Most later steps" 또는 "The steps other than the Cassini rescaling and the pair-channel arithmetic"로 좁히면 표와 일치한다. 재현성을 과소 진술하는 방향이라 해롭지는 않다.
3. **문구 두 곳.** (a) 초록 "no first-order permanent charge in the limit of damped modes"는 본문의 "as t→∞ if all its modes are damped"보다 모호하다; "as t→∞ for damped modes"가 본문과 같다. (b) 125행 "at most about 1.24e−9 a_p²"는 교차항 8.2e−10|a_p||a_o|≤2.9e−12|a_p|를 떨어뜨린 값이다. |a_p|=0.5에서 0.5%라 "about"으로 충분하지만, "for |a_o|≃|α_0|"를 붙이면 근거가 드러난다.

### 4. 직접 확인한 수치

- **순간 응답**: 4×3.08e−10=1.232e−9→1.24e−9(올림) ✓; 4×2.04e−10=8.16e−10→8.2e−10 ✓. 문턱: (2.8e−10/1.24e−9)^½=0.475→"≳0.5" ✓ (4.1e−10 기준 0.575). 교차항 |a_p|=0.5: 1.45e−12 vs 3.1e−10.
- **§4.3 원문**: 275행 고정 동반성 전하, Q_p·m_p 표기; 315행 "feedback stiffness shift −Σ_j C_j/r_pj²", "Neglect requires this small relative to κ" ✓. a_p=Q_p/m_p는 Δ_pj=Q_pQ_j/(m_pm_j)와 일치.
- **𝒬_c/𝒬**: 𝒬=𝒬_c−a_i^(0)δm/m는 4판 식 −(Ψ_out+α₀δm)/m과 a_i^(0)=α₀에서 동일 ✓. 창 식·이력·최대값 3.5e−36이 𝒬_c로 표기 ✓. 모형 q_T=−2.3162e−51→−2.3e−51 ✓; 정적 창 +3.622e−54, 비 639→1/640 ✓.
- **지연 상한**: max_x x/(1+x²)=0.5 ✓. astra 반례 H=s/2+(s/2)/(1+iωτ)의 ωτ=1 직교 성분 s/4 ✓ — 두 가정을 만족하며 지연이 남으므로 정오표 7의 철회가 옳다.
- **구조 척도**: 8.83658e−9×0.781259×(3.07608e−10+2.03106e−10×3.53556e−3)=2.1286e−18→"≤2.13e−18" ✓; ×3.53556e−3=7.5257e−21 ✓. 2차 조석 Cassini 재척도 3.7e−19×0.884=3.27e−19→응답 노트 3.3e−19 ✓; "scaling as φ_∞"는 판독 가중치 α_s(φ₀)∝φ_∞만 들어가므로 타당.
- **변조 식의 차수**: phase274-tides.py 94행 `phi_o*a_in*mp/(mp+mc)/a_out`→"leading order … in a_in/a_out" ✓; e_in=6.92e−4, e_out=0.0354; 질량 1.4378/0.19754/0.4101 M_⊙(§5.9 세트) ✓.
- **자기 결합 정정**: 0.016×5.495e5/(4×8.988e20×1e−3×7.3975e−8)=3.31e−8 ✓(φ₀′ 항 3.3e−8과 같음 — 이중 계산 의심 타당); 0.048→9.92e−8→"order 1e−7" ✓.
- **정오표 9–13**: 2.339e−51/1.415e−54=1653 ✓; 8.8366e−9×(3.076e−10+2.031e−10)=4.51e−18 ✓, ×3.076e−10만=2.72e−18 ✓, 1.7e−9/2.72e−18=6.25e8→6.3e8 ✓; 1.7e−9/(3.53556e−3×3.07608e−10)=1563.1 ✓; 3.53556e−3×4×3.07608e−10=4.350e−12 ✓; 외측 항 제외 7.508e−21 ✓. transit.json 표본 격자: 0.25 ms×80(20 ms까지), 이후 0.0202 s부터 2 ms 간격 ✓; 보간 0.1 s −2.972e−43, 0.23 s −8.04e−41, 0.35 s −2.654e−39 ✓.
- **정오표 14**: 다섯 manifest의 verdict/classification 문자열, `required_kappa_lag_times_alpha_p`=1563.125, `full_beta_relaxing_delta_bound`=4.3503e−12, `freefall_delta_bound_cassini`=7.508e−21, `infinity_tail_bound`=4.22e−6, `j0337_closure`, tidal manifest `assumptions[1]` "white-dwarf spin period longer than about 45 min" — 모두 존재 ✓.
- **정오표 15·유지 문서**: 원장 2563행에 "정적 계수로 붕괴한다(no-go 경계…)" 존재 ✓; 단계285 항목이 6개 문서(adiabatic-limit, dynamic-charge-completion, failure-ledger-dynamic-chi, model-definition, nonadiabatic-regime, observable-targets)에 있음 ✓.
- **재현 경로 표**: 표의 주요 스크립트 이름 전부가 해당 manifest에 해시로 기록됨 ✓; `large_arrays_sha256`가 quad(14항목: source-64·field-source-64·photon seg·balanced-initial-state 등)와 eos(gas.so·levels.so) manifest에 있음 ✓; `def-native-boundary-layer/geometry.npz`는 저장소에 있음; `scan285-inputs.py` 존재·결속 ✓.
- **구조 검사**: §4.6 목록 10개·항목 47개·식 4개·라벨 4개 모두 참조(응답 노트 "확인"과 일치) ✓. 무표지 문단은 수식 뒤 이어지는 줄과 목록 도입문뿐.

### 5. 확인하지 못한 것

- Pandoc 변환과 `verify_unified_paper.py` 실행(통합 후 확인 사항; 내 구조 검사는 텍스트 계수로 대신함).
- phase274-tides.py·phase275·phase279가 단계268 런타임 배경을 읽는다는 표의 진술(스크립트에 경로 문자열이 없어 문자열 검색으로는 미확인; 노트 자체가 이 한계를 적음).
- 되먹임 강성 조건(본문이 미평가로 명시), 열 완화 세기 ≤ 𝒮_struct 가정, 깊은 층 열 상태의 결합(모두 계산되지 않은 가정으로 표기됨).
- 외측 백색왜성 |a_o|≃|α₀|의 수치적 근거(약한 장 표준 결과로 인용만 됨).
- 문헌 원문(Kaplan 2014, Ransom 2014, Bertotti 2003, Damour–Esposito-Farèse 1992)의 웹 대조.
