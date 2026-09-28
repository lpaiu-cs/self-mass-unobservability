## 확인 심사 보고서 — 개정판(4cf06e26f) 반영 확인

읽기 전용으로 작업했습니다. 저장소 파일을 만들거나 고치지 않았고, 커밋하지 않았으며, 2분을 넘는 계산은 하지 않았습니다. 확인 방법: 두 커밋의 diff 전문, Nutimo 릴리스 소스(`Parameters.cpp`, `Utilities.cpp`), 동결 npz에서 근점 각 재계산, 재계산 JSON 5개와 철회 폴더의 대응 파일 대조, `nuisance-intervals.csv`·`sep_limit_curve_10_8e.tsv` 대조, `verify_unified_paper.py` 실행(PASS, exit 0), PDF 3종 텍스트 추출(pdftotext)과 SM 그림 3·4 렌더링(pypdfium2), 제출 소스 zip 내용 대조.

---

### 1. 종합 권고

**경미 수정 후 수락.** 근점 규약 정정(B1)은 코드·동결 매개변수·문헌과 독립적으로 일치하고, 재계산된 27개 새 수치가 모두 근거 파일과 맞으며, 철회된 값은 본문·SM·초록·커버레터·검증기·소스 zip·PDF 어디에도 남아 있지 않습니다. 제가 낸 주요 지적 M1–M4는 해소되었고(M2는 부분 해소이나 수용 가능), 새로 생긴 문제는 모두 한두 문장으로 고칠 수 있는 서술 문제입니다. 차단 사유는 없습니다. 가장 먼저 손볼 것은 §5.5의 2일 지연 단면이 문턱을 0.02만큼 넘는 점을 정량적으로 방어하는 한 문장(N1)입니다.

---

### 2. 이전 지적별 판정표

| 지적 | 판정 | 이유 |
|---|---|---|
| **M1** 새로움·선행 문헌 | **해소** | §1 넷째 문단이 정리 1을 켤레 노드 실계수 보간(Mayo–Antoulas 2007)과 2K차 지속 여기(Ljung 1999), 식 7–8을 절단 Stieltjes 모멘트 문제의 2×2 Hankel 양성(Schmüdgen 2017), 식 9를 완화 스펙트럼 역문제의 악조건성(Honerkamp–Weese 1989)에 놓았고, 초록에 "specialize classical interpolation and moment results"를 넣었습니다. 새로움을 "formulation ... comparator classes ... drive phase gate ... application"으로 한정한 진술은 정확합니다. 커버레터도 같게 고쳤습니다. zip의 `main.bbl`에 17개 키가 모두 있습니다. |
| **M2** 여섯 계수 옴니버스 | **부분 해소, 수용 가능** | 명목 p≈0.012(제가 재계산: 0.01197), 물리 평면 최적 적합 후 통계량(12.84; 14.39–14.59), 후보 원인 네 가지, "the excess is not attributed"가 들어갔습니다. 여섯 계수 추정치·표준오차 표는 넣지 않았습니다(저자: 저장 기록에 없음). 초과가 이제 초록 수준의 진술이 되었으므로 SM에 여섯 계수와 z를 표로 넣는 것을 권하지만(계산은 수 초), 필수는 아닙니다. |
| **M3** 자기 분석 정정 명시 | **해소** | 초록·§1·§4·§5.1·§6에 "our earlier stored analysis", §1에 "unpublished and recorded in the repository". 비율은 "2.10 to 17.41 (Table 1 shows five lags)"로 고쳤고, 65개 지연에서 제가 재계산한 범위(2.1036 at 2 d, 17.415 at 379.4 d)와 일치합니다. K=1 2.34–27.5도 CSV의 여섯 지연에서 2.3404–27.5087로 일치합니다. |
| **M4(a)** 백색왜성 조건구 | 해소 | "Under the stated assumptions, which include linear response, |a_p|≤1 and an outer-companion charge near the Cassini limit", 층별 추정의 성격, 비공통 채널, 남은 열린 상태가 모두 복원되었습니다. |
| **M4(b)** §5.6 "only"·보증 범위 | 해소 | "only" 삭제, "does not certify every lag, large initial amplitudes, or coverage" 복원. 값은 규약 정정으로 0.999998–1.00742로 바뀌었고 JSON과 일치합니다. |
| **M4(c)** §5.5 끝점 | 해소(철회로) | 단면이 모두 비어 끝점 자체가 없어졌고, SM §4.6 척도 목록과 `docs/white-dwarf-*.md`에서도 그 줄이 제거되었습니다. |
| **M4(d)** REML 조건·0/8192 | 해소 | "the same fixed family that generates the stress", "A truncated pointwise K=10 interval also has a 0/8192 stress case". |
| **M4(e)** §4 세 문구 | 해소 | "signed-β fits of Section 5 admit a larger model", "Rescaling the unit-drive interval cannot replace them", "a specified microscopic theory may relate these parameters". |
| **M4(f)** §5.1 주의 두 개 | 해소 | "that period is an analysis domain, not a common period", "(SM Section 5.4)". |
| m1 양의 Debye 부류와 SM §4.6 연결 | **불채택, 수용 가능** | 층별 추정이 층 강도 절댓값의 합이라 양의 스펙트럼을 보장하지 않는다는 이유는 타당합니다. |
| m2 표지 정의 | 해소 | "a state that a static description would miss". |
| m3 "deprojected readout" | 해소(본문) | SM 정리 3에는 남아 있으나 SM은 이전판 유지가 선언된 정책입니다(N5). |
| m4 두 표본 등식 조건 | 해소 | 초록·§1·커버레터에 "recover it when equality holds". |
| m5 Damour–Esposito-Farèse | 해소 | §4 첫 문단에 인용. |
| m6 반송파 근접 정량화 | 부분 해소, 수용 가능 | "differ only by ω_out ... not three equally informative" 문장; 조건수는 §5.4에 있습니다. |
| m7 아홉 진폭·±1 조건 | 해소 | 둘 다 삽입. |
| m8 표 1 K=1 열 | 해소 | 전체 공간 K=1 열이 추가되었고, 절단 K=1 열이 "fails under omitted means"임을 §5.2가 밝힙니다. 새 열이 진짜 4821 원점 포락인지 확인했습니다(아래 5절). |
| m9 AI 절 | 부분 해소, 수용 가능 | 실행 경로("run as a Claude Code subagent"), 전체 원고 심사, "programs named below"로 수정. 도구 버전·기간은 기록이 없어 저자 판단으로 남김—수용. |
| m10 Repository 줄·소속 | 해소(본문) | 본문 제목 블록에서 제거. SM 제목 블록에는 남음(N5). 소속·ORCID는 알려진 미결. |
| m11 절 번호 하드코딩 | 불채택, 수용 가능 | PDF 제출에는 영향 없음. 렌더링에서 "??"·"[?]" 없음을 확인. |

---

### 3. B1 정정 확인 결과

**규약이 코드와 맞는가 — 확인.**
- `Utilities.cpp:160` `long double inversetrigo(long double cosv, long double sinv)`; 분기 `sinv>=0 and cosv<0 → acos(cosv)`이므로 (cos, sin) 순서로 각을 돌려줍니다.
- `Parameters.cpp:1652` `omp = inversetrigo(kappap/ei, etap/ei)`, `:1645` `omB = inversetrigo(kappaB/eo, etaB/eo)`. 따라서 κ = e cos ω, η = e sin ω(ELL1). 주석 처리된 `:1646` `tperio = (tascB + Po*omB/deuxpi)`도 본문의 t_asc = t_peri − Pϖ/2π와 일치합니다.

**문헌과 맞는가 — 확인(기억 기준, 원문 미열람).** 동결값 η_p = 6.937e-4 > 0, κ_p = −8.60e-5 < 0은 η = e sin ω로 읽을 때 ω가 2사분면(97°)이고, 이는 Ransom et al. 2014 표 1의 ω_I = 97.6182°, ω_O = 95.619493°와 같은 사분면입니다. 저자가 옮긴 (e sin ω, e cos ω) = (6.8567e-4, −9.171e-5), (3.5186e-2, −3.4621e-3)은 제 기억의 e_I = 6.9178e-4, e_O = 3.5356e-2와 정확히 재구성됩니다. 이전 읽기(−7.06°, −5.73°)는 문헌과 모순됩니다.

**재계산 값이 옳은가 — 확인.**
- npz에서 atan2(η, κ): ϖ_p = 1.69406454 rad = 97.0627°, ϖ_b = 1.67084424 rad = 95.7323°. 옛 값과의 합이 정확히 π/2(항등식 atan2(κ,η) = π/2 − atan2(η,κ))이므로, 보관 사전이 ϖ = π/2를 썼다는 설명과 위상차 0.12326821 = π/2 − ϖ_p, 0.10004792 = π/2 − ϖ_b, 차이항 π가 모두 맞습니다. 이 세 차이는 시간 원점과 무관합니다.
- 폐합: π + ϖ_p − ϖ_b = 3.16481295. JSON은 (−π, π]로 감은 −3.11837236을 저장하며 검증기가 `% 2π`로 비교합니다. 옛 폐합 3.11837236과 새 폐합이 2π를 기준으로 서로 음수인 것은 규약 교환의 필연적 결과이지 우연이 아닙니다.
- 위상 무관 값은 철회본과 상대 1e-9 이내로 같음: U_* = 1.74999156e-10, 정규화 진폭 (0.24394925, 0.69195682, 0.06409393), 생략 입력 RMS 3.37427 %, 차이항 위상 −2.0434, 조건수 98.187, 옴니버스 16.3525, 문턱 12.8241766, 보정 해시.
- 위상 의존 값: 단면 최소 통계량 12.8446(2 d), 14.3918, 14.5876, 14.5359, 14.4518, 14.4265 → 모두 `empty: true`, `joint_region_abs_beta_upper: null`. 비교기 최소 P1 0.005389, P2 6.286e-5, P3 1.168e-5, P4 4.627e-6(1/√ = 464.9 → "465"), E4 0.0311, P5 잔차 5.38e-30–5.28e-25. 과도 공적합 σ 배율 0.9999985–1.0074182, 반 간격 σ 배율 0.7858–0.8722. 본문·SM 서술과 모두 일치합니다.

**철회된 값이 남았는가 — 없음.** −0.12326821/−0.10004792(음수), 3.11837236, 4.10217, 8.31672, 0.08289, 589, 0.003035, 1.00096, 1.00523, 0.8680, 0.8785, 9.3e-30, 1.4e-24, "2.1--17.4", "Section 7", "deprojected", "listed below"를 `manuscript.md`, `supplement.md`, `main.tex`, 커버레터 tex, 평문 초록, 체크리스트, 검증기, zip의 `main.tex`, PDF 3종 텍스트(본문 13쪽, SM 30쪽, 커버레터 1쪽)에서 찾았습니다. 본문 계열에는 없습니다. SM에 남은 것은 이전판 유지 정책에 해당하는 표현뿐입니다(N5). 검증기는 근점을 npz의 atan2(η,κ)에, 폐합을 3.16481295에, 단면을 `empty`·`minimum_statistic > threshold`에, 새 열을 CSV에 묶었습니다.

**새 결론이 근거를 넘는가 — 대체로 넘지 않음.** "no relaxation is detected", "the excess is not attributed"는 근거 안입니다. "Every physical lag section is therefore empty"는 평가된 여섯 지연에서 성립하되 2일은 한계적입니다(N1). "Every relaxing response to this drive lies in the rejected plane"는 압축이 과합니다(N2).

---

### 4. 새 지적

**차단:** 없음.

**주요:** 없음.

**경미 (우선순위 순):**

- **N1. §5.5·초록 — 2일 단면의 한계적 기각.** 문턱 12.8241766은 네 조건의 95 % 순서통계량(12.2667, 12.8242, 12.6404, 12.4872; `simultaneous-calibration.json`)의 **최댓값**이고, 8192개 표본에서 χ²₆ 밀도 0.017을 쓰면 각 순서통계량의 몬테카를로 표준오차는 약 0.14입니다. 2일 최소 통계량 12.845는 문턱을 0.02(≈0.15 SE)만큼 넘습니다. 본문은 "just above the threshold"라고 쓰지만 정량화하지 않아, 심사자가 "12.84 대 12.82는 기각이 아니다"라고 할 여지가 있습니다. 실제로는 2일 값이 네 조건 각각의 분위수와 명목 12.59를 모두 넘고, 단면 검정(2-매개변수 평면 위 최소 대 6자유도 영역)은 보수적이므로 방어 가능합니다. **수정:** 한 문장 추가 — 네 조건 분위수 범위 12.27–12.82와 SE ≈ 0.14를 밝히고, 2일 기각이 동결 최댓값에 대해서는 한계적이지만 각 보정 분위수와 명목값에 대해서는 성립함을 적을 것. "Every physical lag section is therefore empty"에 "evaluated"를 넣을 것(SM §5.9는 이미 "every evaluated").
- **N2. §5.5 Conjectural·§6.** "Every relaxing response to this drive lies in the rejected plane"는 지연마다 평면이 하나씩 있고 여섯 지연만 평가했으므로 "in a rejected plane at each evaluated lag"로, "the excess is not a relaxation signal"은 "not a relaxation response to the leading drive at those lags"로 좁힐 것. §6 "an excess that the leading physical drive, with or without relaxation, does not reproduce"에 "at the six evaluated lags"를 붙일 것.
- **N3. 초록·커버레터 "factors of 2.1--27.5".** 하한은 K=10(65개 지연), 상한은 K=1(6개 지연)에서 온 혼합 범위입니다. §5.2가 분해해 주므로 오류는 아니나, "2.1--17.4 at K=10 and 2.3--27.5 at K=1"처럼 쓰면 오독이 없습니다.
- **N4. 표 1 새 열과 §5.2.** 새 "Full, K=1" 열은 대각 가중 구성이며, 표 2는 대각 가중이 상관 Fourier 잡음에서 0.5825까지 떨어짐을 보입니다. §5.2는 절단 K=1 열의 실패(생략 평균)만 언급하므로 "both K=1 columns use diagonal weighting, which fails under correlated extra Fourier power (Table 2)"를 덧붙일 것. 그림 1에는 전체 공간 K=1 곡선이 없지만 표에 있으니 선택 사항입니다.
- **N5. SM 일관성(사소).** SM 제목 블록에 "Repository:" 줄이 남았고, 정리 3의 "deprojected readout", §3.5의 "Instantaneous mass at zero time", §5.3의 "factor 2.1--17.4"(65개 지연 한정어 없음)가 남았습니다. SM 초록의 "unchanged except for corrections" 정책으로 설명되지만, Repository 줄은 본문과 맞추어 빼는 편이 낫습니다.
- **N6. 등록 모의실험 행 유지 + 데이터 단면 교체 — 타당, 한 구절 권고.** 6계수 영역·문턱·옴니버스 통계량은 θ 공간에서 정의되어 구동 위상과 무관하므로 단면만 재계산한 처리는 논리적으로 옳습니다. 재실행 적중 수 차이(최대 22/8192)는 표준오차 ≈ 20 안이고, JSON `data_revision`·README·SM 데이터 가용성 절에 기록되어 있습니다. SM §5.9에 "영역·문턱·옴니버스는 위상에 의존하지 않으므로 단면만 재계산했다"는 구절을 명시하면 합성 파일에 대한 의문을 선제합니다. 스레드 수에 따른 고유기저 비유일성 설명은 그럴듯하나 제가 재현하지는 않았습니다.
- **N7. 서론 문헌 문단 — 문제 없음.** 네 배치 모두 정확하고 서지(LAA 425:634–662; Prentice Hall 2판 1999; GTM 277; Macromolecules 22:4372; MNRAS 326:274)는 제 기억과 일치합니다. Ljung의 2K 셈은 이산시간 (0, π) 안의 서로 다른 주파수 조건이 붙지만 "matches the count"라는 표현 수준에서는 문제 없습니다.

---

### 5. 직접 확인한 수치

- **코드:** `Parameters.cpp` 1645·1652행, `Utilities.cpp` 160행과 분기 구조(위 3절).
- **npz:** η_p = 6.9373e-4, κ_p = −8.5951e-5, η_b = 0.0351147, κ_b = −0.0035249, t_asc,p = −575.138 d, t_asc,b = −262.102 d, P_in = 1.629399 d, P_out = 327.2551 d; e_in = 6.990e-4, e_out = 0.035291; ϖ_p = 1.69406454 (97.0627°), ϖ_b = 1.67084424 (95.7323°); 옛 값 −0.12326821, −0.10004792; 합 = π/2.
- **physical-matching.json:** 위상차 −0.12326821/−0.10004792/π, 보관 폐합 3.6e-13, 14개 검사, `common_origin_cannot_fix_closure: true`, `phase_convention_source`가 코드 경로를 가리킴.
- **corrected-physical-drive.json:** U_* = 1.7499915637e-10, 진폭·RMS 철회본과 동일, 위상 (−1.8468, −2.9218, −2.0434) vs 철회본 (−0.0295, −1.1509, −2.0434).
- **simultaneous-validation.json:** 문턱 12.824176622705194(= a=0.25 순서통계량 최댓값), 명목 12.5916, 등록 적중 7796/7818/7805/7816(+형상 스트레스 7858) = 0.951660/0.954346/0.952759/0.954102, 재실행 7799/7818/7803/7804/7836, 옴니버스 16.3525(명목 p 0.01197), 단면 6개 모두 비어 있음(최소 12.8446, 14.3918–14.5876), 철회본은 6개 모두 포함(11.04–11.78; 끝점 4.10217e-10, 8.31672e-9).
- **comparator-audit.json:** 378행, 위상 집합 zero/independent/physical_phase_only, 최소값·P5 잔차·54개 증인(위 3절); 조건수 98.187 불변. SM 그림 3(a)를 렌더링해 P1≈5e-3, P2≈6e-5, P3≈1.2e-5, P4≈4.6e-6, E4≈3e-2 눈금 위치를 확인.
- **runtime12-analysis.json:** 과도 σ 배율 9개 0.9999985(500 d, a=0)–1.0074182(52 d, a=1); 반 간격 배율 0.79207/0.78579/0.84529/0.86472/0.87096/0.87217; 주 사인 0.999913, rank 90, 변화 노름 0.006857, 최소 특이값 1.12207e-6. SM 그림 4(b)의 중간 공분산 값 1.0021/1.0057/1.0000 일치.
- **nuisance-intervals.csv:** 전체 K=1 6.5387e-10/1.5301e-9/5.2645e-9/9.5265e-9/2.7299e-8/6.6420e-8 → 표 1의 6.539e-10 … 2.730e-8 일치; K=1 비율 2.3404–27.5087; K=10 비율 2.1039–17.4150(여섯 지연). 각 (cut, rank, K) 행이 자기 최악 원점(toff)을 따로 갖고(예: 2 d에서 K=1 174.07 d, K=10 174.28 d), 감사 스크립트가 4821 원점(`arange(0, P_out, P_in/24)`) 전부를 평가하므로 새 열은 진짜 등록 격자 포락입니다. 절단 K=1 값은 저장 곡선의 `u95pm_fisher`와 일치.
- **sep_limit_curve_10_8e.tsv:** 65개 지연, K=10 비율 최소 2.1036(2 d), 최대 17.4150(379.379 d).
- **문서·빌드:** 본문 PDF 13쪽, SM 30쪽, 커버레터 1쪽; PDF에 "??" 없음; 절 번호 1–6 유지, AI·데이터 절 무번호; 표 1 여섯 열이 7쪽에 정상 배치; zip의 `main.tex`가 `paper/main.tex`와 내용 동일, `main.bbl` 17개 항목, 그림 1개, README; 본문 그림 1 파일 미변경, SM 그림 3·4 재생성; `failure-ledger-dynamic-chi.md`에 실패 단계와 빠진 최소 가정 기록(AGENTS.md 규칙 준수); 검증기 PASS.

---

### 6. 확인하지 못한 것

- **Ransom et al. 2014 표 1 원문:** 저장소에 PDF가 없고 외부 조회를 하지 않아 기억과 내부 정합성(e sin/e cos → e, ω 재구성)으로만 확인했습니다.
- **Nutimo 시간 원점(treference → 0, t_asc 기준) 주장:** 코드에서 확인하지 않았습니다. 폐합과 세 위상차는 원점과 무관하므로 B1 판정에는 영향이 없지만, 단면·비교기에 쓰인 절대 위상은 이 가정에 의존합니다(B1 이전부터 있던 의존성).
- **재실행:** `simultaneous_inference.py validate` 등은 출력을 덮어쓰므로 실행하지 않았고, 스레드 수 재현성 주장도 시험하지 않았습니다. JSON 간 대조만 했습니다.
- **`compare292.json`(저자 대조 기록)과 다른 심사자 원문:** 열지 않았습니다.
- **새 참고문헌 5편의 DOI 유효성:** 조회하지 않았습니다.
- **소스 zip 단독 컴파일:** 구성과 `main.tex` 동일성만 확인했습니다.
- **공개 스냅숏·소속·ORCID:** 범위 밖으로 두었습니다.