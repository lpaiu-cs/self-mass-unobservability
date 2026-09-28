## 심사 보고서: "Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715"

읽기 전용으로 작업했고 저장소는 수정하지 않았습니다. 확인 방법: 본문·SM·이전 원고(39962a753)의 전문 대조, `verify_unified_paper.py` 실행(PASS, exit 0), 근거 JSON 15개를 읽는 짧은 Python 대조, 본문의 모든 항등식·정리를 sympy로 독립 재검증, PDF 텍스트 추출(pdftotext), 참고문헌·LaTeX 소스·제출 부속 검토.

---

### 1. 종합 권고

**주요 수정(major revision).** 본문에 남은 수학은 전부 옳고(정리 1의 양방향, 식 6–10·12–16을 독립 검증), 본문의 수치는 근거 JSON과 예외 없이 일치하며, SM 본문은 이전 원고와 제목·초록만 다릅니다. 그러나 (i) 서론의 새로움 진술이 시스템 식별·모멘트 문제·유리함수 보간 문헌과의 관계를 전혀 밝히지 않아 정리 1과 식 7–8이 고전 결과의 특수 경우임이 드러나지 않고, (ii) 초록의 J0337 헤드라인 두 개(위상 폐합 위반, 2.1–17.4배)가 데이터의 성질이 아니라 저자 자신의 이전 저장 분석에 대한 정정임을 본문이 명시하지 않으며, (iii) 자체 보정한 95% 문턱을 넘는 여섯 계수 통계량(16.35 > 12.82)을 한 문장으로 지나치고, (iv) 축약 과정에서 의미를 바꾸는 한정어 몇 개가 빠졌습니다. 차단 사유는 없습니다.

---

### 2. 지적 사항

#### 차단
없음.

#### 주요

**M1. 새로움 진술과 선행 문헌(§1 넷째 문단, §3 전체).**
- 문제: "The proofs use only elementary tools ... What is new is the set of identifiability boundaries"라고 하지만, 경계 자체가 고전 결과의 특수 경우입니다. 정리 1은 K개 정현파 입력이 차수 2K의 지속 여기(persistent excitation)라는 시스템 식별의 표준 셈(Ljung, *System Identification*, 13장; Åström–Bohlin 1965)이고, 유리함수의 유한점 보간(Loewner 틀, Antoulas)과 같은 내용입니다. 식 7–8은 Stieltjes 함수의 2×2 Hankel/Pick 양성 조건(Stieltjes 모멘트 문제; Nevanlinna–Pick 보간)이며, 유전 완화·유변학에서 Debye 단일 완화를 Cole–Cole 도에서 판정하는 방식과 같습니다. 식 9의 가까운 극 한계는 역 Laplace/Stieltjes 변환의 부적절성(예: Honerkamp–Weese 1989)이고, §3.4는 유한 주파수 응답이 연산자를 정하지 못한다는 자명한 사실입니다.
- 근거: 본문 참고문헌 10편 중 이 분야 인용이 0편. 이 문헌을 아는 심사자는 §3을 재유도로 볼 것입니다.
- 수정: 위 문헌 3–5편을 인용하고, 새로움을 "알려진 식별 가능성 셈을 nuisance 투영된 자유낙하 타이밍 문제에 적용한 것, EFT 물리(짝수/보존, 양의 상호·속도 간극)로 비교기 부류를 고정한 것, 위상 폐합 게이트, J0337 정보 감사(§5.4)"로 정확히 한정할 것.

**M2. 여섯 계수 옴니버스 통계량(§5.5).**
- 문제: "The recorded data's six-coefficient omnibus statistic, 16.3525, lies above that threshold"—자체 보정 95% 문턱 12.8242를 넘고(명목 χ²₆ p≈0.012), `simultaneous-validation.json`도 `exceeds_threshold: true`입니다. 본문은 영가설이 "all six carrier coefficients vanish"이지 β=0이 아니라는 한 문장으로 넘어갑니다. 물리 지연 단면의 최소 통계량(11.04–11.78)이 문턱 아래여서 "β=0 포함"은 맞지만, 세 반송파에서 nuisance 제거 후 유의한 잔여 신호가 있다는 사실 자체는 설명이 필요합니다(궤도 모형 불일치? 생략 입력 RMS 3.37%? 잡음 모형?).
- 수정: 여섯 계수의 추정값·표준오차(또는 z)를 표나 SM에 제시하고 p값과 후보 원인을 한 문단으로 논의할 것. 초록의 "No relaxation is detected"는 유지 가능하나 이 논의 없이는 설득력이 약합니다.

**M3. 헤드라인이 저자 자신의 이전 분석에 대한 정정임을 명시할 것(초록, §4 마지막, §5.2).**
- 문제: "an archived drive dictionary violates the phase closure"와 "changes by factors of 2.1--17.4 with the retained nuisance span"은 각각 저자의 보조 저장 분석과 저자가 택했던 특이값 절단(10⁻³, rank 71)에 대한 정정입니다. 본문 어디에도 "our earlier"가 없어 독자가 공개 Voisin 자료의 결함으로 오독할 수 있습니다. 71차원 공간은 표준 관행이 아니었으므로 2.1–17.4는 데이터의 성질이 아니라 방법 선택의 효과입니다.
- 수치 근거: 표 1의 비율 최대는 17.35인데 본문은 "2.1--17.4"라고 씁니다. 65개 지연 곡선(`sep_limit_curve_10_8e.tsv`)에서 최대 17.41(τ=379.4 d)이 나오므로 값은 맞지만 표만 보는 독자에게는 불일치입니다.
- 수정: §4·§5.2에 "우리의 이전 보조 분석/절단 선택"임을 명시하고, "over the 65 stored lags (Table 1 shows five)"를 덧붙일 것. 초록에서 2.1–17.4를 방법론적 진술로 낮추는 것을 권합니다.

**M4. 축약에서 잃은 한정어(재구성 충실도).** 개별로는 경미하나 합치면 어조가 바뀝니다.
- (a) §6 백색왜성 문단: SM §6은 "Under the stated assumptions, its structural modulation..."인데 본문은 조건구를 뺐습니다. 가정(|a_p|≤1, |a_o|≃α₀ Cassini 한계, 선형 응답, φ_∞² 재척도)을 한 절로 복원할 것.
- (b) §5.6: "increases beta standard errors only by factors 1.00096--1.00523"의 "only"는 새로 들어간 강조이고, SM의 "does not certify every lag, large initial amplitude, or coverage after adding that coefficient"가 빠졌습니다. §3.4가 과도 응답의 불충분성을 정리로 세운 뒤이므로 이 완화는 특히 눈에 띕니다.
- (c) §5.5 끝점 4.10217e-10, 8.31672e-9: SM의 "smaller endpoints than a pointwise interval are not automatically stronger physical constraints"가 빠져 물리 구동 한계처럼 읽힙니다.
- (d) §5.3 REML: 추정기가 생성 공분산과 "the same fixed extra-Fourier family"를 쓴다는 조건과, 절단 K=10에도 0/8192 사례가 있다는 사실이 빠져 "at least 511/512"가 무조건적 성공처럼 읽힙니다.
- (e) §4: "the generic signed-beta fit admits a larger model", "Rescaling the unit-drive interval cannot do this", 그리고 SM §4.5의 "A specified microscopic theory may relate these parameters"가 빠졌습니다. 마지막은 "Equilibrium information alone therefore does not establish the rate gap" 뒤에 꼭 필요합니다.
- (f) §5.1: 한 외측 주기가 반송파의 공통 주기가 아니라는 SM §5.2의 주의와, β가 진폭이 아니라는 SM §5.4(200 d에서 단위 구동 진폭 배수 0.0013/0.252/0.0013)가 본문에서 참조되지 않습니다. "(SM Section 5.4)" 한 번이면 됩니다.

#### 경미

- **m1. §3.2 비교기 부류의 물리적 동기.** 식 4는 과감쇠 1차 모드만 허용하고(Γ, K 양의 정부호), 관성·진동 모드(f/p/g 모드, 동역학 조석)를 배제합니다. §4는 Khalil 감쇠 행렬이 특이함을 보이므로 이 부류에 속한다고 알려진 J0337 천체는 없습니다. SM §4.6의 층별 열 완화(양의 Debye 스펙트럼)가 정확히 식 5의 부류이므로 그것을 동기로 연결하면 §3.2–3.3이 §4.6과 이어집니다.
- **m2. 표지 "Counterexample candidate".** PRD 독자에게 "무엇에 대한 반례인가"가 본문에 없습니다(§1 정의는 "a proposed physical/model realization"). "Model" 또는 "정적 전용 기술에 대한 반례 후보"로 풀어 쓸 것.
- **m3. §3.1 "deprojected readout"**: 정의되지 않은 용어. "nuisance-projected readout"으로.
- **m4. §3.3 초록 문장** "two exact calibrated complex samples identify one observable relaxation time"은 등식이 성립할 때에 한정됨을 "when the equality holds"로 명시할 것.
- **m5. §4 Q_A=a_A m_A 사용 시 Damour–Esposito-Farèse(bib에 있으나 본문 미인용)를 인용할 것.** 본문 참고문헌 10편은 PRD 기준으로 얇습니다.
- **m6. §5.1 "the two high-frequency carriers are close"**: 2988일 동안 위상차 57 rad로 분해되므로 "close relative to ω_in; condition numbers 98–125"처럼 정량화할 것.
- **m7. §5.3** 아홉 진폭(0, ±2, ±5, ±20, ±50 σ)을 괄호로 제시; **§5.6** "201 assignments"에 "at most two nonzero ±1 steps"를 추가.
- **m8. 표 1의 절단 K=1 열**: §5.3에서 생략 평균 하에 포함률 0인 구성을 표에 그대로 두는 이유를 밝히거나, 전체 공간 K=1 값이 있으면(`nuisance-audit.json`의 intervals) 그것으로 바꿀 것.
- **m9. AI 사용 절**: 모델명·버전은 있으나 도구 버전(Claude Code, Codex CLI)과 사용 기간이 없고, "the verification scripts listed below"인데 아래에는 스크립트 하나만 있습니다. APS 정책의 이름·버전·용도·저자 검증 요구는 그 외에는 충족됩니다.
- **m10. 제목 블록의 "Repository:" 줄**은 데이터 가용성 절로 옮기고, 소속·ORCID를 채울 것(알려진 미결).
- **m11. LaTeX**: "Section 3.1" 등 절 번호가 하드코딩되어 있어 지금은 맞지만 REVTeX 전환 시 깨지기 쉬움; `\ref` 권장. `inputenc` 경고는 XeTeX에서 무해. 미정의 참조·"??" 없음, 인용 10건 모두 bib에 존재, 표·그림·식 번호 본문 언급과 일치.

---

### 3. 직접 확인한 수치

**일치(근거 파일과 대조)**
- 표 1 전 항목 20개와 비율(2.10/5.22/16.02/17.26/17.35) — `sep_phase_marg_10_8e.json`; max|z|=2.2851, p=0.26, detection false, K=10, 4821 origins, 65 lags 2–500 d, 외측 주기 327.26 d.
- rank 90/71, 19차원, 최소 상대 특이값 4.18e-7 — `nuisance-audit.json`.
- 표 2 전 항목(최악 셀 4772/8192=0.5825 포함), 2916행=6×3×18×9 — `coverage-audit.json`; REML 486조건 7749/8192=0.945923 — `estimated-covariance-audit.json`; 완전 격자 K=10 최소 511/512(τ=200 d, extra Fourier 1, full diag).
- 조건수 98.19/111.69/124.82; a=0,0.0831611,1; 정보 최소 0.08289/8.36e-6/3.19e-6/2.88e-6; 1/√(2.88e-6)=589.4; 짝수 0.003035; P5 9.3e-30–1.4e-24; Λ=38.5614=10ω_in; 54개 증인 모두 양 — `comparator-audit.json`.
- ϖ_p=-0.12326821, ϖ_b=-0.10004792, C=3.11837236=π+ϖ_p−ϖ_b, 위상차 1.69406454/1.67084424/π, 보관 폐합 3.6e-13, 14개 검사 — `physical-matching.json`.
- 정규화 진폭 (0.24394925, 0.69195682, 0.06409393): JSON과 일치하고, 질량·반장축·이심률에서 식 13으로 독립 재계산해도 동일; 생략 입력 RMS 3.37427%(128², 256²) — `corrected-physical-drive.json`.
- 문턱 12.8241766, 포함률 0.951660/0.954346/0.952759/0.954102(=95.17–95.43%), 옴니버스 16.3525, 끝점 4.10217e-10/8.31672e-9, 여섯 단면 모두 비어 있지 않음 — `simultaneous-*.json`.
- 과도 0.7984–0.9973, 565 간격, 최약 223.3708 d Δχ² 0.152887(a=1)/27.9406(a=0) — `phase-state-audit.json`; 201 배정, Δχ² 27.9406/11.6183/0.152887 — `gap-pair-audit.json`.
- 살아 있는 과도 τ=2,52,500 d 5% 게이트 통과, SE 비 1.00096–1.00523(9개), 주 사인 0.999913, SE 인자 0.8680–0.8785, 이심률 6.58/13.09/26.11>1 — `runtime12-analysis.json`.
- 3453.88 d=500 ln 1000; 12,474 TOA(ν=12474−90−2=12382 profile_dof로 간접), 2987.9 d(REQUEST10_8H), ω_in=3.856137(증인 간극율/10), ω_out=2π/327.255=0.019200.
- 백색왜성: 4.0e-9(3.978e-9), 3.3e-7(3.348e-7) — `native-thermal-relaxation-manifest.json`; |a_p|≳0.5(0.476), |a_o|≲5e-6(5.5e-6) — REQUEST286; "eight orders"(2.13e-18 대 2.8e-10).
- SM 부록 앵커 U(K=10)=1.67953e-9(β̂=5.51335e-12, σ_F=8.569e-11로 이분법 재계산).
- 재구성: SM은 39962a753 원고와 1행(제목)·9행(초록)만 다름(diff); 본문 11쪽·SM 30쪽; 커버레터의 "Section 7" 일치; 평문 초록 237단어; 소스 zip 구성(main.tex, main.bbl, references.bib, 그림 1개, README).

**수학(sympy 독립 검증, 모두 성립)**: 식 9 항등식; 식 7의 R, D, Q가 3원자 측도 dμ의 0·1·2차 모멘트; 단일 극에서 RQ−D²=0, τ₀=D/R, a₀ 공식; 식 6의 비의 τ 단조성(도함수 인수 2τ(ω_h²−ω_l²)/(1+ω_l²τ²)²); 정리 1의 실계수 구성(K=2 일반값, 차수 3, 여섯 점 일치)과 복소계수 차수 K−1 구성; 식 10의 상수·반송파 소거와 e^{−t/τ} 응답 −∏(1+τ²ω_k²)/τ⁷; 관성 오차 항등식 H/H₀−1=x/(1−x); 단열 나머지; 식 15 템플릿 부호; 식 13의 1/|R+fr| 전개 부호와 A_dif; 폐합 원점 불변성.

**표기 불일치(값은 옳음)**: "2.1--17.4" 대 표 1 최대 17.35(M3).

---

### 4. 초기 반려 위험과 게재 전망

- **편집자 단계 반려 위험: 중간(약 30–40%).** 이유: 단일 저자·소속 미기재(수정 예정), 모든 문단의 표지 형식이 비관행적, 본문 참고문헌 10편, 핵심 결과가 조건부·부정적이며 공개 자료의 저장 재분석에서 비검출, 백색왜성 계산이 30쪽 SM에만 있고 본문은 한 문단 요약, 헤드라인 하나가 저자 자신의 분석 정정. 유리한 점: 범위가 명확하고 주장이 정직하며 PRD는 중력 검증 방법론 논문을 게재하고, AI 공개는 정책 요건을 대체로 충족하며, 커버레터·초록·본문이 일치합니다.
- **심사로 갔을 때**: 시스템 식별이나 펄서 타이밍을 아는 심사자는 M1·M2·M3을 지적할 가능성이 높고 "주요 수정"이 예상됩니다. 수정 후 수락 확률은 25–40%로 봅니다. 수학은 옳지만 규모가 작고, J0337 적용은 내부적으로 일관되나 새 제약을 내지 않습니다. 결론 강도는 근거를 넘지 않습니다(비검출·조건부·척도 비교로 명시). 가장 큰 과학적 공백은 §4 스스로 요구한 "A corrected physical analysis must prescribe ... then validate"가 §5.5의 단면 계산으로만 부분 충족된다는 점입니다.

---

### 5. 확인하지 못한 것

- 그림 1의 실제 곡선: PDF 래스터라이저(pdftoppm, pypdf, fitz)가 없어 캡션과 축 레이블만 텍스트로 확인했습니다. 그림 파일은 4897038c4 이후 변경이 없습니다.
- 저장 응답 열(npz), Nutimo 런타임, 포함률 몬테카를로: 2분 한도 때문에 JSON 요약값만 대조했고 재계산하지 않았습니다.
- SM §4.6 백색왜성 계산 내부: 본문 논의가 인용한 manifest 값만 확인했습니다.
- `check291-numbers.py`의 "139개 수치 모두 이전 원고에 있음": 재실행하지 않고 본문 수치를 SM 텍스트와 수작업으로 대조했습니다(빠진 것 없음).
- 소스 zip의 단독 컴파일: 구성만 확인했습니다.
- APS 2026년 6월 AI 정책 원문: 가져오지 않았고 체크리스트의 요약과 기존 정책 지식으로 판단했습니다.
- 공개 스냅숏 상태: 지시대로 확인하지 않았습니다.