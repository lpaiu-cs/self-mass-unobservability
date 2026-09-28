# 독립 심사 보고서 — fable5.1 (커밋 a322c6462; 원문 보존, 2026-09-27)

## 1. 종합 권고

**주요 수정.** §4.6의 모든 수치는 게시 manifest·JSON과 대조해 정확히 재현된다. 자유낙하 닫힌 식·지연 단극 창 공식·영구 전하 0 논증·J0337 상한의 산술은 모두 성립한다. 그러나 다음 네 가지 문제가 있다.
- (i) 시간 이력 서술이 원천 시각과 판독 시각을 섞어 사실과 다르게 적혀 있다.
- (ii) 기호 β·α·U·a_o가 §3–5의 핵심 기호와 정면 충돌하고, a, B, x, r₀ 등이 정의되지 않았다.
- (iii) 표지 없는 문단과 정의되지 않은 A4가 원고 자체의 규칙(모든 문단에 표지)을 어긴다.
- (iv) 과도 판독값에 "charge"라는 이름과 과도한 정밀도를 부여한 뒤, 나중에 "not a static charge"라고 되돌리는 구성이 독자를 오도한다.

새 계산은 필요 없고 모두 본문 수정이다. 그러나 양이 많고 일부는 진술의 정확성에 관한 것이므로 주요 수정으로 분류한다.

## 2. 지적 사항

**차단: 없음.** 과학적 결론(과도 산란 기억, 영구 전하 0, 궤도 시간척도 붕괴, WD 섹터 no-go)을 무효화하는 오류는 찾지 못했다.

**[주요 1] 시간 이력 첫 항목이 원천 시각과 판독 시각을 혼동한다.**
- 위치: "While the pulse moves inward (to 0.35 s) ... reaching −7.5e−41 when the pulse arrives at the centre."
- 문제:
  - q(t_R)는 R_d에서의 판독 시각 t_R의 함수이고, u = t_R − R_d/c다(`phase278-transit.py` 도입부).
  - 중심 도달은 원천 시각 0.2305 s(`center_arrival_s`)다. 그 신호는 t_R ≈ 2R/c ≈ 0.46 s(`exit_s`=0.4627, `t_max_abs`=0.4626)에야 판독에 닿는다.
  - t_R = 0.23 s에 판독이 보는 것은 깊이 ≲R/2까지 변위된 층이지, 중심 도달이 아니다.
  - 펄스는 0.23 s까지만 안으로 움직이므로 "(to 0.35 s)"는 어느 시각 틀에서도 맞지 않는다.
  - "twelve orders"는 0.35 s 값(−2.7e−39, 12.0자릿수)에, "−7.5e−41"은 0.23 s 값(10.5자릿수)에 대응해 서로 어긋난다.
  - 저장 표본(2e−3 s 간격)의 첫 부호 변화는 0.426 s다.
- 근거: `transit.json`의 `transit/tR,q`, `center_arrival_s`, `exit_s`.
- 수정 문안: "In readout time at R_d the ingoing-pass signal arrives over 0–0.46 s; q stays negative until ≈0.35–0.43 s, growing from −2.3e−51 (T) to −8e−41 (0.23 s) and −2.7e−39 (0.35 s) as denser layers are displaced. The centre-arrival and outgoing-pass signals reach the readout together at 2R/c ≈ 0.46 s, where |q| peaks."

**[주요 2] 기호 충돌과 미정의 기호.**
- (a) β: §3–5에서는 극 세기·정규화 계수인데 §4.6에서는 이론 결합 상수 β=−4다. "|β|=4"와 "U=1.68e−9 of Section 5.3"(β의 envelope)가 한 문단에 나란히 나온다.
- (b) α: 구동 결합 α와 α(φ)=βφ, α₀, α_o, α_p가 겹친다.
- (c) U: §5의 반폭 U, U=rδφ, §4.3의 δU가 겹친다.
- (d) a_o: §4.3의 동반성 전하/질량비 a_o, a_i, a_w가 §4.6에서는 α_o로 바뀐다.
- (e) 그 밖의 충돌: κ=V''(Q₀) 대 κ_lag·κ_struct, φ_k(위상) 대 φ(장), D(가중 차수·모멘트) 대 D(펄스 길이), M(평균 근점이각) 대 M(질량).
- (f) 미정의 기호: a, B, x, x′가 정의되지 않았다. r₀와 R_d는 같은 반지름인데 기호가 다르다.
- 수정: 국소 표기 문장을 둔다. β₀, U_φ(또는 Ψ)로 이름을 바꾸고, ds², x=∫(B/a)dr, r₀≡R_d를 명시한다.

**[주요 3] 표지 없는 문단과 정의되지 않은 A4.**
- "Boundary statement."는 네 표지에 없다. A4도 원고에 정의되지 않았다.
- 수정: Conjectural로 표지하고, "assumption A4 (no orbital-timescale internal state variable in the free-fall sector, in the repository ledger)"로 정의한다.

**[주요 4] 과도 판독값의 "charge" 용어·정밀도·구성 불균형.**
- 문제:
  - q(T)는 프로토콜 시각 T=2D의 지연 단극 판독값인데, 앞 문단들이 이 스냅숏에 수렴 차수, EOS 민감도, 부호 여유, 광구 ×3.2/×6.9를 부여한다.
  - 초록의 "zero permanent charge"는 거의 자명한 명제인데 "shows"로 적었다.
- 수정:
  - 해석을 먼저 쓰고 q(T)를 "endpoint retarded-monopole readout"으로 부른다.
  - 크기 문단은 "for the protocol snapshot"으로 한정한다.
  - 초록은 "indicates"로 낮추고, "minimal free-fall internal state" 대신 "the scalar-driven displacement of the inner white dwarf's own matter"로 쓴다.

**[주요 5] 주장 표지 오류.**
- (a) 5.2e−5 재현은 Counterexample candidate다.
- (b) "corrections are small" 문단은 무표지(=Proven)로 되어 있다. 노트는 고정점을 "부등식 Proven, 적용 Conjectural"로, φ_∞² 상태 부분을 Conjectural로 분류한다.
- (c) 초록·논의 추가문이 Conjectural인 궤도 붕괴 주장을 Counterexample candidate 아래 둔다.
- (d) 영평균 문단은 "in the linear adiabatic radial model"을 명시해야 한다.

**[주요 6] J0337 상한 문단의 제시 방식.**
- (i) U는 단위 구동 계수의 envelope다. §5.4가 "cannot be relabeled as an undifferentiated amplitude of SEP oscillation"이라 경고한다.
- (ii) 척도 선택: K=10 절단값 1.68e−9을 골랐다. §5.5는 전체 랭크 3.53e−9을 기준선으로 삼고, K=1 값 2.79e−10(REML 5.25e−10)이 가장 작다. 요구값 범위는 257–1.4e5이므로 결론은 변하지 않지만, 범위로 제시해야 한다.
- (iii) 선언 φ_∞=1e−3은 Cassini 2σ 한계를 넘는다. 7.5e−21은 φ_∞=8.84e−4에서 재평가한 값인데 원고는 이를 말하지 않는다. 선언값 그대로면 1.1e−20이다.
- (iv) Cassini 값이 2σ라는 점이 빠졌다.
- (v) 논의 추가문이 한정어를 빼고 중성자별로 단정하며, 바깥 백색왜성 지연을 다루지 않는다.
- 수정: 척도를 범위로 제시하고, Cassini 2σ와 재평가를 명시하고, 한정어를 복원하고, 바깥 백색왜성 배제 문장을 더한다.

**[주요 7] 조석 채널의 "collapses to static coefficients" 표현.**
- 진행파 영역의 동적 조석은 연속 흡수, 곧 소산적·위상 지연 응답이다.
- 옳은 진술: l≥2는 선형 차수에서 단극에 들어가지 않고(선택 규칙), 자체 크기는 평형 조석으로 묶인다(상대 조석력 9.1e−12, 근점 이동 1PN의 1.5e−5).

**경미**
- 1. 정적 창 값의 부호가 반대다(+3.6e−54). "of opposite sign and 1/640 the magnitude"로 고친다.
- 2. 흐름 오차는 5.557e−5이므로 5.6e−5로 적는다. 0.25%는 핵 가중값이고 최대 편차는 0.38%다.
- 3. 비회색 보정의 한계를 적는다. 단순 H·He 연속 불투명도의 비를 옮긴 것이고, 그 Rosseland 평균은 표 값의 0.55–0.74배이며, 평행평판이다.
- 4. 빠진 가정: 선형 응답, 전체 별 모형의 뉴턴 자기중력과 평탄 외부 판독, WKB·확산(τ_w는 하한), 단일 완화 직교 상한, 대기 변경 시 펄스 핵 고정, 광구에 단열 연산자를 쓴 점.
- 5. "hot helium-core model"로 적었지만 봉투는 X≈0.89로 수소 우세다.
- 6. §4.3("hold the companion charges fixed", "A responsive companion ...")과 교차 참조하면 κ_struct가 그 가정을 정량적으로 뒷받침한다.
- 7. 조판: "2x-to-4x"를 ×로, "l=2"를 수식으로 바꾼다. 식 라벨은 참조되지 않는다. PDF와 main.bbl은 재빌드가 필요하다.
- 8. 데이터 가용성: 노트가 한국어라는 점, WSL 배열은 해시로만 기록된다는 점을 적는다.
- 9. 참고문헌 세 항목은 기억과 일치하고 키 중복이 없다.

## 3. 수치 대조

수치 대부분이 일치한다(q, 차수, 성분 몫, EOS, 자유낙하, 부호 여유, 창, 이력 최대값, 모드, 보정 항, 광구 배수, 힘 비, k₂, τ_w, κ_struct, δφ_mod, Cassini, 요구값, 상한).

불일치 또는 부정확:
- 시각 틀(0.35 s / 중심 도달)
- 정적 창 부호
- 흐름 오차 5.56e−5
- 0.25%(핵 가중)

## 4. 확인하지 못한 것

- 결합 계산·전체 별 통과의 재실행
- 비회색 풀이기 코드
- 문헌 원문(웹 미조회)
- Pandoc 재현
- 2e−4 s 전 해상도 부호 이력(저장 표본에서는 첫 부호 변화가 0.426 s)
