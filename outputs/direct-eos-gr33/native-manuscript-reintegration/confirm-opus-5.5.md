# 반영 확인 재심: opus5.5, 개정 5판

저장소 파일은 읽기만 했습니다. 스크래치 폴더에 짧은 Python 확인 스크립트 두 개만 만들었습니다. 이 호스트에는 Pandoc이 PATH에 없어서 목록·식 구조는 Python으로 셌습니다. 다른 심사자의 원문은 열지 않았습니다.

## 1. 종합 권고: 경미 수정 후 수락

- **이전 지적:** 4판 재심의 주요 A·B는 해소됐고, C는 대체로 해소됐습니다. 경미 15건은 모두 해소됐습니다.
- **수치:** 새 수치와 정오표 산술은 모두 저장 JSON에서 재현됩니다.
- **결론의 강도:** 차단·주요 사항은 없습니다. §4.6, 요약, 초록, §6은 모두 "scale comparison, not an exclusion"의 강도를 지킵니다. 이는 계산한 범위와 맞습니다.
- **새 경미 지적:** 원고 문장에 닿는 것은 N1–N5입니다. 모두 한 구절 수준이고 새 계산이 필요 없습니다.
- **통합 판단:** N1–N5를 반영하면 추가 재심 없이 단계281 절차(bib 3건 복원 포함)로 통합해도 됩니다. N6은 선택입니다. N7–N8은 노트와 manifest 수준이라 통합을 막지 않습니다. 다만 데이터 가용성 문구(N7)는 원고 문장이므로 함께 고치는 편이 낫습니다.

## 2. 이전 지적별 판정 (4판 재심)

| 지적 | 판정 | 이유 |
|---|---|---|
| A 순간 응답을 §4.3에 귀속 | 해소 | "a zero-lag coefficient, not an internal state"로 적었습니다. 고정 동반성 축약이 두 쌍 인자의 동상 변조를 뺀다고 밝혔습니다. 강성 조건은 별도 문장이고 평가하지 않았다고 적었습니다. 해석 부분은 Conjectural입니다. 정오표 8과 docs의 단계285 항목이 이전 표현을 정정합니다. |
| B 요약·초록의 가정 누락 | 해소 | 초록, 요약, §6 모두 두 가정(완화 세기 ≤𝒮_struct, 단일 완화 시간)을 적습니다. "only a small"은 "computed adiabatic structural term"으로 바뀌었습니다. |
| C 정오표 범위 | 대체로 해소 | 정오표 7–15는 제가 든 노트 행 대부분, manifest 판정 문자열, 분류 문장의 오류를 다룹니다. 데이터 가용성에 대체 문장과 REQUEST284 대응이 들어갔습니다. 남은 누락은 N7에 적었습니다. |
| 1 a_p 정의 | 해소 | a_p=Q_p/m_p |
| 2 변조 식의 차수 | 해소 | "eccentricities and in a_in/a_out" |
| 3 쌍 채널 문장 | 해소 | 구조 척도를 먼저 적었습니다. 남은 문법 문제는 N3입니다. |
| 4 가정 목록 | 해소 | \|a_p\|≤1과 약한 장 전하가 들어갔습니다. |
| 5 𝒬와 순간 응답 | 해소 | 𝒬_c와 𝒬를 나눴고, 긴 파장 극한 문장을 넣었습니다. |
| 6 영구 전하의 시간척도 | 해소 | 본문에 t→∞와 열 모드 잔여를 적었습니다. 초록은 N5입니다. |
| 7 "strongly overdamped" | 해소 | "resonant, two-pole (inertial) or mixed-sign" |
| 8 neither…nor | 해소 | "posits"와 "coupling at the Section 5 scales" |
| 9 보간 | 해소 | "interpolated from" |
| 10 2.5배 추정 | 해소 | 두 곳 모두 추정이라고 밝혔습니다. |
| 11 기호·제목 | 해소 | 𝒫₁·𝒫₂, s, "Structural and thermal"로 바꿨습니다. s는 §5.2의 잔차 척도와 겹치지만 지역 기호라고 선언했으므로 충분합니다. |
| 12 2차 조석 | 해소 | "scaling as φ∞" |
| 13 REQUEST278 불일치 원인 | 해소 | 정오표 13에 넣었습니다. "2 ms 앞선"은 N8에 적었습니다. |
| 14 자기 결합 이중 계산 | 해소 | "Bounding instead the full change…"로 고쳤습니다. REQUEST285의 3.31e−8도 재현됩니다. |
| 15 REQUEST272 19행 | 해소 | 1653배가 재현됩니다. |

## 3. 새 지적

### 차단: 없음

### 주요: 없음

### 경미

**N1. |a_p|≳0.5 문턱이 정적 SEP 한계와 양립하는지 밝히지 않음 (통합 전 반영 권장)**

- **위치:**
  - §4.6 Instantaneous 항목: "reaches the stored Section 5 scales for \|a_p\|≳0.5"
  - §6 추가문: "can reach those scales for \|a_p\|≳0.5 but is not a lag"
- **문제:** 산술은 맞습니다. 그러나 §4.3의 쌍 인자에서 J0337의 정적 SEP 신호는 Δ≈a_o(a_p−a_i)입니다. 그러면 |a_p|≳0.5는 |a_o|가 아주 작을 때만 가능합니다. 이 맥락이 없으면 §6 문장이 열린 채널처럼 읽힙니다.
- **근거:**
  - `request5/REQUEST5_J0337_PHASEA.md`에 기록된 값은 Voisin 2025의 |Δ|&lt;1.5–2.3e−6과 Archibald 2018의 &lt;2.6e−6입니다.
  - 절의 가정대로 |a_o|≃|α₀|=3.54e−3(Cassini)이면 |a_p−a_i|≲6.5e−4이고 |a_p|≲4.2e−3입니다. 이때 순간 항은 ≲2.2e−14입니다.
  - |a_p|≥0.48이 되려면 |α₀|≲4.8e−6, 곧 φ∞≲1.2e−6이어야 합니다.
  - 이 문턱 문구는 제 4판 제안에서 나온 것이므로, 제 제안을 보완하는 지적입니다.
- **수정:** §4.6과 §6에 다음 취지를 넣습니다. "Such \|a_p\| is compatible with the static SEP bound only for \|α₀\|≲5×10⁻⁶; not evaluated." 최소한 "compatibility with static-SEP bounds was not evaluated"만이라도 넣습니다.

**N2. 외측 백색왜성의 같은 지연 없는 응답이 빠짐**

- **위치:**
  - "omits this in-phase modulation of both pair factors"
  - §6의 단수 표현 "The white dwarf's instantaneous zero-lag response"
  - "The outer white dwarf's own internal states were not computed."
- **문제:** 고정 동반성 축약은 외측 백색왜성의 δa_o=β_sδφ도 뺍니다. 그런데 이 절의 용어로 이 항은 내부 상태가 아닙니다. 그래서 마지막 문장에도 포함되지 않고, 어디에서도 다뤄지지 않습니다.
- **근거 (제 추정, 저장값 아님):**
  - 외측 위치에서 펄서장의 변조는 m_p·e_out/a_out=4.25e−10|a_p|이고, 쌍극 항 m_p·f·a_in/a_out²=3.9e−11|a_p|이 더해집니다.
  - 따라서 펄서–외측 쌍 인자에는 ≲1.86e−9 a_p²가 들어갑니다. 내측의 1.24e−9 a_p²와 같은 차수이고, 문턱은 |a_p|≳0.39입니다.
  - N1의 한정이 여기에도 똑같이 적용됩니다.
- **수정:**
  - §4.6에 한 문장을 넣습니다. "The outer white dwarf has the analogous zero-lag response in the pulsar–outer pair factor, of similar size; not evaluated."
  - §6은 복수형으로 고칩니다.

**N3. 구조 변조의 합계와 "its lagged part"**

- **위치:**
  - Cassini 문단: "Both the in-phase structural modulation and … its lagged part"
  - "at least eight orders" (본문과 §6)
- **문제:**
  - (a) "its"가 동상 단열 변조를 가리킵니다. 그러면 정의상 지연이 없는 항의 지연 부분이 됩니다.
  - (b) 완화 부분의 동상 성분이 빠졌습니다. 이 성분은 |ΔS|δφ_mod ≤2.13e−18입니다. 단열 부분과 완화 부분을 합치면 ≤4.26e−18이고, 가장 작은 척도와의 비는 6.6e7(7.8자릿수)입니다. "at least eight"는 각 부분에는 맞지만 합계에는 조금 넘칩니다.
- **수정:**
  - "its lagged part"를 "the lagged thermal part"로 바꿉니다.
  - "each at most 2.13e−18"로 적고, 합계 4.3e−18을 덧붙입니다.
  - "about eight orders"로 씁니다. 본문 136행은 이미 이 표현을 씁니다.

**N4. "reaches the stored Section 5 scales"의 범위**

- **문제:** |a_p|≤1 안에서 1.24e−9 a_p²가 닿는 것은 가장 작은 두 척도뿐입니다.
  - 2.8e−10: |a_p|≥0.48
  - 4.1e−10: |a_p|≥0.58
  - K=10 봉투는 |a_p|≥1.17이 필요합니다.
- **수정:** "reaches the smallest stored Section 5 scales"로 바꿉니다.

**N5. 초록의 영구 전하 표현**

- **위치:** 초록의 "no first-order permanent charge in the limit of damped modes"
- **문제:** 여기서 "limit"이 t→∞를 뜻하는지 분명하지 않습니다. 본문 86행은 t→∞에서만 성립하고, 열 모드 잔여가 관측 기간 내내 남을 수 있다고 적습니다. 초록에는 이 한정이 없습니다.
- **수정:** "as t→∞ if all modes, including slow thermal ones, are damped"로 고칩니다.

**N6. 판독 용어 (선택)**

- 32행은 "the readout"을 𝒬로 정의합니다.
- 그런데 전체 별 모형 문단의 다음 표현은 모두 𝒬_c를 가리킵니다(질량 항 미적용).
  - "The readout rises steeply"
  - "The readout stays negative"
  - "rms readout"
- 끝점에서 두 값의 차이는 2.2e−5라 수치 영향은 없습니다. 첫 언급에만 "𝒬_c"를 적어도 충분합니다.

**N7. 정오표 7–15의 남은 누락 (저장소 수준, 문구 하나는 원고)**

- **정오표 14에서 빠진 필드:**
  - `native-tidal-photosphere`
    - `dissipation.lagged_monopole_delta_bound` 2.718e−18: 외측 항을 뺀 값입니다.
    - `margin` 6.25e8: 철회한 감도 진술입니다.
    - `final_charge_conclusion`: "collapses to static coefficients … radial, tidal and dissipative channels"
  - `native-final-closure`
    - `failure_ledger`: |κ_lag|≳1.6e3/|α_p|와 "static structural response"
    - `closure.freefall_delta_bound_phase274`
  - `native-structure-eft-boundary`
    - verdict의 `SIGN_STRUCTURALLY_ROBUST`
  - 이 값들은 `paper/revision-manifest.json` 47878·47880·47917·47978·48006–48007행에도 그대로 있습니다.
- **데이터 가용성 문구:** 데이터 가용성은 "Verdict fields"만 대체된다고 적습니다. 그래서 수치 필드는 문언상 빠집니다.
- **정오표 11:** REQUEST276 47행 끝의 "선언 이론의 선형 약한장 영역에서는 이 가정이 성립하지 않는다"가 명시적으로 철회되지 않았습니다.
- **정오표 15:** failure-ledger 2578행(단계277–278)의 두 항목이 목록에 없습니다.
  - "중심 도달 −7.5e−41": 정정값은 판독 시각 기준 −8.0e−41입니다.
  - "(Proven)": 2588행에도 있습니다.
- **수정:**
  - 데이터 가용성 문구를 "Verdict, classification and derived-bound fields"로 넓힙니다.
  - 정오표 16을 두고 위 필드와 행을 적습니다.

**N8. REQUEST285의 세부 (노트 수준)**

- **외부 입력 열:** `balanced-initial-state.npz`(268)와 `geometry.npz`(270–273)는 git에 추적되어 있습니다. 곧 이미 공개돼 있습니다.
  - 머리말의 "large_arrays_sha256에만 해시로 있고 저장소에는 없다"는 이 둘과 bank.npz, coupled-128.npz, source-128.npz에 대해 틀립니다.
  - geometry.npz는 native-quad-refined-primary의 large_arrays_sha256에도 없습니다.
  - 다시 실행할 수 없다는 결론은 그대로 유지됩니다. source-64.npz, field-source-64.npz, 광자 seg 파일, captures가 공개되지 않았기 때문입니다.
- **실행 명령:** "실행 명령은 run*.sh에 있다"와 달리, native-final-closure와 native-conditions-closed에는 run*.sh가 없습니다.
- **정오표 13:** "2 ms 앞선"은 실제로는 1.8 ms 앞선 직전 표본입니다. 20.2 ms부터 2 ms 격자이기 때문입니다.
- **정오표 11:** 8.8e−9는 37·43·47행이 아니라 37·39·47행에 있습니다.
- **정오표 12:** 1.0e−27은 REQUEST277 26·33행에 있습니다.
- **정오표 9:** 인용한 문장은 19행의 "첫 문장"이 아니라 둘째 문장입니다. 제 4판 재심의 표현을 그대로 옮긴 것이라 제 실수입니다.
- **Born 값의 출처:** Born 1.7e−10은 `native-incident-response`에 있습니다. 이는 데이터 가용성의 지도와 REQUEST244–285 범위 밖입니다.

## 4. 직접 확인한 수치

| 항목 | 초안·노트 | 재현 | 판정 |
|---|---|---|---|
| 순간 응답 | 1.24e−9\|a_p\|+8.2e−10\|a_o\| | 4×3.07608e−10=1.2304e−9, 4×2.03106e−10=8.124e−10 | 일치 (올림) |
| 문턱 | \|a_p\|≳0.5 | 척도 2.794e−10에서 0.476, 4.10e−10에서 0.577, 1.68e−9에서 1.17 | 근사 일치 (N4) |
| 펄서–내측 순간 항 | ≈1.24e−9 a_p² | 1.2304e−9 a_p² + 2.9e−12\|a_p\| | 일치 |
| Cassini | 3.54e−3, 8.84e−4 | 3.535556e−3, 8.83889e−4 | 일치 |
| 구조 척도 | 2.13e−18, 7.53e−21 | 2.12858e−18, 7.52571e−21 | 일치 |
| 여유 | ≥8자릿수 | 각 부분 1.31e8. 합계 4.26e−18이면 6.6e7 | N3 |
| 지연 상한과 반례 | ≤\|ΔS\|/2, s/4 | max ωτ/(1+ω²τ²)=1/2, Im[(s/2)/(1+i)]=−s/4 | 일치 |
| 정오표 9 | 1653배 | C₊=1.4149e−54, C₋=−2.3387e−51이면 1652.9. 반례는 +4.90e−55 | 일치 |
| 정오표 10 | 4.5e−18 | 8.83658e−9×5.1071e−10=4.513e−18 | 일치 |
| 정오표 13 | 0.0982/0.2282/0.3482 s | 값 −2.672e−43/−7.526e−41/−2.376e−39, 간격 1.8 ms | 값 일치 (N8) |
| 산술 정정 | 3.31e−8, 9.9e−8 | 3.306e−8, 9.92e−8 | 일치 |
| 이력 | −3.0e−43/−8.0e−41/−2.7e−39, 부호 전환 0.43 s | −2.972e−43/−8.037e−41/−2.654e−39, 0.4262–0.4282 s | 일치 |
| 3.5·4.0 ms | −2.8e−51, −8.7e−51 | −2.756e−51, −8.662e−51 | 일치 |
| 구조 | 목록 10, 항목 47, 식 4 | 같음. 라벨 4개 모두 참조, 원고와 충돌 없음, bib 누락 3건 | 일치 |
| 바인딩 | 개정 기록 manifest | 69건 SHA-256 일치, revision_5 블록에 "see REQUEST285 errata 14" | 일치 |
| 정적 SEP (N1, 새) | 없음 | \|Δ\|&lt;2.3e−6과 \|a_o\|=3.54e−3이면 \|a_p\|≲4.2e−3 | 새 지적 |
| 외측 순간 항 (N2, 새, 추정) | 없음 | 4×4.646e−10이므로 ≈1.86e−9 a_p² | 새 지적 |

## 5. 확인하지 못한 것

- Pandoc 변환과 `verify_unified_paper.py`, TeX 빌드. Pandoc이 PATH에 없고, 검증은 파일을 씁니다.
- 정적 SEP 문헌 원문(Voisin 2025, Archibald 2018). 저장소 REQUEST5에 기록된 값만 썼습니다. Δ≈a_o(a_p−a_i)는 §4.3 쌍 인자에서 제가 유도한 대응입니다.
- N2의 e_out. REQUEST274 18행의 ≈0.0354를 썼고, 외측 추정은 저장된 계산이 아닙니다.
- 단계274·279가 저장소만으로 다시 실행되는지. 두 단계는 배경을 `verification` 모듈 사슬에서 읽는데, 끝까지 추적하지 않았습니다.
- §4.3의 강성 조건이 실제로 충족되는지. 이 절에서는 κ가 정해지지 않습니다.
- 여덟 manifest 밖의 결과 manifest들이 모두 스크립트를 해시로 나열하는지.
- 다른 심사자의 원문. 열지 않았고, REQUEST285의 요약만 읽었습니다.
