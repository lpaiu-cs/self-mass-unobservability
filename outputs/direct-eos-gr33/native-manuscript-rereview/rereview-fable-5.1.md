# 재심 보고서 — fable5.1 (개정 3판 `docs/white-dwarf-free-fall-charge-section.md`, 2026-09-27)

## 1. 종합 권고

**경미 수정 후 수락.** 이전 심사의 주요 지적 7건은 모두 해소되었거나(6건) 잔여가 경미 수준으로 줄었고(기호 1건), 경미 지적 9건도 해소되었다. 개정으로 새로 들어간 수치(판독 시각 이력, 3.5·4.0 ms 판독, Born 1.7e−10, 외향 진폭 1e−26, 외측 변조 2.0e−10, 두 쌍 채널 상한, 척도 범위, Cassini 재척도, 자기 결합 1e−7)는 모두 저장 manifest·JSON에서 재현된다. 차단·주요 신규 문제는 없다. 남은 것은 표기·표지·문장 정밀도와 통합 전제(참고문헌 복원 등)이며 새 계산은 필요 없다. 원고 통합은 아래 경미 항목을 반영한 뒤 가능하다.

## 2. 이전 지적별 판정표

| 이전 지적 | 판정 | 이유 |
|---|---|---|
| 주요 1 시간 이력(원천/판독 시각 혼동) | **해소** | 판독 시각 t_R로 통일. 원천 시각(중심 0.2305 s, 출사 0.4627 s)과 판독 값(0.1/0.23/0.35 s), 부호 유지 ≃0.43 s, 최대 2R_d/c≃0.46 s가 `transit.json`과 일치(§4 참조). "twelve orders"와 값의 어긋남도 사라짐. |
| 주요 2 기호 충돌·미정의 | **부분 해소** | β_s, α_s, a_i·a_o·a_p(§4.3과 일치), ds²의 a·B, x, Ψ, R_d, t_p, t_e, m_i, χ_struct가 정의됨. 잔여: q는 §3의 판독 q(t)=c_YF+c_χχ와, B는 §4.3의 B=a_w²/(κm_p)와 여전히 겹친다. α_0, T_0는 미정의. 따라서 "avoids the symbols of Sections 3–5"는 그대로는 부정확. → 새 지적 경미 1. |
| 주요 3 무표지 문단·A4 미정의 | **해소** | "Boundary statement" 삭제. A4를 AGENTS.md 정의("no orbital-timescale internal state variable in the free-fall sector")대로 요약 문단에 명시. |
| 주요 4 과도 판독값의 "charge" 용어·정밀도 | **해소** | "readout"으로 통일, 순간 분극 제외와 Born 1.7e−10·외향 진폭 1e−26 명시, 유효숫자 두 자리, 크기 문단은 "no conclusion below uses them"으로 격리. 초록은 "indicates"와 "scalar-driven displacement of the carrier's own matter". |
| 주요 5 표지 오류 | **해소** | 5.2e−5 재현은 원고 28행 규약(저장소 수치 실험=Imported)에 맞게 Imported로 둔 것이 타당(내 제안보다 규약이 우선). 보정 항목별 표지, 초록·논의 Conjectural, 영평균은 "linear, adiabatic and Newtonian" 모형 문장 뒤에 Proven/Imported/Conjectural로 분리. |
| 주요 6 J0337 상한 제시 | **해소** | 두 쌍 채널 분리, 공통 템플릿과 다름(§5.4 인용), 척도 범위 2.8e−10–1.57e−7, Cassini 2σ·damour1992tensor·선언 φ_∞ 초과·φ_∞² 재척도 명시, 외측 백색왜성 장 변조 포함, 중성자별 단정 삭제, 결론 한정. |
| 주요 7 조석 "collapses to static coefficients" | **해소** | 선택 규칙(Proven)과 평형 조석 크기만 남기고, 동적 조석은 소산적, 감쇠 깊이 8.9e6은 적용 조건 밖이라 공명 부재 미확정으로 적음. 45분 경계 삭제. |
| 경미 1 정적 창 부호 | 해소 | "+3.6e−54 … opposite sign … 1/640". |
| 경미 2 흐름 오차·0.25%/0.38% | 해소 | 5.6e−5(τ_R≤100), 핵 가중 0.25%, 최대 0.38%. |
| 경미 3 비회색 한계 | 해소 | Rosseland 0.55–0.74, 평행평판, LTE H·He 연속. |
| 경미 4 빠진 가정 | 해소 | 선형·뉴턴 반경·평탄 외부·단일 완화·핵 고정·광구 비단열 명시. WKB·확산은 감쇠 깊이 추정 자체를 무효 처리했으므로 별도 표기 불필요. |
| 경미 5 helium-core | 해소 | "helium core and a hydrogen-rich envelope (X≃0.89)". |
| 경미 6 §4.3 교차 참조 | 해소 | "supports holding the white-dwarf charges fixed in Section 4.3". |
| 경미 7 조판·식 라벨·PDF | 대부분 해소 | 2×, l=0 수식화. eq:wd-modulation은 아직 미참조. PDF·zip 재빌드는 통합 단계 사안. |
| 경미 8 데이터 가용성 | 해소 | 한국어 노트, 해시 목록, 검증기 범위 명시. |
| 경미 9 참고문헌 | 통합 전제로 전환 | 현재 `paper/references.bib`에 kaplan2014j0337·ransom2014triple·bertotti2003cassini가 없다(되돌림으로 삭제, a322c6462의 335–357행에 존재). 통합 시 복원 필요. |

## 3. 새 지적

**차단: 없음. 주요: 없음.**

**경미**

1. **기호 잔여 충돌·미정의** (11–15행, 21행, 91행). q는 §3 판독 q(t), B는 §4.3 응답 계수 B와 겹침. α_0(=α_s(φ_∞)), T_0(배경 응력 대각합)는 미정의. 수정: 판독을 δa_i 또는 q_i로, 계량 함수를 다른 글자(예: e^{2Λ})로 바꾸고, α_0·T_0를 한 구절로 정의하며, "avoids the symbols" 문장을 "reuses only a_i, a_o, a_p from Section 4.3"로 고친다.
2. **부호 여유 문장의 과장** (72행). "so no positive density profile reverses the sign". 근거 `phase272`의 m=1/Σ|K_f|=0.9988&lt;1이므로, 기여 층의 밀도를 거의 전부 제거하면(양의 밀도 유지) 뒤집을 수 있다. 이전 판의 "a relative density change of 0.9988 … is needed"가 정확하다. 수정: 그 문장으로 되돌린다.
3. **비회색 인자 범위와 1σ의 불일치** (74–75행). 앞 문장은 Kaplan 1σ(1.6–6.2)인데 "a further 1.7–2.8"은 2σ 사례를 포함한 범위다(`nongray.factors`: 2σ 1.716–2.797, 1σ 1.892–2.412). 수정: 1σ 기준이면 1.9–2.4, 아니면 "over the 1σ and 2σ cases"로 명시.
4. **두 쌍 채널 상한의 반올림** (91행). 저장값으로 재계산하면 χ_struct(Cassini)=6.904e−9, δφ_mod=3.076e−10+2.031e−10×3.536e−3=3.083e−10, δa_i=2.13e−18, ×|a_o|=7.53e−21다. 2.2e−18·7.6e−21은 올림이다. "at most"이므로 허용되지만 재계산 독자가 어긋남을 느낀다. 수정: 2.1e−18·7.5e−21로 쓰거나 "rounded up"을 붙인다.
5. **3.5·4.0 ms 판독의 출처** (30행). "−2.8e−51 at 3.5 ms and −8.7e−51 at 4.0 ms"는 결합 계산(t_e에서 끝남)이 아니라 뒤에 소개되는 전체 별 선형 반경 모형(`transit.json`, 0.25 ms 표본)의 값이다. 결합 판독 문단 안에 출처 없이 놓여 있다. 수정: "in the whole-star model below" 한 구절 추가.
6. **"Deeper layers do have thermal times comparable to the orbital period"** (82행). `tides.json`의 `depth_tau_th_equals_1_over_omega_km`는 모든 구동에서 null이며, 봉투 바닥(624 km)까지 τ_th≤4 s다. 연속성으로 존재는 참이지만 위치는 찾지 않았다. 수정: "by continuity some deeper layer …; its depth was not located".
7. **자기 결합 1e−7의 근거** (80행). "Including the change of the self-coupling raises this to order" 1e−7은 어느 노트에도 없는 새 유도값이다. 저장값 |ε/δφ|≤0.048(phase274)와 g=5.49e5 cm s⁻², 펄스 힘 척도 4c²φ_0×7.40e−8 cm⁻¹로 9.9e−8이 재현된다. 수정: 괄호로 "(effective-gravity change |ε/δφ|≤0.048 under a uniform shift)"처럼 유도 근거를 적는다.
8. **초록·논의의 "small static and structural coefficients"** (109, 112행). 본문 순간 계수는 |β_s|=4로 작지 않다. 변위 상태의 정적 계수(힘 비 ≲1e−7, χ_struct≤8.8e−9)를 뜻하는 것이라면 "small static and structural coefficients of the displacement state"로 한정한다.
9. **표지 구조** (17·19·21·60행, 81행). 불릿 목록 뒤의 무표지 문단 네 개는 원고의 무표지 관행(수식 직후 이어짐)과 다르다. Tidal 불릿은 Proven 표지 아래 Imported 수치(k_2, 9.1e−12, 1.5e−5, 8.9e6)와 Conjectural 판단을 함께 담는다. 수정: 문단을 앞 표지 문단에 합치거나 표지를 붙이고, Tidal 불릿 안에서 Imported/Conjectural을 구분.
10. **A4 문장 문법** (93행). "names the absence of such a state in the free-fall sector assumption A4" → "names … as assumption A4".
11. **부호 유지 0.43 s의 해상도 한정** (56행). 저장 표본은 20 ms까지 0.25 ms, 이후 2 ms 간격이며 2e−4 s 이력은 저장되지 않았다(이전 심사의 미확인 사항 유지). 수정: "in the stored 2 ms samples"를 붙인다.
12. **통합 전제**. (a) bib 세 항목 복원. (b) 서론 문장의 삽입 위치 "at the end of the roadmap": 현재 서론(1–29행)에는 Section 번호 로드맵 문단이 없고 a322c6462도 서론 문장을 넣지 않았다. 위치를 정하거나 문장을 뺀다. (c) README·revision-manifest에 PDF·제출 zip이 §4.6 이전 판이라는 기록(a322c6462 README 13행과 동일). (d) eq:wd-modulation 참조 또는 라벨 삭제.

## 4. 직접 확인한 새 수치

- **판독 시각 이력** (`native-closure-transit/phase278/transit.json`, 820표본): 0.1002 s −3.005e−43(→−3.0e−43 ✓), 0.2302 s −8.094e−41·0.2322 s −8.702e−41, 0.2305 s 선형 보간 −8.19e−41(→−8.2e−41 ✓), 0.3502 s −2.685e−39(→−2.7e−39 ✓). 마지막 음수 표본 0.4262 s, 첫 양수 0.4282 s(→≃0.43 s ✓). |q| 최대 3.489e−36 at 0.4626 s ✓. 중심 도달 0.2305 s, 출사 0.4627 s ✓. 주요 모드 진폭순 주기 50.7/47.9/54.1/44.6/59.2 s(→44–59 s ✓). rms 8.11e−45 ✓. ω_0²=0.035338²=1.249e−3, 주기 177.8 s ✓.
- **3.5·4.0 ms 판독**: 0.0035 s −2.756e−51, 0.0040 s −8.662e−51 ✓ (0.25 ms 표본).
- **Born 1.7e−10**: REQUEST268 "첫 반환/입사 1.74e−10" ✓.
- **외향 진폭**: m_i=GM_g/c²=2.917e4 cm, Ψ_in=ηR_d=6.91e−21 cm. q=2.335e−51 → 9.9e−27(≈1e−26 ✓); q=3.489e−36 → 1.47e−11 ✓.
- **외측 변조 2.0e−10**: `tides.json` 1.2147e−10+8.164e−11=2.031e−10 ✓. 내측 3.076e−10 ✓(m_p/a_in×e_in=4.445e−7×6.92e−4 재현).
- **두 쌍 채널**: 위 경미 4의 계산. 2.13e−18, 7.53e−21(초안 2.2e−18, 7.6e−21은 올림). 2.8e−10 대비 1.3e8·3.7e10배 → "at least eight orders" ✓.
- **척도 범위**: 원고 426행 2.794e−10/1.680e−9/3.534e−9, 589행 4.10217e−10, 406행 1.57e−7 ✓.
- **Cassini**: γ−1 2σ 하한 −2.5e−5 → |α_0|=3.5356e−3, φ_∞=8.839e−4 ✓; χ_struct 재척도 6.904e−9(`closure.json freefall_kappa_at_cassini`) ✓.
- **자기 결합 1e−7**: 0.048×5.49e5/(4×(2.998e10)²×1e−3×7.3975e−8)=9.9e−8 ✓(유도값, 노트 미기재).
- **기타**: 정적 창 +3.622e−54, 비 639.4 ✓; −0.81%, 1.18e−6 ✓; 4.53e−18, 4.22e−6, 1.02e−27, 6.75e−11, 1+2.24e−5(phase260 1×) ✓; 0.9988, 363 km, τ 1.057 ✓; 0.25%/0.38%, 3.245, 1.625–6.172, 2.119, 5.557e−5 ✓; k_2 3.29e−4, 9.10e−12, 1.50e−5, 8.87e6 ✓; 28–47 s, 1.60e−6, 3.32–3.35e−8, 8.84e−9, 2.21e−9 ✓; 0.999996, 9.99e−9, 2.18e−10, 5.23e−5 ✓; ct_p=514.8 km, m_i/R_d=4.22e−6 ✓.

## 5. 확인하지 못한 것

- 0.43 s 이전의 2e−4 s(또는 Newmark 1e−5 s) 해상도 부호 이력: 저장되지 않음.
- 열 시간이 궤도 주기와 같아지는 깊이: 기록에 null.
- 결합 계산·전체 별 통과·대기 재구성·비회색 풀이기의 재실행(2분 한도 밖, 읽기 전용).
- 문헌 원문(Kaplan 2014, Ransom 2014, Bertotti 2003, Damour–Esposito-Farèse 1992)의 웹 대조.
- Pandoc 재생성과 `verify_unified_paper.py` 통과 여부(통합 후 확인 사항).
- 0.23 s 값 −8.2e−41이 0.2305 s 보간인지 다른 읽기인지는 저자 확인 필요(내 보간과는 일치).
