## 재심 보고서 — fable5.1 (개정 4판 `docs/white-dwarf-free-fall-charge-section.md`, 2026-09-27)

### 1. 종합 권고

**경미 수정 후 수락.** 3판 재심에서 낸 지적 12건은 모두 해소되었고(표지 1건은 잔여가 새 지적 1로 이어짐), 4판에서 새로 들어간 수치(순간 응답 상한, 부호 여유 99.87%, |𝒮_struct|≤8.84e−9, 두 쌍 채널 2.13e−18·7.53e−21, 2차 조석 3.7e−19, 자기 결합 ~1e−7, 이력 −8.0e−41, 표본 간격, Cassini 재척도)는 모두 저장 JSON·스크립트에서 재현된다. 차단 문제는 없다. 주요 지적은 한 건뿐이며 문장 수준이다: 순간 응답 항목이 §4.3의 feedback-stiffness 조건을 잘못 끌어다 쓰고, 그 채널의 척도(|a_p|≤1이면 §5 척도 대역 안)를 말하지 않는다. 계산은 필요 없다. 정오표는 정확하지만 REQUEST274의 공명 부재·자전 경계와 REQUEST276 39행의 일반 수동 상한·중성자별 귀속을 빠뜨렸다. 이 두 가지를 고치면 동반 수정(참고문헌 복원, README·manifest 기록, Pandoc 재생성, 검증)과 함께 통합할 수 있다.

### 2. 이전 지적(3판 재심)별 판정표

| 이전 지적 | 판정 | 이유 |
|---|---|---|
| 경미 1 기호 잔여 충돌·미정의 | **해소** | 𝒩·ℒ·𝒬·𝒮_struct·M_wd·m_wd로 교체, α₀·T₀·ρ₀·ρ_max·x′ 정의, "reuse only a_i, a_o and a_p" 문장. 원고에 \mathcal N/L/Q/S/P, α₀, ω₀, ξ, Ψ 및 eq:wd-*·eq:free-fall 라벨 없음을 grep으로 확인. 남는 겹침은 η와 §5.10의 δη_extra1뿐이며 "All other symbols are local" 문장으로 충분. |
| 경미 2 부호 여유 문장 과장 | **해소** | "below 99.87% in the maximum norm"으로 교체. m=1/Σ\|K_f\|=0.998791이므로 최대노름 &lt;m인 섭동은 \|ΣK_fδ\|&lt;1이라 부호 불변(충분조건). 99.87%&lt;99.879%라 진술이 성립. |
| 경미 3 비회색 범위 vs 1σ | **해소** | 1σ 사례 nongray_factor 1.892–2.412 → "1.9–2.4 over those 1σ cases". |
| 경미 4 두 쌍 채널 반올림 | **해소** | 2.13e−18, 7.53e−21(재계산 2.1286e−18, 7.5257e−21). δφ_mod 계수 3.08e−10·2.04e−10은 올림이며 "to leading order in the eccentricities" 명시. |
| 경미 5 3.5·4.0 ms 출처 | **해소** | 전체 별 모형 문단의 불릿으로 이동. |
| 경미 6 깊은 층 열 시간 | **해소** | "by continuity some deeper layer … Its depth was not located." (tides.json의 해당 필드 null과 일치) |
| 경미 7 자기 결합 1e−7 근거 | **해소** | "relative change of the effective gravity of at most 0.048" 삽입, Conjectural 표지. 산술은 응답 노트에 있고 재현됨(§4 참조). |
| 경미 8 초록·논의 "small static" | **해소** | "adds only a small structural term to the static response"; §6은 "scale comparison, not an exclusion". |
| 경미 9 표지 구조 | **대부분 해소** | 목록 뒤 무표지 문단 제거, 조석 항목을 Proven/Imported/Conjectural로 분리. 잔여: Instantaneous 불릿의 Proven 표지가 해석 문장(§4.3 조건)까지 덮음 → 새 지적 1. |
| 경미 10 A4 문장 문법 | **해소(삭제)** | A4 문장을 빼고 "neither establishes … nor shows that none exists"로 대체. 원고에 A4 번호 체계가 없으므로 타당. |
| 경미 11 0.43 s 해상도 한정 | **해소** | "0.25 ms apart up to 20 ms, 2 ms apart after"; transit.json 820표본(20 ms까지 80개)과 일치. 최대값은 0.2 ms 격자라고 명시. |
| 경미 12 통합 전제 (a)–(d) | **(b)(d) 해소, (a)(c) 통합 시 항목** | 서론 문장 삭제, eq:wd-modulation 본문 참조. bib 3항목 복원과 README·manifest의 PDF 이전 판 기록은 동반 수정에 명시되어 통합 단계에서 수행. |
| (초심 주요 2 잔여) 기호 | **해소** | 위 경미 1과 동일. |
| (초심 경미 7) 식 라벨 미참조 | **해소** | 네 라벨 모두 본문에서 참조. |
| (초심 경미 9) 참고문헌 | **통합 전제 유지** | 현재 references.bib에 kaplan2014j0337·ransom2014triple·bertotti2003cassini 없음(a322c6462에서 복원 예정). |

### 3. 새 지적

**차단: 없음.**

**주요 (문장 수준, 새 계산 불필요)**

1. **순간 응답 항목의 §4.3 교차 참조와 척도 누락** (125행, 136행).
   - 문제: "falls under the feedback-stiffness condition that Section 4.3 places on responsive companions"는 대상이 다르다. §4.3의 조건(−Σ_j C_j/r_pj²이 κ보다 작을 것)은 동반성 전하가 펄서 전하 변화 δQ_p에 반응해 펄서 극의 강성 κ를 재정규화하는 되먹임을 다룬다. §4.6의 순간 항은 전하를 고정한 채 궤도 운동이 만드는 장 변조 δφ_mod에 백색왜성 전하가 즉시 따라가 두 비공통 쌍 인자 a_pa_i, a_oa_i를 구동과 같은 위상으로 변조하는 효과다. 감수율은 같지만(단위 장당 C_i=β_s m_wd) 관측량이 다르다. 본문은 되먹임 강성 비율을 수치로 주지 않고, 이 동상 채널의 척도도 §5와 비교하지 않는다. 본문 스스로 택한 \|a_p\|≤1을 쓰면 a_pδa_i≤1.24e−9로 §5 저장 척도 대역(2.8e−10–1.57e−7) 안, K=1 값의 약 4배다. 지연이 없으므로 Section 3의 언어로는 c_Y형 정적 계수이지 상태가 아니지만, 그 말을 해야지 되먹임 조건을 가리켜서는 안 된다. 또 이 해석 문장이 Proven 표지 아래 있다.
   - 근거: 원고 315행 "A responsive companion with susceptibility \(C_j\) produces a leading feedback stiffness shift"; §3.1 q=c_YF+c_χχ; §5.4 "co-fitted instantaneous contribution"; tides.json 변조 계수.
   - 수정 제안: 마지막 문장을 다음 취지로 바꾸고 해석 부분은 Conjectural로 표지한다. "This is a zero-lag coefficient in the sense of Section 3, not a state. Its in-phase modulation of the pulsar–inner pair factor, up to 1.24e−9 a_p², was not projected through the timing model either. In the Section 4.3 model the same susceptibility, C_i=β_s m_wd, enters only as the responsive-companion feedback stiffness." 되먹임 비율을 적고 싶다면 약한 장에서 κ=1/(\|β_s\|m_p)로 두어 C_i/(r_pi²κ)=β_s²m_wd m_p/r_pi²≈4e−13이 나온다(내 추정, 어느 노트에도 없으므로 저자 확인 필요). 136행 "these two non-common channels"에 순간 항도 포함됨을 한 구절로 밝힌다.

**경미**

2. **초록과 요약의 가정 수 불일치** (138행, 155행). "bounded under an assumed single-relaxation form"만 적었는데 상한은 두 가정(완화 세기 ≤ 𝒮_struct, 단일 시간 상수)에 걸린다. 바로 뒤 가정 목록과 §6 동반 문장("Under the stated assumptions")은 둘 다 담는다. "under two assumed relaxation properties" 또는 "under the stated relaxation assumptions"로 통일.
3. **정오표 누락** (REQUEST284 정오표 절). 정확하지만 다음이 빠졌다. (i) REQUEST274 48행·60행과 REQUEST276 37행의 "진행파 영역이라 이산 공명이 없다", 공명 포착 부재, 자전 주기 ≳45분 경계 — 3판부터 본문은 감쇠 깊이 식이 적용 조건 밖이라 공명을 배제하지 못한다고 적는다. (ii) REQUEST276 39행의 일반 수동 상한(정적 감수율 전체 완화 시 \|δΔ\|≤4.4e−12\|α_p\|), "내측 백색왜성의 스칼라 전하로는 도달할 수 없다", "중성자별 전하…에서 와야 한다" — 초심(astra ②) 뒤 본문에서 삭제됐지만 정오표에 없고 closure.json은 full_beta_relaxing_delta_bound를 그대로 저장한다. 데이터 가용성 문장이 "whose last note lists errata to earlier notes"라고 안내하므로 두 항목을 추가한다. 참고: `docs/failure-ledger-dynamic-chi.md` 2563행(단계272–273 항목)은 아직 옛 부호 문장과 "정적 계수로 붕괴한다(no-go 경계)"를 담고 있다. 단계284 항목이 뒤에 붙어 조건부로 정정하므로 append-only 관행상 허용되지만, 정오표에서 이 항목도 가리키면 좋다.
4. **Born 비교의 노름 표현** (34행). "as a maximum over the stored space–time array"에서 저장 필드는 native-refined-primary result.json의 `maximum_first_return_over_incident`=1.738e−10뿐이며 노름은 기록되지 않았다. "the stored maximum of the first-return Born field relative to the incident one" 정도로 낮추면 기록과 정확히 맞는다.
5. **116행 표지**. 쌍 인자가 a_pδa_i, a_oδa_i로 갈라진다는 대수는 Proven 수준인데 Conjectural로 묶였다. 과소 표지라 무해하며 선택 사항.
6. **통합 체크리스트**. 데이터 가용성의 "all are bound in the revision manifest"는 현재 참(현 초안 SHA-256 a5984ba2…와 REQUEST284 해시가 manifest에 있음). 지적 1·2를 반영하면 초안 해시가 바뀌므로 통합 커밋에서 다시 묶는다. 1.24e−9는 1.2304e−9의 이중 올림이며 1.23e−9로도 상한이 성립하지만 "at most"라 허용.

### 4. 직접 확인한 새 수치

- **순간 응답 상한**: 4×3.0761e−10=1.2304e−9 ≤1.24e−9 ✓; 4×2.0311e−10=8.124e−10 → 8.2e−10 ✓.
- **부호 여유**: kernel.json sign_margin_Linf=0.998791 → "99.87%" ✓. 정오표의 반례: C₊=−positive_share·q=1.415e−54, C₋=1.000605·q=−2.339e−51; 1.999C₊+0.001C₋=+4.9e−55 ✓(최대노름 0.999&gt;m).
- **𝒮_struct**: kappa_struct_max=8.83658e−9 → "≤8.84e−9" ✓; /4=2.2091e−9 → "2.21e−9 of \|β_s\|" ✓; eps_per_dphi=4\|α₀\|+2\|α₀β\|=0.048(스크립트 84행) ✓.
- **두 쌍 채널**: 8.83658e−9×(0.883889)²=6.9037e−9(closure.json freefall_kappa_at_cassini와 일치); δφ_mod=3.0761e−10+2.0311e−10×3.5356e−3=3.0833e−10; 펄서–내측 2.1286e−18 → 2.13e−18 ✓; 내측–외측 ×3.5356e−3=7.5257e−21 → 7.53e−21 ✓; log10(2.794e−10/2.13e−18)=8.12 → "at least eight orders" ✓.
- **2차 조석**: 스크립트 101행 \|dq/dε\|·(q(R/a)³)²·6e_in=1.841e−7×(2.204e−5)²×6×6.92e−4=3.713e−19 ✓(이심률 주기 성분).
- **자기 결합**: 0.048×10^5.740/(4c²×1e−3×7.3975e−8 cm⁻¹)=9.92e−8; φ₀′ 항 3.3e−8과 합해 1.3e−7 → "order 1e−7" ✓.
- **Cassini**: 2σ 하한 γ−1=−2.5e−5 → \|α₀\|=3.5356e−3, φ_∞=8.8389e−4 ✓; (2.1±2.3)e−5 표기 ✓.
- **판독 이력**(transit.json 820표본; 0.25 ms×80 → 20 ms, 이후 2 ms ✓): 선형 보간 0.02 s −7.04e−47, 0.1 s −2.97e−43(→−3.0e−43 ✓), 0.23 s −8.04e−41(→−8.0e−41 ✓; 0.2305 s −8.19e−41), 0.35 s −2.65e−39(→−2.7e−39 ✓); 3.5 ms −2.756e−51, 4.0 ms −8.66e−51 ✓; 마지막 음수 0.4262 s·첫 양수 0.4282 s(→≃0.43 s ✓); 저장 표본의 0.35–0.47 s 범위 −3.16e−36…+3.31e−36(정오표 ✓); 최대 3.489e−36 at 0.4626 s는 0.2 ms 격자 ✓; 외향 진폭 9.9e−27≈1e−26, 1.47e−11→1.5e−11 ✓.
- **질량 항**: −2.33502e−51/−2.33497e−51 = 1+2.24e−5 ✓. **고정점**: 2πGρ_max t_e²=4.53e−18 ✓. **열 시간**: 0.41×2.5≈1.0 s, 3.98 s at 624.45 km ✓.
- **광구**: 1σ nongray_factor 1.892–2.412(중심 2.119), gray 1.625–6.172(3.245), 흐름 오차 최대 5.557e−5, 불투명도 비 0.545–0.736 ✓.
- **§5 척도**: Table 1(2 d) 3.534e−9/1.680e−9/2.794e−10, §5.9 4.10217e−10, §5.3 K≈934의 two-day 1.57e−7 ✓. §5.9 질량·장반경이 tides 스크립트 입력과 일치 ✓. ct_p=514.8 km ✓.
- **기록 결속**: 현 초안·REQUEST284 해시와 native-manuscript-rereview 파일이 revision-manifest.json에 있음 ✓. 유지 문서 6개에 단계284 조건부 재분류 항목이 추가됨 ✓. §4.6 내 목록 10개·항목 46개·번호 식 4개(check284 보고와 일치).

### 5. 확인하지 못한 것

- Pandoc 변환과 `verify_unified_paper.py` 통과(통합 후 확인 사항; 실행하지 않음).
- Born 1.74e−10의 노름(저장 필드명만 있음).
- 지적 1의 되먹임 강성 비율 ~4e−13(내 추정, 노트 미기재).
- "열 완화 세기 ≤ 𝒮_struct" 가정의 성립 여부(계산되지 않은 가정으로 명시됨).
- 0.43 s 이전의 0.2 ms 해상도 부호 이력, 열 시간이 궤도 주기와 같아지는 깊이(기록에 없음/null).
- REQUEST278 표 불일치의 원인(첫 실행 산출물 덮어씀, 저자 설명대로 재현 불가).
- 문헌 원문(Kaplan 2014, Ransom 2014, Bertotti 2003, Damour–Esposito-Farèse 1992)의 웹 대조.
