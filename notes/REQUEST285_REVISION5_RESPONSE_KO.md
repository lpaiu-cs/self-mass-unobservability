# 백색왜성 절 5판: 4판 재심 응답, 추가 정오표, 재현 경로

분류: Imported from prior work. **개정 4판의 재심 결과, gpt-6-astra는 차단을 해제했지만 주요 수정을 권고했다. opus5.5와 fable5.1은 경미 수정 후 수락이었다.** 이를 반영해 `docs/white-dwarf-free-fall-charge-section.md`를 5판으로 고쳤다. 새 결합 계산은 없고 저장 기록과 짧은 산술만 썼다. REQUEST284 정오표의 오류와 누락은 아래 추가 정오표로 고친다. 기록=2026-09-27 KST.

## 재심 결과 (개정 4판)

분류: Imported from prior work. 세 심사자는 같은 지시문(`rereview284-prompt.md`)을 받았고 자기 이전 원문 두 편만 열었다. 원문은 `outputs/direct-eos-gr33/native-manuscript-rereview2/`에 보존한다.

| 심사자 | 권고 | 핵심 |
|---|---|---|
| gpt-6-astra (Codex CLI, 읽기 전용) | 주요 수정, 차단 해제 | 정오표의 "두 가정 아래 no-go"와 붕괴 회피 필요조건이 틀렸다(대수적 반례). 경미: 초록의 세기 가정, 𝒬와 compact 판독의 구분, 고정 핵 조건, 재현 경로 |
| opus5.5 | 경미 수정 후 수락 | 주요 A: 순간 응답을 §4.3의 강성 조건에 귀속한 것이 틀렸다. 주요 B: 요약·초록의 세기 가정 누락. 주요 C: 정오표 범위 부족과 manifest 판정 문자열. 경미 15건 |
| fable5.1 | 경미 수정 후 수락 | 주요 1건: 순간 응답의 §4.3 귀속과 척도 누락. 경미: 초록 가정 수, 정오표 누락, Born 노름 표현, 표지, 통합 체크리스트 |

## 지적별 반영

A=gpt-6-astra, O=opus5.5, F=fable5.1이다.

| 지적 | 반영 |
|---|---|
| 순간 응답의 귀속 (O·A, F 주요) | 순간 응답은 Section 3의 지연 없는 계수로 적었다. §4.3의 고정 동반성 축약이 두 쌍 인자의 동상 변조를 빠뜨린다고 적었다. 펄서–내측 쌍에서 약 1.24e−9 a_p² 이하이고, \|a_p\|≳0.5이면 §5 척도에 닿으며, 타이밍 응답은 계산하지 않았다. 같은 감수율이 펄서 전하의 강성도 옮기지만 그 조건은 평가하지 않았다. 해석 문장은 Conjectural로 옮겼다. §6 추가문에도 한 문장 넣었다. |
| 세기 가정 (O·B, A, F2) | 초록과 요약에 두 가정(완화 세기 ≤ 𝒮_struct, 단일 완화 시간)을 모두 적었다. "only a small"은 "계산한 단열 구조 항은 작다"로 바꿨다. |
| 지연 상한의 표현 (A) | 완화 부분을 ΔS/(1+iωτ)로 쓰고, 직교 성분 ≤\|ΔS\|/2 ≤\|𝒮_struct\|/2가 모든 진동수에서 성립한다고 적었다. "상한은 지연의 크기를 묶을 뿐 없애지 않는다"를 넣었다. 강한 과감쇠 대신 "공명, 두 극점(관성), 부호가 섞인 다중 완화"를 상한 밖으로 적었다(O7). |
| 판독 정의 (A, O5) | 𝒬_c=−Ψ_out/m_wd와 𝒬=𝒬_c−a_i⁽⁰⁾δm_wd/m_wd를 나눴다. a_i⁽⁰⁾는 약한 장 값 α₀를 쓴다. 창 식과 전체 별 모형은 𝒬_c이고, 질량 항은 이력에 적용하지 않는다. 𝒬는 분극을 뺀 물질 응답 부분이며, 뺀 분극의 긴 파장 극한이 순간 응답이라고 이었다. |
| 고정 핵 조건 (A) | 부호 문장에 "면 핵(응답 계수와 전파기)을 고정할 때"와 ‖δρ/ρ‖∞를 넣었다. |
| Born 표현 (F4) | "저장된 첫 반환 Born 장과 입사 장의 최대 비"로 고쳤다(`maximum_first_return_over_incident`=1.738e−10). |
| 쌍 인자 대수의 표지 (F5) | Proven으로 올렸다. |
| a_p 정의 (O1) | §4.3에는 a_p가 없으므로 a_p=Q_p/m_p로 정의했다. |
| 변조 식의 차수 (O2) | "이심률과 a_in/a_out의 선도 차수"로 고쳤다. |
| 쌍 채널 문장 (O3) | 구조 척도 \|𝒮_struct\|δφ_mod≤2.13e−18을 먼저 적고, 쌍 인자에는 \|a_p\|, \|a_o\|를 곱한다고 적었다. |
| 가정 목록 (O4) | \|a_p\|≤1과 약한 장 전하 a_i=α_s(φ)를 넣었다. |
| 영구 전하 0의 시간척도 (O6) | t→∞에서만 성립한다고 적었다. 열 모드는 냉각 시간만큼 느리게 감쇠할 수 있어 1차 잔여가 관측 기간 내내 남을 수 있다. |
| 요약 결론 (O8) | "Section 3 기준이 상정하는 §5 척도의 결합을 가진 궤도 시간척도 내부 상태를 입증하지도 배제하지도 않는다"로 한정했다. 연속성으로 느린 열 상태는 있지만 결합은 계산하지 않았다고 적었다. |
| 이력 값 (O9) | "stored samples에서 보간한 값"으로 고쳤다. |
| 이온화 에너지 추정 (O10) | 광구 열 시간은 "추정으로 약 1 s", 봉투 바닥은 "추정으로 최대 약 2.5배"로 적었다. |
| 기호와 제목 (O11) | F₁·F₂→𝒫₁·𝒫₂, 펄스 위상 p→s로 바꿨다. 항목 제목 "Dissipative"는 "Structural and thermal"로 바꿨다. |
| 2차 조석 (O12) | φ∞에 비례한다고 적었다(Cassini 한계에서 3.3e−19). |
| 자기 결합 (O14) | "ε 전체의 상한으로 바꾸면 비가 1e−7 차수"로 고쳤다. 아래 산술 정정을 보라. |
| 데이터 가용성 (O·C, A) | manifest 파일명(-manifest.json)을 적었다. 리뷰 이전의 판정 필드는 정오표로 대체된다고 적었다. 쌍 채널·순간 응답 산술은 REQUEST284라고 적었다. 결과별 재현 경로는 이 노트에 두고, 이후 단계가 끝점 실행의 런타임 배열(해시만 기록, 미공개)을 읽는다고 적었다. |

## 반영하지 않은 것

분류: Conjectural.
- 되먹임 강성 비율의 추정(F, 약 4e−13): 이 절에서 κ가 정해지지 않아 싣지 않았다. 본문은 조건을 평가하지 않았다고만 적는다.
- J0337 정적 SEP 한계로 \|a_p\|를 줄이는 것(O): 하지 않았다. \|a_p\|≤1과 \|a_p\|≳0.5 문턱만 적었다.
- 두 비공통 채널의 타이밍 응답 계산: 하지 않았다. 척도 비교로 한정했다.

## 산술 정정

분류: Proven. REQUEST284의 자기 결합 산술에서 "φ₀′ 항 3.3e−8과 합하면 약 1.3e−7"은 이중 계산일 수 있다.
- ε/δφ의 2|α₀β_s| 항 가운데 |α₀β_s|=0.016 부분만으로 비를 계산하면 0.016×5.495e5/(4×8.988e20×1e−3×7.3975e−8)=3.31e−8이다. 이 값은 φ₀′ 항의 비 3.3e−8과 같다.
- 따라서 ε에 φ₀′ 항이 들어 있을 가능성이 크다. 이 해석은 opus5.5의 추론이고, ε의 항별 유도는 다시 하지 않았다.
- 본문은 "ε 전체의 상한 0.048을 쓰면 9.9e−8, 곧 1e−7 차수"로만 적는다.

## 추가 정오표

분류: Imported from prior work. REQUEST284의 정오표 1–6번에 이어서 적는다. 해시로 묶인 원문은 고치지 않는다.

7. **REQUEST284 정오표 4번과 docs 6개의 단계284 항목**
   - "두 가정 아래에서만 no-go 경계"와 "붕괴를 피하는 최소 추가 조건은 … 완화 세기가 𝒮_struct를 넘거나 단일 완화가 아닌 것"은 틀렸다.
   - 반례(gpt-6-astra): H(ω)=s/2+(s/2)/(1+iωτ), τ=1/ω_orb. 이 응답은 두 가정을 모두 만족한다(세기 s/2≤s, 단일 완화). 그런데도 궤도 진동수에서 직교 성분이 s/4>0이다.
   - 두 가정은 지연의 크기(≤\|𝒮_struct\|/2)를 묶을 뿐 지연을 없애지 않는다. 척도가 작아도 no-go가 되려면 비공통 채널의 감도 조건이 따로 필요하다.
   - 저장소 분류는 "선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress"다. 관측 배제가 아니다.
   - 연속성으로 존재하는 깊은 층의 열 완화 상태는 세기를 계산하지 않은 loophole 후보로 남는다.
8. **REQUEST284 반영표와 docs 6개의 단계284 항목: "순간 응답 … §4.3 조건"**
   - §4.3의 강성 조건은 동반성이 펄서 전하 변화에 반응하는 되먹임을 다룬다. 순간 응답 β_sδφ_mod는 그 되먹임이 아니다. 전하를 고정한 궤도 변조에 대한 지연 없는 응답이고, 두 쌍 인자를 동상으로 바꾼다.
   - 5판 본문의 서술을 따른다.
9. **REQUEST272**
   - 8행: δlnρ₀_f 표기는 유한 상대 변화 δρ/ρ로 읽는다.
   - 19행 첫 문장 "거의 −100%여야 한다": 이 문장도 부정확하다. 반대 부호 면의 밀도를 키워도 부호가 뒤집힌다(반대 부호 면만 약 1.65e3배 넘게 키우면 된다: 2.339e−51/1.415e−54=1653).
10. **REQUEST274**
    - 48행: "진행파 영역 … 이산 g모드 공명과 공명 포착이 없다"와 자전 주기 ≳45분의 경계는 확정되지 않았다. 약감쇠 식을 적용 조건 밖에서 썼다.
    - 54행: 완화 세기를 κ_struct로 "보수적으로 둔다"는 근거 없는 가정이다. 5판에서는 계산하지 않은 가정으로 명시했다.
    - 58행: \|α_o\|≤1이라 하면서 외측 항을 뺐다. 2.7e−18이 아니라 8.8366e−9×(3.076e−10+2.031e−10)=4.5e−18이다. 공통 템플릿의 한계 1.7e−9와 비교한 "6.3e8배 작다"는 감도 진술이 아니다.
    - 60행: "no-go 경계를 조석·소산 채널까지 확장한다"는 철회한다.
11. **REQUEST276**
    - 37·43·47행: "정적 계수로 붕괴한다", "진행파 영역", "J0337 관측량을 만들지 않는다", "theorem progress(no-go)"는 철회한다. 8.8e−9는 상한 8.8366e−9를 아래로 반올림한 값이다.
    - 39행: 다음을 철회한다.
      - 일반 수동 상한(정적 감수율 전체가 완화될 때 4.4e−12\|α_p\|)과 필요 조건 \|κ_lag α_p\|≥1.56e3. 두 값 모두 α_o 채널만 쓴 값이다.
      - "공명 없는 수동 완화의 지연은 절반 이하"를 일반 원리로 쓴 것
      - "도달할 수 없다"
      - "중성자별 전하에서 와야 한다"
    - 순간 계수 β_s는 지연 없이 즉시 응답하므로 완화 세기가 아니다.
12. **REQUEST277 25·31·37행과 REQUEST280 12·18–24·28·37행**
    - 곡률 꼬리 4.2e−6은 상한이 아니라 차수 추정이다.
    - 비선형 1.0e−27은 특정 2차 항의 비교이지 판독 오차의 상한이 아니다.
    - GR 고정점 부등식은 뉴턴 반경 채널에서만 성립한다.
    - 따라서 표의 "닫힘"은 "추정했다"로 읽는다.
    - "영구 전하 0(Proven)"에서 Proven은 이산 선형계의 영평균 정리에만 해당한다. 별에 적용하는 것은 Conjectural이고 t→∞에서만 성립한다.
    - 28행 "모두 결론의 부호·no-go 판정은 바꾸지 않고"와 37행 "theorem progress(no-go)·A4 유지"는 철회한다.
13. **REQUEST278 불일치의 원인 (REQUEST284 정오표 2번 보완)**
    - 옛 값 −2.7e−43, −7.5e−41, −2.4e−39는 표시 시각보다 2 ms 앞선 저장 표본 0.0982, 0.2282, 0.3482 s의 값과 정확히 같다(−2.672e−43, −7.526e−41, −2.376e−39).
    - 원인은 직전 표본을 읽은 것이다. 덮어쓴 첫 실행과는 관계없다(opus5.5 확인, 저장 표본으로 재확인).
14. **manifest 판정 필드**
    - 아래 필드는 독립 심사 이전에 쓴 판정이다. 해시 기록은 그대로 두고, 5판 본문과 이 정오표로 대체한다.
    - `native-structure-eft-boundary`: verdict `…MEMORY_COLLAPSES_AT_ORBITAL_TIMESCALES`, final_charge_conclusion의 "collapses to static coefficients"
    - `native-tidal-photosphere`: verdict `NO_GO_EXTENDED_TO_TIDAL_AND_DISSIPATIVE_CHANNELS…`, assumptions의 45분 자전 경계
    - `native-final-closure`: 다음 값들
      - verdict `…NO_J0337_OBSERVABLE__A4_HOLDS_FOR_THIS_STATE`
      - classification "no-go boundary (theorem progress)"
      - `required_kappa_lag_times_alpha_p`=1563, `full_beta_relaxing_delta_bound`=4.35e−12
      - `freefall_delta_bound_cassini`=7.508e−21(외측 항 제외)
    - `native-closure-transit`: verdict `GR_NONLINEAR_CLOSED…`, `infinity_tail_bound`(차수 추정)
    - `native-conditions-closed`: verdict `…A4_HOLDS`, classification "theorem progress (no-go)", j0337_closure
    - `paper/revision-manifest.json`의 같은 항목들
15. **유지 문서**
    - `docs/failure-ledger-dynamic-chi.md` 2563행(단계272–273 항목)의 부호 문장과 "정적 계수로 붕괴(no-go)"는 대체된다.
    - docs 6개의 단계273·274·276·280·284 항목에 있는 no-go·A4 판정도 대체된다.
    - 단계285 항목이 이를 정정한다.

## 결과별 재현 경로

분류: Imported from prior work. 각 manifest는 스크립트와 결과를 해시로 나열한다. 실행은 WSL(Ubuntu-22.04) 작업 경로 `/home/lpaiu/work/native-refined268-runtime` 등에서 했다. 큰 입력 배열은 `native-quad-refined-primary-manifest.json`의 `large_arrays_sha256`에만 해시로 있고 저장소에는 없다. 따라서 표의 "외부 입력"이 있는 단계는 저장소만으로 다시 실행할 수 없다. 실행 명령은 각 manifest의 `run*.sh`에 있다.

| 결과 | manifest | 주요 스크립트 | 결과 파일 | 외부 입력 |
|---|---|---|---|---|
| 끝점 판독(단계268) | native-quad-refined-primary | phase268-driver.py, phase268-bank.py, compare268.py | compare.json, readout-61–64 | 단계267·268 런타임(광자 seg-61–64.npz, captures), balanced-initial-state.npz |
| EOS 시험(269) | native-eos-sensitivity | phase269-eos.py, phase269-build.py, compare269.py | compare.json, comparison/source-diff-61–64.json | 단계268 런타임, EOS·광자 공유 라이브러리(gas.so, levels.so; 해시는 이 manifest의 large_arrays_sha256) |
| 기작·자유낙하(270–271) | native-charge-mechanism | phase270-components.py, phase270-memory.py, phase271-freefall.py | phase270/components-T.json, phase271/{1x,2x,4x}-charge·compare.json | readout268-quad64-work/gr/source-64.npz, field-source-64.npz, def-native-boundary-layer/geometry.npz |
| 핵·반경 경계(272–273) | native-structure-eft-boundary | phase272-kernel.py, phase273-history.py, phase273-modes.py | phase272/kernel.json, phase273/modes.json, history-*.json | 같은 단계268 배열, geometry.npz |
| 조석·구조 감수율·장 변조(274), 광구 척도(275) | native-tidal-photosphere | phase274-tides.py, atmos274.py, phase275-photosphere.py | phase274/tides.json, atmos.json, phase275/photosphere.json | 단계268 런타임의 배경 |
| Cassini 재척도(276) | native-final-closure | phase276-closure.py | closure.json | 없음(저장 JSON) |
| 보정·질량 정규화·이력(277–278) | native-closure-transit | phase277-bounds.py, phase278-transit.py | phase277/bounds.json, phase278/transit.json | readout268-quad64-work/gr/field-source-64.npz |
| 광구 재구성(279) | native-atmosphere-reconstruction | kappa279.py, phase279-gray-sph.py, phase279-nongray.py, phase279-combine.py | phase279/combined.json 등 | 단계268 런타임의 배경 |
| 쌍 채널·순간 응답 산술(284) | notes/REQUEST284 | (산술) | — | tides.json, closure.json, kernel.json |

목록은 `scan285-inputs.py`로 manifest에서 뽑았다. 외부 입력은 스크립트가 읽는 경로를 문자열로 찾은 것이라 빠진 것이 있을 수 있다.

## 저장소 분류

분류: Conjectural.
- 펄스 시간척도: 끝점 판독은 과도 지연 산란 신호다(Counterexample candidate). 선형·양의 감쇠를 가정하면 t→∞의 영구 전하는 0이다.
- 궤도 시간척도: 세 부분을 나눠 적는다.
  - 순간 응답은 정적(c_Y형) 계수다. §4.3의 고정 동반성 축약이 이를 뺀다.
  - 단열 구조 항은 작다.
  - 열 완화의 지연 부분은 두 가정 아래에서만 진폭이 묶인다. 이는 "선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress"다.
- 관측 배제나 A4 유지는 성립하지 않는다.
- 실패 원장 기준으로, no-go가 무너지는 정확한 단계는 열 완화 세기다. 단열 정적 계산이 이를 묶지 못한다.
- 빠진 최소 계산은 두 가지다. 하나는 깊은 층의 비단열 열 응답(세기와 극점 구조)이고, 다른 하나는 두 비공통 채널의 타이밍 응답이다.

## 확인

분류: Imported from prior work. `check284-pandoc.py`로 5판을 다시 변환했다. 결과는 목록 10개, 항목 47개, 식 4개였고 모든 식 라벨이 참조된다. 빠진 인용 키는 통합할 때 복원할 세 개뿐이다.

## 다음 단계

분류: Conjectural. 같은 세 심사자에게 반영 확인 재심을 맡긴다. 전면 재심이 아니라 4판 재심 지적의 반영 여부와 새로 생긴 문제만 본다. 셋 모두 수락이나 경미 수정 후 수락이면 단계281 절차로 원고에 통합한다.
