# 백색왜성 절의 원고 재통합: 반영 확인 재심과 6판

분류: Imported from prior work. **5판의 반영 확인 재심 결과는 fable5.1 수락, gpt-6-astra와 opus5.5 경미 수정 후 수락이었다. 셋 모두 수락 또는 경미 수정 후 수락이므로, 경미 지적을 반영한 6판을 통합 원고 §4.6으로 다시 넣었다.** 사용자 결정(단계282: 되돌리고 수정 후 재심, 통과하면 재통합)에 따른 것이다. 새 결합 계산은 없고 저장 기록과 짧은 산술만 썼다. 기록=2026-09-28 KST.

## 반영 확인 재심 결과 (5판)

분류: Imported from prior work. 세 심사자는 같은 지시문(`confirm285-prompt.md`)을 받았다. 자기 4판 재심 원문과 4판→5판 차이만 기준으로 삼았다. 원문은 `outputs/direct-eos-gr33/native-manuscript-reintegration/`에 보존한다.

| 심사자 | 권고 | 남은 지적 |
|---|---|---|
| fable5.1 | 수락 | 선택적 문구 3건: −0.81%의 비교 대상, 데이터 가용성의 재실행 불가 범위, 초록의 t→∞ 표현 |
| gpt-6-astra (Codex CLI, 읽기 전용) | 경미 수정 후 수락 | ① 순간 응답의 전체 상한과 0.5 문턱의 범위 ② 요약의 "느린 열 상태가 존재한다" ③ 재현 안내의 예외 ④ 정오표 13의 1.8 ms |
| opus5.5 | 경미 수정 후 수락 | N1 정적 SEP 한계와의 양립, N2 외측 백색왜성의 같은 순간 응답, N3 구조 변조의 합계와 "its lagged part", N4 가장 작은 척도, N5 초록의 영구 전하, N6 판독 기호, N7 정오표 누락, N8 REQUEST285 세부 |

## 6판 반영

A=gpt-6-astra, O=opus5.5, F=fable5.1이다.

| 지적 | 반영 |
|---|---|
| 순간 응답의 상한과 문턱 (A①, O·N4) | 전체 상한 1.24e−9 a_p²+8.2e−10\|a_pa_o\|를 적고, 가장 작은 저장 척도에만 \|a_p\|≳0.5에서 닿는다고 적었다. |
| 정적 SEP 양립 (O·N1) | 원고가 이미 인용하는 J0337 정적 SEP 한계(약 2.6e−6 이하)는 쌍 인자의 선도 차수 차이 a_o(a_p−a_i)를 묶는다. \|a_o\|≃\|α₀\|(Cassini 한계)이면 \|a_p\|≲4e−3이고 변조는 4e−14 아래다. \|a_p\|≳0.5이려면 \|a_o\|≲5e−6이어야 한다. §4.6과 §6에 모두 적었다. |
| 외측 백색왜성의 순간 응답 (O·N2) | 펄서–외측 쌍 인자에 같은 순간 응답이 있고 선도 차수 추정으로 크기가 비슷하다고 적었다. §6은 복수형으로 고쳤다. |
| 구조 변조의 합계 (O·N3) | 단열 부분과 완화 부분이 각각 2.13e−18 이하(합 4.3e−18)라고 적고, 지연 열 부분은 그 상한의 절반 이하라고 적었다. "at least eight"는 "about eight"로 고쳤다. |
| 요약의 열 상태 (A②) | "연속성으로 궤도 주기와 비슷한 국소 열 시간을 가진 깊은 층이 예상되지만, 그 열 응답의 극점 구조와 전하 결합은 계산하지 않았다"로 고쳤다. |
| 초록의 영구 전하 (O·N5, F3) | "t→∞에서, 느린 열 모드까지 모든 모드가 감쇠하면"으로 고쳤다. |
| 판독 기호 (O·N6) | 전체 별 모형의 첫 판독 언급을 𝒬_c로 적었다. |
| 데이터 가용성 (O·N7, A③, F2) | 대체되는 필드를 "verdict, classification and derived-bound fields"로 넓혔다. 노트 범위는 REQUEST286까지로 했다. 재실행 불가는 "끝점 실행의 런타임 배열을 읽는 단계"로 좁혔다. 끝점 항목에 Born 첫 반환 비를 넣었다(phase268-born.log). |
| −0.81%의 비교 대상 (F1) | 본문은 그대로 둔다. 전체 별 모형의 𝒬_c를 질량 정규화된 결합 끝점 𝒬와 비교한 값이다. 𝒬_c끼리 비교하면 −0.80%이고, 차이 2.2e−5는 반올림 밖이다. |

## 새 산술

분류: Proven. 저장값과 인용 한계를 대입한 산술이다.
- 정적 SEP:
  - 한계는 \|Δ\|<2.6e−6(Archibald 2018)이다. Voisin 2025는 1.5e−6(행성 가설)과 2.3e−6(적색 잡음)이다. 여기서는 가장 약한 값을 쓴다.
  - \|a_o\|=3.5356e−3이면 \|a_p−a_i\|<7.354e−4이다. a_i≃α₀이므로 \|a_p\|<4.271e−3이다.
  - 순간 변조는 1.24e−9×(4.271e−3)²+8.2e−10×4.271e−3×3.5356e−3=2.26e−14+1.24e−14=3.5e−14<4e−14다.
  - 가장 작은 척도의 문턱은 \|a_p\|≥0.476이다. 이때 \|a_o\|<2.6e−6/0.476=5.5e−6이다.
- 외측 순간 응답(선도 차수 추정):
  - 입력은 tides.json의 외측 장 3.4312e−9와 §5.9 매개변수 집합(m_p=1.4378 M_⊙, m_i=0.19754 M_⊙, m_o=0.4101 M_⊙), e_out=0.0354(REQUEST274)이다.
  - 펄서가 외측 위치에 만드는 장은 Gm_p/(c²a_out)=3.4312e−9×1.4378/0.4101=1.2030e−8이다.
  - 외측 이심률 변조는 1.2030e−8×0.0354=4.26e−10이다. 쌍극 항(f=m_i/(m_p+m_i)=0.1208, a_in/a_out=0.02706)은 3.93e−11이다.
  - 합은 4.65e−10\|a_p\|이다. 여기에 \|β_s\|=4를 곱하면 1.86e−9 a_p²로, 내측의 1.24e−9 a_p²와 같은 차수다.
- 구조 합계: 2×2.1286e−18=4.26e−18이다. 가장 작은 척도와의 비는 2.794e−10/4.26e−18=6.6e7, 곧 약 8자릿수다.

## 추가 정오표

분류: Imported from prior work. REQUEST284의 1–6번과 REQUEST285의 7–15번에 이어서 적는다. 해시로 묶인 원문은 고치지 않는다.

16. **manifest 수치·분류 필드 (REQUEST285 정오표 14 보완)**
    - `native-tidal-photosphere`: `dissipation.lagged_monopole_delta_bound`=2.718e−18은 외측 항을 뺀 값이다. `margin`=6.25e8은 철회한 감도 진술이다. `final_charge_conclusion`의 "collapses to static coefficients … radial, tidal and dissipative channels"도 대체된다.
    - `native-final-closure`: `failure_ledger`의 \|κ_lag\|≳1.6e3/\|α_p\|와 "static structural response", `closure.freefall_delta_bound_phase274`
    - `native-structure-eft-boundary`: verdict의 `SIGN_STRUCTURALLY_ROBUST`. 부호 여유는 고정 핵에서의 충분조건일 뿐이다.
    - `paper/revision-manifest.json`의 같은 필드
17. **REQUEST276 47행 끝 문장**
    - 원문: "약한장 백색왜성의 정적 감수율은 \|β\|=4와 1e−8 수준의 구조 항뿐이므로, 선언 이론의 선형 약한장 영역에서는 이 가정이 성립하지 않는다."
    - 철회한다. 순간 계수 β_s는 완화 세기가 아니다. 단열 구조 계수는 열 완화의 세기를 묶지 못한다.
18. **`docs/failure-ledger-dynamic-chi.md` 2578·2588행(단계277–278 항목)**
    - "중심 도달 −7.5e−41"은 판독 시각 0.23 s의 −8.0e−41로 읽는다.
    - "(Proven)" 표지는 영평균 정리에만 해당한다. 별에 적용하는 것은 Conjectural이다.
19. **REQUEST285의 세부**
    - 머리말의 "큰 입력 배열은 large_arrays_sha256에만 해시로 있고 저장소에는 없다"는 일부 파일에 대해 틀렸다.
      - balanced-initial-state.npz, def-native-boundary-layer/geometry.npz·bank.npz, def-native-conservative-rates/thermal-refined/bank.npz, native-retained-completion/evolution/coupled-128.npz·source-128.npz는 git에 추적되어 공개돼 있다.
      - geometry.npz는 large_arrays_sha256에 없다.
      - 그래도 저장소만으로 다시 실행할 수 없다는 결론은 유지된다. 공개되지 않은 입력이 남기 때문이다: source-64.npz, field-source-64.npz, 광자 seg 파일, captures.
    - 입력 해시의 위치는 입력마다 다르다. 단계269의 EOS·광자 공유 라이브러리는 native-eos-sensitivity manifest의 large_arrays_sha256에 있다.
    - "실행 명령은 run*.sh에 있다"의 예외: native-final-closure와 native-conditions-closed에는 run*.sh가 없다. 단계276의 실행 형식은 phase276-closure.py 머리말(`python phase276-closure.py <out json>`)에 있다. 저장소 경로는 11행에 Windows 경로로 고정돼 있다.
    - 정오표 13: "2 ms 앞선"이 아니라 1.8 ms 앞선 직전 저장 표본이다. 20 ms 이후 표본은 20.2 ms부터 2 ms 간격이다.
    - 정오표 11: 8.8e−9는 REQUEST276 37·39·47행에 있다(43행이 아니다).
    - 정오표 12: 비선형 1.0e−27은 REQUEST277 26·33행에 있다.
    - 정오표 9: 인용한 "거의 −100%여야 한다"는 REQUEST272 19행의 둘째 문장이다.

## 원고 통합

분류: Imported from prior work.
- `paper/manuscript.md`: CRLF를 보존한 바이트 단위로 고쳤다. 원래 있던 LF 한 줄은 그대로 두었다. 697줄이 850줄이 됐다.
  - 날짜: 28 September 2026
  - 초록: 감쇠 스칼라 전하 문장 뒤에 Conjectural 두 문장
  - §4.5 뒤에 §4.6(6판 본문)
  - §6: 물리적 목표 문단 뒤에 한 문단
  - 데이터 가용성: 재현 명령 문단 앞에 한 문단, 목록, 한 문단
  - 서론은 고치지 않았다.
- `paper/references.bib`: a322c6462의 세 항목(kaplan2014j0337, ransom2014triple, bertotti2003cassini)을 바이트 그대로 복원했다. 새 인용 archibald2018universality·voisin2025planet·damour1992tensor는 이미 있다.
- `paper/main.tex`: Pandoc 3.11로 `paper/build_manuscript.py`를 다시 돌렸다. 고치기 전 원고에서는 기존 main.tex가 바이트 단위로 재현됨을 먼저 확인했다(sha256 29bb4eaf…). 통합 뒤 식은 29개에서 33개로, 목록은 0개에서 11개로 늘었다.
- `paper/README.md`: PDF와 제출 zip이 §4.6 이전 판이라는 줄을 넣었다. 이 호스트에는 TeX가 없다.
- `verification/verify_unified_paper.py`: 바인딩을 갱신한 뒤 통과한다. §4.6 자체는 검사하지 않는다.

## 남은 것

분류: Conjectural.
- PDF와 제출 zip의 재생성: TeX 설치가 필요하다. 사용자 결정 사항이다.
- 계산하지 않은 것: 두 비공통 쌍 채널의 타이밍 응답과 nuisance 투영, 깊은 층의 비단열 열 응답(세기와 극점 구조), 반응하는 동반성의 강성 조건
- 저장소 분류: 선택한 응답족의 지연 진폭 상한에 관한 조건부 theorem progress다. 관측 배제도 A4 유지도 아니다. 깊은 층의 열 완화 상태는 세기를 계산하지 않은 loophole 후보로 남는다.
