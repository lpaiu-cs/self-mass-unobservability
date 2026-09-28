당신은 물리·방법론 논문의 독립 심사자다. 다른 심사자의 의견은 보지 않았고 보지도 않는다. 저장소의 어떤 파일도 수정하지 말고(읽기 전용), 2분 넘게 걸리는 계산은 하지 않는다. JSON을 읽는 짧은 Python 확인은 해도 된다.

저장소 루트: E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672 (git 저장소)

## 심사 대상

통합 원고 `paper/manuscript.md`에 커밋 a322c6462로 들어간 변경이다. `git show a322c6462 -- paper/manuscript.md paper/references.bib`로 차이를 볼 수 있다.
- 새 절 "### 4.6 A worked white-dwarf charge: retarded scattering memory and its orbital-timescale collapse"
- 초록에 더한 한 문장("A worked white-dwarf carrier shows ...")
- 논의(Section 6)에 더한 문단("Section 4.6 works through ...")
- 데이터·코드 가용성에 더한 문장, 날짜
- `paper/references.bib`의 세 항목(kaplan2014j0337, ransom2014triple, bertotti2003cassini)

맥락을 위해 원고 전체, 특히 Section 3–6을 읽는다. 생성된 `paper/main.tex`의 해당 부분도 확인한다.

## 근거 자료

- 연구 노트(한국어): `notes/REQUEST268_*` ~ `notes/REQUEST281_*`. 특히 270–280.
- 게시 manifest와 결과: `outputs/direct-eos-gr33/`의 다음 manifest와 폴더(result.json, 스크립트, 로그).
  - native-quad-refined-primary, native-eos-sensitivity, native-charge-mechanism, native-structure-eft-boundary
  - native-tidal-photosphere, native-final-closure, native-closure-transit, native-atmosphere-reconstruction
  - native-conditions-closed, native-manuscript-integration
- 절 초안: `docs/white-dwarf-free-fall-charge-section.md`
- 저장소 규칙: `AGENTS.md`. 주장 표지 Proven / Imported from prior work / Conjectural / Counterexample candidate의 뜻이 여기 있다.

## 심사 항목

1. 수치 충실도: §4.6과 초록·논의 추가문의 모든 수치를 근거 자료와 대조한다. 불일치가 있으면 위치와 근거(파일·JSON 키)를 적는다.
2. 주장 표지: 각 문단의 표지가 정당한지 본다. Proven은 실제 증명이어야 하고, 계산 결과는 Counterexample candidate, 가정에 기댄 추정은 Conjectural이어야 한다. 잘못 붙은 표지를 지적한다.
3. 물리·논리 검증. 특히 다음을 따진다.
   (a) 자유낙하 닫힌 식
   (b) 끝점 값을 지연 단극장(광속 횡단 창)으로 해석한 논리와 창 공식
   (c) 전체 별 선형 단열 반경 모형의 타당성과 "영구 전하 0"의 조건
   (d) 보정 항: GR 고정점 한계 2πGρT², 평탄 외부 정확성, 곡률 꼬리, 짝함수 논증, φ_∞² 스케일
   (e) 관측 대기 재구성과 비회색 LTE 보정
   (f) 궤도 시간척도 붕괴: 반경, 조석 선택 규칙, 진행파 감쇠 깊이, 구조 감수율
   (g) J0337 상한: δΔ≤|α_o||κ_lag|δφ_mod, Cassini 한계의 사용, Section 5.3의 U=1.68e−9와 K≈934 시나리오 1.57e−7과의 비교, "지연은 중성자별 전하 쪽이어야 한다"는 결론
4. 과장·과소 진술: 근거를 넘는 주장, 한계 누락.
5. 원고와의 정합성: 표기(α, β, φ), 절 참조, Section 3–6과의 모순 여부, 초록·논의 추가문이 §4.6을 정확히 요약하는지, main.tex 변환(목록·수식·인용).
6. 문장과 구성의 명료성.

## 출력 형식 (한국어, 원고 인용은 영어 원문 15단어 이내)

1. 종합 권고: 수락 / 경미 수정 후 수락 / 주요 수정 / 거부 중 하나와 한두 문장의 이유
2. 지적 사항: 심각도 순(차단 / 주요 / 경미). 각 항목에 위치, 문제, 근거, 수정 제안을 적는다.
3. 직접 확인한 수치 목록: 일치 / 불일치
4. 확인하지 못한 것과 그 이유
