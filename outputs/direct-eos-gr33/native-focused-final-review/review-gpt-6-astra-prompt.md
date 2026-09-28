당신은 Physical Review D에 투고하기 직전인 물리·방법론 논문의 독립 심사자다. 다른 심사자의 의견은 보지 않았고 보지도 않는다. 읽기 전용으로 작업하라. 저장소의 어떤 파일도 수정·생성하지 말고 커밋하지 않는다. 2분 넘게 걸리는 계산은 하지 않는다. JSON을 읽는 짧은 Python 확인은 해도 된다. 최종 응답 전체가 심사 보고서다.

저장소 루트: E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672 (git 저장소, 현재 커밋 6e8110603)

## 심사 대상

- 본문: `paper/manuscript.md` (11쪽, 제목 "Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715"). 생성된 LaTeX는 `paper/main.tex`, 렌더링은 `output/pdf/free-fall-identifiability.pdf`다.
- 보충 자료(SM): `paper/supplement.md` (30쪽, `output/pdf/free-fall-identifiability-supplement.pdf`). 이전 전체 원고(`git show 39962a753:paper/manuscript.md`)에서 제목과 초록만 바꾼 것이다. 본문은 "SM Section 4.6"처럼 이 문서의 절 번호로 가리킨다.
- 제출 부속: `output/submission/cover-letter-prd.tex`, `output/submission/abstract-plain.txt`, `output/submission/submission-checklist-prd.md`(한국어).
- 경위: 저자 지시("초점을 좁히고 이미 알려진 사실은 반복하지 말 것")로 30쪽 원고를 11쪽 본문과 SM으로 나눴다. 재구성 기록은 `notes/REQUEST291_FOCUSED_MANUSCRIPT_KO.md`(한국어)와 `outputs/direct-eos-gr33/native-focused-manuscript/result.json`에 있다. 저자 측은 본문의 모든 수치가 이전 원고에 있다고 확인했다. 이 확인도 검증 대상이다.

## 근거 자료

- 저장소 규칙: `AGENTS.md`. 주장 표지(Proven / Imported from prior work / Counterexample candidate / Conjectural)의 뜻은 본문 서론 끝 문단에 있다.
- 표 1과 그림 1: `request10_external/sep_dynamic/sep_phase_marg_10_8e.json`
- 표 2와 감사 수치: `outputs/research-completion/`의 `coverage-audit.json`, `estimated-covariance-audit.json`, `comparator-audit.json`, `physical-matching.json`, `corrected-physical-drive.json`, `simultaneous-calibration.json`, `simultaneous-validation.json`, `phase-state-audit.json`, `phase-refinement.json`, `runtime12-analysis.json`, `gap-pair-audit.json`, `state-identifiability.json`
- 검증 스크립트: `verification/verify_unified_paper.py` (표 행 대조와 해시 확인)
- 백색왜성 계산의 근거: SM §4.6과 그 데이터 가용성 절이 가리키는 `notes/REQUEST244_*`–`REQUEST288_*`, `outputs/direct-eos-gr33/`의 manifest

이전 심사 회차의 보고서(`outputs/direct-eos-gr33/native-manuscript-*/`, `native-final-review/`)는 열지 않는다.

## 심사 항목

1. **수학·물리 정확성.** 본문에 남은 결과를 직접 검증한다.
   - Theorem 1과 그 증명, 첫 장애 차수, 복소 계수의 경우
   - 짝수 비교기 논증, 양의 상호적 스펙트럼과 속도 간극 부등식(식 6)과 증명
   - 모멘트-분산 등식(식 7–8), 단일 극 복원식, 가까운 두 극의 차이(식 9), 양성 없이 반송파에서 사라지는 추가항
   - 과도 응답 연산자(식 10)와 그 결론의 범위
   - 스칼라 전하 축약(식 11–12), 관성 오차 한계, 평형 자료가 완화시간을 정하지 못한다는 명제, 감쇠 행렬의 null 방향
   - 물리 구동 전개(식 13), 위상 폐합(식 14)과 원점 이동 불변성, 템플릿(식 15), 구간 정의(식 16)
2. **재구성 충실도.** 이전 원고·SM과 대조해 본문이 주장·수치를 바꾸거나 의미를 바꾸는 한정어를 잃지 않았는지, 새 주장이 섞이지 않았는지 본다. 표지가 올바른지 보고, 특히 한 문단 안에서 표지를 나눈 곳을 확인한다.
3. **수치.** 본문의 수치를 근거 JSON과 가능한 범위에서 대조한다.
4. **초점과 반복.**
   - 알려진 사실을 여전히 불필요하게 반복하는지 본다.
   - 반대로, SM 없이 본문만으로 논증을 따라가기 어려울 만큼 빠진 정의나 단계가 있는지 본다.
   - SM 절 참조가 맞는지 확인한다.
   - 서론의 새로움 진술("the proofs use only elementary tools ... What is new is the set of identifiability boundaries")이 정확하고 정직한지, 선행 문헌(유리함수 보간, 모멘트 문제, 응답 함수의 극으로 본 내부 모드 EFT)과의 관계를 충분히 밝히는지 본다.
5. **결론 강도.** 초록·서론·논의가 근거를 넘지 않는지 본다(비검출, 조건부 결과, 백색왜성의 척도 비교).
6. **제출 준비.**
   - PRD 범위 적합성, 제목과 초록, 커버레터와 본문의 일치
   - AI 사용 공개가 APS 정책(도구 이름·버전, 사용 방식, 저자의 지시·검증)을 충족하는지
   - 데이터·코드 가용성. 공개 스냅숏을 갱신하는 일은 저자 지시를 기다리는 알려진 항목이므로 지적할 필요는 없다.
   - LaTeX 변환: 수식·표·인용·참조
7. **편집 판단.** 편집자 단계 초기 반려의 위험과 그 이유, 심사로 갔을 때의 전망을 판단한다.

## 출력 형식 (한국어, 원고 인용은 영어 원문 15단어 이내)

1. 종합 권고(PRD 심사자 기준): 수락 / 경미 수정 후 수락 / 주요 수정 / 거부 중 하나와 두세 문장의 이유
2. 지적 사항: 심각도 순(차단 / 주요 / 경미). 각 항목에 위치(절·식·문장), 문제, 근거, 수정 제안을 적는다.
3. 직접 확인한 수치 목록: 일치 / 불일치
4. 초기 반려 위험과 게재 전망: 근거와 함께
5. 확인하지 못한 것과 그 이유
