# Physical Review D 제출 준비: PDF, 투고처, 커버레터

분류: Imported from prior work. **사용자는 2026-09-28에 PDF를 만들고, 투고처를 조사하고, 커버레터까지 만들어 제출 직전까지 진행하라고 지시했다.** 투고처는 Physical Review D로 정했다. PDF, 저널용 소스, 커버레터, 평문 초록, 체크리스트를 준비했다. 실제 제출과 공개 저장소 push는 하지 않았다. 둘 다 저자 계정이 필요하거나 외부 공개에 해당하기 때문이다. 기록=2026-09-28 KST.

## 한 일

분류: Imported from prior work.
1. **TeX 설치:** 이 호스트와 WSL 모두 TeX가 없어서, README가 대안으로 적은 Tectonic 0.17.0을 사용자 영역에 설치했다.
   - 공식 GitHub 릴리스의 windows-msvc zip(21 MB)을 받았고, GitHub API의 SHA-256 digest와 일치함을 확인했다.
   - 설치 경로는 `C:/Users/lpaiu/AppData/Local/Programs/tectonic`이다.
2. **빌드 중 발견한 오류:** §4.6 척도 목록의 `(K\(\approx\)934)`가 LaTeX 오류를 냈다.
   - 원인: Pandoc은 닫는 `$` 바로 뒤에 숫자가 오면 수식으로 인식하지 않아 `\$`로 이스케이프한다.
   - 초안과 원고 모두 `(\(K\approx934\))`로 고쳤다. 같은 유형은 원고에 이것 하나뿐이었다.
   - 앞선 Pandoc 점검은 목록과 식 개수만 셌기 때문에 이 오류를 놓쳤다. 컴파일까지 해야 잡힌다.
3. **제출용 원고 수정:**
   - 머리의 Status 줄을 지웠다. 제목 블록의 "Unified revised manuscript"가 첫 투고에서 재투고처럼 읽히기 때문이다.
   - 데이터 가용성 앞에 "Use of AI tools" 절(7절)을 넣었다. APS 정책에 따른 것이다.
4. **PDF(30쪽):** 경고는 underfull 1건뿐이고 미정의 인용·참조는 없다.
   - pypdfium2를 사용자 영역에 설치해 쪽을 PNG로 렌더링했다. 제목 쪽, §4.6(9–15쪽), AI 절(25쪽)을 눈으로 확인했다.
5. **패키지:**
   - `output/pdf/free-fall-identifiability.pdf`
   - 저널용 소스 zip. 풀어서 따로 컴파일해 같은 크기의 PDF가 나옴을 확인했다.
   - 커버레터(1쪽), 평문 초록, 제출 체크리스트(`output/submission/`)
6. **manifest 보호:** `paper/package_revision.py`는 `revision-manifest.json`을 9월 9일 목록으로 통째로 다시 쓴다. 그대로 돌리면 그 뒤의 바인딩 수천 건이 사라진다. 그래서 쓰지 않고 `phase289-package.py`로 패키징했으며, README에 경고를 남겼다.

## 투고처 조사

분류: Imported from prior work (공식 안내, 2026-09-28 확인).
- **Physical Review D:**
  - 범위와 기준: 중력, 밀집천체, 스칼라-텐서 이론. "significant contribution"을 요구한다([about](https://journals.aps.org/prd/about)).
  - 저자 안내([authors](https://journals.aps.org/prd/authors)):
    - 첫 심사는 PDF만으로 시작한다. REVTeX은 권장이고 필수가 아니다.
    - 연구 논문에는 길이 제한이 없다.
    - 데이터 가용성 진술(DAS)이 필수다.
    - AI의 실질적 사용은 논문 안에 공개해야 한다.
    - 소속과 ORCID가 필요하다.
    - 의무 게재료는 없다. APC와 인쇄 컬러는 선택이다.
- **APS AI 정책(2026년 6월 갱신)** ([정책](https://journals.aps.org/authors/ai-based-writing-tools), [발표](https://www.aps.org/about/news/2026/06/releases-updated-ai-policy-journals)):
  - 과학적 추론, 주장 작성 같은 실질적 사용을 허용한다.
  - 도구 이름·버전, 도움의 내용, 지시·검증 방법을 공개해야 한다.
  - 연구에 쓴 것은 방법이 기술된 곳에, 그 밖의 사용은 Acknowledgment에 적는다.
- **Classical and Quantum Gravity** ([about](https://publishingsupport.iopscience.iop.org/journals/classical-and-quantum-gravity/about-classical-quantum-gravity/)):
  - 게재료가 없고, 어떤 TeX 형식이든 제출할 수 있다.
  - IOP의 생성형 AI 정책([정책](https://publishingsupport.iopscience.iop.org/questions/generative-ai-tools/))은 공개를 Acknowledgements에 요구한다. 또한 심사 답변 작성에 AI가 참여하는 것을 제한한다.
- **판단(Conjectural):** 과학적 범위는 두 저널 모두 맞는다. 이 연구는 AI 에이전트를 계산과 원고 작성, 심사에 폭넓게 썼으므로 APS 정책이 더 잘 맞는다. PRD를 1순위, CQG를 대안으로 둔다.
  - 이전 검토(docs/submission-review-2026-09-09.md)는 당시 Paper B의 큰 수치 주장 때문에 PRD를 늦추라고 권했다. 지금 원고는 그 주장을 조건부로 낮췄다.

## 저자가 해야 할 일

분류: Conjectural. 자세한 절차는 `output/submission/submission-checklist-prd.md`에 있다.
1. **AI 사용 절:** 사실관계를 확인하고, 초기 단계 모델의 이름과 버전을 추가한다.
2. **소속·이메일·ORCID:** 커버레터의 빨간 칸과 원고 Author 줄을 채운다.
3. **공개 저장소 갱신:** 원격 main의 마지막 push는 2026-07-12이고, 로컬이 1,137커밋 앞서 있다. push하지 않으면 데이터 가용성 진술이 가리키는 내용이 공개되지 않은 상태가 된다. 외부 공개라서 사용자 결정 사항이다.
4. **확인과 선택:** 투고 이력과 이해상충을 확인한다. 추천 심사자는 선택이다.
5. **제출:** https://authors.aps.org/Submissions/ 에서 한다.
