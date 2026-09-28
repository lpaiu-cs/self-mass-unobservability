# Physical Review D 제출 체크리스트 (2026-09-28)

제출 직전 상태까지 준비했다. 실제 제출은 저자 계정으로 APS 제출 시스템에서 해야 하므로 하지 않았다. 아래 "제출 전 필수" 항목을 끝낸 뒤 제출한다.

**진행 상황 (2026-09-28 저자 확인 반영)**
- AI 사용 절: 저자가 사실임을 확인했다. 모델 버전(Claude Opus 5.5, GPT-6-Astra)을 넣었다.
- 투고 이력과 이해상충: 문제없음을 저자가 확인했다.
- 공개 저장소: 전체 이력(51.6 GB)은 GitHub 한도를 넘는다. 저자 결정에 따라 대형 배열을 뺀 공개 스냅숏을 push했다(`PUBLIC_SNAPSHOT.md`). 원고의 데이터 가용성 절도 이에 맞게 고쳤다.
- 최종 검토: 끝났다(`notes/REQUEST290_FINAL_REVIEW_KO.md`).
- **남은 일:** 소속·교신 이메일·ORCID(저자가 나중에 채움), 그리고 제출.

## 1. 투고처 결정

**Physical Review D (Regular Article)** 를 1순위로 권한다.
- **범위:** PRD는 중력, 일반상대론과 대안 이론, 밀집천체, 천체물리를 다룬다. 이 원고는 자유낙하·SEP 시험, 스칼라-텐서 이론, 확장 천체의 EFT, 펄서 타이밍 추론을 다룬다.
- **AI 정책:** 이 연구는 AI 에이전트를 계산과 원고 작성에 실질적으로 썼다.
  - APS의 2026년 6월 정책은 과학적 추론, 주장 작성 같은 실질적 사용도 허용한다. 대신 도구 이름·버전, 도움의 내용, 지시·검증 방법을 논문 안에 공개하도록 요구한다.
  - IOP(CQG)의 정책은 심사 답변 작성에 생성형 AI를 쓰는 것을 제한한다. 이 연구의 작업 방식에는 APS가 더 잘 맞는다.
- **형식과 비용:**
  - 첫 제출은 PDF만으로 심사가 시작된다. REVTeX은 권장일 뿐 필수가 아니고, 소스는 수락 뒤에 요구된다.
  - 연구 논문에는 길이 제한이 없다.
  - 의무 게재료는 없다. 오픈 액세스(APC)와 인쇄 컬러는 선택 사항이다.
- **대안:**
  - Classical and Quantum Gravity(Research Paper): 게재료 없음, 어떤 TeX 형식이든 제출 가능. 다만 원고에 "더 넓은 중력 물리 맥락" 요약이 필요하고, AI 공개는 Acknowledgements에 넣어야 한다.
  - General Relativity and Gravitation

## 2. 제출 파일

| 파일 | 용도 |
|---|---|
| `output/pdf/free-fall-identifiability.pdf` | 원고 PDF(30쪽). 첫 제출에 이것만 올려도 된다 |
| `output/submission/cover-letter-prd.pdf` (`.tex`) | 커버레터 1쪽. 빨간 칸을 채운 뒤 다시 컴파일한다 |
| `output/submission/abstract-plain.txt` | 제출 양식에 붙일 평문 초록 |
| `output/submission/free-fall-identifiability-source.zip` | 저널용 소스(main.tex, main.bbl, references.bib, 그림 4개). 선택 업로드이며, 이 zip만으로 컴파일됨을 확인했다 |

## 3. 제출 전 필수 (저자만 할 수 있는 일)

1. **AI 사용 절(원고 7절) 확인:**
   - 문안은 저장소 기록으로 확인되는 사실만 적었다. 저자의 역할("연구 질문·수락 기준·표지를 정하고 검증했다")은 저자가 사실인지 확인한다.
   - APS는 도구 버전을 요구한다. 초기 단계에 쓴 Codex·Claude 모델의 이름과 버전을 알면 추가한다.
   - 고칠 곳은 `paper/manuscript.md`의 "## Use of AI tools"다. 고친 뒤 아래 5절의 빌드 명령을 다시 실행한다.
2. **소속·교신 이메일·ORCID:**
   - 커버레터의 빨간 칸 세 곳을 채운다.
   - APS는 소속과 ORCID를 요구한다. 원고 제목 블록에는 지금 저자 이름만 있다. 소속을 넣으려면 `paper/manuscript.md` 머리의 Author 줄에 적는다(예: "Juneyoung Kim (소속)").
3. **공개 저장소 갱신:**
   - 원고의 데이터·코드 가용성 절과 커버레터는 공개 저장소 github.com/lpaiu-cs/self-mass-unobservability를 가리킨다.
   - 그런데 원격 main의 마지막 push는 2026-07-12이고, 로컬 main이 1,137커밋 앞서 있다. 이대로면 심사자가 §4.6 계산과 manifest를 볼 수 없다.
   - push 전에 공개해도 되는 내용인지 검토한다(심사 기록, CLI 로그, 노트 포함). 원하면 Zenodo 등으로 DOI 스냅샷을 만든다.
   - push는 외부 공개이므로 명시적으로 지시해야 진행한다.
4. **투고 이력과 이해상충:** 커버레터의 "다른 곳에 게재·투고되지 않았다", "경쟁 이익 없음"이 사실인지 확인한다.
5. **(선택) 추천 심사자:**
   - 원고가 인용한 연구의 저자 중에서 고른 후보는 다음과 같다. 친분이나 이해상충이 있으면 뺀다.
     - Gilles Esposito-Farèse (damour1992tensor, 스칼라-텐서 시험)
     - Jan Steinhoff (chakrabarti2013response, steinhoff2016dynamical, khalil2022scalarization, 동적 응답 EFT)
     - Mohammed Khalil (khalil2022scalarization, 비단열 동적 스칼라화)
     - Norbert Wex (voisin2020sep, J0337 SEP 시험)
     - Anne M. Archibald (archibald2018universality, J0337 자유낙하 시험)
   - 넣지 않아도 된다. 제외할 심사자는 없는 것으로 적었다.

## 4. APS 제출 절차

1. https://authors.aps.org/Submissions/ 에 로그인하고 Physical Review D, Regular Article을 고른다.
2. 원고 PDF를 올린다. 소스 zip은 선택이다.
3. 메타데이터를 입력한다.
   - 제목
   - 초록: `abstract-plain.txt`
   - 저자, 소속, ORCID
   - 주제어(PhySH) 후보: Tests of gravity, Equivalence principle, Scalar-tensor theories (Alternative gravity theories), Effective field theory, Pulsars, White dwarfs
4. 커버레터를 붙인다. AI 사용 여부를 묻는 항목이 있으면 원고 7절과 같게 답한다.
5. 확인 화면을 검토한 뒤 제출한다.

## 5. 다시 빌드할 때

저장소 루트에서 실행한다.

```bash
PANDOC=<pandoc 경로> python paper/build_manuscript.py
C:/Users/lpaiu/AppData/Local/Programs/tectonic/tectonic.exe --keep-logs --keep-intermediates --outdir paper/build paper/main.tex
python outputs/direct-eos-gr33/native-submission-prep/scripts/phase289-package.py
```

**주의:** `paper/package_revision.py`를 그대로 실행하지 않는다. 이 스크립트는 `paper/revision-manifest.json`을 2026-09-09의 고정 목록으로 다시 써서, 그 뒤의 바인딩 수천 건을 지운다. 패키징은 위의 `phase289-package.py`로 한다.

## 6. 알아 둘 점

- PDF는 Tectonic 0.17.0(XeTeX)으로 만들었다. 같은 소스는 pdfLaTeX(latexmk)로도 컴파일된다.
- 모든 주장에 네 가지 표지(Proven 등)를 붙이는 방식은 일반적이지 않다. 편집자가 의견을 낼 수 있지만, 원고 서론이 이 방식을 정의하고 있다.
- §4.6의 근거 노트(REQUEST244–288)는 한국어이며, 원고의 데이터 가용성 절이 이를 밝힌다.
