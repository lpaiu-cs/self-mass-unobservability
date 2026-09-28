# Physical Review D 제출 체크리스트 (2026-09-28)

제출 파일은 투고 전 검토의 지적을 반영했다. 실제 저널 제출은 아직 하지 않았다. 업로드할 파일과 공개 스냅숏의 식별자는 아래와 같다.

**진행 상황 (2026-09-28)**
- AI 사용 절: 저자가 사실임을 확인했다. 모델 버전(Claude Opus 5.5, GPT-6-Astra)을 넣었다.
- 투고 이력과 이해상충: 문제없음을 저자가 확인했다.
- 최종 검토: 끝났다(`notes/REQUEST290_FINAL_REVIEW_KO.md`).
- 초점 축소: 저자 지시로 본문을 식별성 경계와 J0337 조건부 적용으로 좁혔다(새 제목). 정적 연산자 목록, 전체 감사와 출처, 백색왜성 계산은 보충 자료(30쪽, 이전 전체 원고)로 옮겼다(`notes/REQUEST291_FOCUSED_MANUSCRIPT_KO.md`).
- 최종 독립 검토: Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra가 모두 "주요 수정"을 권했다. Opus가 물리 구동의 근점 규약 오류를 찾았다. 코드는 η=e sin ϖ를 쓰는데 분석은 반대로 읽었다. 이를 바로잡자 J0337 물리 구동 결론이 바뀌었다. 평가한 여섯 지연에서 여섯 계수 초과를 물리 구동이 재현하지 못한다(2일은 근소하고, β≥0이면 확실하다). 세 심사자의 확인 심사는 모두 "경미 수정 후 수락"이었고, 그 지적도 반영했다(`notes/REQUEST293_CONFIRMATION_REVIEW_KO.md`). 수정 뒤 본문은 13쪽이다(`notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md`).
- 후속 독립 검토: 열 완화 절단·미계산 꼬리, 보충자료 참고문헌, 패키징 경로, 공개본 불일치를 수정했다. 확인 내용과 공개 배포 기록은 `notes/REQUEST295_SUBMISSION_CORRECTIONS_KO.md`에 있다.
- 저자 정보: Juneyoung Kim / Independent researcher / lpaiu.cs@gmail.com. 이전 논문 표기를 재사용하고 이메일은 2026-09-28에 저자가 직접 확인했다. 확인된 ORCID가 없어 빈칸을 제거했다.
- 공개 배포 식별자: `prd-submission-2026-09-28`. 공개 재생 코드와 현재 본문·보충자료를 포함한다. 이 태그는 준비본의 이름이며 저널 제출 완료를 뜻하지 않는다.
- **남은 일:** APS 제출 양식 입력과 최종 제출. 연구의 조건부 한계는 본문에 유지했다.

## 1. 투고처 결정

**Physical Review D (Regular Article)** 를 1순위로 권한다.
- **범위:** PRD는 중력, 일반상대론과 대안 이론, 밀집천체, 천체물리를 다룬다. 이 원고는 자유낙하·SEP 시험, 스칼라-텐서 이론, 확장 천체의 EFT, 펄서 타이밍 추론을 다룬다.
- **AI 정책:** 이 연구는 AI 에이전트를 계산과 원고 작성에 실질적으로 썼다.
  - APS의 2026년 6월 정책은 과학적 추론, 주장 작성 같은 실질적 사용도 허용한다. 대신 도구 이름·버전, 도움의 내용, 지시·검증 방법을 논문 안에 공개하도록 요구한다.
  - IOP(CQG)의 정책은 심사 답변 작성에 생성형 AI를 쓰는 것을 제한한다. 이 연구의 작업 방식에는 APS가 더 잘 맞는다.
- **형식과 비용:**
  - 첫 제출은 PDF만으로 심사가 시작된다. REVTeX은 권장일 뿐 필수가 아니고, 소스는 수락 뒤에 요구된다.
  - 연구 논문에는 길이 제한이 없다. 보충 자료(Supplemental Material)는 별도 파일로 올린다.
  - 의무 게재료는 없다. 오픈 액세스(APC)와 인쇄 컬러는 선택 사항이다.
- **대안:**
  - Classical and Quantum Gravity(Research Paper): 게재료 없음, 어떤 TeX 형식이든 제출 가능. 다만 원고에 "더 넓은 중력 물리 맥락" 요약이 필요하고, AI 공개는 Acknowledgements에 넣어야 한다.
  - General Relativity and Gravitation

## 2. 제출 파일

| 파일 | 용도 |
|---|---|
| `output/pdf/free-fall-identifiability.pdf` | 원고 PDF. 첫 제출에 이것과 보충 자료만 올려도 된다 |
| `output/pdf/free-fall-identifiability-supplement.pdf` | 보충 자료 PDF. Supplemental Material로 올린다 |
| `output/submission/cover-letter-prd.pdf` (`.tex`) | 저자·교신 정보가 반영된 커버레터 |
| `output/submission/abstract-plain.txt` | 제출 양식에 붙일 평문 초록 |
| `output/submission/free-fall-identifiability-source.zip` | 저널용 소스(main.tex, main.bbl, references.bib, 그림 1개). 선택 업로드이며, 이 zip만으로 컴파일됨을 확인했다 |

## 3. 제출 메타데이터와 공개본

1. **저자:** Juneyoung Kim; **소속:** Independent researcher; **교신 이메일:** lpaiu.cs@gmail.com.
   - ORCID는 이번 파일에 기재하지 않았다. 제출 시스템에서 연결하려면 본인의 실제 iD를 사용한다. ORCID가 항상 필수라는 이전 체크리스트의 단정은 삭제했다.
2. **공개 스냅숏:** `https://github.com/lpaiu-cs/self-mass-unobservability/tree/prd-submission-2026-09-28`.
   - `paper/submission-manifest.json`으로 파일 무결성을 확인할 수 있다. `python verification/check_submission_package.py`와 `python verification/replay_public_inference.py`가 공개본에서 동작한다.
   - 추후 수정본은 새 태그로 구분한다. 이미 공개한 제출 태그를 다른 커밋으로 옮기지 않는다.
3. **(선택) 추천 심사자:**
   - 원고가 인용한 연구의 저자 중에서 고른 후보는 다음과 같다. 친분이나 이해상충이 있으면 뺀다.
     - Gilles Esposito-Farèse (damour1992tensor, 스칼라-텐서 시험; 보충 자료에서 인용)
     - Jan Steinhoff (chakrabarti2013response, steinhoff2016dynamical, khalil2022scalarization, 동적 응답 EFT)
     - Mohammed Khalil (khalil2022scalarization, 비단열 동적 스칼라화)
     - Norbert Wex (voisin2020sep, J0337 SEP 시험)
     - Anne M. Archibald (archibald2018universality, J0337 자유낙하 시험)
   - 넣지 않아도 된다. 제외할 심사자는 없는 것으로 적었다.

## 4. APS 제출 절차

1. https://authors.aps.org/Submissions/ 에 로그인하고 Physical Review D, Regular Article을 고른다.
2. 원고 PDF를 올리고, 보충 자료 PDF를 Supplemental Material로 올린다. 소스 zip은 선택이다.
3. 메타데이터를 입력한다.
   - 제목: Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715
   - 초록: `abstract-plain.txt`
   - 저자, 소속, 교신 이메일; ORCID는 본인 iD가 있는 경우 연결
   - 주제어(PhySH) 후보: Tests of gravity, Equivalence principle, Scalar-tensor theories (Alternative gravity theories), Effective field theory, Pulsars
4. 커버레터를 붙인다. AI 사용 여부를 묻는 항목이 있으면 원고의 "Use of AI tools" 절과 같게 답한다.
5. 확인 화면을 검토한 뒤 제출한다.

## 5. 다시 빌드할 때

저장소 루트에서 실행한다.

```bash
PANDOC=<pandoc 경로> python paper/build_manuscript.py
PANDOC=<pandoc 경로> python paper/build_manuscript.py supplement.md supplement.tex
C:/Users/lpaiu/AppData/Local/Programs/tectonic/tectonic.exe --keep-logs --keep-intermediates --outdir paper/build paper/main.tex
C:/Users/lpaiu/AppData/Local/Programs/tectonic/tectonic.exe --keep-logs --keep-intermediates --outdir paper/build paper/supplement.tex
C:/Users/lpaiu/AppData/Local/Programs/tectonic/tectonic.exe --outdir output/submission output/submission/cover-letter-prd.tex
python outputs/direct-eos-gr33/native-focused-manuscript/scripts/phase291-package.py
python verification/check_submission_package.py
python verification/replay_public_inference.py
```

**주의:** `paper/package_revision.py`를 그대로 실행하지 않는다. 이 스크립트는 `paper/revision-manifest.json`을 2026-09-09의 고정 목록으로 다시 써서, 그 뒤의 바인딩 수천 건을 지운다. 수정된 `phase291-package.py`는 자신의 체크아웃을 사용하고, 이전 입력 목록을 유지하면서 제출 파일의 해시를 갱신한다. 과거 phase publisher는 해당 과거 버전의 검증용이며 현재 제출 파일에 다시 적용하지 않는다.

## 6. 알아 둘 점

- PDF는 Tectonic 0.17.0(XeTeX)으로 만들었다. 같은 소스는 pdfLaTeX(latexmk)로도 컴파일된다.
- 모든 주장에 네 가지 표지(Proven 등)를 붙이는 방식은 일반적이지 않다. 편집자가 의견을 낼 수 있지만, 원고 서론이 이 방식을 정의하고 있다.
- 보충 자료는 이전 전체 원고에 근점 규약 정정, 주장 범위·열 완화 꼬리 조건, 공개 재생 안내를 반영했다. 절·식·표 번호가 본문과 다르며, 본문은 "SM Section 4.6"처럼 가리킨다. 보충자료 전용 인용도 본문 참고문헌에 포함했다.
- 초록과 모든 문단의 주장 표지(Proven 등)는 저장소 규칙이다. Opus와 Fable이 편집자에게 낯설 수 있다고 지적했다. 초록에서 표지를 뺄지는 저자가 정한다.
- 보충 자료 §4.6의 근거 노트(REQUEST244–288)는 한국어이며, 보충 자료의 데이터 가용성 절이 이를 밝힌다.
