# 투고 전 지적 해결과 공개 저장소 정리 — 2026-09-28

사용자는 REQUEST294의 지적을 모두 해결하고 공개 저장소의 stale branch도 정리하도록 지시했다. 시작 체크포인트는 `019af22c493fe0dbb490647010009b8d1948d7cb`, 검토된 과학적 기준선은 `6a0279b1276d60d9440f12e211a1f7e38daf7d05`다. 이번 작업은 제출 준비이며 저널에 실제 제출하는 작업은 아니다.

분류: Imported from prior work. 저장된 여섯 계수 추론을 다시 읽어 공개 재생 자료로 내보냈다. 옴니버스는 16.35252153, 임계값은 12.82417662이며 평가한 여섯 물리 지연 단면은 모두 비어 있다. 비음수 beta 제약에서 최솟값은 14.58758030 이상이다. 미분 공간·공분산 모형에 조건부인 기존 결론과 백색왜성의 조건부 과도 신호 해석은 유지한다. 새로운 완화 검출, 경험적 SEP 상한 또는 전체 연구 목표 완료를 주장하지 않는다.

분류: Proven. 층별 Debye 모형의 미계산 꼬리를 `T_tail = sum_tail |s_j|/|S_struct|`라 두면 그 직교 성분은 `T_tail/(omega*tau_cut)` 이하이다. 분류: Imported from prior work. 기존 열 완화 계산은 `tau_cut=1e13 s`까지이며, 이에 대한 인자는 내측 2.24058e-9, 외측 4.50008e-7이다. 분류: Conjectural. `T_tail` 자체는 계산하지 않았다. 전체 지연 1% 이하라는 기준은 추가로 `T_tail <= 2.22e4`를 가정하면 충족하지만, 이 가정은 검증 결과가 아니다. 원고의 4.0e-9 / 3.3e-7은 계산한 층의 기여로만 한정했다.

원 실패 기록은 보존했다. `native-thermal-relaxation`의 당시 `passed`와 전체 상한처럼 읽힐 수 있는 verdict는 계산한 껍질에 대해서만 해석한다. 전체 꼬리를 수치적으로 인증했다는 해석은 이번 노트와 수정된 SM §4.6이 대체한다. 미계산 층을 추가 계산하거나 수락 기준을 완화하지 않았다.

수정 사항:

- 본문과 보충자료에 열 시간 절단·미계산 세기·잔여항 조건을 명시했다.
- 보충자료를 본문 참고문헌 항목으로 인용하고, 보충자료에만 있던 문헌 네 편을 본문 참고문헌에도 포함했다.
- `phase291-package.py`의 Claude 작업폴더 고정을 제거했다. 현재 스크립트의 저장소에서 읽고 쓰며, Markdown·TeX·문헌·그림보다 PDF가 최신인지 확인한다. `submission-manifest.json`을 만들고 기존 master 입력 바인딩을 삭제하지 않고 갱신한다.
- 저자 정보는 과거 causal-spacetime 원고(`docs/paper/paper_a/latex/sections/front.tex:8`)를 확인해 Juneyoung Kim / Independent researcher로 재사용했다. 과거 entropy-arrow의 하이픈 이메일과 달라 사용자에게 확인했고, 이번 사용자 답변인 `lpaiu.cs@gmail.com`을 적용했다. 확인된 ORCID가 없으므로 자리표시자를 제거했다.
- `public-inference.json`과 `replay_public_inference.py`는 공개 JSON과 NumPy만으로 여섯 계수의 옴니버스·물리 단면을 다시 계산한다. `--export`는 사유 배열이 필요한 별도 동작이다. 추론 결과를 원 저장 기록과 비교하는 검사를 내장했고, 전체 엔진·도함수 정확도·보정 재실험의 대체물이 아님을 밝혔다.
- `check_submission_package.py`는 소스 ZIP이 현재 원고와 일치하는지, 모든 보충자료 인용이 본문 `.bbl`에 있는지, 저자 정보와 제출 해시가 맞는지 검사한다. 공개 clone에서도 동작한다.
- README·데이터 가용성 절·CITATION.cff·체크리스트를 현재 패키지와 맞췄다. 과거 phase manifest는 각 과거 커밋의 증거로 보존했다.

공개 배포 대상은 `prd-submission-2026-09-28` 태그다. 해당 태그를 한 번 생성하고 이동하지 않는다. 본문·보충자료·코드·소규모 재생 자료·제출 해시를 같은 스냅숏에 담는다. 제출 파일의 실제 해시는 `paper/submission-manifest.json`을 따른다.

공개 브랜치 조사 시 열린 PR은 없었고, main 외의 세 브랜치는 모두 2026년 4월 이후 갱신되지 않았다. 로컬 작업폴더는 정리 범위에 포함하지 않았다.

| 공개 브랜치 | 기존 tip | 보존 방법 |
|---|---|---|
| clock-timing-dictionary | ccc8f1c63b84a18e89192dbb89a51b444b1f65dc | archive/2026-09-28/clock-timing-dictionary 태그로 고유 이력 보존 |
| collapse-theorem-dynamic-visibility | a6c28705e0ab1bd51fccdec2e65fde6764da3f2a | archive/2026-09-28/collapse-theorem-dynamic-visibility 태그로 고유 이력 보존 |
| lpaiu/dynamic-chi-observable | f3c0a7195f39ac610c0bca1dd73ef4242ae65ba8 | 기존 공개 main과 로컬 main의 조상이므로 이력이 이미 보존됨 |

첫 두 브랜치는 main의 조상이 아니므로 단순히 merged라고 간주하지 않았다. 공개 아카이브 태그의 대상 SHA를 확인하고 원격 브랜치 tip에 lease를 걸어 정리한다. 기본 main과 로컬 작업폴더를 삭제하지 않는다.

검증 결과:

| 검사 | 결과와 범위 |
|---|---|
| 정확한 항등식·비교 모형 | `state_identifiability.py` 14개 항등식과 `comparator_audit.algebra_checks()` 3개 검사 통과. 기존 산출물을 덮어쓰지 않았다 |
| 전체 원고 검증기 | `verify_unified_paper.py` 통과. 해석 경계, 역사적 표, 정정 구동, 저장 공동영역, 실패한 승격 gate, 상태 항등식과 master 해시를 확인했다. 약 17.6GB의 저장 입력을 읽었으며 새 GR·타이밍 계산은 없다 |
| 최종 문구 변경 후 검증 | 전체 검사 후 바뀐 SM 요약·PDF·제출 문서와 패키지의 현재 해시는 29개 제출 바인딩 검사로 다시 확인한다. 나머지 기존 입력은 그대로 보존한다 |
| 빌드·시각 검사 | Tectonic 0.17.0으로 본문 14쪽, 보충자료 31쪽, 커버레터 1쪽을 컴파일하고 전체 페이지를 렌더링해 확인했다. 마지막 SM 요약 변경 뒤 14–16쪽도 다시 확인했다. 잘림·겹침·누락 글리프 없음 |
| 컴파일러 제한 | 앱 내장 LaTeX 미리보기는 환경 디렉터리 오류로 실패했다. 기존 Tectonic의 성공한 출력으로 배포했다. 소스는 보존했다 |
| 독립 소스 ZIP | 임시 폴더에서 컴파일 성공. 14쪽 전체의 추출 텍스트가 배포 PDF와 같고, 다시 생성한 `.bbl`도 ZIP 내용과 일치한다 |
| 교차 폴더 패키징 | 저장소 루트가 아닌 작업 디렉터리에서도 스크립트가 자신의 체크아웃을 묶고, PDF·ZIP·현재 소스의 일치 검사를 통과했다 |
| 공개 파일만의 재생 | 대형 배열이 전혀 없는 40개 파일의 복사본에서 `check_submission_package.py`와 `replay_public_inference.py`가 통과했다. 이후 공개 스냅숏 커밋에서도 같은 두 검사를 수행한다 |
| 해시 보존 | 기존 master의 모든 입력 키를 유지했고 현재 master는 12,323개다. 별도 제출 manifest의 29개 파일과 master의 대응 값이 일치한다. 역사적 phase manifest는 재작성하지 않았다 |

공개 반영은 이 검증을 통과한 로컬 커밋에서 스냅숏을 만들고, 원격 main 갱신·제출 태그·보존 태그·세 브랜치 삭제를 하나의 atomic push로 적용한다. 각 원격 tip에는 조사한 SHA의 lease를 적용해 다른 작업의 갱신을 덮어쓰지 않는다. 실제 공개 상태와 대상 커밋은 위 태그와 원격 refs에서 확인한다.
