# 남은 연구 레버 6--10 및 상태 식별성 결과

기준 revision: 3bc2fce. 승인된 항목을 모두 실행 대상으로 삼아 이론, 저장 배열 분석, 별도 WSL 엔진 검증을 수행했다. 각 연구 항목의 진행과 모든 물리적 결론의 확정은 구분한다.

| 항목 | 확보한 결과 | 완료하지 못한 물리적 주장 |
| --- | --- | --- |
| 6. 수정된 구동 | Imported from prior work: 올바른 위상, 서로 다른 진폭, 명시한 정규화로 조건부 구간을 재계산했다. | Conjectural: 생략 고조파/추가 힘까지 통제한 실제 중력 제약. |
| 7. 동시 추론 | Proven: 여섯 계수 신뢰영역을 통한 위상/완화시간 추론과 알려진 잡음 상계 조건의 보수적 보장. Imported from prior work: 독립 모의 포함률 95.17--95.43%. | Conjectural: 모든 천체물리 잡음, 새 초기상태와 미분 오차까지 포함하는 완전한 likelihood. |
| 8. 초기상태 | Imported from prior work: 2/52/500일의 실제 transient timing 응답을 계산하고 진폭 수렴 검사 통과. | Conjectural: 모든 완화시간과 큰 초기 진폭의 비선형 응답 및 통계 보장. |
| 9. 물리적 매칭 | Proven: 평형 정보만으로 감쇠와 완화율 하한을 정할 수 없다는 불충분성 정리. 문헌의 구체적 매칭 경로를 확인했다. | Conjectural: 실제 J0337의 EOS/중력결합/모드/백색왜성 전하 수치 매칭. 필요한 입력은 임의로 선택하지 않았다. |
| 10. 수치/주기 계수 | Imported from prior work: 28개 미분 재계산, 7개 추가 step, 201개 공백 조합, 국소 비선형 검사. | Conjectural: 엄밀한 미분 오차 상계와 물리적으로 제약된 완전한 pulse 재연결. |
| 추가 이론 | Proven: 양의 응답 클래스의 단일 관측 완화시간 판별식과 유한 정밀도 분리 불가능성. | Conjectural: 실제 별의 내부 자유도 개수 판정. |

## 핵심 판정

Status: Proven. 두 정확한 복소 응답에서 얻는 양의 측도의 moment determinant가 0인 조건으로 단일 관측 완화시간을 판별할 수 있다. 숨은 상태의 개수는 결정하지 못하며, 서로 매우 가까운 두 양의 pole은 하나에 임의로 접근한다.

Status: Imported from prior work. 수정된 구동의 생략항은 공면 Kepler 모형에서 입력 RMS의 약 3.37%였다. 이를 timing 오차 3.37%라고 해석하지 않는다. 독립 검증을 거친 신뢰영역의 실제 데이터 omnibus 통계가 문턱을 넘지만, 이는 모든 carrier 계수가 0이라는 가설의 검사다. 물리적 순간 응답을 허용한 beta=0은 시험한 모든 lag section에서 포함되며, 새 relaxation 검출을 주장하지 않는다.

Status: Imported from prior work. 전용 초기상태 응답의 수렴은 통과했고, 시험한 경우 beta 표준오차 증가는 최대 약 0.523%였다. 미분 열의 최대 변화는 0.665%인데도 약한 nuisance 공간의 최대 principal sine은 0.999913이었다. 고정 행렬의 선형대수 정밀도와 실제 미분의 정확도는 다르다.

Status: Imported from prior work. 가장 약한 223일 공백의 선형 보상은 물리 범위 밖의 이심률을 요구했다. 작은 변위만의 비선형 검사는 통과했으나 전체 보상을 정당화하지 않는다. 다른 비선형 해의 존재 여부는 이 실패로 결정되지 않는다.

현재 연구 분류: theorem progress 및 loophole progress. 통합 방법론 원고를 확장했다. **실제 EOS 기반 제약과 모든 관측계통의 검증이 완료됐다고 판정하지 않는다.**

## 재현 기록

사전 계획: `notes/REQUEST12_COMPLETION_PLAN.md` (abb973e). 구동/통계 문턱 동결: 6d7f3f3. 동적 구현 및 baseline 재현: 150fa43. 단일 상태 정리: d8eb72d. 비선형 실패 후 국소 검증 등록: a08c36a. 실제 응답/미분 결과: 05d7af9.

상세 결과: `notes/REQUEST12_DRIVE_INFERENCE_RESULT.md`, `notes/REQUEST12_THEORY_RESULT.md`, `notes/REQUEST12_RUNTIME_RESULT.md`. 원래 REQUEST10 산출물은 보존하며 새 수치 결과와 전체 응답 벡터는 `outputs/research-completion/` 아래에 둔다. 외부 실행은 기존 WSL 환경의 별도 source/run 디렉터리에서 했다.

## 최종 검수

통합 원고 23페이지 전체를 렌더링해 검수했다. 소스 ZIP을 별도 폴더에서 재빌드한 결과 23페이지 모두 텍스트와 검사 해상도의 렌더링 픽셀이 일치했다. 통합 검증은 수학 경계, 기존 표, 새 구동, 독립 신뢰영역 검증, 실제 응답 벡터와 실패 판정, SHA-256 입력 160개를 확인하고 통과했다. REQUEST10 파일의 기준 revision 이후 diff는 없다. OSK의 기존 프로젝트·Paper B·물리 후보·허브 노드 및 scope를 갱신했다. 저널 제출이나 원격 공개는 수행하지 않았다.


## 단계292 — 최종 독립 검토와 근점 규약 정정

분류: Imported from prior work. Claude Opus 5.5, Claude Fable 5.1, GPT-6-Astra의 독립 심사가 모두 주요 수정을 권했다. Opus가 J0337 물리 구동의 η·κ 규약 오류를 찾았다. 코드는 η=e sin ϖ를 쓰는데 분석은 η=e cos ϖ로 읽었다. 이를 바로잡고 영향받는 계산을 다시 돌렸다(각 수 초). 이전 출력은 `outputs/research-completion/withdrawn-periastron-convention/`에 보존했다.

실패 기록: 실패한 단계는 `symbolic/physical_matching.py`와 `verification/physical_drive_completion.py`의 근점 각 계산이다. 빠진 최소 가정은 런타임 코드의 매개변수 규약 확인이다. 이전 결론 '모든 물리 지연 단면이 β=0을 포함한다'를 철회한다.

분류: Imported from prior work. 수정 뒤 여섯 단면이 모두 비어 있다. 기록된 반송파 초과(옴니버스 16.3525, 명목 p≈0.012)를 물리 구동의 순간+완화 응답이 재현하지 못한다. relaxation 검출은 아니며 원인은 분리하지 못했다. 백색왜성의 최종 전하 결론은 유지된다. 문헌 위치, 두 표본 조건, 한정어, 정의 등 나머지 지적도 반영했다(본문 13쪽). [근거 292](../notes/REQUEST292_FINAL_INDEPENDENT_REVIEW_KO.md).


## 단계293 — 확인 심사와 경미 수정

분류: Imported from prior work. 단계292 개정판을 세 심사자(Opus 5.5, Fable 5.1, GPT-6-Astra)가 확인 심사했다. 모두 경미 수정 후 수락을 권했고, 세 명 모두 근점 규약 정정을 독립적으로 확인했다. J0337 결론의 범위를 평가한 여섯 지연과 순간항+단일 완화 모형으로 한정했다. 2일 단면의 근소한 기각을 정량화했고(보정 순서통계량 12.27–12.82, 표준오차 약 0.13), 보관된 도함수 조건을 명시했다.

분류: Proven. 통계량은 볼록 이차식이고 β=0 직선이 모든 단면에 들어 있다. 제약 없는 β̂은 18일을 빼면 음수다. 따라서 등전하 실현의 β≥0 아래에서는 여섯 단면 모두 최솟값이 14.59 이상이다. 백색왜성의 최종 전하 결론은 유지된다. [근거 293](../notes/REQUEST293_CONFIRMATION_REVIEW_KO.md).
