# GRG submission and decision - updated 30 September 2026

## 2026-11-11: 연구선 종결

사용자 결정으로 이 연구선을 종결했다. 미션 질문은 미결인 채로 닫았다. 판정과 이유는 [종결 기록](research-line-closure-2026-11-11.md)에 있다. Physica Scripta 패키지는 미제출 상태로 동결하며, 제출 여부는 저자의 별도 결정이다. 아래의 개정·거절·접수 기록은 역사적 증거로 보존한다.

## 후속 개정: Physica Scripta 준비 완료, 미제출

거절 원인 분석 뒤 사용자의 “진행해봐”에 따라 방법론 중심의 단일 논문으로 개정했다. 본문 12쪽, 보충자료 35쪽, 커버레터 1쪽과 소스·재현 코드 ZIP을 준비하고 독립 실행·컴파일·전 페이지 조판 검수를 마쳤다. 제목은 *Finite-frequency identifiability of relaxation responses in free-fall timing*이다. [개정 및 검증 기록](physica-scripta-revision-2026-09-30.md)과 [저자 제출용 패키지](../output/submission/physica-scripta/README.md)를 참조한다.

분류: Proven. 비교 모형의 보간 경계, 실제 측정 투영 뒤 식별성, 조건수에 따른 정보량을 하나의 판정 절차와 정확 대조 예제로 연결했다. 분류: Imported from prior work. 저장 J0337 결과는 조건부 사례로 유지한다. 분류: Conjectural. 최종 전하의 물리·관측 폐쇄와 미분 오차·잡음 문제는 이 개정으로 해결되지 않았다. 장기 계산은 재개하지 않았다.

Physica Scripta에는 아직 제출하지 않았다. IOP의 제출 저자 본인 제출 정책에 맞춰 최종 검토·정책 선언·제출은 저자가 수행한다. 아래 거절 분석과 이전 접수 기록, GRG 동결 패키지는 역사적 증거로 보존한다.

## 현재 상태: 편집 단계 거절

2026-09-30에 교신 계정의 결정 메일과 [SNAPP 처리 이력](https://submission.springernature.com/submission-details/adae79fb-fa8e-4e35-87c6-be5672bbcc5a)을 직접 확인했다. 현재 상태는 **Your submission has been rejected**다.

- 9월 28일: 접수 및 기술 검사 시작.
- 9월 29일: 기술 검사 통과, 편집 배정, 편집자 배정, 거절. 표시된 처리 이력에는 외부 심사 단계가 없다.
- 결정 메일: 9월 29일 22:28 KST, 제목 `Decision on your submission to General Relativity and Gravitation`, 서명 Marco de Cesare, Editor.
- 명시된 사유: **scientific novelty, significance, and rigour**에 관한 저널의 출판 기준 미충족. 상세 심사 의견, 특정 수식의 오류 지적, 수정 후 재심사 초청은 없다.
- 이력 화면: [grg-decision-2026-09-29.png](../output/submission/grg-decision-2026-09-29.png). 기존 접수 기록과 동결 제출본은 그대로 보존한다.

이 기록은 편집 단계 거절(desk rejection)에 해당한다. 기술 검사 실패나 동료 심사에서 특정 정리가 반증된 사건으로 해석하지 않는다. 저널 범위 불일치만이 사유라고 단정할 근거도 없다. 아래 세부 원인은 편집자가 밝힌 사실이 아니라 제출 원고에 근거한 우리의 추정이다. 독립 연구자 소속이나 AI 사용이 거절 원인이었다는 증거는 없다.

## 원인 분석

**Imported from prior work.** [제출 원고](../output/submission/grg/manuscript.md)의 1절은 보간 경계가 고전적인 실수 다항식 보간, 단일 극 판정이 양의 모멘트 행렬의 특수 경우, 가까운 극의 구분 한계가 기존 역문제의 악조건성과 연결된다고 명시한다. 새 기여로 내세운 것은 nuisance 투영을 거친 자유낙하 timing에 대한 정식화, 물리적으로 제한한 비교 모형, 구동 위상 검사와 저장 분석에의 적용이다. **Conjectural.** 신규성의 부담은 이 결합이 기존 방법으로 바로 얻기 어려웠던 어떤 판정이나 절차를 가능하게 하는지에 있다. 현재 원고는 익숙한 정리들의 모음과 자체 분석 감사로 읽힐 위험이 크다. 정리의 정확성과 저널이 요구하는 새 기여의 크기는 별개다.

**Imported from prior work.** 원고 4절의 구동 위상 정정 대상은 우리 저장소의 미발표 분석이며, 기존에 발표된 J0337 결과가 아니다. 5.5절은 여섯 지정 지연에서의 조건부 모형 부적합을 보고하고, 6절은 이완 검출도 물리적으로 보정된 SEP 제약도 얻지 못했다고 명시한다. **Conjectural.** 중력 전문 저널의 독자에게 새로 확정된 물리적 결론이 무엇인지가 약하다. 자체 오류 발견은 중요하지만 그 자체를 분야의 기존 결론을 바꾼 성과로 제시할 수는 없다.

**Imported from prior work.** 원고 5.6절에서 유한차분 간격을 절반으로 바꿔 재계산한 28개 timing 미분은 nuisance rank 90을 유지하지만, 두 부분공간의 최대 principal sine은 0.999913이다. 물리 구동 표준오차의 비도 0.7858--0.8722로 변하며 엄밀한 미분 오차 예산은 없다. 5.5절의 2일 지연 부적합 통계량은 약 12.84로 임계값 12.8241766을 근소하게 넘는데, 개별 보정 분위수의 Monte Carlo 표준오차는 약 0.13이다. 양의 결합 부호를 강제한 결과는 저장된 조건부 분석에서 더 강하지만 미분 오차 문제까지 해결한 것은 아니다. **Conjectural.** 이는 데이터 응용의 엄밀성에 대한 실질적인 약점이다. 모의실험에서 지정한 공분산 계열의 포함률이 맞거나 저장 JSON을 재생할 수 있다는 사실만으로 실제 천체 자료의 추론까지 검증되지는 않는다.

**Imported from prior work.** 6절과 보충자료의 백색왜성 계산은 가정한 모형의 과도 신호와 크기 비교를 제공한다. 열 지연 수치는 정해진 열 시간 절단까지의 부분합이며 깊은 층의 기여와 대응하는 timing 응답은 미계산이다. **Conjectural.** 수학적 경계, 관측 추론 감사, 조건부 항성 계산을 모두 강조하면 하나의 완결된 기여가 흐려진다. 한 논문으로 유지하되 중심 주장을 더 좁히는 편이 낫다.

이전 투고 준비에서 확인한 빌드, 파일 해시, 상징 계산, 재생 검사는 제출물의 일관성에 대한 근거였다. 그것을 편집자가 인정할 신규성·중요성이나 전체 관측 추론의 완결성과 동일시해서는 안 된다. 이번 추천은 이 차이를 반영한다.

## 다음 투고처

2026-09-30 공식 범위·출판 비용 안내를 확인했다. 순위는 수정 후 적합성에 대한 판단이며 게재 확률 추정이 아니다.

| 순위 | 저널 | 적합성 및 조건 | 비용 |
| --- | --- | --- | --- |
| 1 | **Physica Scripta**, Research Paper | 수리·통계 물리, 천체물리, 계산·측정 방법을 포함하며 타당한 방법론의 음성 결과도 명시적으로 환영한다. 식별성 방법과 조건부 실패 사례를 중심으로 재구성할 때 가장 맞는다. 원천적 새 기여와 엄밀성은 여전히 요구된다. [공식 범위·논문 유형](https://publishingsupport.iopscience.iop.org/journals/physica-scripta/about-physica-scripta/) | 같은 공식 안내상 구독형 출판 무료, 투고료·기타 출판료 없음. OA는 선택 사항. |
| 2 | **The European Physical Journal Plus**, Original Research | 이론·계산 기법을 포함한 방법론 논문을 명시적으로 받는다. 물리적 활용과 재현 가능한 비교 사례를 선명하게 해야 한다. 폭넓은 범위가 편집 거절 위험을 없애지는 않는다. [공식 범위](https://link.springer.com/journal/13360/aims-and-scope) | 구독형은 APC 없음. [공식 비용 안내](https://link.springer.com/journal/13360/how-to-publish-with-us) |
| 3, 조건부 | **Astrophysics and Space Science** | 천체물리 적용이 분명해야 하며 그 연결이 없는 수리물리·GR 논문은 제외한다. 따라서 J0337 적용의 신뢰성을 강화하는 방향일 때만 적합하다. 현재의 방법론 중심 수정에는 앞의 두 곳이 낫다. [공식 범위](https://link.springer.com/journal/10509/aims-and-scope) | 구독형은 APC 없음. [공식 비용 안내](https://link.springer.com/journal/10509/how-to-publish-with-us) |

Physica Scripta의 짧은 Technical Note도 제도상 가능하지만 거절 분석 당시의 GRG 본문 14쪽과 보충자료 31쪽을 바탕으로 Research Paper를 기준으로 계획했다. 후속 개정판 분량은 위 최신 상태에 기록한다. Foundations of Physics는 개념적·기초적 의미를 요구하므로 이번 방법론 원고의 우선 대안으로 삼지 않는다. [공식 범위](https://link.springer.com/journal/10701/aims-and-scope)

IOP의 [현행 저자 윤리 정책](https://publishingsupport.iopscience.iop.org/ethical-policy-journals/)은 AI의 연구·코드·원고 지원 사용을 허용하되 저자의 검증과 공개를 요구한다. 기존 실질적 AI 사용 공개를 숨기거나 언어 교정만으로 축소하지 않는다. 연구 과정의 사용은 Methods에도 설명하고 저자의 검증 범위를 정확히 적어야 한다. 같은 정책의 Identity fraud and impersonation 항목은 투고 저자 본인의 제출을 요구하므로, Physica Scripta로 진행할 경우 자료 준비와 별개로 실제 제출 절차는 저자가 직접 수행하는 방식으로 맞춘다.

## 재투고 전에 필요한 최소 수정

1. **새 기여를 대조표로 한정한다.** 고전 보간·모멘트 결과와 이번 연구의 추가 판정 절차를 분리하고, abstract와 introduction에서 무엇을 새롭게 할 수 있는지 하나로 연결한다. 제안 제목은 `Finite-frequency identifiability of relaxation responses in free-fall timing`이다. 제목 변경만으로 신규성이 생기는 것은 아니다.
2. **방법의 완결된 사례를 앞세운다.** 이미 저장된 결과를 우선 재사용하여 알려진 구동, 비교 모형 차수, nuisance 투영, 정밀도에 따른 성공·실패 대조를 독자가 재현하게 한다. 합성 예제의 참값 검증과 실제 관측자료의 물리적 검증은 구분한다. 새 계산의 필요성은 이 최소 사례에 빠진 증거가 있는지를 본 뒤 정한다.
3. **J0337은 조건부 사례로 엄격히 제한한다.** 미분에 따른 부분공간 변화, 2일 지연의 근소한 거절, 미지정 잡음·누락 구동을 결과 바로 옆에서 함께 읽게 한다. 이미 원고에 한계가 명시되어 있으므로 숨겨진 오류를 처음 발견했다는 식으로 쓰지 않는다. 검증되지 않은 empirical exclusion을 내세우지 않으며, 실제 관측 제약을 주요 성과로 유지하려면 미분·likelihood 검증이 먼저 필요하다.
4. **중심 주장과 보충자료의 역할을 정리한다.** 상세 실패 이력과 항성 계산은 보존하되 방법론 결과의 증명처럼 사용하지 않는다. SHA와 내부 AI 심사 수가 정확성이나 출판 가치의 대리 지표가 되지 않게 한다. 표·그림·커버레터는 동일한 좁은 기여를 설명해야 한다.

이 수정은 장기 EOS/GR 또는 pulsar 계산을 자동 재개하자는 계획이 아니다. 거절 사유를 더 구체적으로 확인할 근거가 없으므로 현재로서는 항의나 원고 그대로의 재투고보다 위 수정의 가치가 높다. 거절 분석을 마친 시점에는 확인과 투고처 조사까지만 수행했고, 새 투고·전송·원고 변경·장기 계산은 실행하지 않았다. 이후 사용자가 승인한 원고 개정 결과는 위 최신 상태에 기록한다. GRG 처리는 거절로 종료되었고 새 투고는 순차적으로 진행할 수 있다.

## Historical submission receipt - 28 September 2026

The integrated manuscript was submitted to **General Relativity and Gravitation** as a **Research** article. SNAPP confirmed **Submission received** and the tracking page showed **Technical check - in progress**. This is a submission receipt, not editorial acceptance or completed peer review.

- Submission ID: `adae79fb-fa8e-4e35-87c6-be5672bbcc5a`, version `v.1.0`.
- Tracking: https://submission.springernature.com/submission-details/adae79fb-fa8e-4e35-87c6-be5672bbcc5a
- Package: `output/submission/grg/`; frozen public tag `grg-submission-2026-09-28`, commit `2720a4006fe421da63c7ba897b98b621750e2971`.
- Receipt metadata and file hashes: `output/submission/grg-receipt.json`; screenshot: `output/submission/grg-receipt.png`.
- Server-generated reviewer PDF: `C:/Users/lpaiu/Downloads/manuscript-source.pdf` (14 pages).

The author confirmed the Seoul/Republic of Korea address, absence of external funding, and the required author-policy acknowledgement. The corresponding email is the previously confirmed `lpaiu.cs@gmail.com`. The portal records no institutional affiliation, consistent with independent researcher status in the manuscript. Research Square preprint posting was declined. No APC payment or publication licence was executed; publication-route/licence choices occur after acceptance.

The title, 230-word abstract, author and correspondence details, editable source ZIP, Online Resource 1, cover letter, author contribution, funding, competing interests and data availability were checked in the final submission review. The journal converted the ZIP successfully, and its actual reviewer PDF was downloaded and visually checked. All frozen package hashes passed; the independent ZIP compile preserved every page's text and rendered pixels. The existing 14 exact symbolic identities and compact public conditional-inference replay passed. No long numerical calculation was started.

**Imported from prior work.** This journal adaptation preserves Main Sections 2-6 and the supplement's scientific body. Finite-frequency identifiability results and the conditional J0337 application remain the contribution. The derivative uncertainty and uncomputed thermal tail remain stated limitations; submission does not certify a new empirical SEP bound, a detected relaxation, or complete physical/observational closure.

The earlier `prd-submission-2026-09-28` tag is a preparation archive. No PRD submission was made. The pending technical-check status in this historical receipt was superseded by the editorial rejection recorded above; the receipt and frozen submission package remain unchanged.
