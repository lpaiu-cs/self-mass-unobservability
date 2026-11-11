# Physica Scripta 방법론 개정 - 2026-09-30

상태: **수정·검수 완료, 저자 직접 제출용 패키지 준비. 아직 제출하지 않음.** 사용자의 “진행해봐”에 따라 GRG 거절 뒤 제안한 방법론 중심 개정을 수행했다. 원고·보충자료·커버레터·재현 코드가 함께 있는 단일 논문이며 A/B 독립 투고로 되돌리지 않는다.

## 과학적 변화와 보존 경계

**Proven.** 기존 보충자료의 투영·정보량 관계와 유한 주파수 보간 경계를 본문의 하나의 판정 절차로 연결했다. 고정 선형 모형에서 \(A=(I-P_J)LC\), \(I_{\beta,N}=\min_a\|A(w-V_Na)\|^2\)가 진폭 식별성 판정식이다. full column rank 아래에서는 \(\sigma_{\min}^2\delta_N^2\le I_{\beta,N}\le\sigma_{\max}^2\delta_N^2\)다. 표준 보간·최소제곱을 새로운 일반 수학으로 주장하지 않는다. 분류는 **theorem progress: 기존 경계의 관측 적용과 정확 대조를 명료화**한 것이다.

**Proven.** 새 수치는 오직 명시적으로 구성한 무차원 해석 예제에서 나온다. 주파수 (1,2), 시간상수 1, 차수 1 비교 모형에서 계수 공간 거리 제곱은 117/1000이다. 동일한 응답에 지도 I, 비교 모형으로의 정사영 P, \(P+\epsilon(I-P)\)를 적용하면 정보량은 각각 117/1000, 0, \(117\epsilon^2/1000\)이다. 이들은 새 천체 관측값이 아니며 장기 계산이나 Monte Carlo 실험도 아니다. 정확 연산 검사는 별도 스크립트로 제공한다.

**Imported from prior work.** J0337 관련 수치는 저장된 comparator 감사와 공개 여섯 계수 기록에서 재사용했다. comparator 정보량 최솟값, 여섯 물리 지연 단면, 구동 위상 정정, 미분 공간 불안정성, 지정 공분산 안의 검증 범위를 바꾸지 않았다. 새 구동·미분·잡음 모의실험·항성 진화는 실행하지 않았다. 원고 §5.3의 비음수 beta 최솟값은 실제 값 14.587580...의 반올림임을 드러내기 위해 `approximately 14.59`로 표기한다.

**Conjectural.** 최종 전하의 과학적 결론과 남은 물리·관측 오차는 이번 개정으로 달라지거나 닫히지 않는다. 열 지연 꼬리, 실제 타이밍 연결, 미분 오차 예산과 비선형/잡음 추론은 여전히 미해결이다. 공개 재생이나 새 정확 예제를 그 대체 검증으로 쓰지 않는다. 새로 확보한 관측 loophole 또는 empirical SEP 제약은 없다.

## 내용 재구성

| 항목 | 수정 |
| --- | --- |
| 제목·초록 | 방법론에 초점을 맞춘 제목, 215단어 평문 초록. 실데이터 발견을 암시하는 제목 제거. |
| 본문 §2--4 | 고정 모형 가정, 투영 후 rank/정보량, 보간 경계, 물리적 comparator 조건, 정확 대조 예제와 추론 중단 조건. |
| 본문 §5 | J0337을 조건부 실패 사례로 배치. 미분 오차·근소한 2일 거절·noise 및 구동 한계를 결과와 함께 제시. |
| 보충자료 | 과학 Sections 1--6 및 부록은 수정된 GRG 판본과 텍스트 일치. 표지·안내·AI 및 데이터 공개 문구만 새 판본에 맞춤. |
| AI 공개 | 실제 연구·코드·문안 작성·내부 검토의 광범위한 사용을 Methods와 Acknowledgements에 유지. 언어 교정만으로 축소하지 않음. |
| 재현 파일 | 소스 ZIP과 별도 Supplementary Code ZIP. 큰 비공개 배열 없이 실행되는 범위와 불가능한 검증을 구분. |

편집 원본은 `output/submission/physica-scripta/manuscript.md` 및 `supplement.md`다. 이전 `paper/manuscript.md`, 원 연구 manifest, GRG 동결 패키지는 보존한다. 현재 패키지 빌더는 `paper/build_physica_scripta.py`이며 GRG manifest의 19개 바인딩부터 확인한다.

## 검증 기록

- `verification/projected_identifiability.py`: 15개 정확 연산 대조 통과.
- `symbolic/state_identifiability.py`: 기존 14개 정확 항등식 통과.
- `verification/replay_public_inference.py`: 6개 조건부 물리 단면과 옴니버스 통계량 재현.
- 저장 comparator 기록의 최솟값을 본문 표와 대조했다. 새로운 fitting이나 노이즈 보정은 하지 않았다.
- 본문·보충자료·커버레터를 Tectonic 0.17.0으로 컴파일했다. undefined citation/reference 및 overfull hbox는 없다. 일부 underfull 경고는 실제 렌더링에서 확인했으며 겹침·잘림·읽기 어려운 수식은 없었다. 본문은 12쪽, 보충자료는 35쪽, 커버레터는 1쪽이다.
- 내장 독립 LaTeX 컴파일러는 `Unable to find standard directories for platform`으로 실패했다. 소스 오류가 아닌 환경 제한이므로 원고를 보존하고 이미 설치된 Tectonic으로 실제 제출 PDF를 생성했다. 원고 편집기를 여는 요청은 유지한다.
- `supplementary-code.zip`을 독립 디렉터리에 풀어 위 세 명령을 다시 실행했으며 모두 통과했다. 원 저장소의 대형 배열은 사용하지 않았다.
- `manuscript-source.zip`을 독립 디렉터리에서 컴파일한 PDF는 전달 본문 PDF의 12쪽 전체와 추출 텍스트 및 렌더링 픽셀이 일치했다.
- 최종 본문 12쪽, 보충자료 35쪽, 커버레터 1쪽 전체를 렌더링해 시각 검수했다. 본문 마지막 참고문헌 한 항목만 다음 쪽으로 넘어가던 조판은 중복 DOI URL 세 개를 제거해 수정했다. DOI와 문헌 정보는 유지했다.
- 새 패키지의 20개 파일 바인딩과 동결 GRG의 19개 파일 바인딩을 확인했다. 보충자료의 과학 Sections 1--6과 부록은 동결 판본과 텍스트가 일치한다.

검수 기록과 실제 PDF 해시는 [physica-scripta-validation.json](../output/submission/physica-scripta-validation.json), 제출 파일·재현 명령·저자 절차는 [패키지 안내](../output/submission/physica-scripta/README.md)에 있다. 이 검수는 파일과 계산 재현성에 관한 것이며 저널 신규성 판단이나 미해결 물리 가정의 인증이 아니다. 공개 저장소의 이전 GRG 태그는 수정하지 않았고, 이번 개정판을 새 공개 태그로 발행하지 않았다. 새 소스와 코드는 제출용 ZIP에 포함했다.

## 제출 경계

대상은 Physica Scripta Research Paper, 구독형 출판이다. [공식 저자 안내](https://publishingsupport.iopscience.iop.org/journals/physica-scripta/)와 [비용 안내](https://publishingsupport.iopscience.iop.org/journals/physica-scripta/about-physica-scripta/)를 사용한다. 소속·교신 이메일·Seoul/Republic of Korea·외부 연구비 없음은 이 세션에서 사용자가 확인한 값이다.

[IOP 저자 정책](https://publishingsupport.iopscience.iop.org/ethical-policy-journals/)은 투고 저자 본인의 제출을 요구한다. 따라서 접수 사실을 만들거나 저자 선언을 새로 대리 확정하지 않고, 검토 가능한 완성 패키지와 웹 양식 입력 자료를 제공한다. 접수번호·최종 PDF·정책 동의는 저자의 실제 제출 뒤 기록한다. GRG 결정 메일에 답장하거나 항의하지 않았고, 새 자동화나 반복 예약을 만들지 않았다.
