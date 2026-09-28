# 통합 원고 최종 수정 기록 — 2026-09-09

기준: 이전 리뷰 commit `018d255`, 작업 전 checkpoint `1214d06`. 사용자의 요청에 따라 Paper A/B를 **하나의 원고**로 통합했다. 현재 소스는 [paper/manuscript.md](../paper/manuscript.md), 생성 LaTeX는 [paper/main.tex](../paper/main.tex)이다. 이전 A/B 파일은 역사적 초안이라는 표시를 붙였으며 기본 빌드·README·논문 목차를 통합본으로 전환했다.

최종 제목: **Static response and dynamical identifiability in free-fall tests: finite-order boundaries and a pulsar-triple application**.

논문의 중심 질문은 “정해진 차수에서 유한하게 표현되는 국소 응답과 내부 상태의 지연 응답을, 실제 관측에서는 어떤 조건에서 구별할 수 있는가?”이다. 전기 조석 quotient는 구체적인 정적 기준 모형이고, J0337은 구별 조건과 nuisance 의존성을 보여주는 조건부 응용이다. 두 모형 사이에 유도되지 않은 물리적 동일성을 가정하지 않는다.

## 이전 지적의 처리

| 항목 | 수정 내용 | 현재 근거와 한계 |
| --- | --- | --- |
| R1 — A5와 finite jet | Proposition 2에 C^(D+1) 정칙성, 양의 정수 가중치, O(epsilon^(D+1)) 나머지를 명시했다. Lemma 55/58/59, 모델 표, 예제 코드와 호출 테스트도 수정했다. | Proven: smooth-flat은 모든 유한 jet을 가지며 정확한 Taylor 급수 복원만 깨뜨린다. 제곱근 threshold는 별도로 정칙성을 실패한다. |
| R2 — 조석/SEP 동일시 | I2=6(GM)^2/r^6, 조석 가속도 r^-7와 상수 Nordtvedt 가속도 r^-2를 구분하고 LLR 수치의 직접 전환을 삭제했다. | Proven: 계수 간 동일성에는 별도 matching이 필요하다. |
| R3 — 유한성/흡수/유일한 escape | Theorem 3에 실수 degree-N 비교자의 K개 양의 carrier 보간 필요충분조건 N>=2K-1과 구성식을 넣었다. 실제 관측에는 whitened nuisance rank 조건을 추가했다. | Proven: 세 carrier는 자유로운 실수 5차 비교자로 정확히 맞출 수 있다. A4가 유일한 escape라는 문구는 철회했다. |
| R4 — K 및 nuisance | Gaussian U의 정의, K=10의 가정적 성격, 71/90차원 결과를 같은 표와 그림에 표시했다. 전체 공간 결과를 표의 첫 수치 열로 배치했다. | Imported from prior work: 같은 K에서 약 2.10–17.35배 차이가 난다. coverage-calibrated 보편 SEP 상한이나 1000배 개선으로 주장하지 않는다. |
| R5 — beta 진폭/물리 구현 | beta의 drive 정규화, carrier 감쇠, peak/RMS 및 co-fitted c_Y의 차이를 명시했다. 고정 관성질량과 prescribed drive를 갖는 pair potential로 현상론적 힘을 정의했다. | Proven: beta 구간만으로 총 Delta의 peak를 제한할 수 없다. Counterexample candidate: 원래 field-dependent mass와의 물리적 matching은 미완료다. |
| R6 — 독창성/자체 완결성 | 동적 EFT 선행연구 세 편을 인용하고 원고 내부에 quotient 증명, 보간 증명과 구체적 P5 예제를 넣었다. | Imported from prior work: 내부 mode 자체는 기존 EFT에 있다. 새로운 물리 원리나 최초라는 주장은 하지 않는다. |

## 함께 수정한 사항

- Proven: 다섯 대표는 선형 연산자 기저이며, I2와 I2^2가 대수적으로 독립인 좌표라는 표현을 제거했다. Appendix A의 독립성은 주기 궤적의 작용 적분으로 설명하여 점별 value rank와 total-derivative quotient를 구별했다.
- Proven: A9의 harmonic Newtonian potential과 임의의 relativistic electric-Weyl tensor를 구분했다. 선택적 B 등 추가 family census는 별도의 대수적 확장이다.
- Proven: 일반 해의 초기 transient를 명시했다. beta=0 또는 omega=0이라는 이유만으로 모든 초기상태 신호가 사라진다고 하지 않는다.
- Proven: 선형 두 주파수 모델은 sideband를 만들지 않고, 비선형 정적 모델도 sideband를 만들 수 있음을 명시했다.
- Proven: 유한 상태의 rational transfer 논증은 linear time-invariant 모형으로 한정했다.
- Proven: phase marginalization을 유한한 등록 시간원점 grid의 최댓값으로 다시 이름 붙였다. [0,Pout)는 모든 상대위상을 보장하는 영역이 아니다.
- Imported from prior work: 두 기준 원점의 overlap, 국소 injection/recovery, 제한된 turn lattice의 성공을 전체 위상·잡음·pulse-slip 공간의 검증으로 확대하지 않는다. 실패·수정 gate의 이력은 그대로 보존했다.
- verification README의 오래된 E/B/S=33 및 exact high-rank 미검증 문구를 30과 실제 character 검증 상태에 맞췄다.
- 참고문헌을 본문 번호 인용과 BibTeX로 연결했다. 공식 DataCite 기록으로 Zenodo 13899771의 정확한 제목·저자 및 발행연도 **2025**를 확인했다.

공식 확인 출처: [Chakrabarti et al. 2013](https://arxiv.org/abs/1306.5820), [Steinhoff et al. 2016](https://arxiv.org/abs/1608.01907), [Khalil et al. 2022](https://arxiv.org/abs/2206.13233), [DataCite dataset metadata](https://api.datacite.org/dois/10.5281/zenodo.13899771).

## 검증 및 산출물

다음 기존 8개 진입점을 오류 없이 실행했다: symbolic smoke, exact survivor character, 단색 ODE, 두 주파수 ODE, A3/A4 경계, frequency sweep, reduction identities, electric survivors. 정확 character는 8/8 sector 일치, ODE residual과 선형 sum/difference sideband는 0, mixed quartic 항등식의 symbolic residual은 0, electric value-matrix rank는 5였다.

추가된 `verification/verify_unified_paper.py`는 smooth-flat 미분·나머지 극한, threshold 실패, K=1–4의 보간 및 낮은 차수의 불가능성, P5의 여섯 등식, ODE·미분전개 잔차, 조석/SEP r 의존성, 투영 예제, 저장된 추정치로부터의 대표 U 및 표의 모든 행을 검사한다. 이는 수학·산술 검증이며 원시 timing 데이터의 재분석은 아니다.

빌드는 기존 Pandoc builder를 보완해 사용했다. 최종 PDF의 12페이지를 렌더링하여 수식·표·그림·인용·페이지 배치를 검수했다. 최종 LaTeX 로그에는 미해결 인용·참조, missing character, overfull/underfull box 경고가 없다. 로컬 fontconfig 설정 메시지는 있었으나 실제 사용한 Latin Modern 수학·본문 글꼴과 출력은 정상 확인했다.

SHA-256 입력 목록은 [revision-manifest.json](../paper/revision-manifest.json)에 있다. 대표 결과, full/truncated 비교, gate·overlap·Fourier 비교 JSON과 실제 저장된 수치 입력을 결속한다. 기존 REQUEST10 원시 결과와 등록 판정을 수정하지 않았다. 주파수 sweep이 재생성한 analytic JSON/TSV는 기존 내용과 같음을 확인했다.

최종 manifest의 27개 파일 해시와 표의 산술 검사를 통과했다. 수정한 nonanalytic 호출 테스트와 별도 invariant 검증도 통과했다. [LaTeX 소스 ZIP](../output/submission/free-fall-identifiability-source.zip)은 별도 임시 디렉터리에 풀어 다시 컴파일하여 포함된 소스·그림·참고문헌만으로 PDF가 만들어짐을 확인했다.

## 제출 판단과 남은 실제 연구

후속 갱신: 아래 내용은 통합 직후의 상태다. 사용자가 승인한 1/2/3 후속 검증의 완료 결과와 남은 범위는 [연구 완료 기록](research-completion-2026-09-09.md)에 있다. 현재 PDF는 16페이지이며 nuisance/coverage 감사와 조건부 힘 matching을 포함한다.

현재 형태는 하나의 **조건부 이론·방법론 원고**이다. 이전 리뷰의 “A 선행 제출/B 별도 제출” 전략은 사용자의 통합 결정에 따라 철회한다. 통합 원고의 독창성은 제한된 comparator의 정확한 경계와 관측 해석에 있다. 저널 수락이나 선행연구 전체에 대한 우선권은 검증하지 않았다.

Status: Conjectural. 실제 compact-body 이론에서의 coupling/drive matching, 정당화된 공분산·likelihood, nuisance 약한 방향에 대한 prior, 구간 coverage 검증은 여전히 연구 과제다. 이를 완료했다고 가장하지 않고 해당 결과를 요구하지 않는 주장 범위로 원고를 수정했다. 보편적·정밀한 경험적 SEP 제약을 주결론으로 되돌리려면 이 과제들이 필요하다.

원고는 로컬에서 수정·빌드·검수했다. 실제 저널 제출이나 공개 저장소 push, DOI 등록은 수행하지 않았다. 최종 제출 시 저자 소속·연락처, 저널 서식 및 공개 revision 접근성을 실제 제출 정보에 맞추어 확인해야 한다.

작업 분류: **theorem progress 및 loophole progress**. 새 관측 검출이나 새 coverage-calibrated 상한을 추가한 작업으로 분류하지 않는다.
