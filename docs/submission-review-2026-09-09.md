# 연구 성과 및 투고 준비도 리뷰 — 2026-09-09

> 이 문서는 수정 전 A/B 초안에 대한 리뷰 기록이다. 이후 사용자의 통합 결정과 [최종 수정 기록](unified-revision-2026-09-09.md)이 논문 구성·현재 주장·제출 전략을 대체한다. 아래의 미완료 인용/PDF 및 A 우선·B 별도 제출 판단을 현재 상태로 사용하지 않는다.

검토 기준: 로컬 `36740a6`에서 시작한 Paper A/B 원고, 해당 원고가 가리키는 이론 문서와 저장된 결과. 작업 전 checkpoint: `873e88d`. 이 문서는 연구·편집 리뷰이며, 기존 논문과 실험의 판정을 덮어쓰지 않는다.

**심사자 관점의 판단: 논문화할 성과는 있지만, 현재 원고 두 편 모두 그대로 제출하는 것은 권하지 않는다. A는 정리의 의미와 물리적 연결을 수정한 후 우선 제출할 수 있는 후보이고, B는 조건부 분석 방법 논문으로 범위를 좁히거나 통계적 제약의 근거를 보강해야 한다.** 합격 가능성의 수치화는 하지 않는다.

## 1. 현재까지의 성과

| 상태 | 성과 | 확인 범위와 의미 |
| --- | --- | --- |
| Imported from prior work | 내부 자기상호작용의 단순 재가중은 독립적인 COM 단극자 힘을 만들지 않는다는 결과, 자기중력 제거 시 특정 뉴턴 폴리트로프의 평형이 사라지는 결과 | Paper A §1.2와 Requests 1–2에서 가져온 결과. 이번 검토에서 원래 계산을 다시 실행하지 않았다. 모든 별·물질 모형에 대한 무조건적 no-go로 확대하지 않는다. |
| Proven | 양의 연산자 가중치, 고정 차수, cutoff 아래 유한한 primitive 종 수를 전제로 한 연산자 공간의 유한성 | `docs/theorem-package.md`, `docs/power-counting.md`. 유한성 자체와 관측 불가능성은 다른 명제다. |
| Proven | 명시된 전기 조석 블록과 harmonic scalar potential 가정 A9 아래의 5차원 축약 공간 | `E²`, `tr(E³)`, `(E²)²`, `(DτE)²`, `(∇E)²`. 재실행한 정확 character 계산 및 별도 수치 rank 계산이 일치했다. 일반적인 상대론적 Weyl tensor의 완전 기저라는 뜻은 아니다. |
| Proven | 8개 감사 계열의 수정된 survivor 차원 `5,16,30,15,17,23,17,21` | `verification/tier1_survivor_exact.py`에서 8/8 일치. 비교 대상은 전기, E/B, E/B/S, E/V, E/T, E/Q, E/U, E/Z 블록이다. 계산 방법의 일치는 공통 물리 가정의 타당성까지 증명하지 않는다. |
| Proven | 독립적인 새 primitive를 허용하면 기존 전기 기저 밖의 witness가 생기는, 감사 계열에 한정한 최소 sector 유일성 no-go | 자기 norm 등의 witness가 주어진 reduction을 통과한다. 자연이 실제로 어느 primitive를 가지는지 결정한 결과는 아니다. |
| Proven | 한 상태 완화 모형의 단색 응답과 단열 극한 | `τχ χ̇+χ=αF`에서 `G(iω)=cY+β/(1+iωτχ)`, `β=α cχ`. 단색 및 두 주파수 해의 ODE residual이 모두 정확히 0이었다. |
| Proven | 동적 응답의 구별 가능성에 대한 정확한 제한 | 한 주파수는 `{F,Ḟ}`로 흡수 가능하다. 세 양의 주파수는 공통 실수 계수의 미분 차수 N≤4 비교자를 반박하지만 N≥5는 보간 가능하다. 선형 두 주파수 모형은 합·차 sideband를 만들지 않는다. |
| Counterexample candidate | 알려진 구동·제한된 투영·제한된 비교자 아래에서의 공유 완화시간 응답 | 내부 상태를 검출할 물리적 후보이며, 모든 정적 EFT를 배제하는 보편적 관측량은 아니다. |
| Imported from prior work | PSR J0337+1715의 공개 TOA 12,474개를 사용하는 기존 분석과 비검출 결과 | 저장 JSON의 `Z=2.2851`, `p=0.26`, advanced control `p=0.282`, detection=false를 확인했다. 이번에는 pulsar/Nutimo runtime, 원시 데이터 재적합, 신규 관측 분석을 수행하지 않았다. |
| Proven | 저장된 추정치·오차에서 대표 구간 `1.679527523232755e-9`의 산술 재현 | `sep_gateG2wp.json`의 β 추정치와 Fisher 폭으로 `u95(βhat,10σF)`를 재계산하여 대표 JSON과 상대오차 1e-12 이내 일치했다. 실제 데이터로부터의 전체 파이프라인 독립 재현이나 95% coverage 검증은 아니다. |

평가: 가장 방어하기 좋은 기여는 **명시된 조석 연산자 quotient의 계산, 가정별 실패 경계의 구분, 그리고 제한된 주파수·nuisance 조건 아래의 동적 응답 구별 방법**이다. “자기 질량을 관측할 수 없다는 보편 법칙을 증명했다” 또는 “새로운 중력을 발견했다”라는 논문으로 쓰는 것은 현재 근거를 넘어선다.

## 2. 제출 전 해결해야 할 주요 문제

### R1 — A5 반례는 고정 차수 Taylor 근사의 실패를 입증하지 않는다

대상: `paper/paper-A-collapse-theorem.md:365–387`, `lemmas/55-monopole-jet-collapse.md:21–45`.

**Status: Proven.** `f(Y)=exp(-1/Y²)Θ(Y)`는 0에서 모든 차수의 Taylor jet이 존재하며 0이다. 임의의 유한한 양의 N에 대해 `f(Y)/|Y|^N → 0`이다. 따라서 정확한 함수와 그 무한 Taylor 급수의 동일성은 실패하지만, 0을 이용한 유한 차수 점근 근사는 실패하지 않는다. 이번 symbolic 점검에서 N=5의 극한도 0으로 확인했다.

**Status: Proven.** Lemma 55의 실제 결론에는 `O(Δ>4)`가 붙어 있다. 이 의미의 유한 차수 전개에는 충분한 유한 미분 가능성이면 되고, 해석성은 필요조건이 아니다. 반대로 정확한 함수 복원을 뜻한다면 일반적인 해석 함수도 유한 jet 하나와 정확히 같지는 않다.

조치: “A5를 버리면 finite Taylor jet의 존재가 깨진다”를 철회하고, **정확한 analytic germ 복원**과 **고정 차수 점근 동등성**을 나눠 정리를 쓴다. genuinely non-smooth 반례로 교체하려면 cutoff 안에서 어떤 항이 어떤 나머지 추정을 깨는지 새로 보여야 한다. 현재 smooth-flat 예제의 명백한 수학적 성질은 살리되, 이것이 증명하는 sharpness의 범위를 줄인다.

### R2 — E² 계수와 Nordtvedt 계수의 동일시는 성립하지 않는다

대상: `paper/paper-A-collapse-theorem.md:533–544`.

**Status: Proven.** 외부 점질량 M의 뉴턴 조석장에서 `Eij Eij=6(GM)²/r⁶`이다. 위치에너지에 `C Eij Eij`를 추가하면 그 leading 힘은 `∇(Eij Eij)`에 비례하므로 가속도는 `r⁻⁷`로 스케일한다. 반면 보통의 Nordtvedt 항은 물체의 내부 결합에너지 비율을 곱한 외부 단극자 가속도 `∇Φext ∝ r⁻²`를 수정한다. 상수 조석 계수 하나를 Nordtvedt η와 동일시할 수 없다. 조석 norm과 그 미분의 식은 이번 검토에서 SymPy로 확인했다. Nordtvedt 정의의 문헌 근거: [LLR gravity review](https://link.springer.com/article/10.12942/lrr-2010-7).

조치: “LLR이 첫 survivor를 이미 측정한다”는 문장을 삭제하거나, 구체적인 이론에서 질량비·scalar charge·조석 연산자를 관측 가속도로 잇는 별도의 matching을 제시한다. 기존 LLR 수치를 조석 Wilson 계수의 제약으로 옮겨 적어서는 안 된다. 이 문제는 다섯 조석 invariant의 계산을 무효화하는 것이 아니라, 그 계산의 SEP 해석을 제한한다.

### R3 — 유한 차원성에서 fit 흡수와 유일한 탈출구가 따라오지는 않는다

대상: `paper/paper-B-dynamic-sep-limit.md:12–19,73–104,131–156,490–497`.

**Status: Proven.** 신호 템플릿 T가 nuisance Jacobian J에 흡수된다는 조건은 `T∈col(J)`이다. T의 매개변수 수가 유한하다는 사실은 이 포함관계를 주지 않는다. 기존 static-Δ 열의 큰 흡수율도 그 열과 해당 구성에서 측정한 사실이며 모든 static sensitivity의 보편적 흡수를 증명하지 않는다.

**Status: Proven.** Paper A 자체가 A3/A4/A5/A8의 여러 탈출 경로를 열어 둔다. 그 중 A4만 연구 대상으로 선택했다는 사실을 “정리가 허용하는 유일한 관측량”이라고 바꿀 수 없다. 추가적인 모든 상태 모형이 `(β,τχ)` 두 수로 표현되는 것도 아니다.

**Status: Proven.** 전 주파수 대역의 rational transfer와 세 점의 주파수 표본은 다른 문제다. 세 양의 주파수의 켤레쌍을 모두 맞추는 실수 5차 다항식이 존재한다. 예를 들어 `τ=β=1`, `cY=0`, `ω=1,2,3`에서

```text
P5(z) = -z^5/100 + z^4/100 - 3z^3/20 + 3z^2/20 - 16z/25 + 16/25
P5(±iω) = 1/(1±iω),  ω=1,2,3.
```

이 식의 6개 등식과 실수 계수는 정확 symbolic assert로 확인했다. 일반 세 주파수에서도 켤레대칭 보간이 같은 경계를 준다. 레포의 `docs/nonadiabatic-regime.md`와 failure ledger에는 이 경계가 이미 들어 있다.

조치: B를 “선택한 single-pole 응답을 지정한 comparator/nuisance 아래에서 제약한다”로 다시 쓴다. 광범위한 “only”, “any finite comparator”, “any mechanism”을 삭제하고, 주파수 표본 수·미분 차수·구동과 투영 가정을 초록부터 명시한다.

### R4 — 대표 상한은 통계적 가정과 nuisance 절단에 강하게 의존한다

대상: `paper/paper-B-dynamic-sep-limit.md:256–307,414–486`와 `request10_external/sep_dynamic/sep_phase_marg_10_8e.json`.

**Status: Proven.** 저장된 대표 값은 `u95(βhat,KσF)`라는 Gaussian 구간 계산이며 K=10은 측정된 noise parameter가 아니라 선택한 safety floor다. 원고 §4.1과 §7도 coverage-calibrated confidence가 아니라고 명시한다. 그 구간을 계산한 사실과 실제 오차 모형에 대해 보장되는 95% 상한은 구분해야 한다.

**Status: Proven.** 같은 K=10에서 저장된 full-rank와 truncated 결과의 차이는 다음과 같다. 기존 JSON에서 직접 나눗셈하여 확인했다.

| τχ (일) | truncated 구간 | full-rank 구간 | full / truncated |
| --- | ---: | ---: | ---: |
| 2 | 1.680e-9 | 3.534e-9 | 2.10 |
| 5 | 1.852e-9 | 9.673e-9 | 5.22 |
| 18 | 1.976e-9 | 3.166e-8 | 16.02 |
| 52 | 2.600e-9 | 4.488e-8 | 17.26 |
| 200 | 7.155e-9 | 1.241e-7 | 17.35 |

**Status: Proven.** red-noise Fourier 수 30→60의 저장 구간 변화는 약 0.52–0.60%다. 이것은 해당 Fourier-order 비교의 안정성을 보여주지만, 위의 SVD 절단 의존성을 해소하지 않는다. 71차원과 90차원 nuisance 공간을 비교하는 것은 다른 변경이다.

조치: 절단한 19개 방향을 제외할 정당한 prior/물리 가정을 쓰거나, 지원되는 nuisance 공간 전체를 포함한 결과를 주결과로 삼는다. K 보정은 정당화된 likelihood/추가 오차 모형 또는 적절한 반복 검증과 연결해야 한다. 주입·회수의 국소 성공을 모든 오차 모형의 구간 보장으로 해석하지 않는다. 이번 작업은 이 리뷰 범위에서 멈추며 외부 runtime 재실행을 제안된 조치와 혼동하지 않는다.

판단: 현 상태에서는 **“지정된 선형화·잡음·절단 가정 아래의 조건부 민감도/구간”**이 방어 가능한 표현이다. 발표된 static SEP 신뢰한계보다 “1000배 강한 보편적 제약”이라는 표현은 피해야 한다. static benchmark도 모형에 따라 달라진다. 관련 원 논문: [Voisin et al.](https://arxiv.org/abs/2411.10066).

### R5 — β와 실제 Δ(t) 진폭, 현상론적 구현과 질량 모형을 구분해야 한다

대상: `paper/paper-B-dynamic-sep-limit.md:158–175,263–281,416–450`, `request10_external/scripts/sep_phase_marg_10_8e.py:74–92`, `request10_external/scripts/sep_common.py:138–181`.

**Status: Proven.** 실제 계산은 `[TcY,Tβ]`의 두 번째 계수 β에 대한 구간을 만든다. unit-drive single-pole 부분의 각 carrier 진폭은

```text
|δΔχ,k| = |β d_k| / sqrt(1 + ω_k^2 τχ^2).
```

**Status: Proven.** 예를 들어 τχ=200일, unit drive에서 이 배율은 inner 약 0.001297, outer 약 0.2520, difference 약 0.001303이다. 따라서 β, 특정 carrier 진폭, 합성 신호의 peak/RMS는 같은 양이 아니다. 함께 추정하는 `cY` 성분까지 포함한 총 Δ(t)의 최대값은 β 하나의 구간과도 다르다.

**Status: Counterexample candidate.** 초기 `m(Y,χ)` worldline 질량 모형을 적분기의 pairwise coupling `Δ(t)`로 옮기는 것은 추가적인 물리적 모델링 선택이다. 구체적인 작용·운동방정식·구동 정의 없이 두 표현의 일반적 동등성이 확립된 것으로 취급하면 안 된다.

조치: 주 표를 `|β|` 또는 명확히 정의한 normalized template amplitude로 표기하고, carrier별/peak/RMS Δ 진폭은 별도 변환한다. 모든 구동 모델에 독립적인 제약이라는 문장을 제한한다. 기존 적분기 구현을 현상론적 benchmark로 제시하는 선택은 가능하지만, 그 경우 최초 질량 모형에서 자동으로 유도되었다는 주장을 제거해야 한다.

### R6 — 선행연구 대비 새로움과 독립적으로 읽히는 증명이 부족하다

**Status: Imported from prior work.** 동적 multipole, 응답함수의 pole, 비단열적 내부 자유도를 worldline에 추가하는 방법은 기존 연구에 있다. 이번 검색에서 최소한 다음 비교 대상들을 확인했다.

| 원 연구 | 비교해야 할 내용 |
| --- | --- |
| [Chakrabarti–Delsate–Steinhoff 2013](https://arxiv.org/abs/1306.5820) | compact object의 EFT 선형응답, pole과 내부 mode |
| [Steinhoff et al. 2016](https://arxiv.org/abs/1608.01907) | 동적 quadrupole을 포함한 worldline action, 단열 근사의 실패 |
| [Khalil et al. 2022](https://arxiv.org/abs/2206.13233) | 단극자 scalarization의 비단열 동역학; 단순 조석 연구보다 A4와 가까운 비교 대상 |

이 세 arXiv 항목은 현재 `paper/references.bib`에 없다. 이것만으로 새로움 부재를 확정할 수는 없지만, 내부 상태 또는 relaxation pole의 존재 자체를 독창성으로 내세우기에는 선행연구 비교가 부족하다. 이번 검색은 체계적 priority 전수조사가 아니므로 “최초”도 확인하지 않았다.

조치: 기여를 구체적인 quotient, 정확한 comparator 경계, 지정된 J0337 채널의 분석으로 좁히고 위 연구와 차이를 식으로 비교한다. A의 부록은 현재 실제 증명 대신 레포 경로 목록인 부분이 크다(`paper-A-collapse-theorem.md:549–567`). 핵심 정의·power counting·증명과 오차항을 원고/제출용 부록에 넣어 심사자가 저장소를 탐색하지 않고 읽을 수 있게 한다.

### 추가 정리 항목

- **Status: Proven.** A9가 “purely electric”인 기준 sector와 독립적인 B 등 추가 primitive를 허용한 조건부 감사 sector는 같은 물리 가정의 목록이 아니다. 두 경우를 분리해서 서술하고, A9 아래의 전기 gradient 계산을 일반 상대론의 magnetic/gradient 완전성으로 확대하지 않는다.
- **Status: Proven.** B의 `phase marginalization`은 prior 적분이 아니라 저장된 시간원점 grid에서의 최댓값이다. 비정수 inner/outer 주파수비에서 outer period 하나는 모든 common-origin phase의 정확한 공통 주기가 아니다. 현재 `[0,Pout)` grid에 한정된 envelope로 쓰고, 임의 위상 보장을 주장하려면 위상공간 범위와 grid 오차를 별도로 정당화한다.
- **Status: Proven.** B §6.1의 두 기준 위상 template overlap과 양쪽 maximum statistic의 비슷한 값만으로 전체 위상 영역의 near-collinearity는 따라오지 않는다.
- 기존 verification README에는 과거 `E/B/S=33` 문구와 exact high-rank 미검증 문구가 남아 최신 30 및 exact character 결과와 충돌한다. 과학적 결과를 바꾸는 작업보다 제출용 증거 문서를 하나의 상태로 정리하는 작업이 필요하다.
- `paper/`에 완성 PDF가 없고, 두 원고의 인용은 아직 draft 목록이다. 이번 리뷰는 LaTeX 컴파일·시각 검수를 하지 않았으므로 출력 품질을 승인하지 않는다.

## 3. 투고처 검토

다음 우선순위는 공식 범위·논문 유형을 확인한 뒤 내린 편집적 판단이다. 해당 저널의 수락 의사나 acceptance probability를 확인한 것은 아니다. 2026-09-09 확인.

| 투고처 | 권하는 형태와 순위 | 공식 근거와 남은 조건 |
| --- | --- | --- |
| **General Relativity and Gravitation (GRG)** | A의 1차 후보: 범위를 좁힌 이론 연구 논문 | 이론·상대론적 천체물리·중력 검증을 포괄한다. A5, SEP matching, 독창성 및 자체 완결 증명을 먼저 해결한다. [공식 범위](https://link.springer.com/journal/10714/aims-and-scope) |
| **Classical and Quantum Gravity (CQG)** | A를 짧게 정리하면 **Note** 후보; B를 조건부 관측 방법론으로 정리할 경우 우선 검토 | 중력 이론·실험·관련 데이터 분석을 포함한다. 유용한 짧은 새 결과를 위한 Note 유형이 실제로 있다. 완성된 Research Paper는 유의미한 진전을 요구한다. [범위·논문 유형](https://publishingsupport.iopscience.iop.org/journals/classical-and-quantum-gravity/about-classical-quantum-gravity/) |
| **Physical Review D (PRD)** | B의 통계·물리 해석이 정리되었을 때의 도전적 후보 | compact object, 중력 이론 검증, 천체물리 제약과 데이터 분석이 범위에 든다. 공식 수락 기준은 중요하고 실질적인 지식의 추가다. 현재의 K=10 대표 제약만으로 이 기준을 충족한다고 판단하지 않는다. [공식 범위와 기준](https://journals.aps.org/prd/about) |
| **Astronomy & Astrophysics (A&A)** | B를 실제 pulsar timing 분석 논문으로 완성했을 때의 대안 | 바로 같은 시스템의 기준 분석이 A&A 693 A143에 게재되어 독자층이 맞는다. [해당 원 논문·저널 정보](https://arxiv.org/abs/2411.10066). 이번 환경에서는 공식 저자 가이드 전문 접근에 실패했으므로 최신 비용·세부 제출 요건 확인은 미완료다. |

제출 요건에서 현재 구체적으로 준비할 것은 저널 형식의 영문 원고, 자체 완결 부록, 실제 인용 연결, 입력·코드·결과의 고정된 공개 snapshot, Data Availability Statement다. GRG는 원고의 editable source 제출과 원 연구의 data availability statement를 요구한다. [GRG 제출 지침](https://link.springer.com/journal/10714/submission-guidelines)

비용이 제약이면 CQG는 공식 안내상 **subscription 방식 출판과 투고 비용이 없고**, gold OA만 선택 APC가 있다. 이를 OA 출판도 항상 무료라는 뜻으로 해석해서는 안 된다. [CQG 비용 안내](https://publishingsupport.iopscience.iop.org/journals/classical-and-quantum-gravity/about-classical-quantum-gravity/)

## 4. 가장 짧은 제출 경로

1. **A를 먼저 정리한다.** R1·R2를 해결하고, 관측 불가능성 대신 명시된 가정 아래의 연산자 분류와 실패 경계를 논문의 중심에 둔다. B의 미보정 상한을 A의 정리 가치를 보증하는 근거로 쓰지 않는다.
2. **B의 논문 유형을 결정한다.** 기존 결과만 활용한다면 “single-pole SEP template의 조건부 민감도와 식별 가능성”으로 쓰는 것이 정직하다. 보편적·정밀한 관측 제약을 주결론으로 유지하려면 R3–R5에 답해야 한다.
3. **선행연구 비교와 원고를 완결한다.** 새 물리의 발명, 내부 상태가 있는 EFT라는 이미 알려진 틀, 여기서 새로 계산한 quotient/식별 경계를 분리한다. 실제 부록·그림·인용·재현 snapshot을 준비한다.
4. **완결된 A를 GRG 또는 CQG Note에 순차적으로 고려한다.** B의 큰 수치 개선을 앞세워 지금 PRD에 보내는 경로는 권하지 않는다. preprint도 위 핵심 오류를 고친 뒤 공개하는 편이 낫다.

두 편 분리는 조건부로 타당하다. A의 독립적인 수학적 기여와 B의 독립적인 관측 방법이 각각 선행연구 대비 의미를 가질 때 유지한다. A가 기존 EFT의 재진술에 그친다면 긴 독립 논문을 유지하기보다 짧은 Note 또는 B의 이론 부록으로 합치는 편이 낫다.

## 5. 이번 검증의 정확한 범위

Python 실행 환경: `C:/Users/lpaiu/AppData/Local/hermes/hermes-agent/venv/Scripts/python.exe`, NumPy 1.26.4, SymPy 1.14.0.

다음 기존 진입점을 오류 없이 실행했다. 각 스크립트의 내부 assertion 수를 알 수 없는 상태에서 합산 테스트 개수로 포장하지 않는다.

```powershell
rtk proxy python symbolic/checks/test_symbolic_smoke.py
rtk proxy python verification/tier1_survivor_exact.py
rtk proxy python symbolic/chi_relaxation_response.py
rtk proxy python symbolic/chi_two_frequency_response.py
rtk proxy python verification/verify_ce_a3a4.py
rtk proxy python symbolic/frequency_sweep_distinguishability.py
rtk proxy python verification/verify_identities.py
rtk proxy python verification/verify_survivors.py
```

- 정확 character 결과: 8/8 sectors 일치.
- 단색·두 주파수 ODE residual: 각각 0. 선형 두 주파수 sum/difference sideband: 각각 0.
- mixed quartic symbolic identity residual: 0. electric basis numerical rank: 5.
- 별도의 짧은 symbolic 계산: E²의 r⁻⁶ 및 힘의 r⁻⁷ 스케일, smooth-flat/Y⁵ 극한 0, 세 carrier를 맞추는 실수 5차 보간식의 6개 항등식 확인.
- 저장된 headline β·σ에서 구간 산술 재현; 저장 JSON들의 full-rank/truncated 비율 및 Fourier-order 비율 재계산.
- 주파수 스크립트가 재작성한 산출물은 줄바꿈만 달라졌고, 원래 저장본으로 복구했다.

전체 71-lemma 체계의 형식 검증, 전체 monolithic symbolic suite, 원시 TOA부터의 runtime 독립 재현, 모든 사전등록 commit의 시간 선후 전수감사, 출판 novelty 전수검색은 이번 검증 범위에 포함하지 않았다.

**작업 분류: theorem progress 및 loophole progress에 대한 검토와 주장 경계 명료화. 새로운 관측 검출이나 검증된 보편 SEP 상한의 추가가 아니다.**
