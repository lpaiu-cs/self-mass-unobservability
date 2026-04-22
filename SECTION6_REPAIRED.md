# Section 6 Repair: Executable Calculation and Simulation Requests

## 목적

기존 6절 요청서는 좋은 연구 방향을 담고 있지만, 다음 문제가 섞여 있었다.

- "이론 정식화가 먼저 필요한 단계"와 "바로 돌릴 수 있는 수치 적합 단계"가 분리되어 있지 않다.
- 장난감 적분기 수준의 모델로는 실제 LLR/pulsar timing 제약을 재현할 수 없는데, 원문은 둘을 구분하지 않았다.
- strong-field self-gravity 파라미터 `s`를 어떤 수준에서 고정하고 어떤 수준에서 추론할지 명시되지 않았다.
- clock-sector 변형은 관측량으로 내려오는 중간 유도식이 빠져 있다.

아래는 이를 보완한 **실행 가능한 6단계 요청서**다. 핵심은

1. 먼저 이론적 일관성을 닫고,
2. 그 다음 장난감/모의자료로 검증하고,
3. 마지막에 실데이터 적합으로 넘어가는 것이다.

## 공통 정의

전 단계에서 동일하게 다음 EFT를 사용한다.

```math
\frac{m_{P,A}}{m_{I,A}} = 1+\sigma_1 s_A+\sigma_2 s_A^2,
\qquad
s_A := \frac{E_{g,A}}{m_A c^2} > 0.
```

clock sector는

```math
\frac{d\tau_A}{dt}
=
1-\frac{v_A^2}{2c^2}
-\left(1+\zeta_1 s_A+\zeta_2 s_A^2\right)\frac{U_{\rm ext}}{c^2}
```

로 둔다. `s_A`는 천체별 고정 파라미터가 아니라, 필요한 경우 EOS 또는 구조모형 prior를 둔 잠재변수로 취급한다.

## Request 1. COM Decoupling의 일반 증명

**목표**

원문의 핵심 주장인 "self-subtraction만으로는 COM 새 force가 생기지 않는다"를 일반 형태로 닫는다.

**요청**

```text
임의의 compact extended body에 대해

L=\int d^3\xi\,\rho(\xi)\left[\frac12(\dot X+\dot \xi)^2-\Phi_{\rm ext}(X+\xi)-(1-\lambda)\Phi_{\rm self}(\xi)\right]

를 중심질량 좌표계에서 quadrupole order까지 전개하라.

1. 임의의 정적 밀도분포에 대해
   L = (1/2)M\dot X^2 - M\Phi_ext(X) - (1/2)Q^{ij}\partial_i\partial_j\Phi_ext(X) + L_int(\lambda) + O(\ell^3\nabla^3\Phi_ext)
   꼴을 유도하라.
2. Euler-Lagrange 방정식에서 COM force가 `lambda`와 무관함을 보이라.
3. 균일구와 Newtonian polytrope n=1에 대해 `Q^{ij}=0`가 되는 구대칭 평형 예를 써서 결과를 검산하라.
4. 부록으로 self-term이 내부에너지와 정역학 조건에는 들어가지만 COM 운동에는 들어가지 않음을 한 페이지로 요약하라.
```

**산출물**

- 심볼릭 유도 노트북 1개
- 짧은 유도 메모 1개

**성공 기준**

- COM 방정식에 `lambda`가 남지 않는다.
- 구대칭 예시에서는 quadrupole 항도 0이어서 COM 운동이 정확히 test-body 형태로 떨어진다.

## Request 2. 내부구조 일관성 검증

원문 6절에는 빠져 있었지만, 본 연구에서 반드시 먼저 해야 할 계산이다. 그렇지 않으면 "self-mass unobservability를 내부구조까지 literal하게 적용할 수 없다"는 주장 자체가 증명되지 않는다.

**요청**

```text
Newtonian hydrostatic equilibrium과, 가능하면 TOV toy model에서
self-gravity를 내부구조 방정식에 `(1-lambda)` 배로 넣었을 때의 정역학을 조사하라.

1. 균일구 toy model에서 central pressure, binding energy, virial relation을 `lambda`의 함수로 구하라.
2. polytrope n=1에 대해 Lane-Emden 해를 사용해 radius, central density, binding energy의 `lambda` 스케일링을 구하라.
3. `lambda -> 1` 극한에서 bound equilibrium이 사라지거나 비물리적으로 변하는지 보이라.
4. 결론으로, self-unobservability는 내부구조가 아니라 passive-coupling/clock sector EFT로만 써야 함을 명시하라.
```

**산출물**

- `lambda`에 따른 `R(lambda)`, `P_c(lambda)`, `E_g(lambda)` 도표
- 2쪽 이내의 consistency memo

**성공 기준**

- `lambda=1`의 literal 내부구조 적용이 허용되지 않음을 정량적으로 보인다.

## Request 3. LLR: 해석적 민감도 + 모의자료 회수

원문 6절의 LLR 요청은 바로 "실데이터 posterior"를 요구하지만, 단순 1PN 적분기로는 그 수준의 제약을 재현할 수 없다. 먼저 analytic response와 mock-data 회수부터 해야 한다.

**요청**

```text
Sun-Earth-Moon 계의 SEP/Nordtvedt형 구동을 이용해

delta_SEP := sigma_1 (s_E-s_M) + sigma_2 (s_E^2-s_M^2)

가 lunar range에 만드는 synodic 성분을 계산하라.

1. 해석적으로 `cos D` 성분의 leading amplitude를 유도하라.
2. Sun-Earth-Moon barycentric 1PN 적분기에서 passive factor
   1 + sigma_1 s_i + sigma_2 s_i^2
   를 넣어 synthetic normal points를 생성하라.
3. 다음 nuisance를 포함해 injection-recovery를 수행하라:
   - 초기 상태벡터 오차
   - station/reflector range bias
   - solar radiation pressure coefficient
   - thermal expansion coefficient
   - lunar Love number `k2` 또는 등가 tidal parameter
4. 출력은
   (a) range residual,
   (b) synodic `cos D` amplitude,
   (c) recovered `sigma_1, sigma_2` posterior
   로 하라.
5. `sigma_1=sigma_2=0` injected case에서 unbiased recovery가 되는지 확인하라.
```

**산출물**

- analytic note 1개
- mock-data fit notebook 1개

**주의**

- 이 단계의 결과는 "감도 예측"이지, 곧바로 실데이터 경쟁 제약이 아니다.

## Request 4. LLR: 실데이터 적합은 기존 LLR 파이프라인 위에서만 수행

실제 LLR 제약을 주장하려면 toy integrator가 아니라 기존 normal-point/partial-derivative 파이프라인에 EFT 항을 넣어야 한다.

**요청**

```text
공개 LLR normal points 또는 기존 LLR codebase를 사용해
SEP-violation parameter를

delta_SEP = sigma_1 (s_E-s_M) + sigma_2 (s_E^2-s_M^2)

로 재매핑하라.

1. GR baseline fit을 먼저 재현하라.
2. 그 다음 linearized partials 또는 full Bayesian sampling으로 `sigma_1, sigma_2`를 적합하라.
3. nuisance set에는 최소한 다음을 포함하라:
   - Earth-Moon 초기 상태
   - station/reflector bias
   - Earth orientation / troposphere surrogate
   - lunar tidal parameter
   - SRP/thermal surrogate
4. toy integrator 결과와 실데이터 결과를 분리 보고하라.
```

**산출물**

- 실데이터 적합 보고서 1개
- toy-vs-real 비교표 1개

**성공 기준**

- `sigma_1` 한계가 published Nordtvedt/LLR 규모와 같은 order에 도달한다.
- `sigma_2`는 약한 장에서 거의 안 보인다는 점을 posterior로 확인한다.

## Request 5. PSR J0337+1715: 두 단계 분석으로 분리

원문의 J0337 요청도 "공개 자료로 즉시 가능한 것"과 "전용 timing code가 있어야 가능한 것"이 섞여 있다. 먼저 번역 가능한 posterior를 만들고, 그 다음 full TOA fit으로 간다.

**요청**

```text
두 단계로 수행하라.

Phase A. Published SEP bound translation
1. 관측 제약을 `Delta_obs` posterior 또는 상한으로 입력받아
   Delta_SMU = sigma_1 (s_NS-s_WD) + sigma_2 (s_NS^2-s_WD^2)
   로 매핑하라.
2. EOS prior는 `s_NS in [0.1, 0.2]`로 두고 `s_WD`는 고정 또는 매우 좁은 prior로 처리하라.
3. 출력은 `sigma_1, sigma_2` posterior와 EOS-prior sensitivity이다.

Phase B. Full TOA refit
1. 공개 TOA와 검증된 three-body timing code가 있을 때만 수행하라.
2. 외부 WD가 inner binary barycenter에 주는 differential acceleration에 `Delta_SMU`를 곱해 동역학에 삽입하라.
3. white-noise only, white+red-noise, ephemeris-alternative의 세 nuisance model을 비교하라.
4. 출력은 `sigma_1, sigma_2` posterior와 Bayes factor이다.
```

**산출물**

- Phase A posterior notebook
- Phase B timing-fit plan 또는 실행 결과

**성공 기준**

- raw TOA가 없더라도 Phase A만으로 즉시 strong-field posterior를 낼 수 있다.
- Phase B는 가능할 때만 경쟁 제약으로 승격한다.

## Request 6. Clock Sector: 관측식 유도 후 joint pulsar fit

원문의 clock-sector 요청은 그대로는 아직 관측량 수준에 내려오지 않았다. 먼저 timing observable과의 연결식을 닫아야 한다.

**요청**

```text
자유낙하 sector는 GR로 유지하고 clock sector만

d tau / dt = 1 - v^2/(2c^2) - [1 + zeta_1 s + zeta_2 s^2] U_ext/c^2

로 변형한 모형을 관측 가능한 pulsar timing 지연식으로 내리라.

1. eccentric binary에서 Einstein delay parameter `gamma`에 대한 수정식을 유도하라.
2. metric-like tied model (`zeta_i = sigma_i`)와
   decoupled clock-only model (`zeta_i` independent, `sigma_i = 0`)
   을 분리해 정의하라.
3. `gamma`가 측정된 binary pulsar 집합을 골라 joint likelihood를 구성하라.
4. 질량-기하 파라미터와의 degeneracy를 끊기 위해,
   가능한 경우 다른 post-Keplerian parameter 또는 Shapiro 정보도 함께 사용하라.
5. 출력은 `zeta_1, zeta_2` posterior, tied-vs-decoupled Bayes factor,
   그리고 system-by-system pull plot이다.
```

**산출물**

- clock-sector derivation memo
- binary-pulsar joint-fit notebook

**성공 기준**

- `gamma_obs = (1 + zeta_1 s + zeta_2 s^2) gamma_GR`라는 ansatz가 어떤 근사에서 정당한지 명시된다.
- free-fall sector를 건드리지 않고도 clock-sector만으로 남는 제약이 무엇인지 분리된다.

**중단 기준**

- low-side source를 covariance-aware likelihood 또는 TOA-level timing까지
  한 번 더 밀어도 `|kappa_*|_95`가 `10^-2`대로 내려오지 않으면,
  Request 6은 더 이상 주 novelty engine으로 밀지 않는다.
- 그 경우 Request 6은 support/local-audit section으로 고정하고,
  tied-vs-decoupled의 본판정은 joint free-fall-plus-clock analysis로
  넘긴다.

## 실행 순서

실제 작업 순서는 아래가 맞다.

1. Request 1
2. Request 2
3. Request 3
4. Request 5 Phase A
5. Request 6의 유도 부분
6. Request 4와 Request 5 Phase B는 전용 데이터/코드가 확보된 뒤 수행

## 핵심 수정 요약

- 원문 6절의 4개 요청을 **6개 워크패키지**로 쪼갰다.
- 빠져 있던 **내부구조 consistency check**를 별도 Request 2로 추가했다.
- LLR과 J0337은 각각 **mock/inference**와 **real-data refit**을 분리했다.
- clock-sector는 곧바로 fit하지 않고 **관측식 유도 -> joint fit** 순서로 고쳤다.
- 따라서 이제 각 요청이 "지금 당장 가능한 것"과 "전용 파이프라인이 있어야 가능한 것"으로 명확히 나뉜다.
