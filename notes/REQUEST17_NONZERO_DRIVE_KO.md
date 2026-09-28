# 비영 scalar 가지와 동반성 되먹임 검증

실행 전 체크포인트: `043982c`.

분류: Counterexample candidate. DEF coupling `A=exp(−2 phi²)`와 기존 SLy 중성자별을 유지하고, 비영 배경 `phi_infinity=1e−5`를 지정한다. 이는 관측으로 결정한 값이 아니다. 각 별의 영 배경 중력질량으로 기준 별을 찾고 이후에는 바리온 수를 고정해 비영 배경 가지를 계산한다. 두 백색왜성은 `mu_e=2`, 온도 0, 이상 전자 축퇴 압력과 이온 정지 질량을 가진 명시적인 EOS를 사용한다. Coulomb·온도·조성층 보정은 포함하지 않으므로 실제 두 별의 EOS 확정이 아니다.

분류: Conjectural. 0, ±1e−5, 3e−5, 1e−4의 배경에서 전하의 홀수 대칭·작은 장 수렴과 바리온 수 보존을 검사한다. 독립 적분 정밀도·EOS 표본 밀도·중심 반지름 검사를 실시한다. 외부 진공은 `u=R/r`로 압축해 무한대까지 적분한다. 그 뒤 세 천체의 상호 scalar 퍼텐셜을 포함한 정적 선형 응답 행렬을 풀고, 모든 쌍의 전하 변동·상호 힘·되먹임 크기를 계산한다.

분류: Conjectural. 비영 정적 구동의 성공을 궤도 시간 규모의 완화 또는 완전한 타이밍 모형으로 승격하지 않는다. 주파수 응답·복사·관성은 각각 유효 범위와 계수 matching이 필요하다. 이전 실행과 원고의 판정은 동결하고 신규 자료에 별도 결속한다.

## 결과와 판정

분류: Proven. 지정 후보의 비영 배경 평형, 선도 동반성 상호 구동, 정적 축약의 조건부 일·이차 미분 구간을 계산했다. 관성을 생략한 선도 monopole 복사 모형에서는 집단 완화식과 모든 시각의 조건부 추적 오차 상계를 도출했다. 연구 분류는 **지정 후보 계산 및 정리 진전**이다.

분류: Proven. 반면 차가운 안쪽 백색왜성 모형은 발표된 광학 반지름을 재현하지 못한다. 정적 구동의 존재만으로 새 느린 상태가 성립하지도 않는다. 따라서 실제 J0337의 완전한 물리 matching·전체 타이밍 인증·비선형 추론 완료 판정은 유지하지 않고 모두 false로 남겼다.

| 작업 | 분류 | 판정 |
|---|---|---|
| 지정 SLy 및 두 차가운 WD의 유한 배경 평형 | Counterexample candidate | 수치 계산·정밀도 비교 완료 |
| 모든 전하의 선도 상호 응답 | Counterexample candidate | 비영 구동 계산; 고정 동반성 근사의 실패 확인 |
| 명시한 구간의 결합 전하 일·이차 미분 | Proven | 정적 축약 모형에 조건부 인증 |
| 관성을 생략한 선도 복사 모형의 추적 오차 | Proven | 초기 평형과 선언한 궤도·감수율 구간에 조건부 상계 |
| 차가운 WD의 광학 반지름 일치 | Proven | 실패; 실제 별의 열·외피 구조 필요 |
| 느린 실제 내부 상태·완전한 관측 추론 | Conjectural | 미완료 |

## 1. 유한 배경의 별 세 개

분류: Counterexample candidate. 중성자별의 물질 EOS는 Request 13에서 동결한 LAL SLy 압력–에너지밀도 함수를 사용했다. WD EOS는 `x=p_F/(m_e c)`에 대해 다음 압력과 에너지밀도를 사용한다. `mu_e=2`, `n_e=(m_e c/hbar)^3 x³/(3pi²)`, `P0=m_e^4 c^5/(24pi² hbar³)`이다.

```text
p = P0 [x(2x²−3) sqrt(1+x²) + 3 asinh(x)]
epsilon = mu_e m_u c² n_e
        + 3 P0 [x(1+2x²) sqrt(1+x²) − asinh(x)] − m_e c² n_e
```

분류: Proven. 작은 x에서는 차감 소실을 피하는 급수를 사용하고, 독립 전자 운동량 적분과 비교했다. `dp/dx=8 P0 x⁴/sqrt(1+x²)` 및 열역학 항등식을 기호 검증했다. 물질량 보존에는 `d log n=d epsilon/(epsilon+p)`를 사용했다. NS의 `n=(epsilon+p) exp(−h)`, `h=integral dp/(epsilon+p)`에는 공통 정규화 자유도가 있으므로 저장한 바리온 적분의 수치를 실제 바리온 태양질량으로 해석하지 않는다. 같은 EOS에서 이 적분을 고정하는 것은 같은 바리온 수를 고정하는 것이다.

분류: Imported from prior work. 비영 scalar 항성 구조와 바깥 진공 접합은 [Mendes와 Ortiz(2016), 식 31–38](https://arxiv.org/pdf/1604.04175)의 방정식에 근거한다. 아래에서는 `nu=log N`, `f=1−2m/r`, `psi=dphi/dr`, `alpha=−4phi`이며 물질 `epsilon,p`는 Jordan-frame 값이다.

```text
m' = 4pi r² A⁴ epsilon + (r² f/2) psi²
nu' = (m+4pi r³ A⁴ p)/(r² f) + r psi²/2
h' = −nu' − alpha psi
phi' = psi
psi' = −(2/r + nu' − lambda') psi + 4pi A⁴ alpha(epsilon−3p)/f
lambda' = (m'/r − m/r²)/f
Nb' = 4pi r² A³ n / sqrt(f)
```

분류: Proven. 바깥 진공은 `u=R/r`, `z=r² psi`로 바꾸어 `u=1→0`까지 적분했다. 유한 원거리에서 scalar 값을 곧바로 무한대 경계값으로 대체하지 않았다. 이 적분의 ADM 질량·scalar 전하·무한대 배경은 별도의 해석적 진공 접합식과 비교해 상대차 약 2.22e−16 이하, 배경 절대차 2.79e−17 이하로 일치했다. 이는 진공 계산의 교차 검사이며 내부 항성 해의 엄밀한 오차 상계는 아니다.

| 별 | 분류 | 기준 중력질량/태양질량 | 영 배경 반지름 | 영 배경 감수율 chi | phi∞=1e−5의 전하 q |
|---|---|---:|---:|---:|---:|
| SLy 중성자별 | Counterexample candidate | 1.437814408 | 11.69923 km | 41,566.0195 m | 0.415659948 m |
| 안쪽 차가운 WD | Counterexample candidate | 0.197536385 | 14,725.580 km | 1,166.86869 m | 0.0116686869 m |
| 바깥쪽 차가운 WD | Counterexample candidate | 0.410102707 | 10,837.453 km | 2,422.98693 m | 0.0242298692 m |

분류: Proven. EOS 표본 수 4,096→16,384, 적분 상대 허용오차 2e−9→2e−11, 중심 시작 반지름의 절반 변경을 비교했다. 감수율의 최대 상대 변화는 2.75e−8이었다. 이전 Request 13의 native SLy 결과와 새 엔탈피 경로의 감수율 차이는 2.80e−7 상대값이었다. ±1e−5의 전하는 부호만 반전했으며, 정밀 실행의 바리온 적분 상대 편차는 최대 1.32e−9였다. NS의 `q/phi`는 1e−4에서 영 배경 감수율보다 약 2.10e−5 작았다. 선형 응답과 유한 배경 계산을 동일한 정확도로 간주하지 않는다.

분류: Proven. 압력 좌표의 초기 실행 두 개는 내부 적분 단계가 EOS 범위를 벗어나 실패했다. 첫 shooting 방식의 수치 Jacobian도 작은 바리온 차이에 취약하여 명시적 차분 간격과 제한된 shooting 영역으로 바꿨다. 표면 적분은 엔탈피로 변경했다. 첫 엔탈피 실행은 백색왜성에도 NS용 100 m 최대 간격을 사용해 비효율적이어서 중단하고, 동일 방정식에 항성 크기에 맞는 간격을 적용해 다시 실행했다. 실패 소스·로그와 최종 실행 로그를 모두 보존했다. 물리 모형이나 수용 문턱을 결과에 맞춰 바꾼 것은 아니다.

## 2. 광학 관측과의 별도 대조

분류: Imported from prior work. [Kaplan 등(2014)](https://arxiv.org/abs/1402.0407)은 안쪽 WD의 유효온도 `15,800±100 K`, 반지름 `0.091±0.005 R_sun`을 보고했다. 이 수치는 이번 실행에서 새로 측정한 값이 아니다.

분류: Proven. 차가운 EOS의 안쪽 WD 반지름은 약 `0.02117 R_sun`으로 그 중심값의 약 23.3%다. 따라서 차가운 EOS를 실제 안쪽 동반성의 구조 matching 성공으로 채택할 수 없다. [광학 대조 기록](../outputs/nonzero-drive17/optical-domain-check.json)에 이 실패를 별도 판정했다. 질량만 맞춘 항성 모형을 실제 별의 검증으로 간주하지 않는다.

분류: Conjectural. 실제 계로 이어가려면 관측된 질량·반지름·온도를 만족하는 열 구조와 수소 외피를 포함한 WD 모형에서 감수율과 동적 응답을 다시 계산해야 한다. 이번 차가운 모형은 상호 응답을 점검하는 구체적인 대조 후보로 남는다. 그 약한 중력 극한의 유사성만으로 실제 WD의 전체 동적 오차가 작다고 인증하지 않는다.

## 3. 서로 응답하는 세 전하와 힘

분류: Imported from prior work. scalar 전하의 작용, 정적 감수율 matching과 leading monopole 복사는 [Khalil 등(2022)](https://arxiv.org/html/2206.13233v2)의 식 2–9, 18, 23–24를 출발점으로 삼았다. 그 논문의 특정 별 계수를 이번 SLy·WD 계수와 혼합하지 않았다.

분류: Proven. 작은 배경의 선도차수에서 `D=diag(chi_p,chi_i,chi_o)`, `K_AB=1/r_AB (A≠B)`, 대각 `K_AA=0`이면,

```text
L = D^(-1) − K
L q = phi∞ 1
Delta_AB = q_A q_B / (m_A m_B)
V(q,r) = q^T L q / 2 − phi∞ 1^T q
```

분류: Proven. q를 독립 변수로 변분하면 모든 쌍의 힘이 상호적이다. 전하를 평형에서 제거하면 `V_eff=−phi∞² 1^T L^(-1)1/2`이다. 이는 정적인 위치 의존 퍼텐셜이며, 전하가 궤도를 따라 변한다는 사실만으로 새로운 동적 상태가 되지 않는다. 유한 phi의 비선형 항성 해는 별도 계산했지만, 여기의 궤도 격자는 그 해를 정확히 대입한 완전 비선형 다체 모형이 아니라 phi의 선도차수 모형이다.

분류: Counterexample candidate. 동결된 GR 상태의 한 시점에서 얻은 Newtonian 접촉 궤도로 안팎의 타원을 정의하고, 두 평균 근점이각을 독립적으로 변화시켰다. 네 번째 천체, 1PN 전파와 궤도의 scalar 반작용은 이 세 천체 격자에 포함하지 않았다. 비교 격자는 64²와 128²이며 실제 관측 날짜들의 likelihood 계산은 수행하지 않았다.

분류: Proven. 전체 전하에 대한 정적 되먹임의 최대 상대 변화는 NS 2.59e−7, 안쪽 WD 8.72e−6, 바깥쪽 WD 2.52e−7이었다. 그래도 작은 주기 변동에서 WD의 변화는 생략할 수 없다. 아래는 phi∞=1e−5에서 힘의 무차원 결합 `Delta`의 Fourier 진폭이다. TOA 시간 지연이나 측정 상한이 아니다.

| carrier와 쌍 | 분류 | 모든 전하 응답 | WD 전하 고정 | 비율 |
|---|---|---:|---:|---:|
| 내궤도, NS–안쪽 WD | Proven | 4.88372e−17 | 1.33358e−18 | 36.62 |
| 외궤도, NS–바깥쪽 WD | Proven | 7.07282e−17 | 3.79519e−18 | 18.64 |
| 차주파수, NS–안쪽 WD | Proven | 2.20607e−18 | 3.51271e−19 | 6.28 |
| 차주파수, NS–바깥쪽 WD | Proven | 5.14669e−18 | 3.51339e−19 | 14.65 |

분류: Proven. 격자 정밀화에 따른 이 세 carrier·모든 쌍의 최대 상대 변화는 약 6.49e−8이었다. 두 NS 쌍의 신호도 서로 다르므로 기존의 공통 결합 여섯 열을 이 모형의 관측 템플릿으로 그대로 사용할 수 없다. 이 결론은 전체 전하 변화의 작은 비율만 보고 구동 오차를 판단할 수 없다는 구체적인 반례다.

## 4. 조건부 일·이차 미분 인증

분류: Proven. 세 역거리를 `a,b,c`라 쓰면,

```text
det(I−D K) = 1−chi_p chi_i a²−chi_p chi_o b²−chi_i chi_o c²
              −2 chi_p chi_i chi_o a b c
```

분류: Proven. 감수율과 거리를 명시한 유리수 끝점의 구간으로 정하고 바깥 반올림을 사용했다. 그 전 구간에서 행렬식은 약 `[0.9999999999978668,0.9999999999978755]`, `||D K||_infinity<8.955e−6`였다. Neumann 급수로 역행렬이 존재한다. `L`은 대칭이며 `I−sqrt(D) K sqrt(D)`와 합동이므로 이 구간에서 양의 정부호다. 이는 축약된 정적 전하 에너지의 성질이며 전체 항성 PDE의 안정성 증명은 아니다.

분류: Proven. 닫힌 형태의 q를 기호 미분하고 세 역거리에 대한 Jacobian과 Hessian 전체를 구간으로 평가했다. `q_theta=L^(-1) K_theta q`, `q_theta_eta=L^(-1)(K_theta q_eta+K_eta q_theta+K_theta_eta q)`로도 해석할 수 있다. [결합 인증 파일](../outputs/nonzero-drive17/coupled-certificate.json)은 정확한 입력 구간과 도함수 범위를 저장한다.

분류: Conjectural. 이 입력 구간은 **조건부 정리의 정의역**이다. 수치 항성의 참 감수율이나 실제 GR 궤도가 항상 그 안에 있다는 구간 인증은 아직 없다. 따라서 이번 Jacobian/Hessian 인증이 이전의 28개 타이밍 매개변수 초기화·전 기간 흐름·광자 지연 인증을 대신하지 않는다.

## 5. 집단 복사 축약의 빠른 완화와 연속 추적 상계

분류: Proven. 관성을 생략하고 선도 monopole 복사 항만 유지하면 감쇠 행렬은 `1 1^T/c`로 rank one이다. 이를 임의의 양의 대각 감쇠 세 개로 바꾸지 않는다. `S=1^T q`, `Ceff=1^T L^(-1)1`로 놓으면,

```text
L q + 1 Sdot/c = phi∞ 1
tau(t) Sdot + S = Ceff(t) phi∞,   tau(t)=Ceff(t)/c
```

분류: Proven. 앞서 선언한 감수율·거리 구간에서는 `Ceff∈[45153.02,45158.03] m`, `tau∈[0.150614,0.150631] ms`이다. 이 축약 모형의 집단 응답은 일 단위의 느린 상태를 만들지 않는다. 이는 관성·추가 모드·고차 복사가 무시 가능하다는 별도 조건 아래의 결론이다. 실제 항성에 그 조건이 성립한다고 인증한 것은 아니다.

분류: Proven. `e=S−Ceff phi∞`라 두면 `edot+e/tau=−(Ceff phi∞)dot`다. 적분인자를 적용하면 모든 시각에서

```text
|e(t)| ≤ |e(0)| exp(−integral_0^t ds/tau(s))
       + tau_max sup_s |(Ceff phi∞)dot(s)|
```

분류: Proven. 선언한 속도·거리 경계와 구간 Jacobian으로 평형 총전하의 시간 변화율을 `9.091e−12 m/s` 이하로 상계했다. 평형 초기자료이면 각 전하의 추적 오차는 `1.370e−15 m` 이하이고, 세 쌍의 `Delta` 오차 상계는 각각 `9.449e−22`, `4.686e−22`, `2.784e−22`다. 양의 성분을 갖는 `L^(-1)1`에 대해 `q_A=(L^(-1)1)_A S/Ceff`이므로 총전하 오차 상계가 각 전하에도 적용된다. 이것은 유한 시각 표본 비교가 아니라 **선언한 축약 모형과 영역에 대한 연속 상계**다.

분류: Conjectural. 항성의 모드 관성, 높은 다중극 복사, 실제 움직이는 원천의 지연장, 실제 궤도 오차는 이 수치에 포함하지 않았다. 특히 단순히 `1/r(t)`를 `1/r(t−r/c)`로 치환하면 움직이는 원천의 다른 항을 빠뜨릴 수 있으므로 그러한 위상 지연을 물리 신호로 추가하지 않았다.

## 6. 관측 추론의 경계와 남은 순서

분류: Proven. 짝 coupling `A(phi)=A(−phi)`에서 scalar 자료와 전하를 함께 반전하면 물질 타이밍 관측량은 보존된다. 작은 장의 선도 신호는 `phi∞²`에 비례하며, 매끄러운 영 가지에서 phi∞의 score와 Fisher 정보는 영점에서 퇴화한다. 이전의 부호 있는 현상론적 beta에 대한 Gaussian 구간을 phi∞ 구간으로 직접 읽을 수 없다. 부호 퇴화, `lambda=phi∞²≥0`의 경계와 실제 nuisance를 포함한 추론이 필요하다.

분류: Conjectural. 다음 물리 단계는 광학 자료를 만족하는 WD의 열·외피 구조와 결합 모드·힘·광자 전파 도출이다. 그 결과가 정해져야 새 세 쌍 신호를 넣은 실제 비선형 timing/noise·pulse 추론을 실행할 수 있다. 동시에 전체 초기화와 전 기간 변분·보간·역시간 오차 인증이 필요하다. 잘못된 WD 반지름과 미도출된 동적 신호를 기존 데이터에 적합시켜 완료 판정을 만들지 않는다.

## 재현과 보존

WSL Ubuntu-22.04의 기존 Request 13 의존성을 사용한다. 외부 계정·새 서비스는 필요하지 않다.

```powershell
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/nonzero_drive_audit.py check"
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/nonzero_drive_audit.py verify"
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && python3 verification/verify_unified_paper.py"
```

분류: Proven. 원시 항성 결과, 두 궤도 격자 요약, 정밀 격자, 소스·실패 로그·검증 출력은 [`outputs/nonzero-drive17`](../outputs/nonzero-drive17/manifest.json)에 결속한다. Request 16 시점의 문서와 원고 manifest는 별도 바이트 보존본으로 검증한다. 원고 PDF와 ZIP은 Request 12 동결 산출물이며, 이번 계산을 포함한 새 제출본이 아니다.
