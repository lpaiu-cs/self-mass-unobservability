# 열 백색왜성 자료와 scalar 응답 검증

실행 전 체크포인트: `5f99017`.

## 결과 확인 전 선택 규칙

분류: Imported from prior work. Kaplan et al. (2014), https://arxiv.org/abs/1402.0407 의 안쪽 백색왜성 관측치는 Teff=15800±100 K, log g=5.82±0.05 (cgs), R=0.091±0.005 태양반지름이다. 이번 계산의 질량 기준은 동결 IVP의 0.197536385307 태양질량이다. 반지름은 질량과 분광 중력으로부터 유도되므로 log g와 R을 독립 likelihood 항으로 중복 집계하지 않는다.

분류: Conjectural. 공개 진화 모형의 원자료와 내부 구조를 먼저 확보한다. 질량은 목표의 1% 이내를 탐색 후보, Teff와 log g는 각각 발표된 3 sigma 범위 이내를 광학 탐색 후보로 삼는다. 이는 정밀 타이밍 질량을 1% 오차로 간주하거나 확정된 관측 일치를 선언하는 기준이 아니다. 목표 질량에 실제로 맞는 구조와 보간·진화 불확실성이 없으면 실제 EOS matching 완료 판정은 false다. 후보가 없으면 불일치 값을 그대로 보고한다. 모형 선택에 scalar 응답 결과를 사용하지 않는다.

분류: Conjectural. 내부 밀도·압력·조성·온도 분포가 있는 후보에 한해 영 배경의 선형 scalar 응답을 계산한다. 정적 Newtonian 진화 구조를 임의 barotropic EOS 또는 비영 배경의 완전한 relativistic 별로 취급하지 않는다. 공개 구조가 없으면 질량·반지름·온도만으로 내부 응답을 유일하게 결정할 수 있는지 수학적으로 점검하며, 표면 일치를 내부 구조의 검증으로 대체하지 않는다.

분류: Conjectural. Request 17의 코드·원자료·판정은 동결한다. 새 결과는 Request 18로 분리하고, 이전 다섯 유지 문서와 원고 manifest의 현재 바이트를 보존한 뒤 새 절을 덧붙인다. 전체 항성 구간 인증·물리 timing·비선형 추론을 부분 자료 확보로 완료 처리하지 않는다.

## 결과

분류: Proven. 이번 단계는 **열 항성 원자료 확보 및 조건부 정리 진전**이다. 발표된 진화 이력에서 광학 조건에 맞는 시점을 찾았고, 그 앞뒤의 실제 내부 구조를 확보해 scalar 응답을 계산했다. 그러나 저장된 내부 구조 자체는 사전 광학 기준을 통과하지 않으며, 질량도 타이밍 기준과 일치하지 않는다. 따라서 실제 J0337 EOS matching 완료 판정은 false다.

분류: Imported from prior work. 사용 자료는 [Istrate et al. (2016)](https://arxiv.org/abs/1606.04947)의 [CDS 진화 이력](https://cdsarc.cds.unistra.fr/ftp/J/A+A/595/A35/)과 저자들이 공개한 [MESA 7624 단일 실행](https://zenodo.org/record/2634020)이다. 단일 실행의 초기 조건은 donor 1.4, accretor 1.2 태양질량, 주기 3.4일, Z=0.02, 회전과 원소 확산이다. `inlist1`의 초기 diffusion=false만 보면 안 된다. 공개 `run_star_extras.f`는 이후 확산을 켠다. 새 MESA 진화를 실행한 것이 아니라 공개 실행 결과를 읽었다.

분류: Proven. CDS 목록 266개 중 파일명 질량의 반올림 여유를 포함한 사전 필터로 4개 이력을 취득했다. 이는 모든 질량 손실 이력을 망라한 전역 탐색이 아니다. 단일 공개 실행에서 가장 가까운 이력 행은 model 18969로, Teff=15786.15 K, log g=5.82751, 반지름=0.0898423 자료 태양반지름이다. 이 행의 질량은 0.197961063 자료 태양질량이다. 해당 행의 값은 Zenodo history와 CDS에서 일치했다. 사용한 잔차 제곱합은 후보 순위를 정하는 진단량이며 공분산·모형오차를 포함한 posterior가 아니다.

분류: Proven. 공개 `profiles.index`와 이력을 대조하면 이 단일 실행의 저장된 내부 구조 중 질량 1%와 Teff·log g의 각 3 sigma 기준을 함께 통과하는 것은 **0개**다. 최적 이력 행의 바로 전후 저장 구조는 profile 387/model 18950과 profile 388/model 19000이다. 두 구조는 원자료 CRC32 및 원문 SHA256을 검사하고 압축 보존했다. 아래 두 끝점이 중간 진화 구조를 수학적으로 포괄한다고 주장하지 않는다.

| 값 | profile 387 | profile 388 |
|---|---:|---:|
| 분류 | Counterexample candidate | Counterexample candidate |
| Teff (K) | 16409.43 | 14576.54 |
| 반지름 (m, 단위 보정) | 61581483.96 | 64452328.63 |
| 중심 온도 (K) | 2.08608e7 | 2.08006e7 |
| 중심 밀도 (kg/m³) | 1.34134e8 | 1.33517e8 |
| 수소 질량 (자료 태양질량) | 0.00191544 | 0.00191456 |
| 셀 수 | 3088 | 3095 |

## 단위 및 재구성 경계

분류: Proven. `v_rot/omega/radius`로 반지름 단위를, 셀 밀도와 체적·질량 차이로 질량 단위를 복원했다. 회전속도를 km/s로 해석하면 반지름 단위는 695980000 m, 질량 단위는 1.9892e30 kg이다. 두 구조 모두 반지름 상대 잔차 <1e−14, 질량 비율이 전체의 1e−7 이상인 2430/2437개 셀에서 질량 단위 상대 잔차 <2e−9이다. [MESA 출력 열 정의](https://github.com/MESAHub/mesa/blob/main/star/defaults/profile_columns.list)는 속도·각속도·셀 경계의 의미를 확인하는 데 사용했다. 이는 현재 정의이며 옛 상수 소스 파일을 확보했다는 뜻이 아니다.

분류: Counterexample candidate. 이후 계산은 복원한 단위와 명시한 G=6.67430e−11 SI, c=299792458 m/s를 채택했다. GM_sun=1.3271244e20 SI 기준 질량으로 환산하면 0.1980397264 태양질량, 타이밍 기준보다 **0.254809% 높다**. 이는 Newtonian 물질 질량의 환산이며 GR ADM 질량을 별도로 계산한 결과가 아니다. 원자료 숫자를 현재 태양 단위로 곧바로 읽은 초기 계산도 `response-nominal.json`에 남겼다.

분류: Counterexample candidate. 각 셀의 누적 질량과 바깥 반지름을 보존하는 균일 밀도 구각을 정의했다. 내부 중심에는 r=m=0을 추가했다. 압력·열·조성 자료는 물리적 범위 확인에 사용하되, scalar 방정식의 소스에는 정지 질량 밀도만 넣었다. 즉 영 배경에서 고정된 구대칭 밀도 위의 선형·질량 없는 scalar 모형이며, 완전한 열·회전 GR 항성 모형이 아니다. 최대 p/(rho c²)는 약 9.12e−6, 최대 회전/임계 각속도는 약 0.01274다. 이러한 생략 오차는 아래 극소 산술 구간에 들어 있지 않다.

## 고정된 밀도 모형의 응답과 정리

분류: Proven. beta=−4, k=omega/c, rho>=0의 고정 구대칭 구조에서 선형 방정식과 경계는 다음과 같다. 입사장은 l=0 정규파 j0(kr)이고, 시간 규약은 exp(−i omega t)다.

```text
u = r phi,
u'' + [k² + 4 pi |beta| G rho/c²] u = 0,
phi_out = j0(kr) + f(k) exp(ikr)/r,
chi = f(0).
```

분류: Proven. 각 구각에서 위 방정식은 상수 계수이므로 sin/cos 전달행렬을 정확히 곱할 수 있다. U=u(R)/R, V=u'(R), K=kR라 두면 접합 결과는 다음과 같다. 수치 코드는 규격화 u'(0)=1을 쓰며, 이 규격화는 최종 응답에서 소거된다.

```text
chi = R (U−V)/V,
f(k) = R exp(−iK) [U cos K − V sin(K)/K] / [V−iK U].
```

분류: Proven. g_s=|beta|G/c², Q0=g_s M, eta=4 pi g_s integral_0^R rho(r) r dr 로 두자. 유한 지지집합 내부의 적분 연산자는

```text
(T_k h)(x) = g_s integral rho(y) exp(i k |x−y|) h(y)/|x−y| d³y.
```

분류: Proven. 구대칭·비음수 밀도의 Newtonian 퍼텐셜은 중심에서 최대이므로, 실수 k 및 Im k>=0에서 ||T_k||_infinity<=eta다. eta<1이면 Neumann 급수가 수렴해 정적 해가 유일하고 상반평면의 독립 성장 scalar 해도 없다. k=0에서 양성으로 1<=phi<=1/(1−eta), 따라서 Q0<=chi<=Q0/(1−eta)다. 이 증명은 독립적인 유체·회전 모드 또는 하반평면의 모든 pole 부재를 주장하지 않는다.

분류: Proven. 실수 k에 대해 |exp(ikd)−1|<=|k|d 이므로 ||T_k−T_0||<=|k|Q0이고, |j0(kr)−1|<=(kR)²/6이다. Resolvent 항등식으로

```text
||phi_k−phi_0|| <= (kR)²/[6(1−eta)] + |k|Q0/(1−eta)².
```

분류: Proven. 바깥 l=0 산란 진폭은 f(k)=g_s integral rho(r) j0(kr) phi_k(r) d³r이다. 여기에 위 부등식과 chi>=Q0를 적용하면 연속 주파수에 대해 다음 조건부 상계를 얻는다.

```text
|f(k)−chi| / chi <= (kR)²/[3(1−eta)] + |k|Q0/(1−eta)².
```

분류: Proven. 두 보정된 구각 모형에서 eta<0.000166이고 chi는 각각 약 1169.85952 m, 1169.85893 m다. 고정된 이산 입력을 정확한 binary 수로 취급한 구간 연산으로 chi를 감쌌으며, 두 산술 구간 폭은 1.2e−59 m 미만이다. 이는 항성 모형의 물리적 정확도가 59자리라는 뜻이 아니다. 2/4개 셀씩 합친 재구성과의 상대 차이는 5e−10 미만이며, 이것도 원래 진화 모형의 수렴 증명은 아니다.

분류: Proven. 선언한 안쪽 궤도 각주파수 약 4.46311e−5 /s에서는 두 모형 모두 |f−chi|/chi<**2.1e−10**으로 위에서 보장된다. 직접 전달행렬 계산값은 약 1.742e−10이다. 따라서 이 고정 WD scalar 모형은 해당 주파수에서 크기 1 수준의 비단열 증폭이나 지연을 제공하지 않는다. 이것을 실제 삼중계의 타이밍 오차 상계나 중성자별·유체까지 포함한 no-go로 확대하지 않는다.

## 검증·재현 및 남은 작업

분류: Proven. 전달행렬 행렬식·미분방정식의 기호 검증, 균일 구의 독립 닫힌 해 및 DOP853 ODE, 진공 응답 0, 주파수 켤레 대칭, 탄성 unitarity Im f=k|f|² 검사를 통과했다. 원자료 압축파일은 필요한 7z 블록만 가져왔으며 선택한 파일의 CRC와 원문 SHA를 검사했다. 전체 1.4 GB archive의 MD5를 확인했다고 주장하지 않는다. 부분 archive의 파일 길이는 원본 길이와 같아도 전체 원본이 아니다.

실행:

```text
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'cd /mnt/e/lab/self-mass-unobservability && PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/thermal_wd.py response'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'cd /mnt/e/lab/self-mass-unobservability && PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/thermal_wd.py check'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'cd /mnt/e/lab/self-mass-unobservability && PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/thermal_wd.py verify'
```

분류: Conjectural. 다음 실제 matching 단계는 광학 조건에 해당하는 세밀한 시점의 열 구조를 저장하는 진화 재실행과 정확한 질량·외피·조성의 보정이다. 그 뒤에야 비영 배경의 열·유체·metric 응답, 두 동반성의 물리적 구동과 광자 관측식, 전 기간 미분 인증 및 전체 비선형 likelihood를 연결할 수 있다. 두 끝점 사이의 임의 보간이나 질량의 단순 재규격화를 실제 EOS 확정으로 바꾸지 않았다. 원고 PDF/ZIP은 Request 12 동결본이며 이번 결과가 편집 완료된 제출본에 반영되었다고 주장하지 않는다.
