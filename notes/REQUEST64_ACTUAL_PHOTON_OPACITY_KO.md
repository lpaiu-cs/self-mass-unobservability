# Phase64 — 실제 온도 흡수의 직접 계산과 좁은 선 적분 병목 해결

분류: Counterexample candidate. **좁은 원자선을 격자점 사이의 직선으로 이어 생긴 결합 응답의 수렴 실패를 해결했다.** 실제 외곽 온도에서 독립 원자 모형의 흡수를 직접 계산하고, 선 모양을 기존 광자 셀에 적분해 같은 물질·집단 산란·공간 시간식에 연결했다. 수정 결합의 시간·격자·선 구적·보존 대조는 통과했다. **같은 EOS의 실제 흡수 인증은 아직 아니다.** 원자 점유수와 미시 상태의 일치, 누락된 원자 성분과 결합 전자 재분배가 남는다.

분류: Imported from prior work. 직전 단계의 같은 EOS 연속 RPA 집단 산란과 저장된 TOPS 흡수 응답을 재사용했다. 기존 약 0.421초 항성 이력은 재적분하지 않았다. 시작 체크포인트는 `b68e0a10f7e2fa38a205e93ec177c72a91700dce`다. 이번 진전은 **loophole progress 및 선언 모형의 conditional theorem progress**이며 전체 동적 전하 목표는 미완료다.

## 직접 온도 입력과 검증 범위

분류: Counterexample candidate. 입력은 저장 외곽 상태의 T=18827.94717983866K, rho=1.91629314435e−9 g/cm³, ne=1.0445347529258283e15 cm⁻³다. native 단일 모형 입력으로 세 값을 동시에 지정했다. 총 원소비도 기존 EOS 재고에서 가져왔다. 광자는 기존 24,134셀(u≤60), 물질 열용량은 기존 Cm=0.5832057872930577, 초기 물질 섭동 1K·광자 섭동 0, kH=1, 시간 H/c=2.316801ms를 유지했다.

분류: Imported from prior work. TOPS는 선택 가능한 온도 목록을 사용하므로 실제 온도 밖의 직접 조회를 반복하지 않았다. 기존 양 끝 온도 응답 차이는 온도 보간의 엄밀한 오차 상계가 아니다. [TOPS 입력 문서](https://aphysics2.lanl.gov/static/opacdocs/opac-help.html)

분류: Imported from prior work. 직접 계산에는 공개 SYNSPEC 소스를 커밋 `b9149f7208eeca9b4fdd38dd11d9f736c7a050d7`에 고정했다. 원자 모형·선 목록·선폭 자료는 별도 입력이다. 코드의 MIT 라이선스를 원자 자료 전체의 라이선스로 확대하지 않는다. [고정 소스](https://github.com/callendeprieto/synple/tree/b9149f7208eeca9b4fdd38dd11d9f736c7a050d7), [공식 원자 입력 설명](https://tlusty.oca.eu/tlusty/Synspec49/synspec-data.html), [공식 선 목록 설명](https://tlusty.oca.eu/tlusty/Synspec49/synspec-line.html)

분류: Counterexample candidate. 공개 분포의 2,309,497개 선에서 실제 활성 원소 H/He/C/N/O/F/Ne/Mg/Ca의 301,366개를 골랐다. 원 숫자 필드가 붙은 48,543행은 숫자 토큰으로 분리해 native reader가 조용히 건너뛰지 않도록 했다. 최종 선택 파일 SHA256은 `339d463da98221b38695334184c7233acf09c94e68cb9f4ed7e3e44df7b3f33a`다. 격자와 무관하게 추출한 통상 LTE 프로파일은 60,140개다. 같은 중심의 서로 다른 감쇠값과 한 배치 안의 중복 개수는 보존하고, 배치 간 반복 출력만 제거했다.

## 원 실패와 수정

분류: Counterexample candidate. 외부 코드의 단일 상태 실행에서 드러난 밀도 격자 0분모, 고정 배열 범위, F77 가변 길이 배열 선언, 첫/마지막 선과 불균일 주파수 검색 오류를 수정했다. `NATOMS=99`와 비활성 원소 mode=0으로 암묵적인 태양 조성 추가를 막았다. 범위 검사 빌드와 실패 로그를 보존했다. native `mode=-4`도 H/He 선을 포함하므로 이를 완전한 연속 흡수 전용 모드라고 해석하지 않는다.

분류: Proven. 같은 LTE 전이의 양의 비유도 흡수 U_nu에 대해 chi_nu=U_nu(1−exp(−h nu/kT)), j_nu=U_nu(2h nu³/c²)exp(−h nu/kT)이면 j_nu=chi_nu B_nu다. 배치 중심의 유도 방출·Planck 인자를 다른 주파수에 그대로 쓰면 이 항등식을 보존하지 않는다. 주파수별 인자로 바꾸고 연속 흡수도 각 주파수에서 직접 평가했다. 이 항등식은 같은 미시 모형의 LTE를 가정하며 실제 원자 자료의 정확성을 증명하지 않는다.

분류: Counterexample candidate. 첫 직접 흡수 결합은 시간·보존 대조를 통과했지만 **격자 차이 0.01905324320697528 > 0.001**로 실패했다. `response.json`의 `numerical_response_passed=false`를 그대로 보존한다. 이때 새 모형과 저장 TOPS 모형의 차이 0.0367464도 수락된 물리 결과로 채택하지 않는다.

분류: Counterexample candidate. 두 원인을 먼저 수정했다. 상대 원자량에 곱하던 H 원자 질량을 EOS와 같은 1/Avogadro로 바꾸고 핵 질량을 맞췄다. 선의 선택 구간을 균일 격자 인덱스가 아닌 실제 주파수의 이진 검색으로 바꿨다. 전체 전하 재고 차이는 −0.8057%에서 −0.02351%로 줄었지만, 이것만으로 선 적분이나 종별 점유수는 일치하지 않았다. 이 수정만으로는 새 결합 계산을 반복하지 않았다.

분류: Counterexample candidate. 남은 수렴 실패는 **미분해 선의 직선 보간**이었다. C II 1334.532Å 부근에서 두 격자의 중심 흡수는 약 2.39e−5 cm⁻¹로 비슷하지만, 인접 배경점의 간격이 달라 큰 삼각형 면적을 만들었다. u≈5.72537인 같은 광자 셀의 흡수율은 101930.15/s와 180.31/s로 갈라졌다. 질량/주파수 수정 뒤 고정 물질 온도의 독립 감쇠 사전 시험도 0.0235372 차이였다.

분류: Counterexample candidate. 해상도를 전역 확대하는 대신 native 선의 중심·세기·Doppler 폭·감쇠값을 추출했다. 폭은 배치 중심이 아니라 각 선 주파수로 계산한다. 통상 LTE 선은 배경에서 분리하고, 양의 Voigt 모양을 **기존 셀 경계에서 잘라** asinh 좌표의 Gauss4/8로 적분했다. 특수 H/He 프로파일은 native 배경에 남아 별도 오차 경계가 필요하다. native의 흡수율 의존 날개 절단 대신 전체 Voigt 날개를 목표 모형으로 선언하고, 생략한 먼 날개에는 아래 상계를 적용했다. 이 모형 변경을 단순 코드 동등성으로 표현하지 않는다.

## 선 날개의 조건부 응답 상계

분류: Proven. Voigt 함수의 정확한 적분 표현을 H(a,x)=∫exp(−s²)a/[π((x−s)²+a²)] ds로 둔다. a>0, |x|≥R에서 적분을 |s|≤R/2와 바깥으로 나누면

```text
H(a,x) ≤ a / [sqrt(pi)((R/2)^2+a^2)]
         + erfc(R/2) / [sqrt(pi)*a].
```

분류: Proven. 각 선의 비유도 진폭을 곱해 합친 균일 흡수 누락량을 epsilon이라 하자. 양의 셀 적분과 LTE 교환의 대칭 에너지 좌표에서 흡수 연산자 차이는 ||Delta A|| ≤ c epsilon (1+Cgamma/Cm)을 만족한다. 공통 이동 항이 반에르미트이고 충돌이 소산적이면 Duhamel 식으로 ||Delta e(t)||/||e(0)|| ≤ t c epsilon (1+Cgamma/Cm)을 얻는다. 이 진술은 **선택한 유한 셀·Voigt 모형**에 조건부다.

분류: Counterexample candidate. 이번 날개 배분의 상계 식 평가값은 9.976938817e−5다. 수치 구적 오차는 별도 비교했다. 부동소수점 계산을 바깥 반올림한 구간 인증이라고 부르지 않는다. 기호식 검사, 189개 날개 부등식 표본 및 독립 적응 구적의 Voigt 정규화 검사를 통과했다. 마지막 정규화 상대 오차는 최대 1.12e−15였다.

## 최종 결합 판정

분류: Counterexample candidate. 선 적분 뒤 두 배경 격자의 고정 열원 응답 차이는 1.873616134e−5, 선 구적 차이는 7.749856094e−7로 사전 시험을 통과했다. 앞 C II 셀은 두 격자에서 177.127013/s와 177.127016/s가 됐다. 그 다음에만 결합 계산을 실행했다. 아래는 초기 1K의 에너지 노름으로 정규화한 같은 2.316801ms 끝점 비교다.

| 대조 | 결과 | 고정 기준 |
|---|---:|---:|
| 시간 16/32/64 차수 | 2.013557071 | ≥1.8 |
| 마지막 시간 차이 | 1.189233095e−5 | <1e−3 |
| 기존 두 배경 격자 | 2.128295844e−5 | 아래 합계 기준 |
| 배경 격자 차이 + 두 날개 상계 | 2.208217348e−4 | <1e−3 |
| 선 구적 Gauss4/8 | 6.241957556e−7 | <1e−4 |
| 최대 에너지 수지 결함 | 2.392279213e−12 | <1e−9 |
| 최대 원 에너지식 잔차 | 2.820472948e−12 | <1e−9 |
| 최대 선형 방정식 잔차 | 6.626637055e−13 | <1e−11 |
| 최대 산란 에너지·광자 수 영모드 잔차 | 4.334723193e−16 | <1e−10 |
| 에너지 노름 증가 | 0 | <1e−10 |

분류: Counterexample candidate. `profile-integral/response.json`의 `numerical_response_passed=true`다. 최종 물질 온도 섭동은 0.803945291336K이고 저장 TOPS 모형보다 0.001372413740K 크다. 두 모형의 전체 응답 차이는 **0.0348483781970**다. 이는 온도·원자 자료·점유수·선 표현이 다른 후보들의 차이이며 **TOPS 온도 보간 오차나 실제 관측 신호가 아니다**. 같은 배경 두 격자가 공유하는 특수 선/연속 성분의 공통 편향도 이 대조로 제거되지 않는다.

## 남은 물리 병목과 진행 순서

분류: Counterexample candidate. 실제 T/rho/ne와 핵 질량을 맞춘 뒤에도 native H 중성 점유수는 EOS보다 15.06%, N 이중 이온은 −14.78%, O 이중 이온은 −21.79% 다르다. 최대 종별 상대 차이는 희박한 중성 Ca에서 약 154.2이며 이를 지배적 흡수 불확실성이라고 해석하지 않는다. native LTE 방출 잔차도 최대 1.41450e−5로 기존 1e−6 기준에 미달한다. 결합 후보에서 j=chi B를 선언했다고 native 실패가 사라지는 것은 아니다.

분류: Imported from prior work. 현재 분자 EOS는 해당 옵션에서 H/He/분자의 여기·점유 처리를 포함하지만 여기된 금속을 같은 방식으로 포함하지 않는다. SYNSPEC의 독립 원자 준위·분배 함수·임계값과 동일하지 않다. 같은 총 원소비나 작은 총 전하 잔차는 이 차이를 대신 검증하지 못한다.

분류: Conjectural. 다음 우선순위는 **같은 EOS의 종·준위 점유와 원자 흡수/방출을 하나의 미시 규약으로 연결하는 것**이다. 단순히 선 세기에 보정 계수를 곱해 맞추지 않는다. 종별/준위별 불일치를 흡수 성분으로 추적하고, EOS에서 실제 제공하는 점유·점유확률·에너지 기준을 연결하거나 필요한 EOS 미시 확장을 명시해야 한다. F 연속 성분, Mg I/Mg III 및 Ca III 연속 성분과 높은 이온 단계의 누락도 닫아야 한다. 그 뒤 결합 전자·Rayleigh 재분배, 비균일 반경/대기, 유체 운동량·전체 GR 되먹임으로 진행한다.

분류: Conjectural. 기존 전도 속도 시간 수렴, 전하 공간 수렴, 실제 궤도 구동·정적 비교·관측 추론 및 전체 동적 전하 완료 요구는 유지한다. 원 셀 안의 좁은 선을 평균한 흡수는 셀 내 광자 수송의 완전한 해가 아니므로 이 오차도 남는다.

## 자원과 재현

분류: Counterexample candidate. 첫 실패 결합은 155.573초·1.580GB였다. 질량/주파수 수정만으로 실패가 남은 동안 이를 반복하지 않았다. 최종 배경 스펙트럼은 기존 2Å/0.5Å 및 IR 간격을 그대로 사용해 6.986/15.858/1.815초, 선 구적 두 경로는 7.194/10.382초였다. 수정 결합 다섯 경로는 한 산란 행렬을 재사용해 **198.892초·1.580GB**, 사전 480초·4GB 상한 안에서 끝났다. 작업자 CPU 1개, 새 TOPS 조회·새 항성 GR 시간 단계는 0이다. 초기 소규모 선 구적 시험의 배분 비용은 전체 선 수에 따라 달라지므로 단순 선 수 비례 추정을 완료 예상으로 사용하지 않는다.

분류: Imported from prior work. 실행 파일은 `verification/def_photon_actual_opacity.py`, `def_photon_actual_opacity_validate.py`, `def_photon_atomic_repair.py`, `def_photon_profile_integral.py`이며 독립 식 검사는 `def_photon_profile_check.py`다. 모든 결과는 `outputs/direct-eos-gr33/def-photon-actual-opacity/`에 둔다. 원 실패·패치·계획·실행 영수증·단계별 SHA 결속을 보존했다. 개발 중 바뀐 profile 스크립트는 각 계획에 일치하는 원문 스냅샷을 별도 보존한다.

분류: Imported from prior work. 큰 원문 스펙트럼은 원 SHA와 압축 SHA를 확인한 gzip으로 보관한다. `inputs.tar.gz`에는 고정 원 소스·라이선스·원자 입력·정규화된 선택 선 목록이 있어 새 네트워크 조회 없이 재현할 수 있다. 원문 복원과 식 검사는 다음과 같다. 생산 스크립트는 기존 판정 덮어쓰기를 거부한다. 새 생산 재현은 별도 출력/캐시 위치에서 실행하고 기존 수락/실패 파일을 지우지 않는다.

```text
rtk proxy wsl -d Ubuntu-22.04 --cd /mnt/e/lab/self-mass-unobservability -- /usr/bin/env OPENBLAS_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification /usr/bin/python3 verification/def_photon_opacity_artifacts.py unpack
rtk proxy wsl -d Ubuntu-22.04 --cd /mnt/e/lab/self-mass-unobservability -- /usr/bin/env OPENBLAS_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification /usr/bin/python3 verification/def_photon_profile_check.py
```
