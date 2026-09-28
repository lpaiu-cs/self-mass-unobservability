# Phase43 — 비영 배경의 동적 외부장

## 전체 목표에서의 위치

분류: Counterexample candidate. 정적 전하 투영을 넘어 **시간 의존 scalar·계량 외부장의 나가는 파 경계 연산자**를 구현했다. 단계 41의 비영 배경 표면 자료, 선언 궤도의 첫 세 조화 주파수, 강한 비영 진공 배경과 평탄 시공간 대조에서 검증했다. 다음 내부 유체 해는 이 경계에 직접 연결할 수 있다. 내부 유체·자유 표면·목표 구동의 전체 해는 아직 없으므로 전체 동적 전하 목표는 진행 중이다. [전체 완료 판정](../docs/dynamic-charge-completion.md)을 유지한다.

분류: Imported from prior work. [Mendes–Ortiz 2018의 보충자료 식 (16)–(17)](https://arxiv.org/html/1802.07847v2)은 metric을 포함한 구면 scalar 외부 섭동과 무한대의 나가는 파 조건을 제공한다. 이번 독립 진공 유도는 그 식과 일치한다. 원문의 차갑고 barotropic인 내부 물질 가정을 현재의 온도·조성·두 열유속 모형에 이식하지 않는다. 원문은 `exp(+i omega t)`, 코드는 `exp(−i omega t)` 관례이므로 나가는 공간 위상의 부호가 반대다.

## 계량을 포함한 진공 유도

분류: Proven. 이 절은 `G=c=1`이며 `Phi=phi_0'`, `b=1−2m/r`, `F=N sqrt(b)`, `dr*/dr=1/F`로 정의한다. 정적 진공에서

\[
m'=\tfrac12 r^2 b\Phi^2,\qquad
(\ln F)'=\frac{2m}{r^2b}.
\]

분류: Proven. 비영 주파수의 섭동에 대해 Einstein 운동량 제약은 `delta m=r² b Phi delta phi`, `delta ln a=r Phi delta phi`를 준다. 이를 방사 제약에 대입해 일치성을 기호 검산했다. 이 관계의 시간 독립 질량 적분 상수는 별도의 정적 모드이며 방사 해에 임의로 더하지 않는다. 섭동 scalar 식은

\[
\delta\phi''+\left(\frac2r+\frac{2m}{r^2b}\right)\delta\phi'
+\left(\frac{\omega^2}{F^2}+\frac{2\Phi^2}{b}\right)\delta\phi=0.
\]

분류: Proven. `u=r delta phi`로 바꾸면

\[
\frac{d^2u}{dr_*^2}+[\omega^2-V]u=0,
\qquad V=2N^2\left(\frac{m}{r^3}-\Phi^2\right).
\]

계량 섭동을 고정하면 `−2N² Phi²`를 놓친다. 이 식은 선형 진공 섭동의 결과이며 물질의 섭동을 아직 풀지 않는다.

## 내부 해에 연결할 경계와 전하

분류: Proven. `x=R/r`, `Omega=omega R/c`, `r*(R)=0`을 사용하고 나가는 해를 `u_out=exp(i omega r*) h(x)`로 쓴다. `h(0)=1`인 무한대 8차 점근 전개로 수치 적분을 시작한다. 함수 `outgoing(mu,q,Omega)`는

`Z_out=R delta phi'(R)/delta phi(R)`와 `A_out=1/h(1)`

을 반환한다. 따라서 무한대 나가는 파의 `1/r` 계수는 표면장의 `R A_out`배다. 비영 주파수에서 이 복소 전달을 정적 Just 전하로 바꿔 쓰지 않는다.

분류: Proven. 실수 주파수에서는 들어오는 기저가 나가는 기저의 복소 켤레다. `H=h(1)`, 지정된 incoming 진폭을 `a_in`이라 하면 표면에서

`R delta phi' − Z_out delta phi = a_in H* (Z_out*−Z_out)`

인 affine 경계가 된다. 전체 나가는 진폭은 `(delta phi_surface−a_in H*)/H`다. **동일한 배경·incoming 조건에서 유체 응답과 물질 고정 대조의 차이는 `delta phi_surface_difference/H`로 직접 읽을 수 있다.** 큰 incoming/outgoing 값을 각각 계산한 뒤 빼는 방법을 내부 구현의 기본 경로로 삼지 않는다. 이 차분 항등식은 외부 기저의 선형 결합에서 따른다.

분류: Proven. 표면이 `xi`만큼 움직이면 Lagrangian 변화는 `Delta phi=delta phi+xi Phi`, `Delta phi'=delta phi'+xi Phi'`다. 따라서 위 경계의 우변에는 `xi (R Phi'−Z_out Phi)`를 더한다. 같은 incoming 조건의 나가는 진폭 차이는 `(Delta phi_contrast−xi_contrast Phi)/H`다. `def_moving_surface_boundary.py`에 이 어댑터를 구현하고 기호 항등식 및 세 실제 주파수의 비영 이동 대조로 확인했다. 표면 운동을 전하로 오인하는 좌표 항을 제거한다. 이는 scalar 매칭이며 압력·대기의 자유 표면 운동 법칙을 대신하지 않는다.

분류: Conjectural. 동반천체의 실제 near-zone 구동이 어떤 incoming 진폭을 만드는지, 물질 고정 대조와 물리적 표면 운동이 어떤 내부 해를 만드는지는 다음 결합 작업에 남아 있다. 외부 연산자의 존재만으로 scalar 힘이나 궤도 전하를 얻었다고 주장하지 않는다. 정적 극한을 비교할 때도 `delta M_ADM=0`인 방사 섭동의 극한과, 무한대 배경장을 바꾸는 정적 항성 계열을 구분해야 한다.

## 수치 검증

분류: Counterexample candidate. 배경 질량·lapse는 compact 좌표의 실제 진공 제약을 적분한다. 원 Just 해의 ADM 질량과 scalar 유속 상수를 초기값으로 사용하고, 표면의 질량과 lapse로 되돌아오는지 검사한다. 외부에서 별도의 EOS 호출은 없다. 나가는 해의 점근 시작점은 `x0=Omega/40` 및 `Omega/80`, 적분 허용치는 각각 `2e-12`, `2e-13`으로 실행 전에 고정했다. 정적 해는 `x0=1e-5`에서 시작해 정확한 Just 미분과 비교한다.

분류: Counterexample candidate. 총 18개 외부 해가 다음 사전 기준을 모두 통과했다.

| 대조 | 수치 결과 | 기준 |
|---|---:|---:|
| 점근 시작점·적분 허용치 변경에 따른 Z/A 최대 차이 | `2.65e-14` 이하 | `<2e-9` |
| 정확한 정적 Just 미분과 최대 차이 | `4.89e-15` 이하 | `<2e-10` |
| 평탄 시공간 `Z=−1+i Omega`, `A=1` | 저장 연산에서 0 | `<2e-10` |
| 나가는 파의 보존 유속 상대 결함 | `5.0e-15` 이하 | `<2e-8` |

분류: Counterexample candidate. 정적 대조는 `delta M_ADM=0` 조건에서 Just 매칭을 고정밀 미분한 독립 식이다. 강한 진공 대조는 `mu=0.1`, `q=0.03`을 사용하므로 실제 약한 배경의 거의 0인 scalar 계량 항만 검사하는 시험이 아니다. 정적·유속·시작점 대조의 수치 일치는 엄밀한 연속 오차 상계나 내부 유체 검증은 아니다.

분류: Counterexample candidate. 실제 비영 배경 표면의 기본 궤도 주파수 `Omega=9.2070058588e-5`에서는

`Z_out=−1.000004238468154 + i 9.20700588142103e-5`,

`A_out=0.9999957627792986 − i 6.261575847871571e-9`

를 얻었다. 이 허수부는 외부 파동 경계의 일부다. 이를 내부 열 완화의 검출이나 자유로운 미분 nuisance를 벗어난 신호로 세지 않는다.

분류: Counterexample candidate. 내부 실행 시간은 약 `0.64초`, 프로세스 전체는 `2.64초`다. 90초 상한·18해 상한을 지켰고 추가 native 호출과 물질 시간 적분은 0이다. 코드 `verification/def_radiative_exterior.py`, 이동 표면 어댑터 `def_moving_surface_boundary.py`, 출력 `outputs/direct-eos-gr33/def-radiative-exterior/`에 계획·기호 유도·해·대조·SHA를 보존했다. 착수 checkpoint는 `23259c45`, 이전 문서 원본은 `588af618`이다.

## 다음 실행의 직접 의존성

분류: Conjectural. 새 외부 경계를 비영 배경의 물질·계량·scalar 선형 연산자에 연결해야 한다. 먼저 이미 저장된 분자 GR-4/8 정적 구조를 재사용할 수 있는지 확인하고, coarse 격자의 정수압 기준 잔차 제거에 의존하지 않는 배경을 만든다. 표면의 유한 압력·대기와 열수송을 명시하고, 자유 표면 또는 실제 외부 물질과의 접합을 닫아야 한다. 원문 cold barotrope로 내부를 교체하거나 계량을 동결해서 문제를 더 쉽게 바꾸지 않는다. 내부·외부가 연결된 뒤에야 실제 구동, 정적 비교 제거, 지배 오차와 전체 동적 전하 완료를 판정한다.
