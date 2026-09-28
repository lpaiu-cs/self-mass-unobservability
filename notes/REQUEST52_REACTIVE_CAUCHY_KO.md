# 단계 52 — 자유 표면의 반응·물질·scalar·계량·중성미자 시간 결합

분류: Counterexample candidate. 단계 51에서 별도로 계산하던 중성미자 응력을 **같은 자유 표면 배경의 물질 운동·scalar 파동·계량 제약 시간식에 연결했다.** 원 native 초기 핵반응과 화학 에너지 방향을 사용해 좌표시간 `0.23080495568542375 s`의 다섯 선형 경로를 실행했다. 사전 등록한 속도/scalar 시간 간격·복사 집계·외부 경계·선형 잔차 기준을 모두 통과했다. 이번 결과는 loophole progress다. 광자·전도 열유속과 비선형 반응 갱신을 포함한 전체 항성 진화 또는 동적 전하·관측의 완료는 아니다.

## 해결한 연결

분류: Imported from prior work. 단계 49의 자유 표면 단열 연산자, 단계 50의 전 셀 native 반응·화학 에너지 방향, 단계 51의 인과적인 null 광선 수송을 재사용했다. 이전의 단열/준정적 계산과 달리 이번에는 물질 관성과 scalar의 시간 2차 미분을 실제로 전진했다. 배경은 원 5,735셀 바리온·조성을 유지한 비영 기계적 평형과 명시적으로 가정한 희박 기체 외층이다. 열 정상성이 확인된 별은 아니다.

분류: Conjectural. 초기 중성미자장은 0이고 초기 방출은 유체 정지계에서 등방이다. 중성미자를 질량 없는 무충돌 복사로 가정한다. 초기 native 핵반응/열중성미자 방출률과 열역학 미분을 고정해 그 원천 세기에 대한 1차 응답을 구한다. 광자·전도 열유속은 이 부분 문제에서 0으로 둔다. 실제 항성의 전체 열 변화가 작다는 가정이나 인증은 추가하지 않는다.

분류: Proven. 이 원천 전개에서는 복사와 계량 섭동의 곱, 방출과 섭동 유체 속도의 곱이 2차다. 따라서 초기 정적 계량에서 계산한 복사장을 1차 물질·scalar·계량 응답에 넣는 것은 이 차수의 일관된 결합이다. 이는 유한 진폭의 복사·계량 되먹임을 계산했다는 뜻이 아니다.

분류: Counterexample candidate. 방출 집계 격자를 원 native 재고의 바깥 면에 맞췄다. 마지막 집계 셀을 실제 자유 표면까지 넓혀 0밀도 외층에 물질 손실을 부여하지 않는다. 각 집계 셀에서 물질에 빼는 에너지와 복사에 넣는 에너지가 같다. 핵반응의 화학 기여는 보존하고, native 열역학 방향의 중성미자 손실만 이 공유 집계 원천으로 바꿨다. 원 5,735개 방출 위치 모두에서 광선을 개별 발사한 것은 아니다.

## 방정식과 구현

분류: Proven. 아래는 `r,m`을 실제 자유 표면 반경 `Rstar`로, `p,e`를 기하 단위로, 시간을 `tau=c*t/Rstar`로 정규화한 1차 식이다. `A=exp(-2 phi²)`, `alpha=-4 phi`, `b=1-2m/r`, `w=e+p`, `Phi=phi_prime`다. 물질의 Lagrangian 섭동을 사용한다.

```text
xi = r*zeta
eta_ad = Delta p/p + Gamma*rho_ref
Delta ln rho = eta_ad/Gamma
Delta p = p*(eta_ad-Gamma*rho_ref)
Delta e = w*eta_ad/Gamma-loss_J

H = r²*b*Phi*f-(4*pi*r²*A⁴*p+r²*b*Phi²/2)*xi
Delta m = H+J
(N*a*J)_prime = 4*pi*r²*N*a*(E_nu-A⁴*loss_J)
J = -(G/c⁴)*1e-7*integrated_redshifted_face_energy/(Rstar*N*a)

Delta lambda = (Delta m/r-m*xi/r²)/b
xi_prime = -eta_ad/Gamma-2*zeta-Delta lambda-3*alpha*f

Delta g = dg · (xi,Delta m,Delta p,Delta e,f,V)+4*pi*r*P_nu/b
Delta F = dF · (xi,Delta m,Delta p,Delta e,f,V)+4*pi*r*Phi*(E_nu-P_nu)/b

B = g*(2*zeta+Delta lambda+3*alpha*f+eta_ad-Gamma*rho_ref)
    +g*loss_J/w-Delta g
eta_ad_prime = w/p*(B-xi_tautau/(N²*b))
               -g*(eta_ad-Gamma*rho_ref)+(Gamma*rho_ref)_prime
f_prime = V+xi_prime*Phi
V_prime = Delta F+xi_prime*F+(f_tautau-xi_tautau*Phi)/(N²*b)
```

분류: Proven. `g,F`와 그 미분은 기존 정수압·scalar 연산자와 같다. 복사의 trace는 0이지만 에너지와 방사 압력은 계량에 기여하므로 `Delta g,Delta F`에서 사라지지 않는다. `J`의 적분 인자는 `N*a`다. 광선이 운반한 누적 면 에너지로 질량 결함을 직접 구성하면 큰 가역 물질 일을 빼서 작은 손실을 찾는 소거를 피할 수 있다. 이 항등식과 압력 변수 변환을 기호 검산했다.

분류: Counterexample candidate. 공간에는 원 내부 5,863구간과 외부 반경당 128구간을 사용했다. 공간 중점 식을 `K*y-D*y_tautau=f(tau)`로 조립하고 Newmark 평균 가속도법으로 전진했다. `eta_ad`를 사용해 실제 밀도를 두 큰 열역학 항의 차로 계산하지 않으며, `(Gamma*rho_ref)_prime`는 노드 간 정확한 차로 넣었다. 선형계 분해는 경로당 한 번이고 저장 계수와 2회의 확장 정밀도 잔차 보정을 사용한다.

분류: Counterexample candidate. 중심은 정칙 조건, 물질 표면은 동적 압력 정칙 조건이다. 외부의 `zeta`는 좌표 표지이며 진공 물질이 아니다. 외부 끝에서 Eulerian scalar `f-xi*Phi=0`을 두고, 동일 공간 간격으로 끝 반경을 `2Rstar`에서 `3Rstar`로 옮긴 대조를 했다. 한 광행시간의 내부 읽기에 대한 경계 영향이 문턱 아래였다. 이 유한 외부 대조를 무한대 나가는 파 전하의 계산이라고 부르지 않는다.

분류: Counterexample candidate. 광선 구간의 진입·이탈 사건을 미리 정렬하고 누적 기울기를 저장하여 각 시각의 복사 에너지·운동량·압력을 계산한다. 원 구간 적분의 직접 재생과 차이는 공통 방출 에너지 노름으로 `5.52e-16` 이하였다. 복사 에너지와 같은 적색편이 고유 부피 가중치를 사용해 물질 손실과 반경별 질량 결함을 투영했다.

## 사전 기준과 판정

분류: Counterexample candidate. 원 미세 복사 집계에서 16/32/64단계, 성긴 복사 집계에서 64단계, 외부 `3Rstar`에서 64단계를 실행했다. 읽기는 원 native 질량 가중 속도 RMS와 Eulerian scalar RMS다. 다음 차이는 표본 시각 전체의 최대 차이를 가장 미세한 경로의 최대 읽기로 나눈 값이다. RMS 읽기의 수렴이며 모든 공간 성분이나 연속 시각의 균일 오차 보장은 아니다.

| 대조 | 속도 RMS | Eulerian scalar RMS | 사전 문턱 |
|---|---:|---:|---:|
| 16→32단계 차이 | `0.07640%` | `0.20698%` | 감소 여부 |
| 32→64단계 차이 | `0.02192%` | `0.05199%` | 2% |
| 경험적 시간 차수 | `1.80154` | `1.99307` | 1.5 초과 |
| 복사 48×24 / 96×48 차이 | `1.26779%` | `0.42335%` | 3% |
| 외부 2R / 3R 차이 | `2.42e-15` | `1.39e-13` | 상대 0.002 |

분류: Counterexample candidate. 최대 선형 잔차는 `3.47e-16`으로 `1e-9` 기준을 통과했다. 총 방출=복사 재고+외부 탈출 수지의 최대 상대 차이는 `1.25e-15`로 `2e-13` 기준을 통과했다. 새 시간 연산자에 무원천 조화 응답을 대입하여 이전 단열 연산자를 재생한 상대 차이는 `2.10e-16`이다. 기존 단열/준정적 **전하 공간 기준 미달은 그대로 유지한다.**

## 같은 해에서 얻은 물질·계량

분류: Counterexample candidate. 마지막 상태에서 원 native 중심의 밀도·온도·26종 조성 섭동과 물질 속도, Eulerian scalar, 질량 및 radial metric 변화를 복원했다. lapse는 같은 선형 방사 제약을 적분해 후처리했다. 다음 값은 **초기 핵반응/중성미자 부분 응답**이며 실제 별의 관측값이나 전체 열 진화값이 아니다.

| 끝 시각 0.2308049557초의 양 | 값 |
|---|---:|
| 질량 가중 물질 속도 RMS | `6.84161e-8 m/s` |
| 질량 가중 Eulerian scalar RMS | `5.98610e-26` |
| native 최대 `abs(Delta ln rho)` | `3.68066e-12` |
| native 최대 `abs(Delta ln T)` | `3.20268e-11` |
| native 최대 `abs(Delta X)` | `5.11032e-14` |
| 최대 `abs(delta nu)` | `2.15242e-22` |
| 최대 `abs(Delta lambda)` | `1.68535e-19` |

분류: Counterexample candidate. 같은 시각에 물질 표면을 통과한 중성미자 에너지는 `9.20219e32 erg`이고 보존 질량 결함 `Rstar*J`는 `-7.60349e-19 m`다. 독립 광선 교차 이력과의 수지 차이는 총 방출 에너지 노름으로 `9.58e-16`이다. `2Rstar` 밖으로 탈출한 에너지는 아직 0이다. 표면 질량 결함을 무한대 총질량 감소나 scalar 전하로 동일시하지 않는다.

분류: Conjectural. 섭동의 작은 크기만으로 native 도함수 오차·원천 동결 오차·비선형 나머지항을 엄밀히 제한하지 못한다. 열유속을 생략한 부분 응답의 작은 `Delta T`가 전체 별의 열 정상성을 입증하지 않는다. 현재 비교에서 가장 큰 수치 차이는 복사/원천 집계에 의한 속도 차이지만, 빠진 광자·전도 수송의 물리 오차와 비교한 최종 우선순위는 아직 아니다.

## 자원·재현·다음 병목

분류: Counterexample candidate. 8단계 예비 물질 경로는 약 1.11초였고, 미측정 미세 광선 비용을 보수적으로 확대해 생산 97.84초를 예상했다. CPU 1개·BLAS 1개·새 native 호출 0회, 생산 180초 상한과 다섯 경로를 먼저 등록했다. 실제 생산은 16.81초, 전체 프로세스 벽시간은 20.09초, 최대 RSS는 402,348 KiB(약 392.9 MiB)로 예측 600 MiB 이내였다. 생산은 합계 240단계였으며 모든 기존 상태와 실패 기록을 보존했다. 기준 미달을 핑계로 기간·격자·경로를 확대하지 않았다.

분류: Imported from prior work. 구현은 `verification/def_reactive_cauchy.py`, 저장 상태 복원·독립 교차 검사는 `verification/def_reactive_cauchy_readout.py`다. 입력 SHA·계획·예비 측정·광선·다섯 경로·판정·물질/계량 끝점은 `outputs/direct-eos-gr33/def-reactive-cauchy/`에 있다. 생산 실행 중 원천과 계획을 바꾸지 않았다.

분류: Conjectural. 다음 직접 병목은 **같은 자유 표면에서 광자·전도 열유속과 물리적 방사 경계를 넣는 것**이다. 이미 실패가 확인된 열 정상성을 생략하거나 Rosseland 평균을 흡수/Planck 평균으로 대신하지 않는다. 가능한 초기 열 상태와 누락된 물리 입력부터 기존 자료로 구분하고, 원천 갱신이 필요한 범위를 정해야 한다. 그 뒤에도 전하 공간 오차, 실제 쌍성 구동, 무한대 나가는 파 전하, 정적 비교·nuisance 제거와 전체 오차/관측 연결은 남는다. 전체 목표를 완료로 표시하지 않는다.
