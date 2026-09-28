# 단계90 — 실제 EOS로 외향 복사 외피 재구성

분류: Counterexample candidate. 기존 비정상 외곽에 임의의 표면 유속을 붙이지 않고, native132셀의 반경·질량·scalar·scalar 기울기·lapse·압력·온도·조성을 접합값으로 유지한 채 EOS와 상태별 불투명도로 외피를 실제 적분했다. 물질·복사의 압력 수지, GR/DEF 정역학 제약, Tolman 온도 기울기와 일정한 적색편이 광도를 함께 풀었다. 외향 Eddington 경계까지 맞는 회색 외피 해를 얻고 수치 기준 및 독립 적분 검사를 통과했다.

분류: Conjectural. 성과는 loophole progress다. 새 외피 해는 전체 항성의 기존 물질 재고·내부 광도와 동시에 맞춘 해가 아니다. 완전한 시간 의존 Einstein 방정식과 비회색 스펙트럼 대기, 궤도 전하·관측 연결도 미완료다. 이 한계를 해소하기 전 원 결합 진화의 배경을 새 파일로 단순 교체하지 않는다.

## 선언 모형과 실제 식

분류: Imported from prior work. [MESA의 회색 대기 설명](https://docs.mesastar.org/en/24.03.1/atm/t-tau.html)은 온도–광학깊이 관계와 정역학 압력식을 결합하며, 국소 상태마다 EOS·불투명도를 다시 평가하는 선택을 구분한다. 이번 계산은 아래 구면 GR/DEF 제약과 명시적인 Eddington 경계를 사용한 별도의 조건부 구현이다.

분류: Proven. Einstein 반경을 r, g=dln(AN)/dr, 물질 기준 복사 유속을 F로 놓는다. E_rad=3P_rad=a_rad*T^4 및 회색 수송을 가정하면 다음 두 식의 합은 dP_total/dr=-(e_total+P_total)g다. native EOS의 광자 에너지를 다시 더하지 않고 gas 에너지·압력만 분리한다.

```text
F = L_infinity / (4*pi*r^2*N^2*A^4)
dP_gas/dr = -(e_gas + P_gas)*g + rho*kappa*A*a*F/c
dlnT/dr = -g - rho*kappa*A*a*F/(4*P_rad*c)
```

분류: Counterexample candidate. m,phi,phi',lnN과 외피 바리온 질량도 동시에 적분했다. 외곽에서는 P_gas=1 dyn/cm2에서 E_rad=2F/c를 맞추고, P_gas=0.25 dyn/cm2 대조를 추가했다. 이는 유한 압력에서 자른 Eddington 외향 경계다. 실제 무한대 진공 접합이나 스펙트럼 수송을 인증하지 않는다. Rosseland 평균을 Planck 흡수율로 바꾸지 않았으며, LTE와 회색 각도 폐쇄는 별도 가정으로 남는다.

분류: Proven. 본 구면 모형처럼 상쇄 흐름 없이 L_infinity>0이면 정확한 정적 시공간으로 해석할 수 없다. 복사 에너지 수지는 질량 손실을 요구한다. 분류: Counterexample candidate. 본 해는 그 정역학·열적 제약을 푼 준정적 후보이며, 시간–반경 Einstein 방정식은 풀지 않았다. 계산된 방출 에너지의 원 궤도1주기/기존 ADM 에너지 비9.64475e-17과 최대 F/(c*e_total)=2.71660e-7은 크기 진단이다. 이를 동적 전하 오차나 실제 관측 상계로 바꾸지 않는다.

## 결과

| 분류 | 측정량 | 값 |
|---|---|---:|
| Counterexample candidate | 새 경계가 요구하는 L_infinity | 2.1733737400e33 erg/s |
| Counterexample candidate | 기존 내부 확산 광도 | 2.1607503605e33 erg/s |
| Counterexample candidate | 내부 광도 대비 불일치 | +0.584212768% |
| Counterexample candidate | ODE 허용오차 대조의 광도 상대 차이 | 2.50779e-9 |
| Counterexample candidate | ODE 허용오차 대조의 최대 온도 상대 차이 | 2.09389e-9 |
| Counterexample candidate | 외곽 압력 절단 대조의 광도 상대 차이 | 2.30330e-7 |
| Counterexample candidate | 접합점에서 바깥 절단점까지 광학깊이 | 61.4400 |
| Counterexample candidate | 기존 최외곽 압력에서 이전/새 온도 | 18827.9472 / 14326.1249 K |
| Counterexample candidate | 같은 압력의 native 겹침 구간 최대 온도 변화 | 23.9103% |
| Counterexample candidate | 외피 최대 복사 가속도/물질 중력 비 | 1.37924% |

분류: Counterexample candidate. 별도의 적분 검증은 저장된 Tolman 온도 양끝과 반경별 opacity로 광도를 재구성했다. 광도 상대 차이8.72817e-8, 정역학 압력 적분1.89242e-8, 바리온 적분1.56557e-8로 사전1e-3 기준을 통과했다. 압력 좌표 생산과 독립인 native 밀도 좌표5상태에서 압력·단열 기울기를 확인했다. 사용한 저장 상태는 공통 불투명도 표의 R/T 범위 안이며 경계 clipping이 없다. 이것은 표의 물리 정확도를 인증하는 결과가 아니다.

분류: Counterexample candidate. tau>10이고 수송 Kn<0.01인23개 저장 표본에서 nabla-nabla_ad의 최대값은-0.0466664로, 해당 균일 조성·단열 parcel 기준의 불안정 표본은 없었다. 전체 대류·비LTE 안정성을 증명하지 않는다. 기준의 문헌 정의는 [MESA controls](https://docs.mesastar.org/en/stable/reference/controls.html)에 둔다.

분류: Counterexample candidate. 새 외피 바리온 비중은1.78114258e-12다. 기존 native0..131셀과132셀 절반의 비중1.77093573e-12와의 차이는1.02068462e-14다. 이 참조에는 과거 별도로 붙인 등엔트로피 대기 재고가 들어 있지 않으므로 전체 별의 질량 불일치 수치로 해석하지 않는다. 동일 물질 재고의 전체 재접합은 아직 필요하다.

## 이번에 해결한 것과 남은 접합 조건

분류: Counterexample candidate. 실패했던 고정 외곽 온도를 유지하지 않고 외향 방출을 수용하는 실제 native EOS 외피를 찾았다. 새 온도 변화23.9%는 기존5% 선형 온도 창 밖이다. 단계89의 작은 온도·전하 되먹임 보정만으로 이 재구성을 대신할 수 없다.

분류: Proven. 정상 상태에서 접합면에 별도 에너지 저장·열원이 없으면 안팎 광도는 같아야 한다. 분류: Counterexample candidate. 현재 두 광도 차이는 약1.26234e31 erg/s로 수치 대조보다 크다. 그러므로 기존 내부의 압력·온도·확산 광도와 본 외향 경계를 모두 그대로 둔 정상 접합은 성립하지 않는다. 선언 회색 모형에서 확인한 이 불일치를 모든 물리 대기의 no-go로 일반화하지 않는다.

분류: Conjectural. 다음 실제 병목은 새 열적 외피를 전체 물질 재고 및 내부 열유속·구조와 함께 재매칭하는 것이다. 그 뒤 바뀐 EOS/경계의 결합 응답을 다시 계산해야 기존 외부 전하 값의 물리적 해석을 평가할 수 있다. 원 배경의 수치 수렴 성공과 새 물리 배경의 완성은 구분한다.

## 예산과 첫 실행의 수정

분류: Counterexample candidate. 파일럿 경로9.99962초를 측정한 뒤, 독립적인 근 찾기60회에 대한 초기940.49초 예측은600초 예산을 넘어 실행하지 않았다. 파일럿 함수값과 앞 경로의 근을 재사용하고 생산 적분30회 상한을 정했다. 최종 생산은29회,283.745초였다. native 호출은 생산10211회, 파일럿351회, 감사5회로 기록됐으며 첫 실패의 호출 수는 별도 계측되지 않았다.

분류: Counterexample candidate. 첫 shooting 시도는 RK45의 중간 시험점이 물리 경계 사건을 찾기 전에4000K EOS 영역 밖으로 넘어가 중단했다. 해당 시험 스텝을 거부하고 줄여 재시도하도록 수정했다. EOS를 clipping·외삽하거나 물리 조건·수락 기준을 바꾸지 않았다. 첫 소스·계획·실패를 보존했다. 실패에4초를 청구해 파일럿·생산·감사 내부 합298.945초로600초 상한 안이었다. CPU/BLAS1스레드, 메모리 상한4GB이며 전 과정 메모리 최고값을 실측했다고 주장하지 않는다. WSL·import·문서 시간은 내부 타이머와 별도다.

근거: `outputs/direct-eos-gr33/def-native-radiative-envelope/`의 계획·첫 실패·예산 재평가·세 저장 해·결과·독립 감사. 실행은 `verification/def_native_radiative_envelope.py`, 감사는 `verification/verify_native_radiative_envelope.py`다.
