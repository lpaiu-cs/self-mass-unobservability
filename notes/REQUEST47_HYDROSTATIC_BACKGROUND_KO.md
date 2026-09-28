# Phase47 — 같은 재고의 비영 기계적 항성 배경

분류: Counterexample candidate. **저장된 5,735셀 분자 EOS·GR 구조에서 각 셀의 바리온·26종 조성·엔트로피를 유지하면서, 비영 scalar와 물질의 정수압 조건을 함께 풀었다.** 반경·질량·압력·scalar가 조정되고 lapse도 복원됐다. 모든 중간점의 native EOS 대조와 셀별 힘 균형 검사를 통과했다. 이전 24셀 Cauchy snapshot의 힘 불균형을 그대로 정상 배경으로 사용하는 병목을 넘었다.

분류: Conjectural. 이번 상태는 유한 압력을 가진 광구에서 자른 **기계적 초기 상태**다. 물리 대기·자유 표면을 포함한 완전한 항성, 열적 정상성, 궤도 시간 동안 유지되는 배경 또는 동적 전하의 완성은 아니다. 단계 46의 고정 벽 진화와도 아직 같은 초기 상태로 연결하지 않았다.

## 방정식과 유지한 물질

분류: Imported from prior work. 기하 단위에서 `b=1−2m/r`, `A=exp(beta phi²/2)`, `alpha=beta phi`, `N=exp(nu)`, `Phi=phi'`로 쓴다. 물질의 Jordan 값에 대해 `epsilon_E=A⁴ epsilon_J`, `P_E=A⁴ P_J`이다. [Damour–Esposito-Farèse의 식 (3.6a–f)](https://arxiv.org/pdf/gr-qc/9602056)를 다음처럼 사용한다. 그 논문의 `nu`는 여기 `log N`의 두 배다.

```
m'       = 4 pi r² epsilon_E + r² b Phi² / 2
nu'      = m/(r² b) + 4 pi r P_E/b + r Phi²/2
phi'     = Phi
Phi'     = −[2/r + nu' − lambda'] Phi
           + 4 pi alpha (epsilon_E−3 P_E)/b
(log P_J)' = −(epsilon_J+P_J)/P_J * (nu' + alpha Phi)
B'       = 4 pi r² A³ rho_J / sqrt(b)
```

분류: Proven. `lambda'=((m'/r)−m/r²)/b`를 scalar 식에 대입하면 원문 식 (3.6d)가 나온다. `d/dB=(B')⁻¹ d/dr` 변환과 함께 기호 검산을 통과했다. 각 셀에서 엔트로피·조성이 일정하면 기계적 평형의 첫 적분은 `log(h_J A N)=constant`이고 `h_J=(epsilon_J+P_J)/rho_J`다.

분류: Counterexample candidate. 원문의 저온 barotrope를 가져오지 않고 저장소의 유한 온도 분자 EOS를 사용했다. 원 `dm`, `X`, 기준 EOS 엔트로피 배열을 그대로 보존한다. 중앙 압력·광구 반경·광구 질량·중앙 scalar·표면 scalar 기울기의 다섯 미지수를, 중앙과 광구에서 적분한 두 바리온 좌표 가지의 일치로 정했다. 압력 잔차를 빼는 정수압 기준 힘이나 바리온 질량 조정은 없다.

분류: Counterexample candidate. `phi_infinity=0`에서는 정규화 scalar의 선도 극한을 보조적으로 풀고, 물질·계량은 GR 방정식으로 환원한다. `phi_infinity=0.001`, `beta=−4`에서는 `A`, scalar 응력, 물질 scalar 힘을 모두 유지한다. 바깥 scalar와 lapse 정규화는 기존 Just 지도를 재사용했다. 유한 물질 압력이 남아 있으므로 이 지도 자체를 물리적 진공 접합으로 인정하지 않는다.

## 오래 걸린 계산을 재사용한 방법

분류: Imported from prior work. 원 세밀한 GR 재구성에는 EOS 표 조회 약 319만 회와 표 밖 직접 역산 83,512회가 기록돼 있다. 원 상태와 17점 등엔트로피 표는 그대로 재사용했다.

분류: Counterexample candidate. 원 압력 범위가 부족할 수 있는 203셀에 1,902개 압력점만 보충했다. 각 점은 같은 엔트로피의 엄격한 native 역산이며 원 17점은 유지했다. 9회 고정밀 fallback과 그 원 실패 기록도 보존했다. shooting 중 표 범위를 벗어나면 실패시키며 외삽하거나 숨은 native 역산으로 실행 시간을 늘리지 않는다. 표의 연속 오차 인증은 아니다.

분류: Counterexample candidate. 새 식의 한 전체 별 가지 비교는 1.51초였다. 이 속도에 근거해 shooting의 600초 상한, 각 배경의 목적함수 45회 한도와 원 일치 기준 `1e-8`을 유지했다. 비영 해는 26회 평가에서 수렴했다. 새 궤도 적분이나 장기 시간 진화는 실행하지 않았다.

## 결과와 독립 대조

분류: Counterexample candidate. 아래 반경과 질량은 Einstein 좌표/기하 단위의 광구 값이다. 관측된 Jordan 질량·반경으로 읽지 않는다.

| 값 | GR 재생 | 비영 배경 |
|---|---:|---:|
| 광구 반경, m | 68,819,337.45172 | 68,818,637.74977 |
| 광구 질량 `GM/c²`, m | 291.6871724624 | 291.6865890267 |
| 형식적 Just ADM 질량 `GM/c²`, m | 291.6871724624 | 291.6865890366 |
| 형식적 외부 계수 `alpha_A/phi_infinity` | −4.000356685518 | −4.000356686614 |

분류: Counterexample candidate. GR 재생은 기존 subdivision 8 상태와 최대 `7.042e-10`으로 일치했다. 비영 배경의 다섯 경계 일치 잔차는 최대 `3.809e-12`다. 반경은 약 699.70m 줄고 최대 log 밀도·온도 변화는 각각 `4.077e-5`, `2.492e-5`다. 이 변화는 지정한 기계적 배경 사이의 차이며, 시간 진화나 실제 동반성 구동의 응답은 아니다. 표의 외부 계수를 곧바로 고정 재고 질량 미분 또는 검출 가능한 신호로 해석하지 않는다.

분류: Counterexample candidate. **비영 상태의 5,735개 중간점 전부를 원 native EOS로 평가했다.** 배열만 비교한 검사가 아니다. 저장된 압력·온도를 원 EOS에 넣어 표로 구성한 밀도·엔트로피·에너지와 대조했다.

| 검사 | 최대 결함 | 사전 기준 |
|---|---:|---:|
| native log 밀도 차이 | `6.276e-11` | `<1e-7` |
| `T delta s / (du/dlogT)` | `3.051e-9` | `<1e-7` |
| 총 물질 에너지 상대 차이 | `1.738e-15` | 기록값 |
| lapse 복원 시 원 다섯 상태 재생 차이 | `3.809e-12` | `<1e-8` |
| 셀별 `log(h_J A N)` 첫 적분 결함 | `1.888e-15` | `<1e-12` |

분류: Counterexample candidate. 같은 해 위에서 lapse 방정식을 적분해 중심 `log N=−4.17629517e-5`, 광구 `log N=−4.23850043e-6`를 얻었다. lapse 복원은 새 EOS 호출·배경 재맞춤 없이 수행했다. native 검사 문턱은 이 초기 상태의 EOS 대조 기준이며, 선언한 동적 전하 `1e-9`의 전체 오차 예산을 인증하지 않는다.

## 실패를 보존한 복구

분류: Counterexample candidate. 첫 비영 shooting은 실패했다. 상대 유한 차분이 거의 0인 log 반경·질량 보정값에 적용되면서 차분 폭이 `1.377e-17`, `9.526e-18`이 됐고, 두 잔차 열이 정확히 0이었다. 이어 무제약 제안에서 지수 overflow가 났다. 원 실패·소스·경로를 보존하고, 절대 중앙 차분을 쓰는 별도 실행에서 이미 통과한 GR 상태와 EOS 표를 재사용해 비영 해를 수렴시켰다.

분류: Counterexample candidate. 수렴 해를 저장한 뒤에는 결과 메타데이터가 부모 폴더의 EOS 표를 자식 폴더에서 찾는 오류로 실행이 종료됐다. 저장된 해를 버리거나 다시 맞추지 않았다. 원 표 SHA와 저장 해를 확인하고 경계 일치를 재생한 뒤 메타데이터만 완성했다. 원 실행의 종료 코드 1과 복구 근거를 보존한다. 최종 native 대조와 lapse 작업은 종료 코드 0으로 완료됐다.

## 비용·재현·다음 병목

분류: Counterexample candidate. 내부 타이머는 표 보충 123.43초, 비영 shooting 마지막 평가까지 33.56초, 전 셀 native 감사 71.20초, lapse 복원 22.96초를 기록했다. `/usr/bin/time`은 각각 표 보충 wall 118.97초, 실패를 포함한 최초 shooting 27.77초, 비영 복구 실행 40.19초, native 작업 74.76초를 기록했다. 두 타이머를 같은 측정치로 합치지 않는다. CPU 작업이며 네 EOS worker를 사용했다. 각 프로세스의 기록된 최대 RSS는 약 122–155MiB로, worker 전체 합계 메모리 측정은 아니다. 모든 hard timeout 안에 끝났다.

분류: Imported from prior work. 실행 전 checkpoint는 `972195f`다. `verification/def_hydrostatic_background.py`의 `prepare`, `pilot`, `table_build`, `shoot_pilot`, `run`이 원 실패를 포함한 경로다. `def_hydrostatic_absolute_shoot.py`가 비영 해를 저장하고, `def_hydrostatic_native_audit.py`가 메타데이터 복구와 전 셀 native 검사를, `def_hydrostatic_lapse.py`가 lapse 복원을 수행했다. 출력은 `outputs/direct-eos-gr33/def-hydrostatic-background/`에 보존한다. 완료된 결과를 덮어쓰는 재실행은 허용하지 않는다.

분류: Conjectural. 다음 직접 병목은 이 기계적 초기 상태에 **물질 표면·대기와 열유속의 물리 조건**을 붙이고, 동일 내부 해를 단계 43의 동적 외부 및 선언 구동에 연결하는 일이다. 광구의 작은 압력이나 짧은 수치 수렴만으로 자유 표면 오차·열적 정상성을 면제하지 않는다. 같은 재고 정적 비교, 직접/물질 매개 전하 분리와 결론을 좌우하는 오차 상계가 필요하다. [전체 완료 조건](../docs/dynamic-charge-completion.md)은 줄이지 않는다.
