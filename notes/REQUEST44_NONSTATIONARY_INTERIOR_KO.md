# Phase44 — 비영 배경의 비정상 내부 결합 진화

분류: Counterexample candidate. 단계 41의 실제 비영 scalar Cauchy 상태에서 물질·scalar·계량을 함께 진화시켰다. 온도, 두 열유속, 26종 조성 이류를 유지하며, 이전 공간 연산자의 정수압 기준 힘 제거를 해제했다. 동일 초기 상태의 양·음·무구동을 함께 계산한다. **고정 벽 내부의 비정상 결합 문제**이며 자유 표면 항성의 복사 전하 완료는 아니다.

## 선택한 경로와 보존한 제한

분류: Imported from prior work. 저장된 분자 GR-4/8 구조는 바리온·조성·엔트로피를 보존한 정수압 재구성이다. 그러나 비영 scalar, 열유속을 포함한 정상 해와 물리적 대기는 아직 아니다. 단계 41의 비영 상태도 radial constraint를 만족하는 Cauchy 자료이며 정수압 평형으로 인증되지 않았다.

분류: Proven. 구면 polar-areal 계량에서 방사 에너지 운동량 유속이 `S`일 때 `dm/dt = −4 pi (G/c^4) c r² (N/a) S`다. `v=0`에서도 `S=Q`이므로 비영 순 열유속의 자료를 정확히 정지한 계량·물질 배경으로 취급할 수 없다. 기존 유도 및 이번 기호 검사에 따른 제한이며 열·조성을 제거한 cold barotrope로 교체할 근거가 아니다.

분류: Counterexample candidate. 따라서 이번 계산은 정상 전달함수를 가정하지 않는다. 저장된 24셀의 같은 비영 초기 자료와 원래 native EOS를 사용하고, 무구동 배경도 함께 진화한다. scalar 초기 조건은 `phi_infinity=0.001`, `Pi=0`, 단계 41의 공간장과 constraint 계량이다. 이 coarse 상태의 시간 변화에는 비정상 열수송과 압력 불균형이 함께 들어간다. 이를 실제 별의 느린 진화라고 해석하거나 공간 수렴된 해라고 부르지 않는다.

## 제거한 인공 정수압 보정

분류: Proven. 기존 운동량식에서 기준 압력 면 유속, 기준 중력항, 기준 구면 압력항을 모두 0으로 두면 다음 보존 잔차가 남는다.

`Delta(a S) + h c div(F_S) − h [−a_t S − c N nu_r E + c N P (Delta area)/V + c N alpha(phi) trace * phi_r] = 0`.

분류: Proven. 여기서 마지막 항은 `alpha(phi) trace phi_r`의 곱이며 trace의 공간 미분이 아니다. `alpha=beta phi`, `trace=−epsilon+3P`, `beta=−4`다. 원 코드의 면 공유 압력 연산자와 음향 점프 응력은 유지한다. 기준 힘을 빼서 정지시킨 상태를 물리적 평형으로 인정하지 않는다.

분류: Counterexample candidate. 모든 native 잔차 평가는 scalar와 radial metric constraint를 다시 푼다. Newton 반복에 쓰이는 근사 행렬은 최종 물리 잔차의 대용품이 아니다. 물질 에너지의 reciprocal scalar work, 바리온·조성 보존, proper conduction time의 conformal 변환을 유지한다. 초기 상태 파일의 `at=0`은 이전 초기화 관례에서 남은 미사용 자리값이다. 이전 상태의 `at`는 단계 방정식의 입력이 아니며, 별도 재생 검사에서 실제 초기 열유속에 따른 `a_t`를 계산한다. 이 자리값을 초기 정지성 증거로 쓰지 않는다.

## 구동과 판정

분류: Counterexample candidate. 같은 초기 상태와 scalar 표면값에 진폭 `+1e-5`, `−1e-5`, `0`의 `sin^8(pi t/tc)` 펄스를 더한다. `tc=R/c=0.2295566003 s`이며 `tc` 이후 구동은 0이다. 이번 실행은 정확히 구동 종료 시점에서 끝난다. nominal 12/24/48 시간 격자의 앞 6/12/24단계, 세 경로로 총 126단계다. 이는 궤도 구동의 진폭·시간 척도를 대체하지 않는다.

분류: Proven. 읽기는 `odd=(plus−minus)/2`, `even=(plus+minus)/2−undriven`으로 나눈다. Jordan log 온도는 `delta log T_E−log A`로 읽는다. 무구동 상태와 초기 상태의 차이는 배경 변화로 별도 기록한다. 이 조합은 짝·홀 성분의 대수적 분리이며 유한 진폭의 odd 성분을 정확한 선형 미분이라고 인증하지 않는다.

분류: Counterexample candidate. 온도, 방사 속도, 내부 scalar 장의 바리온 가중 RMS와 격자 사이 벡터 차이를 비교한다. 시간 차수 `>=0.7`, 마지막 차이/응답 `<=0.2`, 응답/마지막 차이 `>=5`는 실행 전에 정했다. native 잔차·공유 면 질량-일 항등식·바리온·핵종·표본 특성속도 기준도 유지했다. 시간 차이는 엄밀한 연속 오차 상계가 아니다.

## 예산과 재현

분류: Counterexample candidate. 첫 두 물리 단계는 8.9119초, native EOS 436회로 완료됐다. 직후 pilot 요약 writer가 중복 `classification` 키로 실패했으나 두 상태와 경로 결과는 저장돼 있었다. writer만 수정하고 해당 두 단계를 재사용했다. 실행된 원 소스와 이전 계획을 별도 보존했다.

분류: Counterexample candidate. 실측 속도는 원래 252단계에 약 1123초를 요구해 600초 예산을 넘었다. 최종 응답을 확인하기 전에 실행 범위를 펄스 종료까지 126단계로 줄였다. 필요한 세 시간 간격·세 대조와 수락 기준은 유지하고 후속 꼬리 구간을 생략했다. 이때 예상은 약 561초였으며 600초 hard timeout을 걸었다. 미달 시 추가 세분화·구간 연장을 자동 실행하지 않는다.

분류: Counterexample candidate. 재현 입력·계획·이전 코드·각 단계는 `outputs/direct-eos-gr33/def-nonstationary-interior/`에 있다. 실행은 `verification/def_nonstationary_interior.py run`, 새 EOS 호출 없는 끝점 재생은 `verification/def_nonstationary_interior_audit.py`다. WSL Ubuntu-22.04의 기존 의존성 경로와 `OPENBLAS_NUM_THREADS=1`, worker 4개를 사용한다. 착수 checkpoint는 `7950a546`이다.

## 전체 목표와 다음 연결

분류: Counterexample candidate. **126단계는 완주했지만 응답의 시간 수렴 판정은 실패했다.** native 호출은 22,202회다. 실행 내부 타이머는 428.00초, `/usr/bin/time`의 wall 기록은 6분 52.91초이며 두 기록을 그대로 보존했다. RSS 최대는 130,920 KiB다. 초기 두 단계는 재사용했고 자동 추가 적분은 없다.

| 읽기 | 가장 촘촘한 odd RMS | 시간 차수 | 응답/마지막 시간 차이 | 판정 |
|---|---:|---:|---:|---|
| Jordan log 온도 | `1.32906e-7` | `2.1810` | `0.45854` | 미달 |
| 방사 속도/c | `7.54644e-9` | `−0.51352` | `2.12449` | 미달 |
| 내부 scalar | `2.06068e-6` | `0.80741` | `1.52510` | 미달 |

분류: Counterexample candidate. 이 표는 실패한 수치 경로의 값이다. 수렴된 물리 응답·상한·검출로 인용하지 않는다. 무구동 배경의 Jordan log 온도 RMS 변화는 같은 끝점에서 `0.00215299`, 속도/c는 `3.64139e-5`다. 배경 변화의 원인을 열수송만으로 단정하지 않는다.

분류: Counterexample candidate. 126개의 저장 잔차·재고·열 특성속도와 아홉 끝점의 31개 방정식을 검사했다. native 정규화 잔차 최대 `0.820408`, 총 바리온 상대 결함 최대 `7.45e-16`, 핵종 결함 최대 `7.27e-16`, 표본 특성속도 최대 `0.180993 c`다. 끝점 재생 잔차는 저장 값과 차이 0으로 일치했다. 새 EOS 호출과 추가 물질 적분은 0이다. 이 검사는 구현·이산 수지를 뒷받침하며 위 시간 수렴 실패를 구제하지 않는다.

## 빠른 scalar 전파의 원인 분리와 대체 연산자

분류: Proven. 고정된 초기 물질·계량의 유한 격자 scalar 식은 `y''=L y+f g(t/tc)`인 선형계다. `g=amplitude*(35−56 cos(2 pi tau)+28 cos(4 pi tau)−8 cos(6 pi tau)+cos(8 pi tau))/128`이므로 네 구동 oscillator를 더한 상수 행렬의 지수로 시간 격자 없이 해를 계산할 수 있다. 이는 지정된 유한 행렬계의 해법이며 원 비선형 결합계의 정확해가 아니다.

분류: Counterexample candidate. `verification/def_exact_scalar_carrier.py`에 이 전파기를 구현하고 독립 고유모드의 forced-oscillator 해와 비교했다. 상대 차이는 `1.41e-14`다. 기존 midpoint를 동일 고정 물질 연산자에 적용한 결과는 실제 결합 해의 odd scalar와 각각 상대 `7.04e-9`, `1.07e-8`, `1.23e-8`까지 일치한다. 반면 정확한 시간 전파와의 상대 차이는 `2.2244`, `1.0077`, `0.30996`이다. **관측된 큰 scalar 시간 오차의 주원인은 물질 비선형 결합이 아니라 직접 파동의 midpoint 전파다.** 이 결론은 scalar 읽기에 한정하며 온도·속도 실패의 모든 항을 분리했다고 주장하지 않는다.

분류: Counterexample candidate. 최대 구동 성분의 단계당 위상은 `4.1888`, `2.0944`, `1.0472 rad`다. 가장 거친 경로는 이 성분을 Nyquist 이하로 표본화하지 못한다. 실행 범위를 광행 시간 하나로 줄인 결정은 비용을 줄였으나, 빠른 구동을 충분히 분해하거나 느린 유체 응답을 읽는 설계를 보장하지 못했다. 이 실패를 보존한다. 정확한 직접 전파의 계산은 내부 0.109초, 새 EOS·물질 단계 0이었다.

분류: Conjectural. 다음 수정은 정확한 carrier를 단순히 기존 scalar 값에 덮어쓰는 방식이어서는 안 된다. 그렇게 하면 현재 reciprocal work와 질량-경계 일 항등식이 깨질 수 있다. 먼저 고정된 quadratic scalar Hamiltonian을 정확히 전파하고 나머지 discrete gradient를 결합하는 적분식을 유도한다. 일반적으로 `H(y)=y^T M y/2+U(y)`, `Q^T=−Q`, `A=exp(h Q M)`, `g^T(y1−y0)=U(y1)−U(y0)`이면 `y1=A y0+(A−I) M^−1 g`는 이 에너지를 보존한다. 현재 scalar 변수의 canonical 가중치, 질량 제약의 이산 미분, 물질 매개 일, 외부 구동의 일까지 그 조건을 실제로 만족하는지 확인해야 한다. 이 일반 식을 원 항성 모델의 검증된 수정이라고 세지 않는다.

분류: Conjectural. 실제 항성 응답에는 물질 자유 표면·외부 압력/대기와 단계 43의 동적 외부장을 같은 비정상 내부 해에 연결해야 한다. 이번 내부 구동은 prescribed Dirichlet 표면장이므로 나가는 복사 전하를 제공하지 않는다. 또한 실제 궤도, 같은 재고의 정적 비교 제거, 공간/EOS/경계 오차의 공동 판정이 남아 있다. [전체 완료 조건](../docs/dynamic-charge-completion.md)을 이 내부 실험의 성공으로 축소하지 않는다.
