# 단계 34 — GR 계산에서 실제 구동·응답으로 (2026-09-17)

## 판정과 연구 결정

분류: Counterexample candidate. 승인된 전체 GR 계산의 수렴 통과는 그대로 보존한다. 이번 첫 마일스톤은 새 장기 적분 없이 저장된 초기·끝점에서 읽을 수 있는 질량·전하 응답을 판정하고, 다음 실험의 실제 scalar 구동을 정의하는 것이다. 기존 70/140/280 경로를 다시 계산하지 않았다.

분류: Proven. 고정 반사 경계에서 외곽 에너지 유속이 0인 순수 구면 GR 보존계는 내부 열 이동만으로 총 질량 변조를 만들 수 없다. 공유 면 유속의 합이 경계로 망원경 합을 이루고, 가변 시간 BDF 계수의 합이 0이므로 일정한 초기 총 질량이 이산 방정식에서도 유지된다. 분류: Counterexample candidate. 실제 저장된 세 끝점의 총 질량 변화와 누적 외곽 에너지 유속은 저장 정밀도에서 정확히 0이다. 이를 새로운 질량 신호로 해석하지 않는다.

분류: Proven. 저장된 네 상태에서 정의한 **정확한 binary64 구간별 3차 보간 계수 모형**에서, beta=-4인 정적 scalar 읽기의 끝점 변화는 다음과 같이 제한된다. 표는 `alpha/phi_infinity=-chi/M`의 변화 절댓값에 대한 상계를 바깥쪽으로 올림한 값이다.

| 원 GR 경로 | 초기 대비 끝점 변화의 상계 |
|---|---:|
| 70단계 | 5.626e-13 |
| 140단계 | 9.134e-13 |
| 280단계 | 7.927e-13 |

분류: Proven. 세 상계 모두 이번 읽기 실험 전에 정한 1e-9 척도보다 작다. 이 결과는 해당 **세 끝점**의 판정이다. 중간 시각의 최대 변화, 다른 구동, 더 긴 시간, 보간 계수를 생성한 GR/EOS 오차, 실제 관측 민감도를 제한한 것은 아니다. 더 정밀한 새 관측 목표를 배제하지 않으며 이 상계를 전체 물리 모형의 오차 보장으로 확대하지 않는다.

연구 결정: 기존 자유 열 이력의 끝점 변화를 동적 chi 신호로 승격하거나, 이 결과를 키우려고 같은 GR 적분 기간·격자를 자동 확대하지 않는다. 다음 계산의 대상은 **비영 scalar 구동에 의해 실제로 유발되는 물질·계량·전하 변화**다. 새 분자 EOS 연결은 이 구동 모형에 맞춰 수행한다.

## 보존한 실패와 상계의 근거

분류: Counterexample candidate. 먼저 기존 scalar IVP와 독립 Riccati 식으로 저장된 끝점을 읽었다. 두 ODE 허용오차 1e-10/2e-12와 6/12점 구적을 미리 고정했다. 큰 정적 계수를 각각 계산해서 빼는 방법은 최대 대조 차이 1.2594e-7로 원 1e-9 문턱을 통과하지 못했다. 직접 차분의 부호와 크기를 물리 신호로 사용하지 않는다. 실패 결과는 `gr-driven-readout`에 보존했다.

분류: Proven. 후속 계산은 ODE 허용오차를 줄이는 대신 기존 Wronskian 항등식과 Volterra 상계를 사용했다. 반지름·질량과 scalar 운동 계수 a(x)의 저장 보간식이 네 상태에서 정확히 같은지 먼저 확인했다. b(x)의 차이만 정확한 유리수로 계산하며, Bernstein 계수로 각 구간의 부호·최댓값을 둘러싼다. 여기서 a(x)는 scalar 방정식의 계수 N/sqrt(g_rr)이며 GR 코드의 metric a와 구분한다.

분류: Proven. 식 `(x^2 a phi')'=x^2 b phi`에서

```text
B1 = integral_0^1 x |b| dx,   B2 = integral_0^1 x^2 |b| dx
kappa = (B1-B2)/a_min
L = -log(1-2 mu)/(2 mu)
F = 1/(1-2 kappa-L B2)
```

분류: Proven. Volterra 연산자의 최대노름은 kappa 이하이며, kappa<1일 때 단위 중심 해의 노름은 1/(1-kappa) 이하이다. 외부 정규화의 하한은 `(1-2 kappa-L B2)/(1-kappa)>0`이므로 무한대에서 1로 정규화한 장의 노름은 F 이하이다. 두 해의 Wronskian 항등식에서 동일한 a, mu, R/M를 사용하면

```text
|delta(alpha/phi_infinity)|
 <= (R/M) F_initial F_final integral_0^1 x^2 |delta b| dx.
```

분류: Proven. 실제 최대 kappa는 약 1.500927e-4다. 다항식 모멘트와 Bernstein 계산은 유리수로, 로그·정규화는 60자리 바깥 방향 구간 연산으로 처리했다. 변화가 0인 입력과 해석해를 아는 비영 일정 퍼텐셜 대조를 포함한다. 정확한 모형 계수에 대한 상계이며 반올림 전 계수나 EOS의 물리적 정확도를 포함하지 않는다.

## 실제 구동을 연결할 때 필요한 차수

분류: Imported from prior work. 기존 후보 이론은 massless DEF의 `A(phi)=exp(beta phi^2/2)`다. 영 scalar의 비스칼라화 가지에서 선형 scalar와 유체 섭동이 분리되는 경계는 [Khalil et al.의 IV.2](https://arxiv.org/html/2206.13233v2#S4.SS2)와 기존 단계 13 유도에 따른다. 논문의 특정 중성자별 수치 계수를 현재 백색왜성에 가져오지 않는다.

분류: Proven. Einstein-frame 물질 방정식은 `nabla_mu T^mu_nu=beta phi T d_nu(phi)`이고 열유속과 반지름 방향 속도를 포함해도 `T=-epsilon+3P`다. Eulerian 에너지 E를 comoving epsilon 대신 넣으면 안 된다. `phi=phi0*(psi0+eta delta_psi)`를 대입하면 eta에 선형인 구동은

```text
delta S_nu = beta phi0^2 T
             * (psi0 d_nu(delta_psi) + delta_psi d_nu(psi0)).
```

분류: Proven. 비스칼라화 가지가 매끄럽고 phi0→0에서 전하 feedback과 응답 연산자가 유계일 때, 동반천체가 만드는 `delta_phi_ext=phi0 delta_u`에 대해 열·유체를 거치는 부분은 다음 차수다.

```text
delta U           = O(phi0^2 delta_u)
delta q_thermal   = O(phi0^3 delta_u)
delta force/g     = O(phi0^4 delta_u)
```

분류: Proven. 반면 phi0와 무관하게 외부에서 준 열 변화 delta U의 scalar 힘 읽기는 O(phi0^2 delta U)부터 시작한다. 따라서 자유 열 이력을 동반천체 scalar 구동의 전달함수로 바꾸면 두 차수를 누락한다. 이는 점근 차수이며 크기의 보편 상계는 아니다. `phi0^2/(s+lambda)`에서 `s=0, lambda=phi0^2`이면 억제가 사라지는 대조를 포함했다. 임계 가지·공명·독립적인 조석/가열 구동을 이 결론으로 배제하지 않는다.

분류: Proven. 보존 EOS 연결에는 frame 변환도 필요하다. 보존 바리온 밀도와 온도는 `n_E=A^3 n_J`, `T_E=A T_J`, 자유에너지 밀도는 `f_E=A^4 f_J(n_E/A^3,T_E/A,X)`로 연결된다. 임의의 미분 가능한 f_J에 대해

```text
partial f_E / partial ln A |_(n_E,T_E,X) = epsilon_E - 3P_E.
```

분류: Proven. 따라서 핵 정지에너지를 포함한 바리온당 전체 에너지는 A에 비례해 변한다. 정지에너지를 고정한 채 scalar 힘만 덧붙이는 구현은 이 자유에너지와 일치하지 않는다. A=1 복원, trace가 0인 복사 기체, trace가 비영인 차가운 정지질량 기체를 상징 대조로 확인했다. 수송 구성식 전체의 변환이나 결합 진화 구현을 완료한 것으로 세지 않는다.

## 다음 실행 경계

분류: Conjectural. 다음 마일스톤은 위 EOS frame 변환·물질 교환항·scalar 응력의 계량 항을 같은 차수로 연결하고, 비영 구동과 무구동의 **짝 비교**를 원 보존 변수에서 수행하는 것이다. 아직 진화하는 배경을 정상 배경으로 가정하지 않는다. 정상성이 입증되기 전에는 두 시각에 의존하는 응답 커널을 대상으로 한다. 비교 모형은 사전에 정한 순간 응답 및 물리적으로 정당화한 미분·nuisance 방향이며, 단일 위상 지연만으로 성공 판정을 하지 않는다.

분류: Counterexample candidate. 새 분자 EOS 초기화는 저장 블록 기준 1792/5735셀이고 최종 initial.npz는 없다. 기존 분자 진화 연결은 수정 전 composition-tangent 경로를 사용하므로, 이번 검증을 통과한 호환 압력·정수압 및 donor 보정 경로에 연결해야 한다. 남은 셀·근 기록을 재사용하되 새 EOS에 이전 모델의 수락 시간 상태를 가져오지 않는다.

새 장기 계산 진입 전에는 위 교환항의 에너지 상쇄, EOS frame 복원, 영 배경의 선형 유체 응답 0, 구동 부호/진폭 대조를 작은 대표 상태에서 확인한다. 필요한 응답 정확도와 원 보존 문턱을 동시에 만족할 때만 대표 구간의 실제 비용을 측정해 계산 범위를 정한다. 아직 그 비용을 측정하지 않았으므로 다음 전체 진화의 완료 시각은 제시하지 않는다.

## 실행·재현

단일 CPU, GPU 미사용, native EOS 호출 0. 최초 읽기는 벽시간 8.78초·최대 RSS 191520 KiB, 유리수/구간 상계의 내부 측정은 153.05초다. 각 실행에 300초 외부 timeout을 적용했다. 원 적분은 재실행하지 않았다.

```text
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification OPENBLAS_NUM_THREADS=1
python3 verification/gr_driven_readout.py prepare
timeout 300 python3 verification/gr_driven_readout.py run
# 위 실행은 보존된 numerical_controls_passed=false로 종료한다.
python3 verification/gr_driven_readout_bound.py prepare
timeout 300 python3 verification/gr_driven_readout_bound.py run
python3 verification/def_thermal_drive_selection.py
python3 verification/def_thermal_drive_selection.py --frame-check
```

출력 폴더가 이미 있으면 prepare/run을 덮어쓰지 않는다. 결과·코드 결속은 `outputs/direct-eos-gr33/gr-driven-milestone-manifest.json`에 둔다. 원 GR 결과 및 이전 실패 판정은 변경하지 않는다.

분류: Proven. 이번 성과는 선언된 끝점 읽기의 조건부 상계와 실제 구동의 결합 차수를 확정한 **theorem progress**다. 분류: Conjectural. 실제 구동의 열 모드 잔여 응답, 새 EOS 결합 진화, 물리 EOS·대기/외부 및 비선형 관측 폐쇄는 다음 작업으로 남는다.
