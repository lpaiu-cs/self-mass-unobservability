# 단계 50 — 같은 자유 표면 배경의 반응·화학 에너지와 구조 응답

분류: Counterexample candidate. **새 비영 배경의 전 5,735셀에서 실제 26종 반응·중성미자 손실·조성 의존 내부에너지를 연결하고, 그 초기 원천이 만드는 자유 표면의 준정적 반경·질량 응답까지 계산했다.** 압력 모드 EOS 미분을 잘못 읽는 경로와 에너지 역산의 정밀도 불일치를 수정했다. 반응 방향 대조, 반경 공간 대조와 질량–손실 수지는 통과했다. 전체 비정상 열·유체 진화는 아직 완료되지 않았고, 작은 전하 변화는 공간 기준에 미달한다. 이번 성과는 loophole progress다.

## 무엇을 재사용하고 새로 연결했는가

분류: Imported from prior work. 단계 47의 같은 바리온·26종 조성·기준 엔트로피 비영 배경과 전 셀 native EOS, 단계 49의 수정된 자유 표면·유한 온도 단열 연산자·Just 외부를 재사용했다. 이전 고정 벽의 시간 수렴을 새 자유 표면의 진화 검증으로 옮겨 세지 않는다. [단계 49 기록](REQUEST49_FREE_SURFACE_RESPONSE_KO.md)의 2% 전하 공간 기준 미달도 유지한다.

분류: Proven. 저장 EOS 감사는 압력 모드다. `raw[5:7]`을 밀도 모드의 압력 미분처럼 사용하면 안 된다. `rP=d ln rho/d ln P|T`, `rT=d ln rho/d ln T|P`, `uP=du/d ln P|T`, `uT=du/d ln T|P`라 두면 다음 좌표 변환을 얻는다.

```text
chi_rho = 1/rP
chi_T = -rT/rP
u_lnrho = uP/rP
cv*T = uT - uP*rT/rP
cp*T = uT - (P/rho)*rT
d ln T/d ln rho|s,X = (P/rho-u_lnrho)/(cv*T)
Gamma1 = chi_rho + chi_T*d ln T/d ln rho|s,X
```

분류: Counterexample candidate. 저장 미분을 변환하고 여섯 실제 native 밀도 모드 상태와 비교했다. 최대 상대 차이는 `3.915e-13`, 전 셀 Gamma1 항등식 차이는 `1.537e-10`으로 각각 사전 `1e-7` 기준 안이다. 전 셀 비열이 양수다. 이 과정은 4.66초와 여섯 EOS 호출만 사용했다. 물리 EOS 전역 인증을 뜻하지 않는다.

분류: Counterexample candidate. 같은 저장 `rho,T,X`, 실제 자유전자량과 퇴화도를 native 반응 코드에 넣었다. 값 계산에 쓰지 않는 퇴화도 미분 인자를 `(0,0)`에서 `(+1,-1)`로 바꾼 두 전 셀 실행의 반응·열·중성미자 값이 비트 단위로 같다. 독립 profile 내보내기에서도 실제 입력이 복원됐다. 바리온 원천 상대 잔여는 `1.6704e-16`이다. 반환된 전체 반응 Jacobian을 인증하거나 사용한 것은 아니다.

## 화학 에너지를 빠뜨리지 않는 반응 방향

분류: Proven. 원자 질량 초과의 조성 변화가 만드는 정지에너지 전환과 native 내부에너지의 조성 의존성을 함께 넣어야 한다. 고정 순간 압력에서 같은 바리온 유체 원소의 원천 방향은 다음 엔탈피식으로 정의했다. `dt_J=A*N*dt`는 정지 배경의 국소 고유시간이며 두 중성미자 손실을 포함한다.

```text
X_trial = X + dt_J * R(rho,T,X)
Delta e_rest = sum_i [(W_i/A_i-1)c^2 * Delta X_i]
h_native(P,T_trial,X_trial)
  = h_native(P,T,X) - Delta e_rest - dt_J*(epsilon_nu,nuc+epsilon_nu,thermal)
```

분류: Proven. 이 식에 별도의 핵반응 Q 가열을 다시 더하면 정지에너지 전환을 중복 계상한다. 열 생성량만으로 온도를 바꾸는 근사는 `u_X Delta X`를 누락한다. 이번 native 엔탈피 역산은 그 조성 항을 포함한다.

분류: Counterexample candidate. 8초·4초는 초기 반응 벡터의 EOS 방향 대조 간격이다. 그동안 반응률·유체·열유속을 함께 갱신한 실제 시간 진화가 아니다. 각 셀의 에너지 잔여 기준은 `max(2 erg/g, 32 ulp(h_target), 1e-8*|Delta e_rest+loss|)`를 유지했다. 전 5,735셀의 최대 기준 점수는 `0.991425`, 질량 가중 밀도·에너지 방향 차이는 각각 `0.00171548`, `0.00171386`으로 사전 1% 기준을 통과했다. 반응 원천과 화학 내부에너지의 연결 병목은 이 범위에서 해소됐다.

## 반응 원천에서 자유 표면의 반경·질량·전하까지

분류: Proven. 엔트로피·조성 변화가 있을 때 단계 49의 단열 질량 대수식을 그대로 사용하면 일반적인 질량 변화를 놓친다. `zeta=xi/r`, `eta=Delta p/p`, `f=Delta phi`, `V=Delta Phi`, `w=e+p`와 고정 압력 반응 방향 `rho_ref`, `e_ref`를 쓰면 다음 준정적 Lagrangian 변분을 얻는다. 소스는 좌표시간당 값이고 공간 좌표·질량은 단계 49와 같은 기하 단위다.

```text
Delta ln rho = eta/Gamma1 + rho_ref
Delta e = w*eta/Gamma1 + e_ref
e_ref = w*rho_ref - loss              [loss = rho*A*N*epsilon_nu]
Delta lambda = (Delta m/r - m*xi/r^2)/b
xi' = -eta/Gamma1-rho_ref-2*zeta-Delta lambda-3*alpha*f
(Delta p)' = -g*(Delta e+Delta p+w*xi') - w*Delta g
f' = V+xi'*Phi
V' = Delta F+xi'*F
(Delta m)' = Delta M+xi'*M           [M=m']
```

분류: Proven. 위 질량식에서 큰 가역 일 항을 직접 빼면 작은 실제 손실을 잃을 수 있다. 다음 독립 질량 결함 `J`를 쓰면 소거를 해석적으로 끝낼 수 있다. 기호 대조로 원 변분과의 항등식을 확인했다.

```text
H = r^2*b*Phi*f - (4*pi*r^2*A^4*p+r^2*b*Phi^2/2)*xi
J = Delta m-H
J' = -(4*pi*r*A^4*w/b+r*Phi^2)*J - 4*pi*r^2*A^4*loss
(N*a*J)' = -4*pi*r^2*A^4*N*a*loss
```

분류: Proven. 고정 `phi_infinity`, 초기 열유속 0의 준정적 접선에서는 이 적분식이 `dot M_geom=-(G/c^4)*1e-7*sum(dm*A^2*N^2*epsilon_nu)`를 준다. 여기서 `dm`은 g, `epsilon_nu`는 erg/g/s, 결과는 m/s다. 핵 정지에너지가 내부에너지로 바뀌는 부분은 별도의 총질량 생성이 아니다. 이 선형 접선의 항등식은 비선형 시간 경로나 그 오차 상계가 아니다.

분류: Counterexample candidate. 중심 정칙성, 자유 표면 압력 정칙성, 같은 `phi_infinity`의 정확한 Just 외부를 적용했다. 무원천·단위 외부장 대조는 기존 단열 정적 변위와 `7.24e-14` 상대 차이로 일치한다. 동차 질량 결함은 이론적으로 정확히 0이므로 저장된 반올림 잔여만 제거하고 **원 행렬의 모든 행을 같은 기준으로 재검사**했다. 그 대조의 잔차는 `2.901e-16`, 반응 경로 최대 잔차는 `4.114e-16`이다.

분류: Counterexample candidate. 다음 수치는 원 5,863개 응답 구간과 그 두 배 격자, 동일 배경 보간에서 얻은 **초기 준정적 변화율**이다. 시간 진화에서 실제로 측정한 표면 속도나 보편적인 물리 오차로 해석하지 않는다.

| 항목 | 값 또는 차이 | 사전 판정 |
|---|---:|---|
| 표면 반경 변화율 | `1.86074e-4 m/s` | 준정적 접선 |
| 반경 변화율의 격자 상대 차이 | `2.96378e-5` | 2% 기준 통과 |
| 반경 변화율의 8초/4초 원천 차이 | `5.21091e-5` | 1% 기준 통과 |
| ADM 기하 질량 변화율 | `-4.25315e-17 m/s` | 준정적 접선 |
| 중성미자 제1법칙 예상값 | `-4.25377e-17 m/s` | 원 재고 합으로 독립 평가 |
| 질량 수지 상대 차이 | `1.45341e-4` | 2% 기준 통과 |
| 정규화 전하 변화율 | `-1.19891e-19 /s` | 미수락 후보값 |
| 전하 변화율의 격자 상대 차이 | **`22.6268%`** | **2% 기준 미달** |

분류: Counterexample candidate. 정규화 전하에는 Just 외부의 scalar 전하와 ADM 질량 분모의 변화를 모두 넣었다. 작은 값이 나왔다는 이유로 전하를 검출하거나 no-go로 확정하지 않는다. 원 광학 경계 밖의 반응 원천은 기준 계산에서 0으로 선언했고, 원천의 별도 연장 대조도 수행했다. 이번 반경·전하 차이는 표시 정밀도에서 0이지만 두 처방의 일치는 물리적 외층 반응 상계가 아니다.

## 실패와 자원 사용을 보존한 방식

분류: Counterexample candidate. 최초 32개 비연속 셀/광도 0 carrier는 native profile 복원 전에 실패했다. 양의 광도를 요구한 후속 guard도 원 저장 광도 435개가 음수여서 실패했다. 음의 저장 광도 자체는 허용된다. 전 5,735셀 재고와 저장된 부호 있는 광도를 함께 복원한 carrier가 통과했다. 두 변경을 같이 했으므로 최초 실패의 단일 원인을 광도 0으로 확정하지 않는다.

분류: Counterexample candidate. 최초 엔탈피 역산, multiprocessing 함수 식별과 200초 예상 예산 기준의 실패를 보존했다. 따뜻한 worker의 실측 뒤 두 반응 방향·CPU 8개·BLAS 1개·GPU 없음·최대 600초 한 번으로 예산을 재평가했다. 해당 실행은 셀 3968에서 에너지 잔여 `-60.3110 erg/g`가 원 `32 erg/g` 기준을 넘어 중단됐다. 완료 블록은 5,483셀을 보존했다.

분류: Proven. 기존 역산은 후보 선택에 binary64 엔탈피와 반올림 target을 쓰고 최종 수락에는 확장 정밀도 엔탈피와 원 target을 썼다. 따라서 두 검사가 다른 후보를 최선으로 판단할 수 있다. 반복·선택·수락의 잔차식을 같게 하고 제한된 인접 표현값 탐색을 사용했다. 에너지 수락 문턱은 바꾸지 않았다.

분류: Counterexample candidate. 실패 블록의 남은 252셀만 90초 한도로 재개했고 9.43초·940회 native EOS 호출로 완료했다. 수락된 방향 레코드 전체에 기록된 호출 합은 24,909회다. 실패 셀의 저장되지 않은 호출과 이전 pilot 호출은 그 합에 포함되지 않으므로 전체 소비량으로 부르지 않는다. 두 전 셀 source 실행은 약 20.05초, 각 구조 대조 묶음은 약 3초였다. 장기 항성 궤도 적분, 격자 자동 확대 또는 GPU 전환은 하지 않았다.

분류: Counterexample candidate. 처음 독립 질량식의 직접 이산화는 가역 일 소거 때문에 선형 잔차 약 `5e-8`을 남겼다. 질량 수지 보고에서도 에너지 밀도 변환 `0.1`을 luminosity 변환에 사용한 단위 오류가 있어 `1e-7`로 정정했다. 보존 질량 결함식과 확장 정밀도 잔차 보정이 후속 후보를 구성한다. 무원천 `J=0`의 반올림 수치에 상대 잔차를 정의해 1이 나온 원 기록, species 목록 JSON을 객체로 가정한 감사기 오류도 보존한다. 후속 감사기는 자료형을 구별하고, 물리 계산은 재실행하지 않았다.

## 현재 닫힌 연결과 남은 직접 병목

분류: Counterexample candidate. **닫힌 부분은 같은 비영 배경의 native 반응 → 조성·화학 내부에너지 → 바리온 보존 준정적 구조 → 자유 표면 반경·ADM 질량이다.** 164개 실행 계획 결속을 확인했고 에너지 방향·반경 공간·질량 손실 수지 대조를 통과했다. 작은 정규화 전하까지 통과했다고 표시하지 않는다.

분류: Proven. LTE 내부 확산을 쓰면 `Theta=A*N*T_J`, `K_E=A^2*K_J`, `L_infinity=-4*pi*r_E^2*N*A^2*K_J/a_E*dTheta/dr_E`다. 같은 면 유속을 쓰는 저장 구적은 내부 에너지를 망원경 합으로 상쇄한다. 그러나 이는 선언한 내부 수송식이며 희박한 자유 표면의 실제 복사 조건을 결정하지 않는다.

분류: Conjectural. 다음 직접 병목은 **같은 자유 표면 배경에서 관성·반응·열유속·물리 복사 경계를 갖춘 유한시간 경로를 완결하는 것**이다. 현재 초기 열유속은 0이므로 초기 에너지 발산은 0이고 수송은 열유속의 시간 미분으로 시작한다. 확산 정상 유속을 초기 Cauchy 유속에 대입해서는 안 된다. Rosseland 평균을 Planck 흡수율로 바꾸거나 외층에 LTE를 자동 연장하지 않는다. 현재 외곽 경계에 의존하는 큰 냉각률은 물리 냉각시간으로 사용하지 않는다.

분류: Conjectural. 전하 공간 기준 실패에는 배경·섭동의 보존형 이산 일치와 작은 잔여의 읽기를 먼저 점검한다. 기준 미달을 이유로 격자나 적분 기간을 자동 확대하지 않는다. 실제 조화 구동, 동일 재고 정적 비교, 미분 nuisance, EOS·수치·경계 결합 오차 및 관측 연결은 계속 미완료다.

분류: Imported from prior work. 재현 진입점은 `verification/def_free_surface_thermal.py`, `verification/def_thermal_native_inventory.py`, `verification/def_reactive_energy.py`, `verification/def_reactive_structure_balanced.py`, `verification/def_thermal_reactive_audit.py`다. 실패한 선행 구현도 보존한다. 실행 입력·결과·중단 블록·방향 레코드는 `outputs/direct-eos-gr33/def-free-surface-thermal/`에 있으며 최종 상태는 `completion.json`, 결속은 단계 50 milestone manifest를 따른다. 이미 존재하는 결과 디렉터리는 덮어쓰지 않는다.
