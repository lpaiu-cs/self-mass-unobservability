# Request 37 — 회전점에서 정칙한 scalar–계량 보존 적분

작성: 2026-09-17. 선행 단계: [Request 36](REQUEST36_SPHERICAL_COUPLING_KO.md).

## 판정과 범위

분류: Proven. **운동량으로 나누지 않는 구면 진공 scalar–계량 적분식을 얻었다.** 이산 질량 제약식의 두 상태 차분과 방사 방향 수반을 사용한다. 양의 정칙 가지에서 운동량이 0이거나 상태 변화량이 0이어도 계수가 매끄럽고, 정확히 푼 갱신은 외곽 질량을 보존한다. 고정 유한 격자의 매끄러운 Hamiltonian 계에 대한 시간 2차 합치성을 아래 조건으로 증명한다.

분류: Counterexample candidate. 수정 공간 구적을 쓰는 8개 작은 파동 경로가 원 수락 기준을 통과했다. 비선형 자기중력 파동의 시간·공간 차수는 두 scalar 성분 모두 약 2이며, 정확히 모든 운동량이 0인 시작 상태와 이후 부호 전환을 포함한다. 저장 끝점으로 질량·오차·차수를 재계산해 같은 결과를 얻었다. 이는 **theorem progress와 진공 scalar 부문의 loophole progress**다.

분류: Conjectural. 새 적분기를 기존 5,735셀 열유체 항성에 적용한 결과는 아니다. 비평형 열유속의 정준 상태와 엔트로피 폐쇄, 새 EOS 전체 시간 수렴, 균일한 비선형 PDE 오차 상계, 전체 Einstein 제약의 독립 수렴 검사, 외부 구동과 관측 연결은 남는다. 과거 native 결합 첫 단계의 `1/Pi_mid` 보정이 소급해 증명된 것도 아니다.

## 정의와 정확한 에너지 차분

분류: Proven. 다음은 선언한 이산 모형의 정의다. `G=c=1`, DEF 정규화 `R-2(∂phi)^2`, `ds²=-N²dt²+a²dr²+r²dΩ²`, `b=a^-2`, `Pi=(a/N)phi_t`를 사용한다. 중앙 질량은 0이고 외곽 `phi=0`을 고정한다. 보통 물질과 외부 에너지 유입이 없는 닫힌 구면 영역이다. 셀 중심 질량은 `m_i=(1-f_i)m_L+f_i m_R`이다.

```text
K_i = w_i Pi_i²/2
e_f = face_weight_f Phi_f²/(8π),  Phi_f = 차분(phi)/거리
L_i = e_left/2,  R_i = e_right/2  (외곽 면만 R=e_outer 전부)
m_R - m_L = K_i b_i + L_i b_L + R_i b_R
C_i = 1 - 2 K_i(1-f_i)/r_i - 2 L_i/r_L
D_i = 1 + 2 K_i f_i/r_i + 2 R_i/r_R
m_R = (C_i m_L + K_i+L_i+R_i)/D_i
```

분류: Proven. 중앙 면의 `L=0` 항은 정확히 0으로 정의하여 `0/0`을 만들지 않는다. 두 시간 상태의 산술 평균을 위줄로 표시하면 곱셈 차분 항등식에 의해 정확히

```text
Dbar_i Δm_R - Cbar_i Δm_L
  = bbar_i ΔK_i + bbar_L ΔL_i + bbar_R ΔR_i
lambda_outer=1
H_i=lambda_(i+1)/Dbar_i,  lambda_i=lambda_(i+1) Cbar_i/Dbar_i
k_i=H_i bbar_i
kf= bbar_f (H_left+H_right)/2  (외곽은 bbar_outer H_last)
```

가 성립한다. 방사 방향 소거와 공유 면의 이산 부분적분을 적용하면 외곽 질량 `M=m_outer`의 정확한 두 상태 gradient는

```text
G_phi = -w div(kf Phi_mid)
G_Pi  =  w k Pi_mid
ΔM = Σ(G_phi Δphi + G_Pi ΔPi)
Δphi = dt k Pi_mid
ΔPi  = dt div(kf Phi_mid)
```

이다. 마지막 두 식을 대입하면 `ΔM=0`이다. 보존을 위해 운동량, 장의 변화량, 에너지 결함으로 나누는 과정이 없다. 비선형 반복의 유한 종료 오차는 별도의 잔차와 보존 문턱으로 검사한다.

분류: Proven. 고정 격자에서 `C>0`, `b>0`, `bf>0`인 열린 가지를 택한다. `D>=1`이고 모든 반지름과 가중치는 고정이므로 위 gradient는 이 가지 안에서 매끄럽다. 두 상태를 바꾸어도 동일하고, 두 상태가 같으면 질량 제약을 직접 미분한 실제 gradient와 일치한다. 갱신 잔차의 새 상태 Jacobian은 `dt=0`에서 항등행렬이므로 암시함수정리는 충분히 작은 `dt`의 유일한 근접 해를 준다. 이 해의 시간 대칭성과 대각선 합치성으로 국소 오차는 `O(dt³)`이며, 매끄러운 해가 정칙 가지의 compact subset에 머무르면 고정 격자의 전역 시간 오차는 `O(dt²)`이다. 큰 시간 간격에서 반복 수렴을 보장하거나 격자 크기에 균일한 PDE 오차 상수를 주는 정리가 아니다.

## 공간 경계의 실제 수정

분류: Counterexample candidate. 첫 후보는 구각 체적과 `r_face²`를 gradient 가중치로 썼다. 질량 보존과 회전점은 통과했지만 flat 정확해의 운동량 차수는 1.55–1.58, 비선형 시간 운동량 차수는 0.051로 실패했다. 평탄 시공간의 구면 sine 모드에 대한 공간 연산자 결함은 0.51537이었다. 이 실패 판정과 원 소스·계획·결과를 보존한다.

분류: Proven. 새 후보는 `u=r phi`, `p=r Pi`의 중점 운동 에너지 구적과, 두 중심 사이에서 선형인 `u`의 정확한 gradient 에너지를 사용한다.

```text
w_i = r_i² Δr,   f_i=1/2
area_f = 4π r_left r_right  (외곽은 4π r_last R)
face_weight_f = area_f * distance_f
∫[r_left,r_right] r² phi_r² dr
  = r_left r_right (phi_right-phi_left)²/(r_right-r_left)
```

분류: Proven. 중앙부터 첫 중심까지 `u`가 원점을 지나는 선형이면 `phi`는 상수이고 gradient 에너지는 0이다. 균일 격자에서 이 정의의 flat 연산자는 `u`의 셀 중심 2차 차분을 `r_i`로 나눈 것이며, 양 끝의 홀대칭 ghost가 각각 정칙 중앙과 외곽 Dirichlet 조건을 구현한다. 따라서 `phi_i=sin(πr_i)/(πr_i)`는 양 경계까지 포함해 고유값 `-4sin²(πΔr/2)/Δr²`를 갖는다.

분류: Proven. 매끄러운 내부 구간에서는 수정 질량 제약의 잔차를 `Δr`로 전개하면

```text
F = Δr [m_r - r² b (Pi²+Phi²)/2] + O(Δr³)
log(C/D) = -r Δr (Pi²+Phi²) + O(Δr³)
```

이다. 수반 극한은 `H=N a`, `H(R)=1`, `(log H)_r=r(Pi²+Phi²)`이며, 장의 극한 식은 `phi_t=Hb Pi`, `Pi_t=r^-2 ∂r(r²Hb Phi)`다. 내부 전개와 dual-cell gradient 에너지 항등식을 기호 검산했고, flat 경계 고유함수 관계는 위 유도와 수치 대조로 확인했다. 비선형 경계 전체에 대한 균일 오차 증명은 별개다.

분류: Imported from prior work. 연속 구면 scalar 파동식과 fully constrained radial evolution의 기준은 [Gundlach–Martín-García, §3.1](https://link.springer.com/article/10.12942/lrr-2007-5)이다. 해당 문헌의 canonical scalar와 여기의 DEF scalar는 `phi_DEF=sqrt(4πG) phi_canonical`로 정규화가 다르며, lapse의 시간 게이지도 구분한다. [Gonzalez의 conserving time integration 자료](https://web.ma.utexas.edu/users/og/numerics.html)는 discrete-gradient 방법의 선행 맥락이다. 이 보고서의 구체적인 질량 수반 구성과 검산을 외부 문헌의 정리로 대신하지 않는다.

## 동결한 대조와 결과

분류: Counterexample candidate. 초기 `phi=.05 sinc(r)`, `Pi=0`, 영역 `[0,1]`, 종료 시각 1.7을 고정했다. flat 정확해 공간 격자 32/64/128 및 각 4배 단계 수, 비선형 시간 대조 64셀·128/256/512단계, 비선형 공간 대조 32/64/128셀·각 4배 단계 수를 사용했다. 겹치는 한 경로를 재사용해 총 8개다. 첫 실패 뒤에도 진폭·기간·격자·차수·오차 문턱을 바꾸지 않았다.

분류: Counterexample candidate. 아래 값은 수정 공간식의 측정 결과다.

| 항목 | 측정 | 원 수락 기준 |
|---|---:|---:|
| flat 정확해 phi 차수 | 2.00076 / 2.00019 | >=1.7 |
| flat 정확해 Pi 차수 | 1.99630 / 1.99908 | >=1.7 |
| 비선형 시간 차수 phi / Pi | 2.00002 / 1.99882 | >=1.8 |
| 비선형 공간 차수 phi / Pi | 2.00094 / 1.99527 | >=1.6 |
| 최대 끝점 상대 질량 차이 | 1.10203e-18 | <2e-12 |
| 최대 단계 상대 질량 증가, 안정 차분 | 2.59054e-20 | <5e-14 |
| 최대 장 방정식 잔차 / 진폭 | 3.32037e-17 | <2e-15 |
| 경로별 운동량 부호 전환 최소 횟수 | 32 | >0 |
| 128셀 비선형 경로 최대 compactness | 약 0.00152 | >1e-5 |

분류: Counterexample candidate. 모든 운동량이 0인 시작 상태에서 첫 갱신, 정확한 영 상태 유지, 정방향 뒤 역방향 시간 갱신도 통과했다. 두 상태 gradient 교환 대칭은 배열 수준에서 같았고, secant 상대 결함은 3.264e-20이었다. flat 고유함수 결함은 8.102e-14로 감소했다. 이 마지막 수치는 binary64에서 가져온 π와 sinc 인자의 반올림을 공간 2차 차분이 증폭한 결과를 포함한다.

분류: Counterexample candidate. 수정 후보의 첫 실행기는 복제한 함수의 기본 인자를 전달하지 못해 첫 경로 이전에 실패했다. 실패 파일을 보존하고 `__defaults__`만 복원하는 실행기를 별도로 결속했다. 성공 실행의 파동 경로 합계 벽시간은 2.666초, 프로세스 전체는 3.69초, 최대 RSS는 82,036 KiB였다. 실패한 첫 공간식 경로도 2.815초였으며 CPU 1개·BLAS 1스레드를 사용했다. 새 native EOS 호출·장기 항성 이력 재계산은 0회다.

분류: Counterexample candidate. 소스·계획·실패와 저장 끝점의 SHA를 확인하고, 8개 끝점에서 질량 제약, flat 정확해 오차 및 시간·공간 차수를 다시 계산하여 저장 결과와 정확히 일치했다. 중간 궤적 전체를 독립 재생한 것은 아니다. 중간 단계 잔차·부호 전환은 실행 중 집계한 기록에 근거한다.

## 실제 열유체로 넘어가기 위한 정확한 조건

분류: Proven. 보통 물질에서 `W²=(1-v²)^-1`, `w=epsilon+P`, proper 열유속 `Q`를 두고

```text
E = (epsilon+P v²+2Qv) W²
S = (w v+Q(1+v²)) W²
R = (epsilon v²+P+2Qv) W²
B = a rho W,  J = a² S
```

라고 하자. 고정 조성·평형 비엔트로피에서 native first law `d epsilon=w d ln rho-alpha T_trace dphi`를 가정한다. `B,J`를 고정하고 `ell=d ln a`라 하면 직접 미분으로

```text
dE + (E+R) ell + alpha T_trace dphi
  = v dQ + 2Q W² dv + 2Qv ell
```

를 얻는다. 따라서 `Q=0,dQ=0`에서는 물질 질량 제약의 계량 미분이 lapse의 `E+R` 항을 공급하고, scalar 미분은 reciprocal trace source를 공급한다. 여기서 고정할 정준 운동량은 `a²S`이며 기존 native 저장 변수 `aS`를 그대로 고정하는 것과 다르다. 이 결과는 **완전유체의 조건부 해석적 연결**이고 물질까지 포함한 유한 단계 적분기의 구현·보존 증명은 아니다.

분류: Proven. `Q≠0`인 기존 열수송 상태에서 평형 entropy와 proper `Q`를 무조건 고정하는 확장은 실패한다. `v=0`, `dQ=0`, `J` 고정이면 `dv=-2Q ell/w`여서 위 우변은 `-4Q² ell/w`다. 따라서 필요한 `dE/d ln a=-(E+R)`가 성립하지 않는다. 이는 해당 변수 고정 가정의 정확한 반례이며 모든 열유체 Hamiltonian에 대한 no-go가 아니다.

분류: Conjectural. 다음 목표는 **기존 quadratic-entropy 열수송의 비평형 엔트로피·열 상태를 포함해 위 추가 항을 처리하는 정준 결합**이다. 먼저 연속 수준의 계량·scalar 교환과 정칙한 `v=0,Q≠0` 변분을 맞춘 뒤 이산 에너지 교환을 구성해야 한다. 이 유도와 작은 국소 대조가 끝나야 저장한 native 초기 상태의 새 결합 첫 단계를 검증할 근거가 생긴다. 단순히 `Q`를 고정한 채 진공 kernel을 장기 항성 적분에 넣지 않는다.

## 실행 및 증거

검산: `OPENBLAS_NUM_THREADS=1 PYTHONPATH=verification:<기존 request13_deps> python3 verification/def_scalar_hamiltonian_proof.py`. 이미 결과가 있으면 `algebra()`, `model.check()`, `replay()`를 호출하여 읽기 전용으로 다시 검사한다. `prepare/run`은 동결된 산출물을 덮어쓰지 않는다.

소스: `verification/def_scalar_hamiltonian.py`, `def_scalar_regular_radius.py`, `def_scalar_regular_radius_run.py`, `def_scalar_hamiltonian_proof.py`. 첫 실패는 `outputs/direct-eos-gr33/def-scalar-hamiltonian/`, 수정 계획·실행기 실패·성공 결과는 `def-scalar-regular-radius/` 및 그 아래 `execution-fixed/`에 있다. 결속은 `gr-regular-scalar-milestone-manifest.json`에 기록한다. Phase36 문서 snapshot은 commit `05dd4e57`에 남으며 과거 판정을 수정하지 않는다.
