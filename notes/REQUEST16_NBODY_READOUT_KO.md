# 4체 관측 지연과 연속 보간 오차 연결

분류: Conjectural. 이번 단계는 실제 4체 운동에 대응하는 추가 천체의 Einstein·Shapiro 항을 공통 관측 경로에 연결하고, 자연 3차 스플라인의 값·미분·적분 오차를 잔차 기반으로 제한하는 것이다. 초기화 전체 미분과 전 관측 기간 변분 인증은 별도 미충족 조건으로 유지한다.

분류: Imported from prior work. Request 15는 초기화의 두 결함을 고쳤고 원래 물리 초기 상태·스핀을 재현하는 좌표를 얻었다. 이 수정본과 보정 매개변수를 기준으로 사용한다. Request 10·12·13·14·15의 데이터·코드·판정은 동결하며 신규 결과와 구분한다.

분류: Counterexample candidate. GR의 같은 1PN 근사 아래 추가하는 항은 `dE_extra/dt = GM_extra/(c² r_p,extra)` 및 `S_extra = -2 T_sun M_extra log[(r_p,extra − n·r_p,extra)/c]`이다. 코드의 시간 단위는 일이며 Shapiro 계수는 초에서 일로 변환한다. 이 확장은 GR 관측 기준선의 누락 항을 다루며 동적 scalar 힘·광자 전파를 도출한 모형은 아니다. 수차 스핀 축은 기존의 세 천체 각운동량 정렬 규약으로 명시적으로 유지한다.

분류: Conjectural. 실행 전 검증 기준은 (1) 추가 질량 0의 동일성, (2) 합산·회전 공변성·로그 영역과 충돌 거부, (3) 실제 입력에서 독립 계산의 일치, (4) 원 엔진의 고정 펄스 번호와 비교, (5) 천구 좌표만 바꾸는 경로 및 모의 관측 호출의 일관성이다. 관측 잔차가 작다는 이유로 연속 오차 인증을 통과시키지 않는다.

## 판정

분류: Proven. 추가 천체의 Einstein·Shapiro 지연을 공통 호출 경로에 구현하고, 실제 12,474개 TOA와 128개 모의 관측에서 실행했다. 비균일 자연 스플라인의 조건부 연속 오차식을 도출했으며, 지정된 영 scalar 가지의 정확한 무구동 경계를 분리했다. 이번 단계는 **정리 진전**이다. 새 scalar 검출, EOS에 대한 관측 제약, 완전한 비선형 관측 추론을 확립한 결과는 아니다.

| 항목 | 분류 | 결과 |
|---|---|---|
| 추가 천체 GR Einstein·Shapiro 항 | Proven | 실제·초기화·모의 관측의 공통 경로로 연결, native 대조 검사 통과 |
| 천구 좌표 변경 시 저장 상태 회전 | Proven | 모든 천체의 위치·속도에 적용; 영 이동 복귀 검사 통과 |
| 비균일 스플라인의 값·시간 미분·셀 적분 상계 | Proven | 끝점 오차와 곡률 잔차가 주어졌을 때의 조건부 정리 |
| 곡률과 운동·시선 경계의 연결 | Proven | 기호식 검증; 실제 전 기간 입력 경계는 미확보 |
| 지정된 영 scalar 가지의 무구동 경계 | Proven | 영 초기·입사 자료 및 고전 해의 유일성에 조건부 |
| 전 기간 및 28개 적합 매개변수의 결합 미분 인증 | Conjectural | 미완료 |
| 전체 동반성·힘·광자 전파 matching 및 전역 추론 | Conjectural | 미완료 |

## 1. 실제 호출 경로의 수정과 측정

분류: Proven. 기존 공통 함수는 추가 천체를 받지 않았다. 신규 `Nbody_nogeometric_delays` 멤버가 보존된 전체 상태와 추가 질량을 전달하게 하여 실제 관측, 초기화, 모의 관측의 호출을 함께 수정했다. 클래스 데이터 멤버는 추가하지 않았고 기존 Cython 인터페이스의 바이너리 해시를 확인했다. RA/DEC만 변하는 캐시 경로에서는 기존 `sp/si/so`뿐 아니라 추가 천체를 포함한 원시 상태의 위치·속도도 회전시킨다.

분류: Proven. 각 추가 천체에 대해 다음 항을 합산한다. `d=x_p−x_extra`, `r=|d|`, `n`은 코드의 SSB→PSB 방향이며, 질량은 태양질량 단위다. 내부 시간은 일이다.

```text
U_extra = (G M_sun / c²) M_extra / r
dE_extra/dt = U_extra
S_extra = −2 (T_sun / 86400) M_extra log[(r−n·d)/c]
```

분류: Imported from prior work. 원래 엔진의 Einstein·Shapiro 지연 규약과 1PN 근사 수준을 그대로 확장했다. 상대론적 타이밍 모형의 배경은 [Voisin 등(2020)](https://arxiv.org/pdf/2005.01388)을 참조한다. 이 문헌을 새 네 번째 천체나 scalar 효과의 경험적 증거로 사용하지 않는다.

분류: Proven. `TIMING16_EXTRA_SCALE=0,1,2`는 지연 항에만 적용한 진단 배율이다. 세 실행에서 운동, 매개변수와 펄스 번호는 고정했다. 수차의 스핀 축은 기존의 세 천체 각운동량 정렬 규약을 유지했다. 따라서 추가 GR 지연의 작동 검증과 전체 4체 상대론·scalar 관측 모형의 완성을 구별한다.

| 비교 | 분류 | RMS | 최대 절댓값 |
|---|---|---:|---:|
| 추가 지연 0 − Request 15 저장 잔차 | Proven | 6.7549 ps | 19.7481 ps |
| 추가 지연 1 − 추가 지연 0 | Proven | 1.3273033 ns | 2.6056848 ns |
| 추가 지연 2의 효과 − 추가 지연 1 효과의 두 배 | Proven | 1.8461e−7 ns | 1.4804e−6 ns |

분류: Proven. 이 표는 재적합 전 잔차의 차이다. 기준선의 작은 변경이므로 이전 동결 likelihood·한계값을 새 모형의 결과로 재표기하지 않았다. 배율 검사도 연속 오차 보장이나 매개변수 전체에서의 선형성 증명은 아니다. 원시 근거는 [`live.npz`](../outputs/nbody-readout16/live.npz)와 [`live.json`](../outputs/nbody-readout16/live.json)에 있다.

분류: Proven. native 검사에는 추가 질량 0, 두 추가 천체의 합과 순서 교환, 시선과 위치를 함께 회전한 공변성, 상수 퍼텐셜의 Einstein 적분, 닫힌 형태의 Shapiro 값이 포함된다. 양의 질량 충돌·로그 특이점·음의 질량·잘못된 진단 배율·잘못된 상태 차원을 거부한다. [`native-control.log`](../outputs/nbody-readout16/native-control.log)는 실제 공유 라이브러리에 링크한 검사 결과다.

분류: Proven. 실제 각 실행에는 보간점 458,938개가 있다. Einstein/Shapiro 각각에서 간격을 두고 추출한 65개 입력을 80비트 hexadecimal로 저장하고 정확한 유리수로 읽었다. 좌표·시선 값 주위에 명시적으로 선언한 작은 구간을 놓고 독립적인 바깥 반올림 대수 계산과 native 값을 비교했다. 3개 배율×2개 성분의 검사가 통과했다. 이 구간은 참 궤도를 포괄한다는 인증이 아니며, 전 격자 검사도 아니다. 근거는 [`native-sample-audit.json`](../outputs/nbody-readout16/native-sample-audit.json)에 있다.

분류: Proven. 모의 관측 API의 BAT는 일, 네 지연 성분은 초다. 감사 코드가 지연에 일→초 변환을 중복 적용했던 오류를 수정했다. 원시 실행 파일은 바꾸지 않았다. 수정 후 Einstein 차이의 최대값은 4.4865232 ns, Shapiro는 0.7652230 ps다. 방출시각도 달라지므로 같은 시각의 항을 빼는 대수 검사와 다르다. [`fake-call-audit.json`](../outputs/nbody-readout16/fake-call-audit.json)에 단위와 수정 이력을 남겼다.

분류: Proven. 첫 회전 검사의 사전 절대 문턱 1 cm는 통과하지 못했다. 출력 위치가 약 10¹² m 규모의 binary64이며 측정된 거리 차이는 0.01611328125 m였다. 이를 숨기지 않고, 후속 회귀 검사에는 `128*eps_float64*max(|x_extra|+|x_p|)=0.02865695162 m`를 사용했다. 복귀 후 출력된 추가 상태의 최대 차이는 0이었다. 후속 문턱은 경험적인 수치 회귀 허용치이며 엄밀한 회전 오차 인증이 아니다. 실패 로그 `live.log`와 성공 로그 `live-rerun.log`를 모두 보존했다.

## 2. 연속 보간 오차 정리와 적용 경계

분류: Proven. 각 비균일 셀 `[a,b]`에서 `h=b−a>0`, `g∈C²`라 하자. 저장 계수를 정확한 실수로 읽은 cubic을 `S`라 한다. `e=g−S`의 두 끝점 오차가 `eps` 이하이고 셀 전체에서 `|e''|≤rho`이면,

```text
sup |g−S|  ≤ eps + rho h²/8
sup |g'−S'| ≤ 2 eps/h + rho h/2
|integral_a^b (g−S) dt| ≤ eps h + rho h³/12
```

분류: Proven. `e`에서 끝점 선형 보간을 빼면 영 Dirichlet 경계의 `e''` 문제다. Green kernel의 절댓값 적분은 `x(h−x)/2`, 그 x 미분 kernel의 L1 노름은 `[x²+(h−x)²]/(2h)`이고, 첫 식을 셀 전체에 적분하면 `h³/12`가 된다. 선형 끝점 항을 더하면 위 세 부등식을 얻는다. `e=x(x−1)/2`는 단위 셀에서 값과 적분의 상수를 달성한다. 누적 적분은 통과한 셀별 상계를 더해 제한할 수 있다.

분류: Proven. native 자연 스플라인으로 비균일 격자의 `g(t)=t⁴`를 보간하고 저장한 `y2` 및 값·적분 출력을 검사했다. 셀 전체의 곡률 잔차 구간과 도함수의 극점을 함께 검사했다. 자연 경계에서 일반적으로 성립하지 않는 전역 C4 보간 공식을 사용하지 않았다. 근거는 [`spline-error-audit.json`](../outputs/nbody-readout16/spline-error-audit.json)에 있다.

분류: Proven. 실제 관측식에 필요한 `rho`도 명시할 수 있다. 한 셀에서 `r≥rmin>0`, `z=r−d·n≥zmin>0`, `|d|≤R`, `|d'|≤V`, `|d''|≤A`, `|n|≤N`, `|n'|≤N1`, `|n''|≤N2`라 하자. 모든 시간 미분은 같은 시간 단위를 써야 한다. 그러면 다음 상계가 성립한다.

```text
|U_extra''| ≤ (G M_sun/c²) M_extra (A/rmin² + 2V²/rmin³)
Z1 = V(1+N) + R N1
Z2 = A(1+N) + V²/rmin + 2V N1 + R N2
|S_extra''| ≤ 2(T_sun/86400) M_extra (Z2/zmin + Z1²/zmin²)
```

분류: Proven. `r''=(v²+d·a)/r−(d·v)²/r³`와 `(log z)''=z''/z−(z'/z)²`를 사용한다. `1/r` Hessian의 방사·접선 고유값은 `2/r³,−1/r³,−1/r³`이다. Einstein 운동 항 `|v_p|²/(2c²)`의 이차 미분은 `(a_p²+v_p·j_p)/c²`이다. 여기서 `v_p`는 물리 속도, `a_p,j_p`는 선택한 시간 단위에 대한 그 첫째·둘째 미분이다. 일 단위 미분에 초 단위 가속도를 그대로 넣을 수 없다. `S''`는 인접 `y2` 값의 볼록 보간이므로 `rho≤상계(|g''|)+max(|y2_left|,|y2_right|)`로 연결된다. 기호 검증은 [`curvature-transfer.json`](../outputs/nbody-readout16/curvature-transfer.json)에 기록했다.

분류: Conjectural. 실제 전 기간에서 쓸 만한 위치·가속도·jerk·시선 경계, 끝점 오차, 저장 계수로 cubic을 평가·적분하는 산술 반올림이 아직 모두 연결되지 않았다. 위 정리는 **시간 미분**의 정리이며, 매개변수 변화에 따른 초기 상태·이동 격자·Tempo2·역시간 미분을 자동으로 보장하지 않는다. 따라서 전체 timing D2와 전체 28개 적합 매개변수 미분 인증은 여전히 false다.

## 3. 지정된 EOS 후보의 영 구동 경계

분류: Imported from prior work. Request 13은 `A(phi)=exp(beta phi²/2)`, `beta=−4`, 배경 scalar 0의 비회전 SLy 별에 대해 외부 입사 scalar에 대한 감수율과 외향파 pole을 계산했다. 감쇠시간 약 0.1995 ms는 그 지정 문제의 수치 결과다. 이것은 동반성들이 실제로 만드는 비영 외력을 계산한 결과가 아니다.

분류: Proven. 모든 천체를 비스칼라화 가지로 두고, scalar 값과 법선 시간 미분의 초기자료가 영이며 외부·입사 scalar 자료도 영이라고 가정한다. 고전 초기값 문제가 유일한 영역에서는 `phi=0`인 GR 해가 정확히 유지된다. 이때 `alpha(phi)=d log A/dphi=beta phi`이고 물질원 `−4pi alpha T`는 영이다. `phi=epsilon H`, `T=T0+epsilon deltaT`로 전개하면 선형 원천은 `−4pi beta T0 H`이며 `deltaT`만의 독립 구동 항은 없다. 움직이는 동반성이 계수를 시간에 따라 바꾸어도 영 자료에서 비영 해를 생성하지 않는다. scalar 방정식과 coupling 규약은 [Damour와 Esposito-Farèse(1996)](https://arxiv.org/pdf/gr-qc/9602056)의 식 (1.2b)와 quadratic coupling을 따른다.

분류: Proven. 선형 감수율 관계 `q=chi F`에서 `F=0`이면 `q=0`이다. 영 전하 주위에서 scalar 힘 `q*grad(phi)`는 두 scalar 섭동의 곱이므로 일차에서 영이다. 또한 `(A²)'(0)=0`이다. 따라서 이 정확한 영 가지에서 산란 pole을 구했다는 사실만으로 궤도 시간 규모의 강제 응답이나 비영 `beta` 타이밍 신호를 도출할 수 없다. 기호 검사는 [`zero-scalar-branch.json`](../outputs/nbody-readout16/zero-scalar-branch.json)에 있다.

분류: Proven. 이 주장은 영 해의 안정성을 증명하거나 다른 scalarized 별·쌍성 해를 배제하지 않는다. 무한대 경계가 0이라는 조건만으로 모든 정적 해가 0이라고 주장하지도 않는다. 동반성의 비영 전하, 비영 우주론적 배경, 비영 초기·입사 scalar, 불안정성과 작은 섭동은 별도의 문제다. 선택한 가지의 정확한 영 초기값 해와 그런 가지의 존재·안정성을 혼동하지 않는다.

분류: Conjectural. 실제 비영 관측 후보로 이어가려면 이 경계를 벗어나는 물리 자료를 먼저 정하고, 같은 가정 아래 별·동반성 전하, 되먹임, 상호 힘, 운동과 광자 전파를 함께 도출해야 한다. 임의의 외력을 추가한 뒤 기존 EOS matching이 그 외력까지 입증했다고 설명할 수 없다. 기존 J0337의 조건부 계수 상한을 DEF의 `beta` 제약으로 변환하지 않는다.

## 4. 남은 순서와 재현

분류: Conjectural. 남은 실행 순서는 (1) 영 구동 경계를 벗어나는 특정 물리 가지와 초기·경계 자료를 정하고 비영 동반성 구동을 도출, (2) 전체 초기화와 전 기간 변분 경계를 확보하고 위 셀별 오차 및 역시간 연쇄를 연결, (3) 그 물리 관측식으로 nuisance/noise와 pulse 할당을 포함한 비선형 추론 및 포함률을 검증하는 것이다. 이 세 조건이 성립하기 전에는 전체 물리 논문의 완료 판정을 내리지 않는다.

다음 명령은 WSL Ubuntu-22.04에서 기존 의존성을 사용한다. 원시 실행을 덮어쓰지 않는 검증 명령이다.

```powershell
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && python3 verification/nbody_readout_audit.py check"
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && python3 verification/nbody_readout.py verify"
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && python3 verification/verify_unified_paper.py"
```

분류: Proven. 코드·입력·실행 로그·실패 기록·소스 연결은 신규 [`manifest.json`](../outputs/nbody-readout16/manifest.json)에 결속한다. Request 15가 묶은 문서와 원고 manifest의 기존 바이트는 `request15-notes`에 보존하고, 기존 문서에는 이번 결과를 덧붙인다. Request 14의 과거 문서는 Request 15의 보존본으로 검증한다. 원고 PDF와 소스 ZIP은 Request 12 동결 산출물이며 이번 후속 보고를 포함한 새 제출본이 아니다.
