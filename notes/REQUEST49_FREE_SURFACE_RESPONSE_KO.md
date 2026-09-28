# 단계 49 — 자유 표면의 내부 유체·scalar·계량과 동적 외부 연결

분류: Counterexample candidate. 비영 분자 EOS 배경에서 유한 온도 **단열 유체·scalar·계량과 이동 자유 표면·동적 진공 외부를 같은 연립 방정식으로 풀었다.** 원 5,735셀 정보를 재사용하고 정적 및 세 조화 주파수, 두 응답 격자를 계산했다. 변위의 공간 대조는 통과했다. 물질 운동의 outgoing scalar 기여도 직접 분리했지만 전하의 공간 기준은 미달이다. 열·조성 변화와 실제 복사 대기를 포함한 전체 동적 전하는 완료되지 않았다.

## 먼저 정정하는 단계 48의 입력

분류: Counterexample candidate. 저온 기체식을 native 출력과 연결하자 압력·단열 미분은 맞지만 열에너지에 일정한 1,599,868 erg/g 차이가 나타났다. 원 `def_free_surface_response.py`의 build는 사전 에너지 기준에서 실패했다. 원 출력과 코드는 그대로 보존했다.

분류: Proven. `direct_ion_eos.EOS.__call__`의 실제 단위 변환에서 분자 바닥 에너지 출력 `raw[11]`의 중성 원소 질량 기준과 구조 계산의 동위원소 정지 질량 기준은 다르다. 바리온 질량당 바닥 엔탈피는 `h0=CX_rest-CX_EOS*E_bind/c²`다. 단계 48에서 두 CX를 같게 취급한 식은 정정한다. 실제 변환 경로는 `ym=(X/A)@mapping`, `CX_EOS=ym@weights`다.

분류: Counterexample candidate. 이번 조성에서 `CX_rest=1.0070231612220188`, `CX_EOS=1.007022322061226`이다. 올바른 계수를 사용한 다섯 native 저온 상태 127–248 K의 압력·열에너지·Gamma1 최대 상대 차이는 각각 1.665e-15, 1.998e-15, 1.110e-16이다. 이는 유한한 겹침 구간의 일치이며 전체 물리 EOS 인증이 아니다.

분류: Counterexample candidate. 수정된 조건부 자유 표면은 **69,193,584.98351425 m**로, 단계 48의 값보다 0.02901402 m 작다. 수정 구간은 `[69193584.98351128, 69193584.98351721] m`다. **단계 48의 이전 표면 위치·마이크로미터 구간은 철회한다.** 반폭 2.954e-6 m는 동결 입력과 기체 가지 가정 아래 생략 물질원에 대한 상계일 뿐 전체 EOS·배경 오차가 아니다. 상수 바닥 에너지 이동이 소거되는 warm native 엔탈피 차이와 제1법칙 보간 수정은 유지한다.

분류: Conjectural. 127 K 아래에서는 같은 조성의 이상 중성 원자·H2 기체, 같은 302개 분자 준위 및 광자를 명시적으로 사용했다. 마지막 native 상태에서 엔트로피 상수를 맞추고 기체 가지를 진공까지 잇는다. 응축·상전이·외부 복사욕을 제외한다는 가정과 저온 native 가지의 전역 존재 미인증은 유지한다.

## 실제로 연결한 방정식

분류: Imported from prior work. scalar·metric의 방사형 섭동 체계와 진공 영역은 [Mendes–Ortiz 2018](https://arxiv.org/html/1802.07847v2), 정적 DEF·Just 해는 [Damour–Esposito-Farèse 1996](https://arxiv.org/pdf/gr-qc/9602056)을 참고했다. 기존 외부 연산자를 재사용한다. 차가운 barotropic 내부 가정은 가져오지 않고, 아래 유한 온도 단열 내부식을 별도로 정리했다.

분류: Proven. `r=r_phys/R0`, `m=m_geom/R0`, 압력·에너지는 같은 기하 단위, `Phi=phi'`, `b=1-2m/r`, `A⁴=exp(2 beta phi²)`, `alpha=beta phi`, `w=e+p`로 둔다. 주파수는 `omega=Omega R0/c`, 시간 관례는 `exp(-i Omega t)`다. 배경 함수는 다음과 같다.

```text
g = (ln N)' + alpha Phi
  = m/(r²b) + 4 pi r A⁴ p/b + r Phi²/2 + alpha Phi
F = Phi'
  = 4 pi A⁴/b [alpha(e-3p)+r Phi(e-p)] - 2(r-m)Phi/(r²b)
```

분류: Proven. 물질을 따라가는 변위 `xi`, `zeta=xi/r`, `eta=Delta p/p`, `f=Delta phi`, `V=Delta Phi`를 쓰고 `Delta s=Delta X=0`를 가정한다. 각 셀의 native `Gamma1`을 사용하면 `Delta ln rho=eta/Gamma1`, `Delta e=w eta/Gamma1`이다. 엔트로피·조성의 배경 기울기를 따로 차분할 필요 없이 다음 식을 얻는다. 여기서 `Delta g`, `Delta F`는 각각 배경 함수의 여섯 변수 `(r,m,p,e,phi,Phi)`에 대한 전미분이다.

```text
Delta m = r² b Phi f - (4 pi r² A⁴ p + r² b Phi²/2) xi
Delta lambda = (Delta m/r - m xi/r²)/b
xi' = -eta/Gamma1 - 2 zeta - Delta lambda - 3 alpha f
zeta' = (xi'-zeta)/r
eta' = (w/p)[omega² xi/(N²b) + g(2zeta+Delta lambda+3alpha f)
                          - Delta g + g eta] - g eta
f' = V + xi' Phi
V' = Delta F + xi' F - omega²(f-xi Phi)/(N²b)
```

분류: Proven. 중심에서 `eta+3 Gamma1(zeta+alpha f)=0`, `V=0`를 적용한다. 표면에서는 압력 방정식의 발산 계수 앞 대괄호를 0으로 놓는 정칙 조건과 이동 scalar 접합을 함께 사용한다. outgoing 진공 연산자의 `Z=R delta phi'/delta phi`, 전달계수 `h`에 대해 경계식은 다음과 같다.

```text
R V - R² zeta F - Z f + Z R zeta Phi = drive
drive = exp(-i omega R)/(N_surface sqrt(b_surface) h)
```

분류: Proven. 구동 정규화는 단위 regular incident monopole의 Wronskian 조건이며 영주파수는 그 정칙 극한이다. 물질 자유 표면과 외부 outgoing 파동은 같은 선형 해에서 결정된다. 핵 반응·열 섭동을 생략한 기계적 단열 블록이라는 범위는 식 자체의 가정이다.

## 결과와 실패를 구분한 판정

분류: Counterexample candidate. 원 내부 5,735셀에 native 외층과 명시적 저온 기체 구간을 연결했다. 배경 절점은 5,864개이며 응답 격자는 5,863/11,726셀이다. 기준 쌍성의 기존 주기 15,665.751496628753초를 유지했다. 표면 반경이 바뀌었다고 주기를 다시 정의하지 않았다. 각 주파수의 진폭은 단위 입사파이며 실제 쌍성의 개별 Fourier 계수는 아직 적용하지 않았다.

분류: Counterexample candidate. 아래는 해당 동결 배경 위의 계산값이다. 변위 차이는 전 격자의 최대 변위 차이를 미세 격자 최대 변위로 나눈 값이다. 전하 계수는 outgoing scalar monopole 차이를 **고정한 배경 ADM 질량**으로 나눈 값이다. 전체 동적 질량 정규화나 관측 전하로 해석하지 않는다.

| 주파수 / 기본 주파수 | 표면 xi/R, 단위 구동의 실수부 | 변위 격자 상대 차이 | 물질 운동 scalar 계수의 실수부 | 해당 계수의 격자 차이 |
|---:|---:|---:|---:|---:|
| 0 | -0.020578093 | 9.42845e-5 | -6.25084594e-9 | 2.03932% |
| 1 | -0.020587672 | 9.42842e-5 | -6.25096160e-9 | 2.03975% |
| 2 | -0.020616442 | 9.42834e-5 | -6.25130867e-9 | 2.04104% |
| 3 | -0.020664495 | 9.42821e-5 | -6.25188743e-9 | 2.04320% |

분류: Counterexample candidate. 여덟 전체 선형 해의 잔차는 최대 7.35e-16이며 변위 공간 기준 2%를 통과했다. 두 해의 O(1) scalar 값을 그대로 빼면 차이가 약 2.6e-14에 불과해 binary64 뺄셈의 손실이 나타났다. 최초 전하 비교는 네 주파수 중 세 개가 2% 기준에 미달했다. 그 실패를 보존했다.

분류: Proven. 이산 유체 고정 비교는 전체 행렬 `M`에서 운동량 행만 `xi=0`으로 바꾼 `M_c`로 정의한다. 같은 우변 `b`에 대해 `M_c y_c=b`라면 `M Delta y=(M_c-M)y_c`이다. 우변은 교체한 운동량 행에만 존재한다. 공통 scalar·바리온 행의 잔차를 수치적으로 빼지 않으므로 O(1) 해의 뺄셈이 필요 없다. 비교를 유지하는 외력이 필요하며 별도의 자유 진화 항성이라는 뜻은 아니다.

분류: Counterexample candidate. 이 차이 방정식으로 물질 운동의 outgoing 기여를 직접 얻었다. 같은 격자에서 확장 정밀도 잔차 보정을 두 번 적용하자 유체 고정·차이 해의 잔차는 최대 3.24e-16, 기존 전체 해의 재구성 상대 차이는 최대 2.26e-13이었다. 전하 계수 변화는 보정 전후 상대 약 2e-9 이하였지만 **공간 차이 2.039–2.043%는 남았다. 사전 2% 기준은 미달**이다. 추가 격자나 주파수를 실행하지 않았다. 정밀도 손실의 해결과 공간 정확도 수락을 구분한다.

분류: Counterexample candidate. 선언된 최대 scalar excursion을 공통 척도로 곱하면 `delta alpha/phi0` 크기는 약 1.178e-15, 정적 값을 뺀 차이는 첫 세 주파수에서 1.111e-19, 2.348e-19, 3.814e-19다. 이 값은 실제 조화 진폭을 적용한 결과도, nuisance 제거 후 잔여도, 인증된 상계도 아니다. 작은 수치만으로 관측 no-go를 선언하지 않는다.

## 열 정상성의 판정과 다음 병목

분류: Proven. 정지 물질·열유속 0을 가진 정적 배경이 양의 전도율의 수송 법칙까지 정상적으로 만족하려면 Jordan 온도에 대해 `d ln(T_J A N)/dr=0`이어야 한다. 초기 열유속을 0으로 지정하는 것은 이 정상 조건의 증명이 아니다.

분류: Counterexample candidate. 저장된 비영 배경에서 `ln(T_J A N)` 범위는 7.48641839, 최대/최소 비는 1,783.6523이다. 해당 무열유속 정상 조건은 실패한다. 온도 차이 자체로 냉각 시간이나 궤도 동안의 변화량을 정할 수는 없다. 이번 주파수 해는 단열 기계적 기준 해이며 약 4.35시간 궤도 전체의 정상 응답으로 승격할 수 없다.

분류: Conjectural. 다음 우선 작업은 같은 배경·자유 표면 위에 열원·열유속·조성의 비정상 변화를 연결하는 것이다. 기존 수송·반응 계수와 저장 전하 미분을 재사용해 필요한 대표 시간 구간과 관측 오차 예산을 먼저 정한다. 배경 변화의 영향이 목표보다 작다고 제한할 수 없으면 단열 주파수 해를 최종 결과로 사용하지 않고 비정상 보존 접선 진화로 넘어간다. 물리적 복사 경계, 전하 공간 오차, 실제 조화 구동과 동일 재고 정적 비교·미분 nuisance·결합 오차는 모두 남아 있다.

## 실행·재현·보존

분류: Counterexample candidate. 새 native EOS 호출과 새 시간 적분 단계는 각각 0이다. 저장된 전 셀 native 출력·분자 준위·배경·lapse·외층을 재사용했다. 전체 응답 여덟 sparse solve의 기록된 합은 약 3.70초, 직접 차이 및 두 번 잔차 보정 묶음은 각각 4.57초와 4.65초다. 이는 내부 기록 시간으로 WSL 시작·import·문서 작업을 포함하지 않는다. 단일 BLAS 스레드 CPU를 사용했고 GPU·장기 궤도 적분을 실행하지 않았다. 최대 각 명령 예산은 build/pilot 30초, 원 묶음 90초, 차이 묶음 60초였다.

분류: Imported from prior work. 실행 전 체크포인트는 `b499ce84`, 직전 단계는 `0156fc60`이다. `outputs/direct-eos-gr33/def-free-surface-response`의 원 실패와 `normalized`의 수정 배경·응답, `direct-contrast` 및 `direct-contrast-refined`를 모두 보존한다. 각 plan에 실행 소스와 입력 SHA가 있다. 같은 경로의 `manifest.json` 및 `gr-free-surface-response-milestone-manifest.json`이 새 근거를 묶는다. 기존 단계 48 packet의 17개 SHA를 문서 추가 전에 대조했다.

분류: Imported from prior work. 재현 환경은 WSL Ubuntu-22.04, `OPENBLAS_NUM_THREADS=1`, `PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification`이다. 순서는 `def_free_surface_response_normalized.py build`, `correct_surface`, `pilot`, `run`, 이어서 `def_free_surface_charge.py`, `def_free_surface_contrast.py`, `def_free_surface_contrast_refined.py`다. 산출물 덮어쓰기를 차단하므로 재실행은 별도 깨끗한 작업 사본에서 수행한다. 원 실패 build는 정정 전 에너지 기준 미달을 재현하는 별도 근거다.

분류: Counterexample candidate. 이번 분류는 **loophole progress**다. 자유 표면의 단열 내부·외부 연결은 실제 계산으로 진전됐고, 전하 읽기의 소거 오차와 표면 입력 오류를 고쳤다. 전하 공간 기준 미달과 전체 열·복사·관측 미완료를 유지한다.
