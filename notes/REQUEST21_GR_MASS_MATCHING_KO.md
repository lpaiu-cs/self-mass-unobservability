# Request 21 — 실제 열 EOS를 이용한 GR 질량 매칭

분류: Counterexample candidate. **명시한 열·조성 조건의 구대칭 GR 모형에서는 질량 매칭을 달성했다.** 원래 온도 분포를 고정한 모형은 광학 조건을 실패했다. 이 실패와 별도로 등록한 온도 배율 조정에서는 질량·조건부 유효온도를 함께 맞추는 후보를 얻었다. 원래 MESA 진화 모형 전체를 GR로 완성했다는 결론은 아니다.

분류: Proven. 이번 진전은 내부에너지 영점과 정지질량의 일관된 연결, TOV 구조 재계산, 고유 바리온 질량 적분, 영압 진공 경계 및 독립 적분 검증이다. 단순히 기존 질량에 U/c²와 Newtonian 결합에너지를 더한 계산이 아니다.

## 1. 사용한 EOS와 질량 정의

분류: Imported from prior work. [FreeEOS 3.0.0 공식 배포](https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/)의 EOS1 옵션 `(3, 1, -2)`을 사용했다. 이 버전은 20개 원소의 295개 이온화 단계를 다루는 열 EOS를 제공한다. 내려받은 소스 원본은 수정하지 않았으며 C ABI 연결만 추가했다. 배포 압축 파일과 실제 사용 소스의 바이트 일치를 확인했다. 배포 서명 검증은 수행하지 않았다.

분류: Proven. FreeEOS의 실제 구현은 비수소 원소에 중성 원자의 전자 바닥상태를, 수소에는 H₂ 바닥상태를 내부에너지 영점으로 사용한다. 다음 원본 소스에 정의와 계산이 함께 있다.

- [FreeEOS API와 조성 정의](../outputs/gr-mass21/sources/src/mod_free_eos.f90)
- [내부에너지 계산과 H₂ 영점 이동](../outputs/gr-mass21/sources/src/free_eos_detailed.f90): 원본 3487–3492행 및 3693행 이하
- [이온화 에너지 합산](../outputs/gr-mass21/sources/src/ionize.f90)
- [물리 상수](../outputs/gr-mass21/sources/src/mod_free_eos_constants.f90), [H₂ 해리 에너지](../outputs/gr-mass21/sources/src/mod_ionization_data.f90)
- [연결 코드](../verification/gr_eos_bridge.f90), [실행·검증 코드](../verification/gr_mass.py)

분류: Proven. 중성 원자 질량을 정지질량으로 사용하기 위해 연결 코드는 원본의 수소 영점 이동을 정확히 뺀다. 원본 상수로

\[
u_{\rm atom}=u_{\rm FreeEOS}
-\frac12(c_2\mathcal R)D_{H_2}\epsilon_H
\]

를 계산한다. 여기서 `c2*cr = h c N_A`, `h2diss = 36118.3 cm^{-1}`이며, \(\epsilon_H\)는 FreeEOS가 정의한 단위 질량당 수소 원자 수 계수다. 중성 원자의 질량에 포함된 전자 정지질량을 다시 더하지 않는다.

분류: Proven. 바리온 분율 \(X_i\), 질량수 \(A_i\), 배포 원자량 \(W_i\)에 대해

\[
C_X=\sum_i X_i\frac{W_i}{A_i},\qquad
\epsilon_Z^{\rm EOS}=\frac{1}{C_X}\sum_{i\in Z}\frac{X_i}{A_i},\qquad
\rho_B=\frac{\rho_{\rm atom}}{C_X}
\]

로 입력 조성과 바리온 밀도를 연결했다. FreeEOS에 전달하고 돌려받는 밀도는 중성 원자 질량으로 정의한 밀도다. 따라서 총 에너지 밀도는

\[
\mathcal E=\rho_{\rm atom}(c^2+u_{\rm atom})
\]

이며, 이 밀도에 \(C_X\)를 한 번 더 곱하면 이중 보정이 된다. 원자량과 질량 단위는 결속한 수치 EOS의 규약을 따른다. 이 규약의 수치 일치를 계측 상수의 불확실성 인증으로 해석하지 않는다.

분류: Proven. \(u\to u+\delta\), 정지 에너지 계수 \(c^2\to c^2-\delta\)의 동시 변화는 총 에너지 밀도를 보존한다. 반대로 정지질량 규약을 고정한 채 내부에너지 상수만 바꾸면 다른 GR 모형이 된다. 이를 기호 계산으로 확인했다.

## 2. 실제로 계산한 항성 모형의 범위

분류: Counterexample candidate. 시작점은 Request 20의 공간 격자 개선 후보 `mesh`, 모델 19041, 5735개 구역이다. 원래 모형에서 \(T(P)\)와 조성의 압력 의존성을 가져왔다. 온도는 log P에서 PCHIP, 조성은 정규화와 원자 수의 선형 관계를 보존하는 선형 보간을 사용했다. 중심의 마지막 표본을 넘어서는 압력에서는 중심 온도와 조성을 고정했다. 매 압력에서 실제 EOS를 다시 평가해 밀도와 내부에너지를 얻었다. 원래 밀도·에너지 곡선을 그대로 GR EOS로 간주하지 않았다.

분류: Counterexample candidate. 이는 **온도와 조성의 압력 경로를 지정한 평형 구성**이다. 질량 좌표에 붙은 엔트로피·조성을 보존하는 GR 진화 계산은 아니다. 원래 수소 외피 질량, 생성 이력, 열수송, 광도 방정식, 회전 평형을 함께 고정하거나 재계산하지 않았다. 광학 진단에는 원래 모형의 광도를 유지했다.

분류: Proven. FreeEOS의 지원 원소에는 불소가 없다. 이번 계산에서는 f17·f18·f19를 제외하고 나머지 바리온 분율을 다시 정규화했다. 제외 분율은 원래 입력 프로필에서 최대 \(1.8143\times10^{-7}\), 원래 구역 질량으로 가중한 평균 \(2.6945\times10^{-8}\)이다. 이 값은 새 GR 별 전체의 물리 오차 상계가 아니다. 별도의 불소 대체 민감도 시험은 실행하지 않았다.

분류: Counterexample candidate. 따라서 조정 후보는 원래 22개 동위원소 모형과 동일한 별이라는 인증을 갖지 않는다. FreeEOS로 교체한 미시물리, 지정한 압력 경로, 극미량 조성 생략 및 외곽 대기가 모두 명시된 가정이다. EOS의 동위원소별 분배함수·열역학 근사까지 정밀 핵물리 인증으로 승격하지 않는다.

## 3. 완전한 구대칭 TOV와 영압 외곽

분류: Proven. 계산에서는 \(G=c=1\)의 기하 단위를 사용한다. \(e\)는 총 에너지 밀도, \(b\)는 바리온 질량 밀도의 기하 단위다.

\[
\frac{dm}{dr}=4\pi r^2e,\qquad
\frac{dP}{dr}=-\frac{(e+P)(m+4\pi r^3P)}{r^2(1-2m/r)},
\]
\[
\frac{dM_B}{dr}=\frac{4\pi r^2b}{\sqrt{1-2m/r}},\qquad
\frac{d\nu}{dr}=\frac{m+4\pi r^3P}{r^2(1-2m/r)}.
\]

분류: Proven. 주 계산은 \(\ln P\)를 독립변수로 바꿔 풀고, 별도의 검증은 반지름을 독립변수로 직접 풀었다. 정지질량, 내부에너지, 압력의 관성·능동 중력 항, 고유 부피를 포함한다. 질량은 영압 경계에서 Schwarzschild 외부 해의 \(m(R)\)로 정의했다.

분류: Counterexample candidate. 원래 첫 구역의 압력이 0이 아니므로 그 지점을 정확한 진공 경계로 취급하지 않았다. 이 지점을 광학 반지름으로 정하고, 아래의 별도 등록한 \(\gamma=5/3\) 외곽을 이어 붙였다.

\[
P=P_s z^{5/2},\quad \rho=\rho_s z^{3/2},\quad
u=u_0+\frac{3P}{2\rho},\quad
u_0=u_s-\frac{3P_s}{2\rho_s},\quad 1\ge z\ge0.
\]

분류: Proven. 이 외곽은 연결점의 밀도·에너지를 연속으로 맞추고 \(du=P\,d\rho/\rho^2\)를 만족하며, \(z=0\)에서 압력·밀도가 0이 된다. 제1법칙과 압력 좌표 변환은 [기호 검증](../outputs/gr-mass21/symbolic-audit.json)을 통과했다.

분류: Counterexample candidate. 최종 조정 후보에서 이 수학적 외곽의 질량 비율은 약 \(1.1\times10^{-14}\), 두께는 약 124.3 km다. 이는 광학 대기 수송·분광 계산을 대신하지 않는다. 광학 반지름과 진공 반지름을 구분해 보존했다. 표면 중력 진단에서 이 외곽의 미소 질량은 최종 총질량에 포함했다.

## 4. 수치 결과와 유지한 실패 판정

분류: Counterexample candidate. 아래 질량은 \(GM/(GM_\odot)\) 단위이며 `GM_sun=1.3271244e20 m³/s²`를 사용한다. 반지름 단위는 앞선 MESA 자료와 동일한 `R_sun=695980000 m`다. 유효온도는 유지한 모형 광도와 재계산한 광학 반지름에서 얻은 조건부 값이다.

| 모형 | 중력원 질량 | 고유 바리온 질량 | 광학 반지름/R_sun | 조건부 Teff/K | logg/cgs | 광학 컷 |
|---|---:|---:|---:|---:|---:|---|
| 동일 EOS Newtonian, 온도 고정 | 0.197536385286 | 0.197398808525 | 0.0903804223 | 16562.332 | 5.82121648 | 실패 |
| GR, 온도 고정 | 0.197536385365 | 0.197400320165 | 0.0903451153 | 16565.568 | 5.82155788 | 실패 |
| GR, 별도 온도 배율 조정 | 0.197536385339 | 0.197399166639 | 0.0993123229 | 15799.999991 | 5.73936077 | 통과 |

분류: Counterexample candidate. Newtonian 행의 중력원 질량은 정지질량 밀도의 적분이고 GR 행은 총 에너지의 TOV 질량이다. 각 이론의 같은 수치 목표에 맞춘 결과다. 모든 행의 질량 잔차가 작다는 사실을 이들 질량 정의가 같다는 주장으로 사용하지 않는다.

분류: Counterexample candidate. 최종 조정 후보의 온도 배율은 **1.01996801252538**, 중심 압력은 \(1.15785222688\times10^{21}\) dyne/cm²다. 중심 압력과 온도 배율 두 개로 질량과 조건부 Teff를 맞췄다. 원래 고정 온도 시험의 실패 후 별도 계획을 등록했으며 그 실패 판정을 보존한다. **맞추는 데 사용한 Teff·반지름은 독립 예측이나 관측 증거가 아니다.**

분류: Counterexample candidate. 같은 EOS·고정 온도 경로의 Newtonian–GR 비교에서는 반지름 차이가 −24.573 km, 조건부 Teff 차이가 +3.236 K다. 원래 MESA 진화 후보와의 훨씬 큰 차이에는 EOS 교체와 구성 조건 변화가 포함된다. 전체 차이를 GR 효과로 돌릴 수 없다.

분류: Proven. 초기 직접 호출·낮은 정밀도 시험은 더 엄격한 재적분에서 질량 허용오차를 통과하지 못했다. 그 결과는 `gr-match.json`, `newtonian-match.json`에 그대로 남겼다. 이후 EOS 표본을 고정하고 정밀도를 높인 결과는 `*-refined-*`, 추가 조정 결과는 `calibrated-*`로 분리했다. 실패한 시험을 성공한 수치로 덮어쓰지 않았다.

## 5. 검증 범위

분류: Proven. [EOS 검증](../outputs/gr-mass21/eos-checks.json)은 지정한 31개 상태에서 압력→밀도→압력 역변환 최대 상대차 \(5.56\times10^{-15}\), 밀도 상대차 \(1.78\times10^{-15}\), Maxwell 관계 상대차 \(2.04\times10^{-13}\)를 확인했다. 희박한 중성 He의 내부에너지는 광자 항을 포함한 독립 이상기체 계산과 해당 부동소수점 정밀도에서 일치했다. 이들은 표본 검증이다.

분류: Proven. 균일 에너지 밀도 별의 해석적 Schwarzschild 내부 해를 별도의 양성 대조로 사용했다. 압력 곡선의 정규화된 최대 오차는 \(1.28\times10^{-15}\)였다. 주 TOV 적분과 반지름 좌표 적분의 차이도 별도로 계산했다.

분류: Proven. EOS 표본은 원래 압력 구역을 1·2·4분할해 5835·11670·23340개로 늘렸다. 온도를 고정한 GR 후보의 반지름은 약 62878.209, 62878.365, 62878.393 km였다. 표본마다 중심 압력을 다시 맞췄으므로 각각의 작은 질량 잔차를 표본 오차 자체로 해석하지 않는다.

| 검증량 | 온도 고정 GR | 별도 조정 GR |
|---|---:|---:|
| 최종 주 계산 질량 상대 잔차 | 2.92×10⁻¹⁰ | 1.64×10⁻¹⁰ |
| 실제 EOS 직접 호출·반지름 적분의 질량 상대 잔차 | 1.60×10⁻⁹ | 2.39×10⁻⁹ |
| 독립 감사에서 직접 적분과 압력 적분의 반지름 상대차 | 7.00×10⁻⁸ | 1.02×10⁻⁷ |
| 67개 확인점의 최대 EOS 표–직접 밀도 상대차 | 8.36×10⁻⁵ | 8.35×10⁻⁵ |

분류: Proven. 위 표는 [고정 온도 독립 감사](../outputs/gr-mass21/primary-audit.json)와 [조정 후보 독립 감사](../outputs/gr-mass21/calibrated-audit.json)의 수치 결과다. 두 직접 EOS 질량 잔차는 등록한 \(10^{-8}\) 기준을 통과했다. 이 허용오차는 **지정한 수치 모형의 목표값 일치 기준**이며 관측 질량의 불확실성이 아니다.

분류: Conjectural. 표본 차이와 적분기 tolerance는 엄밀한 전구간 오차 상계가 아니다. 특히 국소 EOS 표–직접 밀도 차이가 남으므로 EOS 자체의 전구간 정밀도, 물리적 정확도, 진화 경로의 오차를 인증하지 않는다. 이번 조정 반지름의 1e-8 수치 목표도 독립 직접 적분의 전구간 반지름 인증으로 승격하지 않는다.

## 6. 완료 경계와 다음 연구

분류: Proven. 선택한 EOS에 대해서는 절대 에너지 기준을 결속하고 고유 부피 및 구조 반작용을 포함해 GR 질량을 계산하는 경로가 구현됐다. 앞선 “에너지 진단만 존재한다”는 경계는 이 지정 모형에서는 넘어섰다. 결과 유형은 에너지·TOV 정의에 관한 정리 진전과 조건부 GR 평형 후보에 관한 허점 연구 진전이다.

분류: Conjectural. 다음 경계는 조정한 후보를 자체 열수송·진화 해로 연결하는 일이다. 원래 22개 동위원소 전체, 회전, 광도와 대기, 형성 이력 및 실제 궤도 조건을 일관되게 맞추는 검증이 남았다. 정적 평형의 동역학적 안정성도 인증하지 않았다.

분류: Conjectural. 유체·metric·scalar 결합 동역학, 실제 자유낙하 관측 전달함수, 전구간 미분 오차, 완전한 비선형 관측 추론은 미완료다. 기존 영 배경·평탄·고정 구각 scalar 상계를 이 GR 유체 별의 오차 보장으로 옮기지 않는다. 이번 정적 질량 매칭 자체는 새 동적 관측량의 증명이 아니다.

분류: Proven. Request 20의 고정 나이 창 실패 및 별도 시점 후보 회복, 이전 raw verdict와 원고 동결 경계는 보존했다. 다섯 동적 연구 문서와 `paper/revision-manifest.json`에는 추가 근거만 연결했다. 원고 PDF·ZIP은 Request 12 동결본이며 이번 결과를 반영한 최종 투고본은 아니다.

## 7. 재현과 보존

실행 스크립트: [gr_mass.py](../verification/gr_mass.py). 결과 전체와 등록 계획·실패 기록: [gr-mass21](../outputs/gr-mass21). 런타임은 `/home/lpaiu/work/gr-mass21`에 분리했다. 기존 MESA 실행 환경은 변경하지 않았다.

소스 빌드에는 이미 설치된 gfortran·CMake·BLAS/LAPACK을 사용했다. 배포판 CMake에서 `BUILD_TEST=OFF`는 문서 유틸리티의 누락된 타깃 오류를 냈으므로 `BUILD_TEST=ON`으로 구성하고 `free_eos` 라이브러리 타깃만 빌드했다. 실패 로그도 보존했다. 새 호스트에서는 보존한 압축 파일을 풀고 같은 구성으로 라이브러리를 만든 뒤 `build` 단계로 연결 코드를 컴파일해야 한다. 런타임과 소스 SHA는 provenance에 결속했다.

현재 호스트의 핵심 재검증 명령:

```powershell
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/gr_mass.py verify'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/gr_mass.py symbolic'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/gr_mass.py eos_checks'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/gr_mass.py audit calibrated'
rtk proxy wsl -d Ubuntu-22.04 -- bash -lc 'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps python3 verification/verify_unified_paper.py'
```

계산을 처음부터 반복할 때는 별도 복제에서 `prepare → build → probe → atmosphere_plan → match → refinement_plan → refined → calibration_plan → calibrate → polish → audit` 순서로 실행한다. 등록 계획은 최초 생성본을 유지한다. `maintain`·`seal`은 문서·manifest를 갱신하는 제작 단계이므로 이미 동결한 산출물을 검증할 때 다시 실행하지 않는다.
