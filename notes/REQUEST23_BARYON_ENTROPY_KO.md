# Request 23 — 바리온 좌표의 물질·엔트로피 보존 GR 재구성

분류: Counterexample candidate. 압력에 조성을 고정하던 경계를 바꿔, 원래 질량 구역의 조성과 선택한 EOS 엔트로피를 보존하는 구대칭 TOV 해를 계산했다. **원래 내부 물질을 보존하는 해**와 **목표 중력질량에 맞춰 모든 구역 질량을 균일하게 바꾼 별도 모형족**을 구분했다. 후자는 추가 온도·엔트로피 조정 없이 목표 GR 질량과 조건부 광학 컷을 통과했다. 실제 열·핵반응 진화까지 해결한 것은 아니다.

분류: Proven. 물질 좌표 TOV 항등식, 엔트로피 정규화, 질량 배율의 보존 경계를 기호 검산했다. 해석적 균일 밀도 GR 해와 비교하고, 직접 FreeEOS 및 다른 적분법으로 최종 연결 조건을 재검증했다. 이번 결과는 theorem progress와 물리 후보의 조건을 개선한 loophole progress다.

## 1. 보존한 물질과 엔트로피의 정의

분류: Imported from prior work. 출발 자료는 Request22에서 원래 37개 구조 열을 정확히 재현한 MESA 모델 19041의 5,735개 구역이다. 질량은 원래 바리온 질량 요소 `dm`, 조성은 원래 22개 동위원소 중 FreeEOS가 지원하지 않는 불소 3종을 제외하고 재규격화한 값이다. FreeEOS의 중성 원자 에너지 기준, \(C_X=\sum_i X_iW_i/A_i\), \(\rho_B=\rho_{\rm atom}/C_X\)는 Request21 정의를 사용한다.

분류: Counterexample candidate. 각 원래 구역 \(j\)의 바리온 질량 \(\Delta B_j\), 재규격화한 \(X_{ij}\), 기준 엔트로피 \(s_{Bj}\)를 물질에 붙인다. 구역 안에서는 조성과 엔트로피가 상수인 모형을 명시적으로 정의했다. 각 구역의 기준 엔트로피는 원래 중심 \(P_j,T_j,X_j\)에서 FreeEOS를 평가한 뒤 \(s_{Bj}=C_{Xj}s_{{\rm atom},j}\)로 정한다. **원래 MESA 전체 EOS의 절대 엔트로피가 인증됐다는 뜻은 아니다.**

분류: Proven. 원래 내부 물질을 보존하는 모형에서는

\[
B=\sum_j\Delta B_j,\qquad
M_i=\sum_jX_{ij}\Delta B_j,\qquad
S=\sum_js_{Bj}\Delta B_j
\]

가 정의상 보존된다. 원소·동위원소의 수 역시 이 재고와 고정한 질량수로 결정된다. 압력·밀도·온도·반지름은 이 조건 아래에서 바뀔 수 있다. 구역별 상수 템플릿을 연속적 실제 항성의 엄밀한 근사 오차까지 인증한 것은 아니다.

분류: Counterexample candidate. 외곽은 원래 비영압 광학 경계에서 끝내고, Request21의 \(\gamma=5/3\) 수학적 대기를 붙여 영압·진공 경계를 만든다. 내부의 정확한 보존과 대기가 추가하는 극소량의 물질은 따로 보고한다. 실제 대기 수송·스펙트럼 모형은 아니다.

## 2. 물질 좌표의 GR 방정식과 수치 구성

분류: Proven. \(G=c=1\), \(f=1-2m/r>0\), 총 에너지 밀도 \(e\), 바리온 질량 밀도 \(b\)에서

\[
\frac{dr}{dB}=\frac{\sqrt f}{4\pi r^2b},\quad
\frac{dm}{dB}=\frac e b\sqrt f,\quad
\frac{d\ln P}{dB}=
-\frac{(e+P)(m+4\pi r^3P)}{4\pi r^4b\sqrt f\,P}.
\]

분류: Proven. 위 식은 반지름 좌표 TOV를 \(dB/dr=4\pi r^2b/\sqrt f\)로 나눈 결과다. 각 고정 조성 구역에서 \(s_{\rm atom}(P,T,X)=s_{Bj}/C_{Xj}\)를 풀어 온도를 정한다. Newton 역산에 필요한 \(\partial s/\partial\ln T|_{P,X}=c_P\)는 EOS의 엔탈피 미분으로 계산하고 엔트로피의 독립 유한차분과 비교했다.

분류: Counterexample candidate. 중심과 표면에서 각각 질량 좌표로 적분하고, 원래 구역 면 하나에서 반지름·중력질량·압력을 연결했다. 중심에서는 \(\ln(B/B_{\rm tot})\), 표면에서는 \(\ln[(B_{\rm tot}-B)/B_{\rm tot}]\)를 사용한다. 각 끝에서 질량을 따로 누적해 매우 가벼운 표면 구역을 큰 전체 질량의 차로 계산하지 않게 했다. 구역 경계마다 적분을 나누므로 조성·엔트로피의 계단을 건너뛰지 않는다.

분류: Counterexample candidate. 실제 EOS 역산 결과를 구역마다 기준 압력의 \(\Delta\ln P\in[-1,1]\)에서 17점, 33점으로 표화했다. 이 범위 밖으로 나가면 계산을 중단한다. 초기 해와 각 단계 결과는 보존했고, 최종 검사에서는 표를 거치지 않는 직접 EOS 역산을 사용했다.

분류: Imported from prior work. 보존하는 rest-mass 좌표와 변화 가능한 중력질량을 구분하는 실제 GR 항성 진화의 예로 [Althaus et al. (2022)](https://www.aanda.org/articles/aa/pdf/2022/12/aa44604-22.pdf)를 확인했다. 그 논문의 고질량 백색왜성 수치나 진화 결과를 이 저질량 후보에 전용하지 않았다.

## 3. 원래 내부 물질을 보존한 해

분류: Counterexample candidate. 원래 바리온 질량을 고정하고 중심 압력·광학 반지름·중력질량을 연결 조건으로 결정한 결과다. 광도는 기존 모형 값을 유지한 진단 조건이다.

| 항목 | 원래 내부 물질 보존 해 |
|---|---:|
| 내부 바리온 질량 | \(3.927844753722042\times10^{32}\) g |
| \(\ln P_c\), \(P_c\)는 cgs | 48.49569697001607 |
| 광학 반지름 | 68,803,179.234 m |
| 대기 포함 ADM 질량 / 기준 \(GM_\odot\) | 0.197674444822770 |
| 지정 중력질량 목표 대비 | +0.0698906764% |
| 조건부 \(T_{\rm eff}\) | 15,836.266 K |
| \(\log g\), cgs | 5.743647 |
| 대기의 추가 바리온 비율 | 약 \(1.01\times10^{-14}\) |

분류: Counterexample candidate. 따라서 이 보존 해는 목표값 0.197536385307을 수치 허용오차 \(10^{-8}\)로 맞추지는 않는다. 그 허용오차는 **모형 매칭의 수치 조건**이며 실제 관측 질량의 오차막대가 아니다. 관측적으로 배제됐다고 판단하지 않는다.

분류: Counterexample candidate. 원래 조성의 중성 원자 정지질량 합은 약 0.197675936084420 \(GM_\odot\) 단위다. 실제 ADM 질량과의 순 차이는 약 \(-1.49126\times10^{-6}\)이며 내부에너지와 중력 부피 보정이 함께 기여한다. 바리온 질량과 중성 원자 정지질량의 기준 차이를 구분해야 위의 0.0699% 차이를 해석할 수 있다.

## 4. 별도로 등록한 바리온 질량 조정 모형족

분류: Counterexample candidate. 보존 해의 목표 질량 불일치 후 `mass-family-plan.json`을 등록했다. 모든 \(\Delta B_j\)에 같은 양의 배율 \(\lambda\)를 곱하고, 각 구역의 조성과 단위 바리온 질량당 엔트로피는 유지한다. 중심 압력·반지름·\(\lambda\)로 연결 조건을 풀며 중력질량 목표를 부과했다. 온도나 엔트로피를 조정하는 인자, 광학량의 적합 조건은 추가하지 않았다.

분류: Proven. 이 변환은 \(X_i\)와 \(s_B\)를 보존하지만 총 바리온·각 동위원소 재고·총 엔트로피는 모두 \(\lambda\)배가 된다. 원래 물질을 보존한 동일 항성의 변환과 구분해야 한다.

분류: Counterexample candidate. 최종 결과는 다음과 같다.

| 항목 | 별도 바리온 질량 조정 후보 |
|---|---:|
| \(\lambda\) | 0.9993015744453465 |
| 원래 총 재고 대비 변화 | −0.0698425555% |
| 내부 바리온 질량 | \(3.925101446571331\times10^{32}\) g |
| \(\ln P_c\), \(P_c\)는 cgs | 48.49328832066844 |
| 광학 반지름 | 68,819,648.589 m |
| 대기 포함 ADM 질량 / 기준 \(GM_\odot\) | 0.197536385307002 |
| 조건부 \(T_{\rm eff}\) | 15,834.371 K |
| \(\log g\), cgs | 5.743136 |
| 대기의 추가 바리온 비율 | 약 \(1.01\times10^{-14}\) |

분류: Counterexample candidate. 조건부 광학 컷 15,500–16,100 K 및 \(\log g=5.67\)–5.97을 통과한다. 이번 단계에서 광학량을 적합하지는 않았지만, 이미 선택한 물질 템플릿과 기존 광도를 사용했으므로 독립적인 실제 항성 예측은 아니다. Request21의 온도 배율 후보와 Request22의 압력 고정 조성 실패는 그대로 보존한다.

분류: Counterexample candidate. 새 질량 조정 후보에서는 수소도 다른 재고와 동일하게 약 0.06984% 줄어든다. 이전 압력 고정 후보에서 수소만 5.42% 감소하던 비보존적 재배치 대신, 조정한 물질의 양과 조성 분포를 명확히 지정할 수 있게 됐다. 이 차이를 새로운 동적 관측량으로 세지 않는다.

## 5. 독립 검사와 수치 한계

분류: Proven. 82개 비격자 상태에서 직접 엔트로피 역산과 독립 유한차분 \(c_P\) 검사가 통과했다. 유한차분 \(c_P\)의 최대 상대 차이는 \(7.40\times10^{-8}\)이었다. 전체 33점 표의 최대 엔트로피 역산 잔차는 상대 \(2.00\times10^{-12}\) 이내였다.

분류: Proven. 해석적 균일 에너지·바리온 밀도의 Schwarzschild 내부해와 해석적 고유 질량 적분으로 새 좌표의 RHS를 비교했다. 19개 지점의 최대 상대 차이는 \(5.33\times10^{-15}\)였다.

분류: Counterexample candidate. 원래 보존 해에서 RK4 세분화 2→4에 따른 반지름 상대 변화는 약 \(9.44\times10^{-12}\), 중력질량 변화는 \(1.35\times10^{-11}\)이었다. EOS 17→33점에 따른 반지름 변화는 \(1.33\times10^{-8}\), 중력질량 변화는 \(9.89\times10^{-13}\)이었다. 질량 조정 후보의 세분화 2→4에서 바리온 배율 변화는 약 \(1.35\times10^{-11}\)이었다. 구역 수 5,735 자체를 바꾼 실제 항성 격자 수렴 검사는 아니다.

분류: Proven. 최종 매개변수를 고정하고 **직접 FreeEOS 엔트로피 역산 + DOP853**로 양쪽 가지를 다시 적분했다. 아래 잔차는 각각 기준 반지름, 내부 바리온 질량, \(\ln P\)로 정규화한 연결 차이다. 질량을 경계에서 지정한 것만으로 검증 완료라고 하지 않고 내부 연결을 별도로 확인했다.

| 직접 EOS 연결 잔차 | 원래 물질 보존 해 | 질량 조정 후보 |
|---|---:|---:|
| \(\Delta r/R_{\rm ref}\) | \(-1.431\times10^{-10}\) | \(3.020\times10^{-11}\) |
| \(\Delta m/B_{\rm tot}\) | \(-8.853\times10^{-13}\) | \(-1.196\times10^{-12}\) |
| \(\Delta\ln P\) | \(2.155\times10^{-9}\) | \(-8.303\times10^{-9}\) |

분류: Proven. 모두 사전에 정한 연결 허용오차 \(10^{-8}\)를 통과했다. 직접 역산의 엔트로피 잔차 역시 상대 \(2.00\times10^{-12}\) 이내였다. 이 수치 일치는 전체 EOS·연속체 오차의 엄밀한 상계나 항성 안정성 증명이 아니다. 표와 적분 알고리즘은 독립적으로 비교했지만 EOS 자체는 같은 FreeEOS다.

## 6. 해결 범위와 다음 단계

분류: Counterexample candidate. 이번에 해결한 범위는 **지정한 물질·기준 엔트로피를 보존하는 GR 정수압 재구성**과 **그 프로필을 유지한 별도 질량 조정 모형족의 수치 GR 매칭**이다. 후자는 추가 온도 조정 없이 조건부 광학 컷도 통과했다. 원래 총 물질의 보존 해와 목표 중력질량 해를 하나의 결과로 섞지 않는다.

분류: Conjectural. 다음 단계에서는 조성이 변할 때 정지질량·핵반응 Q값·내부에너지·엔트로피의 기준을 중복 없이 연결해야 한다. 이후 새 상태의 반응·중성미자 손실·복사·전도·대류를 재계산하고 GR 열·조성 진화를 풀어야 한다. 고정 조성의 \(Tds\) 관계를 반응하는 조성에 그대로 적용하면 필요한 화학·기준 변화 항을 빠뜨릴 수 있다.

분류: Conjectural. 원래 22개 동위원소 전체, MESA EOS 절대 엔트로피, 회전·형성 경로·열수송·안정성 및 유체·metric·scalar 응답은 아직 연결하지 않았다. 기존 평탄 고정 구각 응답 오차 상계도 이 물질 좌표 GR 모형으로 자동 이전되지 않는다. 완전한 관측 비선형 추론과 전구간 미분 인증은 미완료다. 원고 PDF와 ZIP은 Request12 동결본이다.

## 재현 자료

분류: Proven. 시작 체크포인트는 `729aca0`, 구현은 `verification/baryon_entropy.py`, 결과는 `outputs/baryon-entropy23/`에 저장했다. 과거 문서는 바이트 접두부와 별도 스냅숏으로 보존하고, 원고 revision manifest에 이번 보조 자료만 추가한다. 기존 EOS·Python 실행 환경을 재사용했다.

분류: Proven. 저장된 결과의 검증과 실제 EOS 재검사는 다음 명령으로 실행한다. `verify`와 `recheck_direct`는 동결 결과를 덮어쓰지 않는다. 기호·EOS 대조 결과를 생성하는 명령은 동결 후 별도 Git export에서 실행한다. 실행 시간은 기록값이므로 바이트 재현 대상이 아니다.

```sh
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/baryon_entropy.py verify
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/baryon_entropy.py symbolic
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/baryon_entropy.py uniform_control
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/baryon_entropy.py eos_checks
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 python3 verification/baryon_entropy.py recheck_direct
python3 verification/verify_unified_paper.py
```
