# Request 24 — 조성이 변할 때의 정지질량·핵반응·내부에너지 연결

분류: Proven. **에너지 기준의 연결과 화학적 조성 항을 검증했다.** Request23의 `e_B=C_X c²+u_B`를 유지하면서 MESA 표준 핵반응 Q값을 쓸 수 있는 내부에너지 기준 변환을 구성했다. 기존 GR 해의 질량, 조성, 엔트로피와 과거 실패 판정은 바꾸지 않았다. 이번 성과는 에너지 보존식에 관한 theorem progress이다.

분류: Counterexample candidate. 실제 FreeEOS를 사용한 지정 CNO 반응량 시험에서 유한 에너지 역산과 독립 미분 적분을 대조했다. 기준 보정 누락은 약 **6.21335 ppm**, 화학적 조성 항 누락은 약 **0.481113%**의 열수지 오차를 만들었다. 이를 제거한 계산은 아래 수치 검증을 통과했다. 이 시험은 실제 반응률을 적분한 항성 진화가 아니다.

분류: Conjectural. 실제 반응률·중성미자·열수송을 새 GR 상태에서 계산하는 단계, 완전한 열·조성 진화, 안정성·유체/scalar 응답과 관측 추론은 남아 있다. 이 보고서는 연구 전체의 완료나 투고 준비 완료 판정이 아니다. 원고 PDF/ZIP은 Request12 동결본이다.

## 1. 배포된 질량표와 Q값을 직접 대조한 결과

분류: Imported from prior work. MESA r7624는 `X_i`를 바리온 분율로 정의하고, `Y_i=X_i/A_i`, `C_X=ΣW_iY_i`를 사용한다. 배포 `isotopes.data`에는 원자량 `W_i`와 질량초과가 별도 열에 있다. 표준 Q값은 원자량을 직접 빼서 만드는 대신 `binding_energy`, `del_Mp`, `del_Mn`으로 계산한다. 특히 A≤1에서 결합에너지를 0으로 지정하므로 H1에는 표의 `7.288970947 MeV` 대신 `del_Mp=7.288969 MeV`가 쓰인다. 차이는 1.947 eV이다. [배포 소스](../outputs/reactive-energy24/sources/chem/private/chem_isos_io.f90), [Q 계산](../outputs/reactive-energy24/sources/chem/public/chem_lib.f90).

분류: Proven. 원래 `basic.net → add_hot_cno → add_cno_extras` 및 Ca40 추가를 해석해 **22종·72항목**을 재구성했고 저장된 원래 네트워크 진단과 대조했다. `reactions.list`에 없는 11개 자동 등록 반응은 배포 Reaclib의 실제 입출력 기록으로 확인했다. Q는 Reaclib 피팅 헤더의 Q값을 임의로 가져오지 않고 `set_reaction_info → get_Qtotal` 경로로 재구성했다.

분류: Proven. 72항목 중 **64개의 닫힌 반응식**에서 바리온 수 보존과 에너지 기준 변환을 Decimal로 검산했다. 개별 반응의 원자량 기반 Q와 표준 Q의 차이는 최대 **151.07056 eV**이다. 8개 보조율은 닫힌 독립 반응식이 아니므로 이 인증에 포함하지 않는다: `rbe7ec_li7_aux`, `rbe7pg_b8_aux`, `rc12ap_aux`, `rn14pg_aux`, `rna23pa_aux`, `rna23pg_aux`, `rne20ap_aux`, `ro16gp_aux`. [전체 대조 기록](../outputs/reactive-energy24/reaction-audit.json).

분류: Proven. CNO 여섯 반응을 더하면 중간 C/N/O는 소거되고 `4 H1 → He4`가 된다. 같은 소스 기준으로 얻은 값은 다음과 같다.

| 양 | CNO 순환 1회당 값 |
|---|---:|
| 표준 총 Q | 26.730960209 MeV |
| Request23 원자량 기준 총 Q | 26.730804698711122 MeV |
| 지정 평균 반응 중성미자 에너지 합 | 1.7024 MeV |
| 표준 열 침적 에너지 | 25.028560209 MeV |
| 원자량 기준 열 침적 에너지 | 25.02840469871112 MeV |
| 원자량 기준 − 표준 총 Q | −155.510288879 eV |

분류: Imported from prior work. 원자 질량초과 기반 Q는 전자 정지질량 장부를 포함한다. 소스는 양전자 소멸 에너지를 다시 더하지 않도록 설명하며, 네트워크는 `Q−Q_neutrino`를 `eps_nuc`에 더한다. 따라서 `eps_nuc`에서 반응 중성미자를 다시 빼면 중복이다. 열적 중성미자는 별도 손실이다. [실제 누적 경로](../outputs/reactive-energy24/sources/net/private/net_derivs_support.f90).

분류: Conjectural. 이 대조는 표준 Q와 지정 평균 중성미자 값에 대한 것이다. 실제 `net_eval`은 weak-rate 계산에서 Q와 중성미자 Q를 교체할 수 있다. 이번에 그 상태 의존 계산이나 반응 유량을 실행하지 않았으므로 64개 표준 반응의 일치를 전체 런타임 네트워크의 에너지 인증으로 확대하지 않는다. 반올림된 원자량의 물리적 정확도도 인증하지 않는다.

## 2. 기존 GR 질량을 유지하는 기준 변환

분류: Proven. Request23의 고정된 정의를 W 기준이라 쓰면

\[
 \rho_B=\rho_{\rm atom}/C_X,\qquad
 u_W=C_Xu_{\rm atom},\qquad e_B=C_Xc^2+u_W.
\]

분류: Proven. `K=MeV→erg × N_A`를 배포 상수 그대로 사용하고, 표준 Q를 생성하는 질량초과를 `Δ_i`라 쓰자. 바리온 분율 합이 1인 물질에서

\[
 g(X)=\sum_iY_i[(W_i-A_i)c^2-K\Delta_i],\quad
 u_Q=u_W+g(X),\quad e_{0,Q}=c^2+K\sum_iY_i\Delta_i
\]

이면 **`e_B=e_{0,Q}+u_Q`**이다. 총 에너지는 같고 `u`와 정지에너지 사이의 구분만 바뀐다. 새로운 TOV 질량 적합은 필요하지 않다. 64개 표준 반응의 Q 차이는 이 `g(X)`의 변화와 대조했다. 수소 분율이 바뀔 때 H2 바닥상태 보정도 바뀌므로 기존 FreeEOS bridge를 매 상태에서 그대로 호출했다.

분류: Proven. 단위 바리온 질량의 국소 보존식은, 중성미자를 `q_ν`, 그 외 순 유입을 `q_ext`라 할 때

\[
 \dot u_W+P\dot v_B=-c^2\dot C_X-q_\nu+q_{\rm ext},\qquad
 \dot e_B+P\dot v_B=-q_\nu+q_{\rm ext}.
\]

분류: Proven. Q 기준으로 바꾸면 가열률도 `q_Q=q_W+ġ`로 함께 바뀐다. 기존 EOS의 `u_W`를 유지하면서 표준 Q 가열만 직접 넣는 식은 `ġ`를 빠뜨린다. 총 에너지 식에 정지질량 감소를 이미 넣고 핵 가열을 다시 넣는 식은 중복 계산이다. 반응의 실제 Q가 선택한 질량 퍼텐셜과 다르면 그 차이를 별도로 추적해야 하며, 아무 차이나 하나의 기준 이동으로 없앨 수 있는 것은 아니다.

## 3. 조성이 변하면 엔트로피 항만으로는 부족하다

분류: Proven. 핵 정지에너지를 뺀 열역학적 화학 퍼텐셜을 `μ_i^th`라 쓰면, 단위 바리온 질량의 제1법칙은

\[
 du_B+Pdv_B=Tds_B+\sum_i\mu_i^{\rm th}dY_i.
\]

분류: Proven. 고정 조성에서는 마지막 항이 0이지만 반응 중에는 일반적으로 그렇지 않다. 따라서 Request22의 준정적 GR 광도 식을 조성이 변하는 물질에 적용할 때 다음 항을 명시해야 한다.

\[
 \frac{dL_\infty}{dB}=e^{2\nu}\left[
 q_{\rm nuc}-q_{\rm thermal\,\nu}
 -T\frac{ds_B}{d\tau}-\sum_i\mu_i^{\rm th}\frac{dY_i}{d\tau}
 \right].
\]

분류: Proven. 여기서 `q_nuc`는 해당 내부에너지 기준에서 반응 중성미자를 이미 뺀 가열률이다. 조건은 보존되는 동반 바리온 좌표, 전기적으로 중성인 평형 EOS, 빠져나가는 중성미자와 종별 확산 유속 없음이다. 확산을 포함하면 화학 에너지 유속도 포함해야 한다. 정확히 정적인 계량에서 순 열유속이 0이어야 한다는 Request22 결과와 모순되지 않는다. 기준을 바꾸면 `μ^th`도 `∂g/∂Y`만큼 바뀐다.

분류: Imported from prior work. 조성 변화의 내부에너지 항을 빠뜨리면 에너지가 보존되지 않는다는 점은 이후 MESA 문서에도 명시되어 있다. 이는 원리의 교차 확인이며 이후 버전의 구현을 r7624가 실행했다고 주장하는 근거가 아니다. [MESA 조성 항 설명](https://docs.mesastar.org/en/22.11.1/reference/controls.html#include-composition-in-eps-grav).

분류: Counterexample candidate. Request23의 별도 바리온 배율 후보에서 네 물질 구역의 중간 상태를 얻고, 수소·헬륨이 충분한 세 상태에서 직접 FreeEOS 조성 차분을 했다. 중심의 수소 고갈 상태에는 양방향 수소 연소 차분을 적용하지 않았다. 기준 밀도는 `ρ_B`로 고정하며, 조성 변화에 따라 EOS에 넣는 `ρ_atom=C_Xρ_B`도 바꿨다. 19종의 분율을 원소별 수밀도로 집계하는 기존 FreeEOS 사상을 사용하며 별도의 동위원소 혼합 엔트로피는 추가하지 않는다. 실제 중간 상태를 계산했지만 구역 내 보간 및 원래 구역 수의 연속체 오차는 인증하지 않는다.

분류: Counterexample candidate. 원래 핵 가열 정점에 대응하는 GR 구역 2592에서는 `ρ_B=273.114678 g/cm³`, `T=31,501,931.496 K`이다. H 분율을 줄이고 He4 분율을 같은 양 늘리는 방향에서 `du_B/dx≈−4.97581×10^15`, `T ds_B/dx≈−3.19633×10^16`, 두 값의 차이는 `2.69875×10^16 erg/g`이다. 차분 간격을 두 번 줄였으며 세 상태의 최대 마지막 상대 변화는 `2.212×10^-7`이다. 고정 조성 온도 차분은 EOS 열용량 및 `du=Tds`와 대조했다. **유한 차분 수렴은 엄밀한 전역 미분 오차 상계가 아니다.**

## 4. 유한 반응량과 독립 적분 시험

분류: Counterexample candidate. 위 구역의 바리온 밀도를 고정하고 H 분율을 최대 `10^-4`만큼 줄이는 CNO 순반응을 지정했다. 32개 반응량에서 FreeEOS 내부에너지를 역산했다. 이 반응량에는 시간, 반응률, 연소 확률을 부여하지 않았다. 지정 중성미자 손실 외에는 닫힌 정적 부피 시험이며 팽창·전도·복사·대류는 풀지 않는다.

분류: Counterexample candidate. 최대 반응량에서 W 기준 침적 열은 `6.03718532834×10^14 erg/g`, 중성미자 손실은 `4.10641606075×10^13 erg/g`, 최종 온도는 **35,506,279.219698 K**이다. 두 기준의 유한 에너지 역산은 같은 온도를 냈다. 다음 비교의 분모는 모두 올바른 W 기준 침적 열이다.

| 검사 | 최종 상대 에너지 잔차 |
|---|---:|
| 보정한 유한 에너지 역산 | −2.61×10^-15 |
| 보정한 독립 DOP853 미분 적분 | −9.71×10^-14 |
| 기준 보정 `g(X)` 누락 | +6.21335×10^-6 |
| 화학 항을 빼고 `T ds`만 적분 | +0.00481113 |
| 총 에너지에 핵 가열을 중복 계상 | +1.00000000 |

분류: Counterexample candidate. 독립 계산은 반응량 좌표에서 EOS 열용량과 고정 `ρ_B,T` 조성 편미분으로 온도를 적분했다. 두 차분 간격에서 결과를 비교했으며 유한 에너지 역산 대비 온도 상대 차이는 약 `10^-14`였다. 화학 항 누락 시 온도는 `35,525,420.048 K`로 약 **19,140.828 K** 높아졌다. 이는 지정 시험의 결과이며 전체 항성의 광도나 관측 온도 오차가 아니다. [유한 반응 시험](../outputs/reactive-energy24/finite-burn-controls.json), [독립 미분 적분](../outputs/reactive-energy24/differential-control.json).

## 5. 완료 경계와 재현

분류: Proven. 완료한 것은 ① 표준 반응 Q와 기존 정지질량 기준의 연결, ② 조성 의존 에너지·엔트로피 식 정리, ③ 직접 EOS 유한 조성 변화 및 독립 적분 검증이다. 기호 검산, 원래 기록 보존, SHA 연결과 재계산은 `verification/reactive_energy.py`에 있다. 기존 FreeEOS runtime을 재사용하며 MESA 환경을 다시 만들지 않는다.

```bash
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 \
  python3 verification/reactive_energy.py verify
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 \
  python3 verification/reactive_energy.py recheck
```

분류: Proven. `recheck`는 임시 출력 폴더에서 일곱 소스·수치·기호 결과를 다시 계산해 동결 결과와 비교한다. 과거 문서 접두부와 원고 manifest의 기존 항목을 보존한다. 처음 구현할 때 자동 등록 반응을 `reactions.list`만으로 찾지 못한 문제와 중복 가열 음성 대조의 온도 탐색 구간 부족은 검증 전에 수정했으며, 물리적 판정이나 통과 문턱을 바꾸지 않았다. 전체 재실행에서 작은 부동소수점 차이가 있어 JSON 실수의 비트 일치는 요구하지 않는다. 수치 재현 기준은 `|새 값−기록 값|/max(1,|기록 값|)≤10^-12`이며 실제 최대값은 `2.121×10^-13`이었다. 각 시험의 원래 에너지·온도·차분 통과 기준도 재실행하며 원래 결과 파일을 덮어쓰지 않는다.

분류: Conjectural. 다음 순서는 **새 상태의 실제 반응률과 중성미자 → 열수송 → 보존형 GR 열·조성 시간 적분 → 안정성 및 동적 관측 응답**이다. 원래 22종 중 불소를 생략한 FreeEOS의 경계와 전체 EOS의 물리적 오차도 남는다. 이번에 정지질량/Q 연결을 해결한 것을 실제 반응률, 엄밀한 전역 미분 오차 보장, 완전한 비선형 관측 추론의 완료로 세지 않는다.
