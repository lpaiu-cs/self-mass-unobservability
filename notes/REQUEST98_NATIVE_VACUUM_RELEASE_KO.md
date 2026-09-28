# 단계98 — 유한 압력 절단의 native 비선형 진공 팽창

분류: Counterexample candidate. 저장된 항성 표면의 기체 압력 지지를 제거하고, **복사를 뺀 동일 native EOS로 비선형 희박파를 계산했다.** 기체가 기존 반경 밖으로 실제 이동하는 국소 해와 보존 응력 입력을 얻었다. 단계97의 압력 유지 경로를 고립된 자유 표면 해로 간주할 수 없다는 점이 구체화됐다. 이번 해는 국소 평면·단열·화학평형 조건부이며 구면 GR, 광자 교환과 최종 전하의 완결 해는 아니다.

## 바뀐 물리 경계와 유지한 입력

분류: Proven. 진공과 접한 무응력 기체 물질면은 기체 압력0을 요구한다. 양의 압력을 외부 성분이 유지하는 정지 경계와, 그 지지를 없앤 진공 팽창은 다른 초기·경계 문제다. 후자는 기존 작은 표면 변위에 경계 계수만 보정해서 얻을 수 없다.

분류: Counterexample candidate. 단계92 최종 외피의 반경·밀도·온도·조성·에너지 기준을 재사용했다. 표면 기체 압력은1 dyn/cm², 반경은약6.90985e9 cm, 온도는약13368.84 K다. FreeEOS의 기존 기체 전용 분기를 사용해 LTE 복사 성분을 제거했다. 다항식 EOS나 동결 이온식으로 교체하지 않았다. 밀도비1에서exp(−18)까지 native 기체 엔트로피를 보존하는73/145상태를 계산했다. 이 온도·밀도 범위에서의 화학평형을 실제3.434 ms 안에 달성하는지는 인증되지 않았다. 압력 지지 제거는 명시한 실험 조건이며 실제 천체의 주기적 구동 발견으로 해석하지 않는다.

분류: Proven. `c_s`와 `v`를 광속 단위로 쓸 때, 동일 엔트로피 기체의 상대론적 희박파는 다음 관계를 따른다. 이는 국소 이상화에 대한 관계이며 중력·구면·가열 항을 이미 푼 방정식은 아니다.

```text
d ln T/d ln rho = (P/rho − u_lnrho)/cvT
c_s² = Gamma1 P/(e+P)
Y(rho) = integral_rho^rho0 c_s d ln rho
v = tanh Y,   xi = (v−c_s)/(1−v c_s)
```

분류: Imported from prior work. 상대론적 Riemann 문제의 희박파를 자기유사 파동으로 연결하는 방법은 [Martí–Müller의 원 연구 설명](https://www.uv.es/astrorela/simulacionnumerica/node38.html)을 따른다. 해당 문헌의 이상기체 결과가 이번 native EOS나 항성 조건을 인증하는 것은 아니다.

## 실제 진공 팽창과 검증

| 분류 | 항목 | 결과 |
|---|---|---:|
| Counterexample candidate | 안쪽으로 전파하는 희박파 머리의 속도 | 17.8943 km/s |
| Counterexample candidate | 계산된 바깥쪽 최대 속도 | 152.3853 km/s |
| Counterexample candidate | 3.434431 ms 뒤 머리의 원 표면 대비 위치 | −61.4564 m |
| Counterexample candidate | 원 절단 반경 밖으로 이동한 기체 | 5.71763e11 g |
| Counterexample candidate | 마지막 native 온도 | 3735.325 K |
| Counterexample candidate | 속도73/145상태 상대 차이 | 0.01948% |
| Counterexample candidate | 바리온 적분 상대 잔여 | 0.001182% |
| Counterexample candidate | 비정지 에너지 적분 상대 잔여 | 0.001426% |
| Counterexample candidate | 운동량 적분 상대 잔여 | 0.03258% |
| Counterexample candidate | 엔트로피 정규화 잔여 | 1.60e−12 |

분류: Counterexample candidate. 원 수락 기준과 독립적인 정확한 상대론적 감마 법칙 대조를 통과했다. 감마 법칙은 수치 도구의 대조에만 쓰였으며 생산 EOS를 대체하지 않았다. 압력 유지 경로의 거의 정지한 표면과는 질적으로 다른 유동이다. 원 EOS 수렴 실패나 내부 속도 공간 실패를 이번 국소 성공으로 변경하지 않는다.

분류: Conjectural. 미계산 저밀도 꼬리에서 `c_s(rho)<=c_s,cut*(rho/rho_cut)^0.05`를 추가로 가정하면 진공 전면은 원 표면에서523.35–935.58 m 범위이며 꼬리 물질은1.944e5 g 이하로 제한된다. 이 수치는 **조건부 꼬리 가정의 결과**이고 native 진공 끝점 인증이 아니다. 실제로 계산된 바깥쪽 상태의 위치는 약502.74 m이며,523.35 m는 그 상태의 유속을 이용한 완전 진공 전면 위치의 하한이다.

## 큰 정지질량 상쇄를 보존한 GR 원천 연결

분류: Proven. `E=(e+P)W²−P`, `S=(e+P)W²v`, `W=(1−v²)^−1/2`일 때 Eulerian 응력의 trace는 `e−3P`다. 자기유사 에너지 보존식의 원시함수와 전체 적분 관계는 다음과 같다. `H`는 초기 기체가 존재한 반공간을 표시한다.

```text
K_E = xi*[E−E0 H(−xi)]−S
dK_E/dxi = E−E0 H(−xi)

integral delta(e−3P) dxi
    = −integral [S v+3(P−P0 H(−xi))] dxi
```

분류: Proven. 마지막 식은 전체 에너지 보존 및 완전한 진공 꼬리를 포함하는 조건에서 성립한다. 정지질량의 큰 두 적분을 따로 빼는 대신 이 식으로 trace를 계산할 수 있다. 가변 가중치에는 `integral w delta(trace)=[w K_E]−integral w' K_E−integral w[Sv+3 delta P]`가 필요하다. 기호 검산을 실행했다.

분류: Counterexample candidate. 현재 표에서 정지질량의 구적 잔여를 그대로 전하에 넣으면 약−2.02e−27이라는 가짜 원천 변화가 생긴다. 보존식으로 계산한 국소 trace 에너지는−1.92458e24 erg이고, 실제 표면의 `alpha A N` 및 면적을 고정한 **정규화 scalar 원천 모멘트**는+2.18076e−32다. 73/145상태 차이는0.17413%다. 이는 시간 의존 파동을 무한대까지 보낸 전하가 아니며, 연속 오차의 엄밀한 상계도 아니다.

분류: Conjectural. 꼬리에 위 음속 가정과 `0<=u<=u_cut`, `0<=P/rho<=P_cut/rho_cut`를 추가하면, 그 원천 모멘트의 꼬리 기여 절댓값은1.65e−36 이하로 제한된다. 다른 근사·복사·역반응의 오차를 이 꼬리 상계로 대신하지 않는다.

분류: Counterexample candidate. 실제 반경에 희박파의 에너지·운동량·법선/접선 압력을 배치하고 얇은 층의 선도 GR 제약에 연결했다. 저장된 질량 섭동 최대값은3.99847e−17 cm, lapse 섭동 최대값은약5.44e−33이다. 이 원천 묶음은 반경별 계수 변화, 미계산 꼬리와 양방향 유체·광자 되먹임을 닫지 않았으며 독립 구면 수렴을 주장하지 않는다. 기존 결합 진화에 그대로 더하면 바깥 기체를 이중 계산할 수 있으므로, 다음 결합은 해당 층의 물질 응력을 교체·접합해야 한다.

## 실제 남은 경계

분류: Counterexample candidate. 현 구간에서 중력 충격/초기 음속은0.001054, 최대 조건부 두께/반경은1.443e−5, 원 압력 기울기로 평가한 머리 이동/압력 척도는0.001670이다. 코드의 역사적 키 `initial_density_scale_height_fraction`은 실제로 **압력** 척도를 계산한다. 이 수치는 국소 근사의 규모 지표이며 해 오차의 증명은 아니다. Rosseland 광학 척도1.98e−8 역시 흡수·가열 또는 반응률 상계가 아니다.

분류: Conjectural. 다음 결정적 작업은 이 자유 팽창 층을 원 압력 유지 경계 대신 결합 진화에 접합하고, 구면 보존·광자 에너지 교환·화학 반응 시간의 효과를 제어하는 것이다. 그 뒤 같은 모델의 무한대 전하와 정적 비교가 필요하다. 이번 국소 모델을 전체 목표로 축소하지 않는다. `full_goal_complete`, `final_charge_solved`, `full_spherical_GR_fan_coupled`는 모두false다.

## 실행 예산과 원 실패 보존

첫 실행은1600회 native 호출 한도에서 종료됐다. EOS나 수락 기준 실패가 아니라 호출 예산 실패다.73개 거친 상태는 저장됐지만, 세밀한 표의993회 호출 결과는 개별 저장되지 않았다. 원 소스·계획·실패를 보존했다. 이후 원 fine 격자·영역·물리 기준을 그대로 두고73상태를 재사용하며 엔트로피의 정확한 미분 `dS/dlnT=cvT/T`로 빠진72상태만 풀었다. 상태마다 재시작 캐시를 저장했다.

추가 예산은400회/40초였고 실제226회/4.976초였다. 누적 native 호출은1826회다. 별도 보존/GR 원천 검산은 native0회/0.223초로30초 한도 안이었다. 최초 실패 실행의 정확한 전체 벽시간은 저장되지 않아 총 소요시간을 꾸며내지 않는다. 알려진 최초 coarse 내부 시간은10.838초다. 새 전체 별 적분·추가 밀도 영역·고차 경로·장기 계산은 실행하지 않았다.

근거: [진공 팽창 구현](../verification/def_native_vacuum_release.py), [보존식·GR 원천 검산](../verification/verify_native_vacuum_release.py), [실제 결과](../outputs/direct-eos-gr33/def-native-vacuum-release/result.json), [원천 검산 결과](../outputs/direct-eos-gr33/def-native-vacuum-release/audit.json), [원 호출 예산 실패](../outputs/direct-eos-gr33/def-native-vacuum-release/budget-stop.json).
