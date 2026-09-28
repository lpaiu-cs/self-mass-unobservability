# 단계103 — 바리온 부피 되먹임과 퍼텐셜 응답의 조건부 상계

분류: Counterexample candidate. 단계102의 실제 저장 유출을 재사용해 보존 부피의 계량 되먹임을 적용하고 안팎으로 전파하는 scalar 파형을 계산했다. 지정한1792셀 원천·선형 배경 모형에서, scalar 퍼텐셜을 모든 차수로 포함한 끝점 정규화 계수는 다음 범위에 있다.

```text
[-2.5925617092e-30, -2.5869839773e-30]
```

분류: Counterexample candidate. 이 구간은0을 포함하지 않는다. 따라서 **지정한 퍼텐셜 항만으로 기존 성분이 소거되지는 않는다**. 실제 물질 변위의 되먹임, 화학·흡수, 희박 꼬리의 실제 응력, 산란 광자의 지연 원천 및 Bondi 질량 정규화까지 포함한 전체 물리적 전하 구간은 아니다. 전체 목표는 계속 미완료다.

## 보존 제약을 반영한 식

분류: Proven. 기하 단위, 정적 등방 배경, 선형 계량 섭동에서 f=δφ, b=1−2m/r, Phi=phi_prime이라 놓는다. 좌표 셀의 바리온 재고·엔트로피·조성을 고정한 단열 계량 반응은 다음과 같다. 추가 유체 변위와 가열은 eF,pF에 속하며 이 식만으로 풀리지 않는다.

```text
delta_m = r^2 b Phi f + J
delta_lnV = (3 alpha+r Phi) f + J/(r b)
delta_E = eF - (E+P) delta_lnV
delta_P = pF - Gamma1 P delta_lnV

J_prime + (nu_prime+lambda_prime) J = 4 pi r^2 A^4 eF
nu_prime+lambda_prime = r Phi^2 + 4 pi r A^4(E+P)/b
```

분류: Proven. 따라서 적분 인자는 N/sqrt(b)이고 J는 저장된 강제 Killing 에너지의 누적에 sqrt(b)/N을 곱해 얻는다. 계량을 바꾸면서 Eulerian 물질 밀도를 그대로 두면 나타나던 별도 f 의존 원천은 바리온 부피 반응을 포함해야 상쇄된다. 원 단계102의 고정 밀도 성분 계산을 전체 보존 계량 반응으로 해석하지 않는다.

분류: Proven. U=r f, dx=dr/(N sqrt(b))에서 파동식은 U_tt/c²−U_xx+Veff U=S_forced이다. 부피 반응에 따라 Veff=Vfixed−Dv(3alpha+rPhi)/r, KJ_eff=KJ+Dv/(rb)이며 Dv=4pi N²A⁴r[alpha(E+P−3Gamma1 P)+rPhi(E+P−Gamma1 P)]다. 기호 대입으로 질량식·부피 응력식을 확인했다. 이 선형화가 비선형 유출 전체의 물리적 폐쇄를 보증하지는 않는다.

## 실제 적용과 퍼텐셜 오차

분류: Counterexample candidate. 이전 두 격자의 보존 이력, 같은 native EOS와 같은 시간 구간을 사용했다. 새 기체 진화나 native EOS 호출은 없다. 기존 trace·계량 응력 원천의 outgoing/ingoing Green 적분을 계산하고 위 보존 J를 적용했다. 단계102와 보존 J의 차이는 저장 산술 수준이다.

분류: Proven. 지연 Green 연산자를 K, 자유 원천 해를 U0라 하고 ||U0||≤M0라 하면, eta=(c Delta_u/2)∫|Veff|dx<1 아래에서 다음이 성립한다.

```text
||U-U0|| <= eta M0/(1-eta)
||U-U0-K U0|| <= eta^2 M0/(1-eta)
```

분류: Counterexample candidate. 내부 인과 영역에는 저장된 선형 보간 배경의 극값과 N≤1, A≤1, m≤M을 적용했다. 외부에는 양의 질량 진공의 해석 상계를 사용했고 무한 반경 적분에는 임의의 외부 절단을 두지 않았다. 노름은 큰 정지 에너지의 관측 상쇄를 이용하지 않는 절댓값 상계다. 부동소수점 입력 축약의 여유를 더한 뒤 상계 산술은 바깥 방향 구간 연산으로 계산했다.

| 분류 | 1792셀 결과 | 값 |
|---|---|---:|
| Counterexample candidate | 연산자 수축 상계 eta | 3.19338e−8 |
| Counterexample candidate | 자유 장의 절댓값 노름 상계 M0 | 2.54738e−21cm |
| Counterexample candidate | 퍼텐셜 모든 차수의 정규화 변화 상계 | 2.78887e−33 |
| Counterexample candidate | 위 상계/자유 파형 최대값 | 0.107688% |
| Counterexample candidate | 계산한 첫 Born 보정의 끝점 | +2.70792e−38 |
| Counterexample candidate | 첫 보정을 더한 끝점 추정 | −2.5897728162e−30 |
| Counterexample candidate | 두 원천 격자의 보정 파형 차이 | 0.109583% |
| Counterexample candidate | 두 퍼텐셜 구적의 상대 차이 | 9.16e−11 |

분류: Counterexample candidate. 첫 Born 보정은 원천의 좁은 지지 영역 밖에서 직접 적분했다. 생략한 중앙 영역의 첫 보정은 정규화1.20346e−37 이하, 그 다음 모든 Born 항은8.90590e−41 이하라는 조건부 상계를 따로 남겼다. 두 구적의 작은 차이를 엄밀한 구적 오차라고 부르지 않는다. 위 최상단의 최종 구간은 **구적 추정에 의존하지 않는 전체 퍼텐셜 노름 상계**로 얻었다.

## 확인한 범위와 남은 실제 물리

분류: Counterexample candidate. 독립 감사에서 저장 계수의 양성·lapse·질량 전제를 확인했다. 촘촘한 계수 적분은 해석 상계 안에 있었고, 부호가 양수·음수인 제조 점 퍼텐셜의 정확한 해는 Neumann 나머지 상계를 만족했다. 두 저장 격자 모두 원2% 파형·퍼텐셜 영향 기준을 통과했다. 격자 차이는 연속 해 오차의 엄밀한 인증이 아니다.

분류: Counterexample candidate. 기존 보존 중심화는 제거한 희박 물질의 기준 질량을 x=0에 둔 수학적 점 원천과 동등하다. 이번 노름에는 이 점 원천도 포함했다. 실제 꼬리의 위치·압력·화학 상태는 여전히 별도 문제다. 단계102 NPZ의 scalar_source_coefficient_cm3 이름은 차원 표기가 잘못됐으며 값의 실제 단위는cm⁻²다. 수치식은 맞았고 원 파일은 보존했다. 새 파일은 coefficient_cm_minus2를 사용한다.

분류: Counterexample candidate. 이번 결과로 지정한 퍼텐셜의 고차 반복이 기존 성분을 없앨 가능성은 제한했다. 반면 기체의 실제 계량 유도 변위, 내부의 전체 에너지·운동량 접합, 희박 꼬리, 산란 광자와 질량 분모, EOS의 순간 화학평형·흡수 가정은 이 구간 밖이다. 원 전체 GR 공간 실패도 해소됐다고 표시하지 않는다. 이것은 보존 계량 식의 theorem progress와 실제 파형 적용의 loophole progress다.

분류: Conjectural. 다음 우선순위는 저장된 원천을 다시 계산하는 것이 아니라 남은 물리적 원천의 크기와 접합을 닫는 것이다. 기존 유한 이온 점유·광자 교환 코드를 새 대기 궤적에 재사용하고, 산란 광자에는 계산 시작 전에 이미 존재한 광자까지 포함한 인과적 재고·질량 수지를 적용해야 한다. 단순히 새로 방출한 광자만 세는 상계로 전체 외부장을 대신하지 않는다. 동일 퍼텐셜에 더 큰 격자나 장기 적분을 반복할 근거는 이번에 얻지 못했다.

생산은 등록한120초 중44.37초, 별도 독립 감사는4.00초였다. 첫1792셀 작업20.56초로 다음896셀 작업을31.73초 이내로 예상했고 실제19.82초에 끝났다. 생산 전 절댓값 노름에서 기준 꼬리 점 원천을 빠뜨리지 않도록 수정했으며 실행되지 않은 원 계획·소스를 해시 일치 상태로 보존했다. 생산 실패나 계산 범위 확대는 없었다.

근거: [보존 파동 계산](../verification/def_native_conserved_wave.py), [독립 감사](../verification/verify_native_conserved_wave.py), [파형 결과](../outputs/direct-eos-gr33/def-native-conserved-wave/result.json), [조건부 상계](../outputs/direct-eos-gr33/def-native-conserved-wave/bound-1792.json), [구간의 적용 범위](../outputs/direct-eos-gr33/def-native-conserved-wave/audit.json).
