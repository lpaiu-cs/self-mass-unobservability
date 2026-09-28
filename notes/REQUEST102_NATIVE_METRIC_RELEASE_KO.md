# 단계102 — 안쪽 경계 오염 수정과 실제 계량 원천 적용

분류: Counterexample candidate. 같은3.434431ms·896/1792셀 대기 유출을 수정한 안쪽 경계로 완료했다. 실제 저장 응력·에너지·바리온 이력을 직접 지연 전하, 계량 응력 원천과 층 내부의 GR 질량 제약 원천에 적용했다. 아래 성분 합은 전체 GR 전하가 아니며 전체 목표는 미완료다.

## 실제로 고친 병목

분류: Counterexample candidate. 이전 ghost가 첫 셀의 배경 밀도·온도를 그대로 복사하면서 첫 셀의 압력 기울기를 평탄화했다. 원1792셀 초기 두 셀 가속도는−411955, +136358cm/s²였지만 공통 배경 면을 적용하면−384, −719cm/s²가 된다. 초기 질량 유속은 원래도0이었으므로 이를 초기 HLL 질량 확산으로 설명하지 않는다. 중력·구면 기하·광자 힘을 없애거나 초기 잔여를 강제로 빼지도 않았다.

분류: Counterexample candidate. 모든 면을 배경 섭동 방식으로 바꾼 첫 시도는 희박 유출면의 EOS 온도 겹침 영역을 벗어나 중단했다. 소스·원 계획·마지막 수락 상태를 보존했다. 문제의 ghost가 영향을 주는 안쪽 두 면만 고쳤으며 나머지는 단계101에서 통과한 재구성을 유지했다. EOS 표·시간 구간·셀 수·수락 기준을 확대하지 않았다.

| 분류 | 항목 | 수정 후1792셀 |
|---|---|---:|
| Counterexample candidate | 실제 적분 | 1586단계, 원 구간 완료 |
| Counterexample candidate | 안쪽으로 넘어간 바리온 | 353764.40g |
| Counterexample candidate | 안쪽 바리온 이송의896/1792셀 차이 | 약0.2221% |
| Counterexample candidate | 바리온 장부 상대 잔여 | 2.92e−16 |
| Counterexample candidate | 에너지 장부 잔여/응답 척도 | 2.24e−15 |
| Counterexample candidate | 유출 질량의 격자 차이 | 0.061086% |
| Counterexample candidate | trace 적분의 격자 차이 | 1.705901% |
| Counterexample candidate | 실제 내부 에너지 이송과 국소 음향 추정의 차이 | 2.36536e15erg |

분류: Counterexample candidate. 이전 안쪽 바리온 이송은1792셀에서5.91942e6g였고896/1792 차이가 약88.55%였다. 이를 수정 후 값으로 대체하되 이전 결과는 보존한다. 내부 에너지의 음향 불일치도2.58817e18erg에서 크게 줄었다. 그 차이가0은 아니며 국소 음향 모형을 실제 내부 전체 해로 수락한 것은 아니다. 수정한 경계 두 셀과 희박 유출을 포함한8개 실제 상태의 native 압력·보존 에너지 대조는 최대1.21e−6 상대 차이로 원0.2% 기준을 통과했다.

## 계량 연결에 사용한 식과 실제 결과

분류: Proven. 정적 등방 배경의 선형 제약에서, 기하 단위·Einstein 반경을 쓰고 f=δφ라 두면 다음이 성립한다. 이는 비선형 실제 항성 전체 해의 존재 증명이 아니다.

```text
b = 1 - 2m/r, Phi = phi_prime, alpha = dlnA/dphi
delta_m = r^2 b Phi f + J
J_prime + r Phi^2 J
  = 4 pi r^2 A^4 [delta_E + (E+P)(3 alpha+r Phi)f]
(ln(mu sqrt(b)/N))_prime = -4 pi r A^4(E+P)/b
mu_prime/mu = r Phi^2
```

분류: Proven. 방사형 boost에서 물질의 E_lab−Pr_lab=ρ(c_X c²+u)−p이므로 계량 응력의 직접 Green 가중치는 A N Phi이고 기존 trace 가중치는 alpha A N/r이다. 일반 원천의 진공 환원은 기존 V=2N²(m/r³−Phi²)와 일치한다. 기호 검사와 제조한 시간 선형 원천의 지연 적분 대조를 통과했다.

분류: Counterexample candidate. 새 이력을 이 식에 실제 적용했다. 안쪽 절단면의 J 입력에는 저장된 바리온·Killing 에너지 차감을 사용하고, 층 안에서는 실제 lab-frame 에너지를 적분했다. 큰 정지 에너지 차감과 적분 인자의 작은 변화를 별도로 계산했다. 다음 값은 u=3.430428ms 끝점의 정규화 전하 **성분**이다.

| 분류 | 성분 | 끝점 값 |896/1792 파형 차이|
|---|---|---:|---:|
| Counterexample candidate | 직접 trace 및 기존 국소 내부 음향 응답 | −2.58976186e−30 |0.109583%|
| Counterexample candidate | 물질의 계량 응력 원천 | −1.08782202e−35 |0.108815%|
| Counterexample candidate | 층 내부의 강제 질량 제약 원천 | −1.09146725e−37 |0.273070%|
| Counterexample candidate | 위 세 성분의 합 | −2.58977284e−30 |전체 GR 오차 아님|

분류: Counterexample candidate. 경계 수정이 직접 파형에 준 최대 변화는0.0009892%다. 이번에 추가한 두 계량 성분은 끝점 직접 성분의 약0.0004243%이며, 계산한 범위에서는 직접 성분을 상쇄하지 않는다. 내부 이송·층·꼬리의 보존 에너지를 합한 순에너지는 실제 기체가 받은 산란 일과2.39e−7 상대 차이로 일치한다. 이 값은 작은 산란 일을 분모로 쓴 별도 대조이며 연속 GR 오차 보증이 아니다.

## 완료 경계와 다음 결정

분류: Counterexample candidate. 위 J 계산은 f=0일 때의 물질 강제 원천이다. f 의존 질량항·scalar 퍼텐셜의 반복, 기체의 계량 되먹임, 실제 내부의 응력 분포, 제거된 꼬리의 위치·응력, 산란 광자의 인과적 계량 원천을 아직 더해야 한다. 질량 제약의 전체 바리온·에너지 재고 확인을 이 누락 항들의 소거 증명으로 바꾸지 않는다. 화학·흡수 폐쇄와 단계95의 원 전체 GR 공간 실패도 남아 있다.

분류: Conjectural. 다음 결정적 작업은 이미 저장된 이 원천을 scalar 퍼텐셜·f 의존 질량항 및 인과적 산란 광자 원천과 함께 풀어, 생략 성분을 포함한 전하 판정을 얻는 것이다. 동일한 기체 이력을 다시 진화하기 전에 저장 결과와 짧은 응답 계산을 재사용한다. 직접 성분의 작은 격자 차이만으로 최종 전하나 정적 EFT 탈출을 선언하지 않는다.

생산 최초100초 중 실패 비용15초를 보수적으로 예약하고 남은85초 안에서 두 경로를51.32초에 완료했다.896셀의10.39초 실측으로1792셀59.52초를 예상한 뒤 실제33.36초로 끝났다. 전하 후처리24.07초는 별도45초, native·기호·제조 원천 대조15.47초는 별도30초 예산 안이다. 생산·후처리의 새 native 호출은0이며 감사에서8회 호출했다. 새 장기 항성 계산은 시작하지 않았다. 보존 에너지·경계·원천 연결의 theorem/loophole progress로 기록하며 전체 목표는 유지한다.

근거: [경계 수정 진화](../verification/def_native_metric_release.py), [실제 계량 원천 적용](../verification/def_native_metric_charge.py), [native·기호·제조 원천 대조](../verification/verify_native_metric_release.py), [진화 결과](../outputs/direct-eos-gr33/def-native-metric-release/result.json), [전하 성분 결과](../outputs/direct-eos-gr33/def-native-metric-release/charge/result.json).
