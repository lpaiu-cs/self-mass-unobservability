# 단계97 — 같은 해의 계량·열유속·표면 광도·이동 광자 되먹임

분류: Counterexample candidate. 새 EOS 내부·외피의 결합 경로에 **반경별 lapse와 열유속 기하, 표면 면적·적색편이 광도 및 이동 광자 되먹임을 매 단계 함께 연결했다.** 실제32/64단계 시간 대조, 독립 lapse 적분과 이동 에너지 장부를 통과했다. 최종 전하와 고립된 항성의 전체 표면 접합은 아직 미완료다.

## 같은 물질·scalar 해로 lapse와 열유속을 복원

분류: Proven. `b=1−2m/r`, `A4=A⁴`, `Phi=phi′`, `psi=delta phi`인 반경 좌표에서 lapse 제약의 선형 변화는 다음과 같다. delta p는 Eulerian 압력이며 같은 EOS의 물리 Lagrangian 압력에서 배경 이동을 빼서 계산한다.

```text
delta nu′ = (1+8 pi r² A4 p) delta m / (r² b²)
            +4 pi r A4 (delta p+4 alpha p psi)/b +r Phi psi′

delta nu_s = Phi_s psi_s −2 integral_s^infinity(Phi psi/b) dr
             −integral_s^infinity[J/(r² b²)+4 pi r Pr/b] dr
```

분류: Proven. 고정된 수송 계수에서 적색편이 온도 `A N T`와 열유속 prefactor의 변화를 함께 미분하면 내부에는 `op*(delta lnT+alpha*psi+delta nu)` 및 `L0*(delta nu+2alpha*psi−delta m/(r b))`가 들어간다. 표면에는 다음 추가 항이 필요하다.

```text
delta L/L0 = 4 Delta lnT +(2+4 alpha Phi_s+2 nu′_s) zeta
             +4 alpha psi_s+2 delta nu_s
```

분류: Counterexample candidate. 같은 압력·질량·scalar 해에서 이 lapse를 전 반경에 복원하고, 내부 유속과 표면 방출에 함께 넣었다. 매 시간 단계의 현재 방출량·유속·표면 변위·속도로 지연 광자를 다시 평가하고 모든 추가 항을 반복 결합했다. 기준 열유속 분리와 셀별 변화분 열수지는 단계96의 수정식을 유지한다. 이 결과는 원 opacity·열전도도·조성·EOS 계수가 고정된 접선 모형이다. 수송 계수의 물질 상태 의존성 전체나 비선형 배경 재진화를 인증한 것은 아니다.

## 실제 결과와 독립 대조

분류: Counterexample candidate. 반경·물질 재고·원천 보간·p4 일관 질량·3.434431ms 구간을 유지했다. 추가 결합은 매 단계 두 번의 풀이로 등록된1e-8 상대 문턱을 통과했다. 최대 수렴 잔여는3.33e-22이며 연산자 전체의 원 선형/국소 열수지 검사도 통과했다.

| 분류 | 실제 계산량 | 결과 |
|---|---|---:|
| Counterexample candidate | 온도32/64 상대 차이 | 0.77245% |
| Counterexample candidate | 속도32/64 상대 차이 | 0.14800% |
| Counterexample candidate | scalar32/64 상대 차이 | 0.03096% |
| Counterexample candidate | 끝점 표면 delta nu | −4.42571e-29 |
| Counterexample candidate | 끝점 내부 최대 abs(delta nu) | 1.19660e-26 |
| Counterexample candidate | 내부 lapse 독립8점 구적 차이 | 7.21e-14 상대 |
| Counterexample candidate | 표면 lapse 독립 구적·96각도 차이 | 0.05566% |
| Counterexample candidate | 전체 lapse 최대값 정규화 차이 | 0.0002060% |

분류: Counterexample candidate. 끝점 표면 광도 섭동 `delta L/L0`는7.76075e-29에서2.74300e-29로 변했다. 새 변화분은 이전 변화분의35.3446%다. 기준 광도에 대한 섭동 자체는약10^-29 수준이다. 반면 내부 온도·속도·scalar의 전체 최대값으로 정규화한 전후 차이는6e-16 이하였으며, 이를 별도의 물리 신호 검출이나 정확한0으로 해석하지 않는다. 단계95의6.912% 내부 속도 공간 실패를 이 시간 대조로 수락하지 않는다.

## 이동 에너지·운동량 장부와 남은 물리 경계

분류: Proven. 국소 comoving 광자 유속F, 압력Pr와 외부 지지 압력Pg가 있을 때, 움직이는 물질 표면과 광자의 에너지 유속 차이는 `v Pg`, 법선 운동량 유속 차이는 `Pg`다. 따라서 압력을 유지하는 외부 성분은 이 두 양을 공급해야 한다. 단순히 Pg를 고정했다는 조건은 그 성분의 중력 응력 텐서를 정하지 않는다.

분류: Proven. 반구 광자에는 `Egamma=3Pr`이고, 이동 광자 누적질량의 추가 표면 항은 `delta J_gamma=−4 pi A4 (Egamma+Pr) zeta`다. 지지 성분의 에너지 저장 변화 `delta m_support=4 pi A4 Pg zeta`를 명시하면 물질 쪽 total pressure와의 Lagrangian 질량 접합 장부가 맞는다. 이는 필요한 장부 항등식이며 외부 성분의 상태방정식이나 scalar 응력을 정한 것이 아니다.

분류: Counterexample candidate. 현재 실제 끝점 변위에 대해 독립적인 표면 광선 모멘트가 이 질량 항을 상대1.49e-14로 재현했다. 지지 성분에 필요한 무차원 질량 변화는1.63972e-61이다. 이 작은 변화량을 외부 성분의 전체 정적 질량·응력 또는 전하 오차 상계로 바꾸지 않는다.

분류: Conjectural. 다음 물리 병목은 유한 압력 절단 밖의 실제 대기 또는 그 영향의 통제다. 현재의 유지 압력만으로는 외부 중력·scalar 응력이 정해지지 않으며, 복사 응력을 포함한 scalar 접합 역시 완전히 닫히지 않았다. 최종 전하에는 이 경계, 내부 공간 오차의 전하 영향, 무한대 outgoing 읽기 및 같은 물리 모형의 정적 비교가 연결돼야 한다. 이번 결과는 그중 계량–수송 되먹임을 실제로 연결한 loophole progress이며 목표 전체를 완료로 바꾸지 않는다.

실행 기록: CPU1스레드, 새 native EOS·배경 근0회, 최대 측정 메모리약0.750GB. 실제 진화87.59초, 독립 감사9.72초로 내부 수치 계상은97.31초다. 등록 한도240초를 확대하지 않았으며 Python import·프로세스 시작·문서 작업은 이 내부 계상과 구분한다. 원3.434ms 이외의 장기 진화나 고차/격자 확대는 실행하지 않았다.

기억 기록: 단계94–96에서 확인한 OSK MCP 엔진 불일치가 해소됐다는 근거가 없어, 이번 통합 내용도 저장소에 대기시킨다. 장기 기억 갱신 완료로 보고하지 않는다.

근거: [실제 결합 구현](../verification/def_native_metric_transport.py), [독립 검산](../verification/verify_native_metric_transport.py), [생산 결과](../outputs/direct-eos-gr33/def-native-metric-transport/result.json), [독립 감사](../outputs/direct-eos-gr33/def-native-metric-transport/audit.json).
