# 단계131 — 실제 대기 끝점의 native 보존 역산과 충돌 연결

분류: Counterexample candidate. 단계130에서 새 열·반응 입력으로 완주한 fine 경로의 활성 대기223셀 모두를 native EOS로 독립 역산했다. 밀도만 고정하는 대신 저장된 바리온 D, 운동량 S, Killing 에너지 K, 중성수소 재고를 함께 유지했다. native 압력 차이는최대1.4359e-8, 온도의 로그 차이는최대2.9092e-8이다. 실제 끝점에서 대기 보간을 다시 수정해야 할 큰 오차는 발견되지 않았다.

분류: Proven. 같은 D,S,tau에서 Q=cx*D+tau, v=S/(Q+p)다. 따라서 압력 교체의 정확한 trace 변화는(v_old*v_new-3)*delta_p이고 반경 응력 변화는(v_old*v_new-1)*delta_p다. 전체 에너지와 바리온 변화는0이다. 속도 변화항을 빠뜨리거나 초기 압력 차이를 임의로 빼지 않았다. 이 대수 항등식은 symbolic 검사를 통과했다.

| 분류: Counterexample candidate — 검사 | 결과 |
|---|---:|
| fine 끝점 전수 역산 | 223/223셀 |
| 추가 보존된 pilot 상태 | 초기2셀,coarse1셀 |
| 에너지 역산 최대 상대 잔차 | 6.4948e-13 |
| 운동량 최대 상대 잔차 | 1.9218e-16 |
| native 밀도 최대 상대 잔차 | 1.4211e-14 |
| 대기 trace 원천 보정L1 | 1.2595e17erg |
| 기존 대기 끝점 응답L1 | 9.9014e24erg |
| trace 원천 보정/기존 응답 | 1.2720e-8 |

분류: Counterexample candidate. 저장한 native 준위 분율·affinity·전자 밀도를 실제 끝점 광자와 이동 충돌 코드에 넣었다. 보존 역산 후의 압력·온도·속도와 전자 산란 변화도 포함했다. 양의 채널은최대1.7679e-6 상대 차이, 순H교환은7.6563e-5(0.0076563%), 순 Killing 에너지 교환은6.4746e-6(0.00064746%), 반경 운동량 교환은2.8120e-8이었다. 원 양의 채널0.2percent·순 교환2percent 기준을 통과했다. 실제 충돌의 물질/광자/주파수 밖 에너지 교환 합은최대7.08e-15 상대 잔차로 닫혔다. 유체 수송은 두 판독에서 함께 제외하여 국소 충돌 비교임을 명시한다.

초기128셀과 두 끝점각223셀, 총574상태의 원 계획은 실측 상한420.53s로180s 예산에 실패했다. 원pilot의eligible=false와 생산자를 보존하고 fine 끝점223셀로 범위를 줄였다. 완료한fine3셀을 재사용하여220셀을 추가했고 등록 축소 예상은약169s, 실제 생산68.17s/native2387호출이었다. pilot10.19s/151호출, 실제 충돌 판독8.16s였다. 충돌 판독에는 추가 native 상태 계산이나 새 시간 적분이 없었다.

분류: Counterexample candidate. 생산 결과의forecast/eligible 필드는 미실행인 원574상태 전체 범위를 외삽한 값이다. 축소된223셀 실행의 수락은execution-plan.json의eligible=true와production.json의passed=true다. 예산 실패를 성공으로 바꾸거나180s 한도를 늘리지 않았다.

분류: Counterexample candidate. 이것은 실제 새 fine 끝점의 EOS·반응 원천 확인이다. 초기와coarse의 전 셀, 새 경로의 중간 보존 상태·광자 분포, 연속 EOS/미분 오차는 모두 인증하지 않았다. 끝점 자료로 임의 시간 이력을 만들어 새 retarded 전하를 계산하지 않았다. 단계130의 양의 조건부 구간은 그대로이며 이번 숫자를 그 구간의 전체 시간 오차 상계로 넣지 않는다.

분류: Conjectural. 이번 끝점 차이는 대기 보간 수정이나 같은 물리 궤적 재실행의 근거가 되지 않는다. 다음 우선순위는 단계130의 실제 새 원천으로 GR 귀환을 물질·광자에 다시 연결하여 결합 고정점과 삭제 물질의 후속 효과를 판정하는 일이다. 대기 전체 이력의 오차와 작은 물질 포트 공간 수렴은 열린 조건으로 유지한다. 미시적 모델 자체의 정확도·완전한 비선형GR·정적 모형으로의 흡수 여부·관측 연결도 미완료다.

근거: verification/def_native_atmosphere_inverse.py, verification/verify_native_atmosphere_inverse.py, outputs/direct-eos-gr33/def-native-atmosphere-inverse/, native-atmosphere-inverse-manifest.json.
