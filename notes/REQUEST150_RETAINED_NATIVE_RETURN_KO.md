# 단계150 — 현재 native EOS의 실제 광자·물질·GR·무한원 반환

분류: Counterexample candidate. **현재 저밀도 물질 경로에서 native EOS 압력 힘과 충돌 차이를 실제 광자·열·H 및 유한 물질 진화에 넣고, 그 응답을 compact GR와 무한원 정규화 전하까지 반환했다.** source-only 교체에서 빠졌던 실제 기계적 압력 힘을 적용한 loophole progress다. 최종 정규화 전하 변화의 명목 값은 `5.2381168299432462e-27`이며 양수를 유지한다. 전체 연구의 완료 판정은 아니다.

분류: Counterexample candidate. 기존 531셀(내부19·대기512), 17개 저장 배경, 64/128 응답 시간 계산과 3.4344311179287023ms를 그대로 사용했다. 기존 배경을 재진화하지 않았고, 일치하는 보존 셀 역산을 재사용했다. 재구성된 활성 경계면에는 실제 native 압력·에너지를 평가했다. 독립 검산에서 경계면 상태 6,392개의 밀도 일치를 확인했다. native 음속 미분을 새로 인증한 것은 아니며 기존 표의 gamma와 내부 face_K를 명시적으로 유지했다.

분류: Counterexample candidate. 물질 방정식은 저장 배경 Q에 대한 `F_table(Q+deltaQ)-F_table(Q)+[F_native(Q)-F_table(Q)]`와 실제 광자 충돌 전달을 더한 유한 결함 방정식이다. native 압력의 초기 오프셋을 빼지 않았다. 같은 재구성·HLL·기하 원천·공유 경계면·SSP 시간 규칙을 유지했다. 이 구성식 교정은 저장 배경에서 평가하며 새 물질 운동이 다시 광자와 계량을 바꾸는 반복은 이번 계산에 포함되지 않았다.

다음 수치의 분류는 모두 Counterexample candidate다.

| 판독 또는 검증 | 결과 | 원 기준과 범위 |
|---|---:|---|
| compact GR 추가 교정 | `8.26273103e-34` | 실제 native 힘·광자·물질 응답 |
| compact 교정의 시간 차이 | `0.00013907147` | 0.02 미만 |
| compact 교정의 구적 차이 | `5.17555699e-16` | 0.002 미만 |
| 독립 직접 GR 판독 차이 | `1.03511294e-15` | 1e-9 미만 |
| 원 전체 시점의 배경을 유지한 compact 끝점 | `1.0211270550483363e-27` | 단계149 배경+이번 교정 |
| 무한원 정규화 추가 교정 | `3.40613247e-34` | 실제 signed 방출 패킷 포함 |
| 무한원 정규화 최종 명목 값 | `5.2381168299432462e-27` | 이전 `5.2381164893299992e-27` |
| 정규화 추가 교정의 시간 차이 | `0.000270354082` | 0.02 미만 |
| 작은 외부 scalar 교정의 각도 차이 | `0.000333131116` | 0.002 미만 |
| native 힘 유량의 최대 산술 지표 | `0.000137467278` | 0.002 미만, 균일 오차 정리 아님 |
| 물질·광자 원천 독립 항등식 차이 | `9.79759196e-17` | 1e-12 미만 |
| 무한원 정규화 독립 재구성 차이 | `2.73939981e-16` | 1e-12 미만 |

분류: Counterexample candidate. 광자 64/128 경로는 모두 끝났고 최대 시간 비교는 `0.00173749117`다. 물질은 각각 397/775개 실제 CFL 하위 단계를 마쳤다. 에너지·H 이력의 광자-물질 잔차는 최종 경로에서 `[4.4034428351672085e-08, 1.2291976643559308e-09]`다. 이를 결합 연산자의 수축 상계로 해석하지 않는다. 원 infinitesimal 선형 상태 기준은 통과하지 않으므로 실제 유한 진폭 물질 방정식을 사용했다.

분류: Counterexample candidate. compact 적용 대조의 `previous_same_cadence_endpoint=1.0209754059927525e-27`는 17시점으로 재표현한 배경이다. 이 수치를 단계149의 전체129시점 끝점으로 바꿔 쓰지 않았다. 최종 무한원에서는 단계149의 전체 시점 compact 파형을 정확한 공통 시각에 대응시키고 이번 교정만 더했다. 두 응답 시간 계산에도 동일한 fine 배경 방출을 사용했다.

## 실패를 해결한 위치

분류: Counterexample candidate. 첫 native 압력 힘 계산은 큰 float64 유량의 차이를 취할 때 후기 neutral-H의 산술 지표가 약28.15%여서 실패했다. 뺄셈 직전에만 긴 정밀도로 바꾸는 것으로는 해결되지 않았다. 동일한 conserved/HLL/중력/내부 유량 연산자의 대수 자체를 longdouble로 계산해 원 0.2% 기준을 통과했다. EOS·원시 복원·재구성 좌표는 binary64 그대로이므로 이 산술 수정은 EOS나 복원 오차의 엄밀한 인증이 아니다.

분류: Counterexample candidate. 기존17개 배경 지연 구간에 signed SDIRK 패킷을 그대로 적분한 무한원 결과는 추가 교정의 각도 차이1.28%, 도착 에너지1.43%, 외부 항의 반경 차이1.50%로 원 기준에 미달했다. 각 패킷의 실제 도착 경계로 적분 영역을 제한하고, 도착 각도 적분은 정확한 면적식으로 계산했다. 그래도 선형 mu 좌표의 외부 항 각도 차이0.22155%가 남아 실패했다. 방사 방향 가까이의 `1-mu` 척도를 `z=log(1-mu)`와 정확한 Jacobian으로 처리하자 원4/8차·물리4빈을 유지한 채 모든 원 기준을 통과했다. 과거 세 실패 결과는 별도 원본으로 보존했으며 수락 기준을 완화하지 않았다.

분류: Proven. 하나의 공유 유량은 양옆 셀 보존식에서 정확히 상쇄된다. 각도 빈에서 도착 경계를 mu0로 자르면 광자 에너지 각도 가중치는 `(hi^2-max(lo,mu0)^2)/2`의 양의 부분이다. `q0=(s0+alpha*epsilon0)/(1-epsilon0)`일 때 새 정규화 증분은 `(ds+(alpha+q0)*delta_epsilon)/(1-epsilon0-delta_epsilon)`이다. 이 항등식들을 기호 검산했다. 이것은 구성식·공간·연속 시간 또는 비선형 GR 오차 상계가 아니다.

## 재사용과 계산 비용

운영 기록: 최초 충돌·압력 pilot·계수 bank의 시간 초과를 보존했다. WSL 공유 경로의 파일 읽기 비용을 줄이기 위해 SHA가 같은 입력과 중간 산출물을 ext4에 배치했고 원 배경의 중복 사본은 발행에서 제외했다. 최초603개 경계면 상태의 raw native 호출 총수는 기록되지 않았으며 추정값으로 채우지 않는다. 이후 기록된 native 호출은 15,328회다.

운영 기록: native 생산의 uncached 실측 비용을 바탕으로 한 두 배 벽시간 예측은 508.324초, CPU 예측은 1269.642초였다. 저장된 상태로 산술만 다시 조립한 빠른 시간을 새 native 호출의 예측에 사용하지 않았다. 이후13시점 생산은 87.499초였다. 광자 생산 524.113초, 물질 생산 63.952초였고 각각 승인된 prefix를 재사용했다. 물질의 두 배 생산 예측 451.909초는 원780초 안이었다.

운영 기록: 마지막 도착 경계·각도 좌표 계산도 대표3개 age의 실측 후 나머지 예측 36.883초가 미사용분 안인 것을 확인했다. 최종 독립 검산을 포함한 원천/GR/무한원 판독 그룹은 74.235초로 원180초 안이다. 모든 보존된 실패·재개·진단·staging을 포함해 영수증에 청구한 과학 계산 벽시간은 1316.422초, CPU는 1409.852초로 각각 원2880초 이내다. 이는 코드 작성·WSL 시작 지연·문서·Git까지 포함한 세션 벽시간이 아니다. 최대3개 native worker와 각1개 BLAS 스레드, 프로세스당3GiB 주소 공간 제한을 사용했다.

## 판정과 남은 직접 연결

분류: Counterexample candidate. **현재 native 구성식 교정은 실제 광자·기계 운동을 거쳐 최종 판독까지 적용해도 명목 양의 잔여를 지우지 못했다.** 그러나 source-only에서 한 단계 더 진행한 유한 반환이지 완전 결합 해의 고정점은 아니다. 새 물질 운동의 광자 재귀환과 새 GR의 진화 재적용, 균일 native 음속/충돌 미분 오차 및 연속 상태 범위, 더 작은 floor의 전체 피드백, 완전 비선형 GR, 같은 재고 정적 모형과 자유 미분 nuisance를 제거한 관측 식별성은 남아 있다. 기존 단계148의 작은 floor 교정 자체가 시간 기준에 미달한 판정도 유지한다.

분류: Conjectural. 다음 직접 레버는 저장한 유한 물질 운동을 같은 광자 충돌·수송에 반환해 전하 변화가 얼마나 달라지는지 판정하는 것이다. 새 전체 배경·격자·장시간 경로보다 현재 보존 이력과 경계 패킷을 재사용한다. 작은 현재 귀환으로 전체 수축 또는 관측 검출을 선언하지 않는다.

분류: Counterexample candidate. `native_collision_applied`, `native_pressure_force_applied`, `free_material_evolved`, `represented_GR_applied`, `actual_new_angular_emission_at_null_infinity`, `actual_emission_mass_normalization`은 true다. `actual_material_motion_returned_to_photons`, `native_acoustic_derivatives_certified`, `coupled_fixed_point_verified`, `final_charge_solved`, `full_goal_complete`는 false다. 현재 명목 값에 이전 단계의 조건부 구간을 자동 이전하지 않았다.

재현 자료: `outputs/direct-eos-gr33/retained-native-return/completed/`가 최종 산출물을 소유한다. `result.json`은 compact 판독, `infinity/result.json`은 최종 패킷 판독, `verification.json`은 독립 대조다. 원래 폴더의 실패 영수증과 `completed/infinity-uncut/`, `completed/infinity-packets-linear-mu/`는 실패 당시 결과다. `<f16` 보존 배열은 WSL에서 읽는다. 검산 프로그램은 `verification/verify_retained_native_return.py`이며 `RETAINED_NATIVE_OUTPUT`을 completed 디렉터리에 지정한다. OSK에는 재사용할 결론과 계산 규율을 보존했고 세부 실행 기록은 이 저장소에 남겼다.
