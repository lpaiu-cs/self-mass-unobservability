# 단계133 — 수정된 광자 전달량의 자유 유체·GR 귀환

분류: Counterexample candidate. 단계132의 실제 충돌 에너지·중성수소·운동량 전달을 수정 EOS의 공유 물질 유속에 적용했다. 원531셀,3.434ms,64/128 시간 경로를 유지하고 세 경로를 완주했다. 추가 바리온·운동량·기준 에너지·중성수소를 진화한 뒤 직접 보존 primitive 역산으로 압력·비등방 응력·trace를 구해 실제 compact GR에 적용했다. 이번 결과는 loophole progress이며 최종 전하 인증은 아니다.

| 분류: Counterexample candidate — 검사 | 결과 |
|---|---:|
| 물질 시간 차이 최대 | 0.089710% |
| 물질 배경 경로 차이 최대 | 0.219070% |
| 응력 배경 차이 최대 | 0.242779% |
| 전 이력 보존 잔차 최대 | 1.167673e-15 |
| 압력 직접 probe 차이 최대 | 5.497006e-11 |
| GR 시간/배경 차이 | 0.000874% / 0.137854% |
| 독립 GR 직접 적분 일치 | 4.662937e-15 |
| 추가 compact 전하 | 1.645511424011e-44 |
| 단계130 기존 성분 대비 | 1.611453229913e-17 |

분류: Counterexample candidate. 같은 경로의 원 유체 유속 대조,1e-8 보존,0.2% directional/압력 및2% 시간/배경 기준을 통과했다. 기존 작은 비영 donor 분기를 유지했고 실제 변화량/분기 유량 비는최대약2.87e-14다. reference64는 앞서 확인한 힘 차분 소실에 대한8배 probe를 처음부터 등록하고4/8 및8/16 대조를 유지했다. 새로운 상태에서8/16 비교최대약5.78e-5도0.002 안이었다. 원 forward probe 지표는최대2.44로 작지 않다. 이를 통과로 바꿔 쓰지 않으며, 수락은 Richardson 결과의 두 probe 크기 대조다. 연속 EOS/모든 미분의 균일 보장은 아니다.

분류: Proven. 유지하는 내부 비상대론적 운동에너지 T=유효관성질량*v²/2에서 반경 운동 응력은2T다. 앞서 누락 항목으로 따로 측정했던 이 응력을 실제 radial metric-work 항에 포함했다. 대수 항등식의 symbolic 검사를 통과했으며 상대론적 보정 자체의 오차 인증은 아니다.

분류: Counterexample candidate. 원천 export에서 바리온·에너지·trace·계량 응력뿐 아니라 nonrest 반경 응력과 접선 압력도 새 값으로 채웠다. 이전 배경 값이 템플릿에 남아 새 원천처럼 보이지 않도록 했다. 초기 canonical GR 응답을 중복 계산하지 않고 실제 광자 압력/에너지·반경 포트를 함께 사용했다. 원천 항등식 잔차는1.444593e-16, 입력 결속22건을 확인했다.

분류: Counterexample candidate. 추가 compact 전하가 작아도 현재 결합은 닫히지 않았다. 직전 광자 풀이가 사용한 물질 기준 에너지/H 이력과 실제 자유 유체 반환 이력의 상대 불일치는1.022374/0.439882다. 이는 유체 수치 오차나 관측 오차율이 아니라, 아직 서로 재적용하지 않은 두 하위계 이력의 차이다. 광자 풀이에 새로운 질량·운동량·재고와 비충돌 기계 수송을 넣는 다음 단계가 필요하다. 작은 추가 전하 한 번으로 전체 잔차를 상계하거나 수렴률/수축성을 선언하지 않는다.

생산은 독립3CPU,각1스레드·가상 메모리2GiB 상한으로 실행했다. 후반 CFL까지 포함한 이전 실제 raw-call 수와 새 동시 pilot의 처리율을 사용했고, 상한179.90s가350s 안인 뒤 실행했다. 실제 물질 생산74.96s/350s, CPU 합계175.79s, 최대 RSS 합0.786GiB였다. 원천 변환24.60s/75s, GR25.16s/120s다. 물리 배경·광자 궤적·native 상태는 다시 계산하지 않았다.

분류: Conjectural. 다음은 이 새 유체 운동을 실제 이동 광자 충돌에 적용해 에너지/H와 동시 풀이하고, 반환한 교환량을 다시 물질·GR에 넣어 잔차를 줄이는 것이다. 그 후에도 새 GR장 재적용·외부 scalar·삭제 물질 수송·전체 EOS 이력/미분·비선형 GR·정적 비흡수성·관측 연결을 판정해야 한다. coupled_fixed_point_verified와final_charge_solved는false이며 원 목표는 계속 활성 상태다.

근거: verification/def_native_updated_material_return.py, verification/verify_native_updated_material_return.py, outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/material-return/ 및 native-updated-material-return-manifest.json.
