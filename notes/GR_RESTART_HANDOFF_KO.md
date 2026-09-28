# 최종 GR 인계 — 2026-09-17 09:15 KST, 승인된 세 경로 및 최종 판정 완료

분류: Counterexample candidate. 70/140/280 세 경로가 전체 0.42117120910640804초를 완주했고 원 끝점 시간 수렴 기준을 통과했다. 밀도/온도/속도/총 열유속/복사 열유속의 측정 차수는 1.974655/2.035766/2.003874/2.062330/2.062330이다. 원 native 잔차·정확 BDF 이력·보존·특성속도와 24개 반복 기록 한도를 유지했다. 마지막 경로 적분은 08:30:17 KST, 전체 재생·비교 및 프로세스 종료는 08:34:03 KST다. 정상 종료 코드 0이며 원 결과는 c31098f0에 보존했다.

분류: Counterexample candidate. 종료 후 manifest 결속 513개 파일과 원 계획·소스 결속을 확인했고, 저장 끝점의 차이·차수를 다시 계산해 동일 판정을 얻었다. symbolic 검사도 통과했다. 결과 분류는 조건부 유한 GR 모형의 전체 기간 끝점 수렴 시험을 통과한 loophole progress다. 과거의 수치 실패와 운영 중단은 그대로 보존한다.

분류: Conjectural. 엄밀한 연속시간 오차 상계, 물리 EOS 및 새 EOS 보존 초기화·실제 진화 연결, 물리 대기·외부 경계, 반응·실제 구동·비선형 관측 연결은 미완료다. 유한 격자 수렴을 실제 항성 또는 관측 완성으로 해석하지 않는다.

운영 상태: 실행 중인 GR 계산은 없다. SelfMassGR-Consoleless-20260917 작업은 정상 종료 확인 후 정리했고 gr-156-2 정기 점검을 삭제했다. 승인된 70/140/280 비교는 종료하며 더 촘촘한 격자·더 긴 기간·새 실험을 자동 실행하지 않는다.

근거: outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-resume271-consoleless의 time-refinement.json, execution-manifest.json, final-review.json, closure-check.json. 상세 결과와 해석은 [GR 실험 최종 판정](GR_EXPERIMENT_REDESIGN_KO.md)을 따른다.

이하 실행 중·미완료·예상시각 기록은 모두 역사이며 현재 상태는 위 완료 판정이다.

# 최신 GR 인계 — 2026-09-17 06:32 KST, 사용자 터미널 분리 및 마지막 9단계 재개

운영 기록: production-resume238은 05:54:59 KST에 마지막 경로 271/280단계를 저장한 뒤 실행기와 WSL이 종료됐다. 첫 70/140 경로는 완료됐고 원 출력은 71eaeb83에 보존했다. production-resume271도 06:16:39 KST에 원 이력 123단계 재생 중 종료되어 a9bfbb6d에 보존했다. 후자는 새 적분 0단계이며 원본 271단계를 잃은 것이 아니다. 두 시도 모두 native failure.json과 실행기 종료 파일은 없고 Windows 부팅 시각은 유지됐다.

사용자가 해당 시각에 터미널 창을 닫았다고 직접 확인했다. 숨김 PowerShell 예약 작업은 터미널 종료와 함께 중단됐으며 Task Scheduler 소유권과 WindowStyle Hidden만으로 충분한 분리가 아니었다. 정확한 OS 제어 신호 전달 경로를 별도 추적한 것은 아니다. 실행기만 GUI Python과 CREATE_NO_WINDOW 방식으로 바꾸었다. 콘솔 검사에서 pythonw 실행기는 콘솔이 없고 WSL은 창 없는 전용 콘솔에 해당 계산의 WSL 프로세스 두 개만 포함한 것을 확인했다. CREATE_NO_WINDOW를 콘솔 객체 자체가 없다는 뜻으로 해석하지 않는다. 사용자의 터미널을 실제로 다시 닫는 시험은 하지 않았다.

현재 실행은 verification/gr_resume271_consoleless.py chain --workers 15, 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-resume271-consoleless다. 06:25:48 KST에 시작했다. 완료된 70/140 경로와 수락된 271단계를 원 파일 hardlink로 검증 재생했고 원 272단계 실제 계산 진입을 확인했다. 새 적분은 원 272~280의 9단계와 최종 비교뿐이며 이미 수락된 이력은 재적분하지 않는다.

분류: Counterexample candidate. bce2851f의 계획·실행기는 70/140/271 연결점 native 잔차 0.002705579 / 0.895389962 / 0.000900528, 정확 BDF 상태·누적 유속 복원, 변조 유속 거부와 symbolic 검사를 통과했다. 원 방정식·donor 해법·70/140/280 시각·24개 반복 기록 한도·보존·특성속도·최종 시간 차수 문턱은 동일하다. 이는 실행 수명 문제에 대한 운영 수정이며 새로운 물리적 성과가 아니다.

- 실행 식별: SelfMassGR-Consoleless-20260917. Schedule 3356 → pythonw.exe 52888 → WSL 클라이언트 47680. Linux PID 313, boot_id 2c3e281e-7a3c-4c5e-8b2e-5f342a4722be, start_ticks 10189. 자동 트리거 없는 수동 작업이다. 실행 중 소스·계획은 변경하지 않는다.
- 계획 SHA: eb7ca8d7fb57efb548aedbd7dfcda36450ba8de2695c48c828fe2e404804779c. plan.json, restart-check.json, windows-launcher.json, windows-wsl-child.json, launcher-check.json, console-isolation-check.json, launch.json, replay-launch-check.json을 근거로 한다.
- 예산: 최근 20/10/5단계 평균은 9.11/8.16/8.09분, 최근 10단계 범위는 2.78~20.30분이다. 향후 8~15분/단계 가정으로 새 적분 1.2~2.25시간이며 재생 여유를 포함한 예상 종료는 오늘 07:45~09:00 KST다. 후반 비선형 비용·추가 중단은 미보장이며 최종 재생·비교 시간은 별도다. CPU 0~15, 15 worker, BLAS 1 thread, GPU 미사용을 유지한다. 새 실행 RSS 합은 3.42 GiB로 공유 페이지 중복을 포함한 관측값이며 최대치나 고유 물리 메모리가 아니다.
- 판정: 같은 조건부 GR 모형의 전체 기간 시간 격자 의존성을 확인한다. 다섯 변수의 끝점 차이가 감소하고 lnT/Qtotal/Qrad 차수가 1.5 이상이어야 한다. 미달하면 원 실패를 보존하며 수락 기준을 완화하지 않는다. native·보존·특성속도 실패 시 중단한다. 어느 경우에도 더 촘촘한 격자·긴 기간·새 실험을 자동 확대하지 않는다.
- 사용자 터미널 분리 확인은 Windows 재부팅·로그오프·명시적 프로세스 종료·WSL shutdown 내성을 뜻하지 않는다. 기존 gr-156-2의 조용한 1시간 점검은 새 실행을 추적하며 최종 판정 후 삭제한다. 중단된 이전 수동 작업만 정리한다. 재사용할 실행 수명 교훈은 기존 OSK 결합 GR 상세 노드에 정정했다.

분류: Conjectural. 마지막 경로 완주·세 경로 시간 수렴과 물리 EOS·외부·반응·관측 연결은 아직 미완료다.

이하 기록은 역사이며 최신 실행은 위 경로를 따른다.

# 최신 GR 인계 — 2026-09-17 00:15 KST, 238단계 복원 및 239단계 재개

운영 기록: 직전 production-reboot-recovered는 마지막 경로 238/280단계를 9월 16일 23:43:44 KST에 저장한 뒤 중단됐다. 첫 70/140 경로는 완료 상태다. Windows 실행기는 사라졌고 작업은 Ready/3221225786이며 WSL 부팅 식별자가 바뀌었다. native failure.json과 실행기 종료 파일은 없다. Windows 부팅 시각은 9월 16일 07:06 그대로이고 최근 재부팅 이벤트도 없어, 앞선 OS 업그레이드 재부팅과 구분한다. 종료 주체는 미확정이며 Task Scheduler 상세 이력은 비활성화되어 있었다. 원 파일·운영 근거는 3866ef51의 interruption-20260917.json과 windows-runtime-20260917.json에 보존했다.

현재 실행은 verification/gr_resume238.py chain --workers 15, 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-resume238다. 00:09:32 KST에 시작했다. 완료된 두 경로와 수락된 238단계를 각각 약 38/71/116초에 검증 재생했고, 모든 재사용 저장 상태의 원 파일 hardlink 일치를 확인했다. 00:14 KST에 원 239단계의 실제 계산이 시작됐다. 원 239~280단계만 적분한다. 새 적분은 42단계이며, 원 시각과 두 이전 BDF 상태·누적 유속을 유지한다. 미완료 상태를 수락 이력에 넣지 않는다.

분류: Counterexample candidate. 연결점 70/140/238의 native 잔차 점수는 0.002705579 / 0.895389962 / 0.019836321이다. 정확 BDF 상태·누적 유속 복원, 변조 유속 거부와 symbolic 검사를 통과했고 복구 코드·계획을 557be051에 고정했다. 방정식·donor 해법·원 70/140/280 시각·24개 반복 기록·보존·특성속도·최종 시간 차수 문턱은 바꾸지 않았다. 이번 조치는 운영 복구이며 새로운 물리적 검증 결과가 아니다.

- 실행 식별: 수동 작업 SelfMassGR-Resume238-20260917. Schedule 3356 → PowerShell 19872 → foreground wsl.exe 27144. Linux PID 312, boot_id 9c560e03-cc70-4aed-bb84-a7dbba8f8714, start_ticks 9364. 실제 명령·경로·식별을 대조했고 자동 트리거는 없다. 종료된 이전 수동 작업은 정리했다.
- 계획 SHA: 7949f293526c357aeee9adda72496ee7f624b00e806c231a79a695a7465f599c. 근거는 plan.json, restart-check.json, launcher-check.json, launch.json, replay-launch-check.json이다. 기존 계획·출력은 보존하며 실행 중 소스·계획을 수정하지 않는다.
- 재개 예산: 남은 42단계. 최근 미세 경로 20/10/5단계 평균은 8.77/10.02/9.65분, 최근 10단계 범위는 4.36~15.97분이다. 향후 10~15분/단계 시나리오로 7~10.5시간, 9월 17일 07~11시 KST를 예상한다. 후반 비선형 비용·추가 중단은 미보장이고 최종 재생·비교 시간은 별도다. CPU 0~15, 15 worker, BLAS 1 thread, GPU 미사용을 유지한다. 재개 직후 RSS 합은 약 3.37 GiB이며 공유 페이지 중복을 포함한 관측값으로 최대치나 고유 물리 메모리가 아니다.
- 승인 범위는 기존 세 경로와 동결 최종 판정까지다. native·보존·특성속도 실패 시 그대로 중단·보존하며 추가 격자·기간·실험을 자동 확대하지 않는다. 최근 중단 주체를 모르는 상태에서 런처 소유권 확인을 원인 해결이나 WSL 종료 내성으로 표현하지 않는다. 기존 gr-156-2의 조용한 1시간 점검은 새 실행을 추적한다.

분류: Conjectural. 전체 기간 세 경로 시간 수렴과 물리 EOS·외부·반응·관측 폐쇄는 여전히 미완료다.

이하 기록은 역사이며 최신 실행은 위 경로를 따른다.

# 최신 GR 인계 — 2026-09-16 10:56 KST, 두 번째 전체 기간 경로 완료

운영 기록: Windows System 로그에서 07:03/07:06 종료·재시작과 TrustedInstaller의 계획된 OS 업그레이드 재부팅을 확인했다. 직전 production-budget-repaired는 첫 경로 70/70, 두 번째 134/140단계까지 저장했고 미완료 135단계 반복 21에서 끊겼다. native failure.json과 Windows 실행기 종료 파일은 없으며, 작업은 Ready/3221225786이고 기존 프로세스가 사라졌다. 첫 재시작의 개별 요청 이벤트까지 특정한 것은 아니다. 원 파일·시스템 이벤트·중단 기록은 f3380ae8에 보존했다. 이 중단을 수치 실패로 분류하지 않는다.

현재 실행은 verification/gr_reboot_recovery.py chain --workers 15이며 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-reboot-recovered다. 08:04:29 KST에 시작했다. 완료된 첫 경로는 수치 적분 없이 약 30.4초에 검증 재생했고 두 번째 수락 이력 134단계는 약 53.4초에 재생했다. 두 경로의 모든 저장 상태가 원 파일과 동일한 hardlink임을 확인했다. 미완료 135단계부터 다시 계산하며 복구 시작 시 추가 적분은 두 번째 6단계와 세 번째 152단계, 총 158단계였다. 저장되지 않은 반복은 수락 이력에 넣지 않는다.

분류: Counterexample candidate. 두 번째 140단계 경로가 2026-09-16 10:49:37 KST에 전체 0.42117120910640804초를 완료했다. 완료 manifest의 145개 파일 SHA·계획 결속과 모든 단계의 원 24개 반복 기록 한도를 확인해 aefec58f에 보존했다. 마지막 native 잔차 점수는 0.8953899624로 원 수락 기준 이내다. 첫 두 경로의 완주와 최종 세 경로 시간 수렴은 구분하며, 같은 chain이 마지막 280단계 경로를 계산 중이다. 완료된 두 경로를 다시 적분하지 않는다.

분류: Counterexample candidate. 연결점 70/134/128의 31개 native 잔차 점수는 0.002705579 / 0.001932835 / 0.092894838이다. 원 BDF 상태·누적 유속의 정확한 복원, 변조 유속 거부와 symbolic 검사를 통과했다. 새 계획·복구 코드는 eac15a4f에 고정했다. 원 70/140/280 시간 격자, donor 수정, 24개 반복 기록 한도, 잔차·보존·특성속도·최종 시간 차수 문턱은 동일하다. 이번 검사는 운영 복구이며 새로운 물리적 결과로 세지 않는다.

- 실행 식별: Windows 수동 작업 SelfMassGR-Resume134-20260916. Schedule 서비스 3356 → PowerShell 6224 → foreground wsl.exe 13232. Linux PID 311, boot_id bbe896ba-03b3-4628-931b-bb0e0e01cd8d, start_ticks 8050. 실제 명령·경로·식별과 계획 SHA를 대조했다. 자동 트리거는 없으며 종료된 이전 수동 작업은 정리했다.
- 복구 계획 SHA: c438da9a98b672727d841f75e8fa8a28a40fe9b0dc50a2821a60d90495430a9e. 근거: restart-check.json, launcher-check-20260916.json, linux-launch-check.json, replay-launch-check.json. 실행 중인 소스·계획은 수정하지 않는다.
- 원 검증 목적은 같은 조건부 유한 GR 모형에서 전체 기간의 시간 격자 의존성을 판정하는 것이다. 다섯 변수의 끝점 차이가 감소하고 lnT/Qtotal/Qrad의 관측 차수가 1.5 이상이어야 한다. 통과하면 이 모형의 수치 근거를 보완하고, 미달하면 실패를 보존한다. 어느 경우에도 자동으로 더 촘촘한 격자나 새 장기 실험을 시작하지 않는다.
- 복구 시에는 남은 158단계와 미측정 미세 경로의 15~20분/단계 가정으로 9월 18일 00~17시를 예상했다. 2026-09-16 13:22 KST 실측 갱신: 첫 두 경로는 완료, 마지막은 175/280 수락·176 계산 중이며 105단계가 남았다. 최근 20/10/5단계 평균은 각각 5.66/6.21/7.14분이다. 최근 둔화를 고려해 향후 8~12분/단계 시나리오를 쓰면 약 14~21시간, 9월 17일 03~10시 KST다. 후반이 15~20분/단계로 느려지는 시나리오는 17일 16시~18일 00시다. 후반 비선형 비용은 미측정이며 통계적 신뢰구간·보장 시각이 아니다. 최종 재생·비교와 새 중단 시간은 별도다. 실행 소스·계획·문턱·계산 범위는 바꾸지 않는다.
- CPU 0~15, 15 worker, BLAS 1 thread를 유지하고 GPU는 사용하지 않는다. 08:07의 부모+worker RSS 합은 약 3.35 GiB이며 공유 페이지 중복을 포함한 시작 구간 관측값이다. 실제 고유 메모리·최대치로 해석하지 않는다. 실패 시 원 문턱에서 중단하고, 예상보다 크게 느려지면 남은 승인 범위의 비용·가치를 재평가하며 반복·기간·해상도를 자동 확대하지 않는다.
- 하네스 외부 실행과 Windows 재부팅 내성은 별개다. 이번 복구는 OS 업데이트 정책을 바꾸지 않는다. 기존 gr-156-2의 조용한 1시간 점검을 새 실행으로 갱신했다. 최종 세 경로 판정 뒤 예약을 삭제한다.

분류: Conjectural. 전체 기간 세 경로 시간 수렴과 물리 EOS·외부·반응·관측 폐쇄는 아직 남는다. 장기 계산의 주장·정확도·저비용 대안·실측 예산·중단 기준 선행 원칙은 AGENTS.md와 기존 OSK 결합 GR 진화 노드에 보존돼 있다.

이하 기록은 역사이며 최신 실행은 위 경로를 따른다.

# 최신 GR 인계 — 2026-09-15 16:44 KST, 첫 전체 기간 경로 완료

운영 기록: production-harness-recovered는 09:01:35 KST에 native 54단계 반복 한도 소진으로 종료됐다. 저장된 53단계는 유효하다. failure.json과 windows-runner-exit.json의 exit_code=1이 있으므로 앞선 WSL 수명 중단과 구분한다. 실패와 마지막 iterate를 9390050b에 보존했다.

분류: Counterexample candidate. 기존 donor 보정은 선 탐색 실패 때만 호출된다. 54단계에서는 잔차가 조금씩 감소해 그 조건 없이 0..23의 24개 기록을 소진했으며 마지막 잔차는 5.241030758002932였다. 같은 저장 상태의 donor 보정 대조는 0.0027873065998553315를 얻었다. 그러나 한도 밖에서 수행한 대조이므로 이 상태는 수락·재사용하지 않는다.

현재 실행은 verification/gr_budget_repaired_full_duration.py chain --workers 15이며 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-budget-repaired다. 09:41:04 KST에 시작해 53단계까지 정확 재생한 뒤 54단계를 원 한도 안에서 수락했다. 첫 70단계 경로는 15:47:19 KST에 전체 0.42117120910640804초를 완료했다. 16:42 점검에서는 두 번째 경로 79/140단계가 저장됐고 80단계를 계산 중이다. 실제 Windows 계보와 Linux 실행 식별은 동일하며 새 실패 파일은 없다.

분류: Counterexample candidate. 완료된 path-1/manifest.json의 75개 파일 SHA와 원 계획 결속을 확인했고 모든 수락 단계가 원 24개 반복 기록 이내임을 점검했다. 완료 경로를 2ad5c00e에 보존했다. path-1/result.json의 completed=true는 단일 경로의 완주이며 세 경로 시간 수렴 판정은 아니다. 같은 chain이 140/280단계 후속 경로를 계속 계산하므로 첫 경로를 재실행하지 않는다.

분류: Counterexample candidate. 54단계의 처음 23개 기록은 원 실패 실행과 같고 마지막 보정만 변경되어 기록 23에서 잔차 0.0025248383846417357로 통과했다. 저장 상태의 31개 native 잔차, 정확 BDF 누적 유속, 국소 에너지·바리온·핵종 수지와 특성속도도 재생했다. 근거는 production-budget-repaired/step54-accepted-replay.json과 190181ee에 보존했다.

- 새 해법은 미수렴한 반복 기록 22 직후 마지막 허용 보정을 기존 donor 해법에 배정하고 결과를 기록 23으로 판정한다. 초기 기록 0과 최대 23회 보정으로 원 24개 기록 한도를 유지한다. 24번째 추가 보정은 하지 않는다.
- 원 70/140/280 시간 격자와 53/64/128 수락 이력, 31개 native 잔차·보존·특성속도·최종 시간 차수 문턱은 유지한다. 새 물리식이나 counterfactual 상태를 넣지 않는다.
- 수정 소스·메커니즘 대조: 52baf888. 실행 계획·연결점 native/BDF 재생·변조 거부·symbolic 검사: b4850370. 새 계산은 총 245단계다.
- Windows 수동 실행 작업 SelfMassGR-Budget24-20260915의 Schedule 서비스 → PowerShell → foreground WSL 계보를 확인했다. Windows launcher PID 67308, wsl.exe PID 47724, Linux PID 390, boot_id 2f3bc8e2-f2ed-4b2d-a640-6d9a7285a11e, start_ticks 25345이며 launch.json과 실제 값을 대조한다.
- launch-check.json에 실행 소유권·원 파일 해시·53단계 복원 확인을 보존했다. 종료된 이전 수동 작업은 정리했고 gr-156-2의 1시간 점검은 새 경로로 갱신했다. 실행 중인 소스와 계획을 바꾸거나 중복 실행하지 않는다.
- 한 경로 완료 및 세 경로 시간 수렴은 별도 판정이다. 새 실패 시 원 상태를 보존한다. 현재 실행은 WSL 강제 종료·Windows 재부팅·로그오프·전원 종료를 견디도록 설계한 것이 아니다.

분류: Conjectural. 전체 시간 수렴과 물리 EOS·외부·반응·관측 폐쇄는 아직 남는다. 현재 결과는 원 한도 안의 실제 실패 단계 복구를 확인한 loophole progress다.

이하 기록은 역사이며 최신 실행은 위 경로를 따른다.

# 최신 GR 인계 — 2026-09-15 08:45 KST, WSL 중단 후 복구

운영 정정: 아래 08:14의 하네스 재시작 가능 판단은 충분하지 않았다. 재접속 뒤 WSL 부팅 식별자가 달라졌고 기존 계산이 사라졌다. 원 실행은 53단계까지 수락·저장했고 54단계 반복 15에서 중단됐으며 native 실패 파일은 없다. WSL 종료를 유발한 주체는 확인되지 않았다. Linux session 분리와 파일 로그만으로 WSL 자체의 수명을 보장할 수 없다는 경계를 확인했다.

현재 실행은 verification/gr_harness_recovery.py chain --workers 15이며 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-harness-recovered다. 08:38:46 KST에 시작했고 53단계까지 정확 재생한 뒤 54단계를 다시 계산하고 있다. Windows 작업 SelfMassGR-Resume53-20260915는 트리거 없는 수동 실행 전용이다. Schedule 서비스 → PowerShell → foreground wsl.exe 계보를 실제 확인했고 하네스 프로세스에 속하지 않는다. Windows launcher PID 65236, wsl.exe PID 52252, Linux 계산 PID 316이며 PID만으로 식별하지 않는다.

- 원 실행 중단 근거와 닫힌 파일은 production-donor-repaired/interruption-20260915.json 및 fdd7fbef에 보존했다. 완료되지 않은 54단계 반복은 재사용하지 않는다.
- 복구 소스는 229e5898, 원 시각의 53/64/128 수락 이력·계획·복원 검사는 636bd2cd에 고정했다. 원 70/140/280 격자에서 추가 적분은 17/76/152, 총 245단계다.
- 세 연결점 native 잔차는 0.0031088194 / 0.4212937549 / 0.0928948380이며 정확 BDF 상태·누적 유속 복원, 변조 유속 거부, symbolic 검사를 통과했다. 방정식·검증된 donor 수정·수락 오차·24회 반복·보존·특성속도·최종 시간 차수 문턱은 그대로다.
- 실제 실행은 windows-launcher.json, launch.json의 boot_id·start_ticks·명령·경로와 대조한다. launcher-check-20260915.json에는 Windows 프로세스 계보와 Schedule 서비스 소유 확인을 보존했다. Linux boot_id는 f477aae0-085f-4ef1-868e-26f98a64bbb0, start_ticks는 7817이다.
- 하네스/앱만의 종료가 계산을 끊지 않도록 실행 소유권을 옮겼다. 실제 하네스 재시작을 통한 지속 검증은 아직 하지 않았다. WSL 강제 종료·Windows 재부팅·로그오프·전원 종료는 보호 범위가 아니다. 앱이 닫힌 동안 Codex의 시간별 점검도 실행되지 않는다.
- 기존 gr-156-2의 1시간 점검을 새 경로로 갱신했다. 실행 중 소스·계획을 바꾸거나 chain을 중복 시작하지 않는다. 세 경로 완료 뒤 원 native 재생과 동결 시간 차수로 판정하며 계산 종료 후에만 수동 Windows 작업을 정리한다.

분류: Counterexample candidate. 수락된 이력 복원과 중단 지점의 실제 재계산을 확인했다. 분류: Conjectural. 전체 기간 수렴과 물리 EOS·외부·반응·관측 폐쇄는 여전히 미완료다. 이번 조치는 운영 복구이며 새 물리적 검증 통과로 세지 않는다.

이하 기록은 역사이며, 이전 재시작 가능 판단은 위 정정과 최신 실행 경계를 우선한다.

## 하네스 재시작 직전 점검 — 2026-09-15 08:14 KST

운영 기록: 원 첫 경로 53/70단계가 저장되어 있고 54단계 계산 중이다. 새 실패 파일은 없다. PID 387의 boot_id·start_ticks·명령·작업 경로를 launch.json과 대조했다. 계산의 session/process group은 자체 PID이고 제어 터미널은 없다. stdin은 /dev/null, stdout/stderr는 run.log다. 부모는 WSL의 /init Relay이며 하네스 셸 프로세스가 아니다. 마지막 두 수락 상태의 ZIP 무결성과 SHA-256도 확인했다.

운영 판단: 이 분리 구조에서 하네스/앱만 재시작하면 계산이 계속될 것으로 판단한다. WSL 종료·배포판 terminate·Windows 재부팅은 포함하지 않는다. 실제 하네스 종료 실험을 실행한 것은 아니다. 1시간 heartbeat 설정은 ACTIVE로 저장되어 있으며 로컬 예약 점검에는 앱 실행이 필요하다. 재시작 뒤 같은 launch.json의 실제 프로세스를 먼저 확인하고 기존 chain을 중복 실행하지 않는다.

근거: outputs/gr-harness-restart-check-20260915.json. 예약 동작 참고: https://learn.chatgpt.com/docs/automations?surface=app

# 최신 GR 인계 — 2026-09-15 05:21 KST 재개

운영 기록: 현재 실행은 verification/gr_donor_repaired_full_duration.py chain --workers 15다. 출력 디렉터리는 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-donor-repaired다. 실행 직후 42단계까지의 정확한 재생을 마치고 새 43단계 계산에 진입했다. 현재 진행률은 각 path-1/2/4/progress.json과 iterations.jsonl, run.log에서 읽는다.

- 실행 소스/복구 결과 체크포인트: 95d46b7f. 재개 계획/수락 이력/복원 검사: 2008f420.
- 실제 식별은 launch.json의 pid, boot_id, start_ticks, command, cwd로 대조한다. PID 387은 과거 실행에서도 쓰인 값이므로 PID만으로 같은 실행이라 판단하지 않는다.
- CPU 0–15, 15작업자, OPENBLAS_NUM_THREADS=1. 백그라운드 계산을 중복 실행하지 않는다.
- 계획은 원 70/140/280단계와 원 시각을 유지한다. 42/64/128 수락 이력을 재사용하고 28/76/152단계를 새로 계산한다.
- 세 경로 완료 뒤 native 저장 상태/누적 유속/31개 잔차/보존/특성속도 재생과 동결 시간 차수 문턱을 적용한다. 최종 time-refinement.json과 execution-manifest.json으로 판정한다.
- 활성 heartbeat gr-156-2는 1시간 점검이다. 변화 없으면 반복 보고하지 않고 의미 있는 실패/수정/최종 판정을 한국어로 알린다.
- 실행 중인 코드·계획을 수정하지 않는다. 실패하면 보존된 실제 상태에서 수정 대조와 새 실행을 준비하며 원 파일을 덮어쓰지 않는다.
- gr_step42_subdivision.py는 두 번째 절반에서 실패해 종료됐다. gr_subdivided_full_duration.py의 99/198/396 진입 조건은 충족되지 않았다. 이를 실행하지 않는다.

분류: Counterexample candidate. 원 42단계 실패는 donor 전환과 반복 접선의 불일치였다. donor·조성 미분의 분기를 일치시키되 원 native 조성 기준점의 affine 항을 유지한 수정이 정규화 잔차 18.768935482897138을 0.002270952319166062로 낮추고 실제 수지·특성속도를 통과했다. 세 연결점 복원/변조 거부/symbolic 검사를 통과했다. 실패 진단을 다시 확장하기보다 실제 남은 결합 진화를 완료한다.

분류: Conjectural. 전체 시간 수렴, 물리 EOS·외부·반응·실제 구동·관측 폐쇄는 미완료다. 이전 실패와 짧은 예비실험 통과는 별도 판정으로 보존한다.

이하 기록은 이전 실행의 역사다. 현재 명령과 활성 경로는 위 인계를 우선한다.

# GR 계산 재시작 인계 — 2026-09-14 19:32 KST

하네스 재시작 전에 실제 WSL 프로세스와 저장 파일을 확인했다. 현재 실행 중인 해당 GR/분자 Python 계산은 없다. 코드 수정이나 재계산을 시작하지 않은 상태에서 사용자가 재시작을 요청했다.

## 저장 상태

- 작업 디렉터리: `E:\lab\self-mass-unobservability`
- 결과 디렉터리: `outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/pilot-cached-16-32-64`
- 기존 16/32단계와 새 64단계는 모두 `result.json`의 `completed=true`이다.
- 새 64단계 결과 저장 시각: **2026-09-14 19:20:04 KST**. 계산 프로세스 PID 400은 종료됐다.
- 19:32 확인에서 세 경로의 manifest에 기록된 21/37/69개 파일 모두 SHA-256이 일치했다. 새 path-4 결과는 디스크에 보존되어 있으나 아직 Git에 추가하지 않았다.
- `time-refinement.json`과 `time-refinement-manifest.json`은 아직 없다. **최종 수렴 판정은 미완료**다.

## 중단 지점

`verification/gr_cached_material_evolution.py`의 `chain`은 `run(4)`를 마친 뒤 `compare()`에서 중단됐다. 로그의 예외는 다음과 같다.

```text
verification/gr_conservative_two_carrier.py:489, completed
assert cp['time_seconds'] == t
AssertionError
```

저장 시각과 비교기가 재구성한 시각이 정확히 일치하지 않았다. 어느 경로/단계에서 얼마나 다른지, 이력 재사용 시 시간 격자 재구성의 반올림 때문인지는 아직 조사하지 않았다. 원인 확정이나 수렴 실패로 해석하지 않는다.

## 재개 순서

1. 현재 Git, 세 경로 manifest, 계산 프로세스를 다시 확인한다. 기존 원39/78/156 실패 및 수정8/16/32의 열 차수1.35 미달 판정은 보존한다.
2. 첫 시각 불일치의 경로/단계를 찾아 원본 pilot 격자, 복사한 16/32 이력, 새 계획의 격자를 비교한다. 기록된 이력에 실제 사용된 시간 간격과 BDF 계수를 확인한다.
3. 원인을 고친 최소 후처리 경로로 저장 이력을 재검증하고 16/32/64 비교를 완료한다. 정확 일치 검사를 무조건 생략하거나 수렴 문턱을 완화하지 않는다. 소스 해시로 묶인 기존 실행 코드를 수정해 과거 바인딩을 무효화하지 않는다.
4. 비교 결과와 근거를 관련 연구 문서에 반영하고 체크포인트를 만든다. 전체 기간 GR 실행 여부는 이 판정을 근거로 결정한다.

이번 상태는 적분 완료와 판정 미완료를 구분해야 한다. 현재 새 장기 계산이나 자동 재실행은 시작하지 않는다. 보고는 한글로 하고, 이미 승인된 후속 작업의 허락을 반복해서 요청하지 않는다.

## 확인한 manifest SHA-256

| 경로 | SHA-256 |
| --- | --- |
| path-1 (16단계) | `47e397eec5057a4de7cb36924343aad389cc74d8efd4c4616dc88a8611e5b375` |
| path-2 (32단계) | `3467457f8a25dd0a6995e5278a0f1a7f2a6f30c4f1fce1fb29d32209891facbd` |
| path-4 (64단계) | `90ce33ed4f369007b2d3e5378e35141844fde22e84c1798c3fccd627f4d7c143` |


## 재개 후 갱신

후처리 원인을 원 계획과 재분할 격자의 1 ULP 차이로 확인하고 복구했다. 세 경로의 엄격한 재생이 완료됐고 time-refinement.json에 온도 차수 1.10825 미달을 기록했다. 저장 자료 재검증과 후처리 수정은 완료된 작업이다. 다음 병목은 외곽 셀 5733의 온도 시간 수렴이며 세부 상태는 GR_EXPERIMENT_REDESIGN_KO.md를 따른다.


## 추가 128단계 실행

운영 기록: 2026-09-14 19:48:48 KST에 pilot-cached-32-64-128을 시작했다. 실제 부모 PID는 launch.json에 있고, CPU 0–15·15작업자·정확한 native 캐시를 유지한다. 기존 32/64 경로는 파일을 중복 복사하지 않고 원 계획 및 manifest로 참조하며 새 128단계만 원 초기 상태에서 계산한다. 새 시간 격자는 앞선 16구간 base를 그대로 사용해 8분할하므로 그 2단계마다 기존 64단계 시각과 정확히 일치함을 준비 단계에서 확인했다. 실행 소스와 새 계획은 체크포인트 4ab523e7로 고정했다.

분류: Conjectural. 시간 해상도를 한 번 더 높이면 외곽 온도 차수 미달을 해소할 수 있는지 검증한다. 저장 이력에서 이 셀은 진동하며, 다음 세분화의 성공은 아직 확인되지 않았다. 전체 초기 상태·EOS·경계·공간식·정밀도·수락 오차 및 차수 1.5 기준을 유지한다. 완료 뒤 32/64/128의 native 재생과 같은 다섯 변수 문턱을 자동 판정하며, 실패는 보존한다. 전체 기간 본실험이나 추가 세분화를 자동 시작하지 않는다.


## Request 33 예비실험 32/64/128 수렴 통과

분류: Counterexample candidate. 128단계 완료 뒤 32/64/128의 엄격한 native 재생·보존·특성속도 및 시간 문턱을 모두 통과했다. 다섯 끝점 차이가 모두 감소했고 밀도/온도/속도/총 열유속/복사 열유속 차수는 1.97680/1.73409/2.05741/2.03979/2.03979다. 온도와 두 열유속의 기준은 그대로 1.5다. 검증한 기간은 0.0016452000355719064초의 전체 격자 예비실험이며 원 세 실패는 보존한다. 실제 예비실험의 온도 시간 수렴 병목을 넘어선 loophole progress다.

분류: Conjectural. 전체 0.42117120910640804초의 수정 GR 시간 수렴과 연속 오차 상계·물리 EOS·외부·반응·관측 폐쇄는 남는다.

세부 근거: [GR 예비실험 최종 통과](GR_EXPERIMENT_REDESIGN_KO.md).


## 전체 기간 70/140/280 실행 인계

운영 기록: 2026-09-14 22:09:59 KST에 verification/gr_compatible_full_duration.py chain --workers 15를 시작했다. 실행 소스·시간 격자·복원 검사는 3b33ed68159903516e6950b3357ac015d2bfb8a8에 고정되어 있다. 실제 PID·명령은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-cached-prefix/launch.json, 진행은 각 path-1/2/4/progress.json 및 run.log를 확인한다. 기존 32/64/128단계는 재생이며 새로운 계산량에 포함하지 않는다. 각 경로는 38/76/152단계를 추가 계산한다. 같은 명령을 중복 실행하거나 실행 중인 소스·계획을 수정하지 않는다.

분류: Counterexample candidate. 세 연결점 복원과 기존 symbolic 검사를 통과했다. 세 경로의 0.4211712091064080384초 도달 및 최종 판정은 진행 중이며, 마지막 time-refinement.json과 execution-manifest.json이 저장되어야 최종 판정을 확인할 수 있다. 실패 또는 중단이 있으면 원 경로를 보존하고 저장 이력부터 이어 간다. 전체 시간 수렴은 아직 주장하지 않는다.


## 2026-09-15 최신 상태: 전체 기간 실패 후 구간 분할 대조

운영 기록: production-cached-prefix의 첫 경로는 41단계까지만 수락됐고 42단계에서 종료됐다. 140/280 경로는 미시작이다. 현재 실행 대상은 verification/gr_step42_subdivision.py이며 출력은 outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/step42-subdivision이다. 이 명령의 실제 프로세스, iterations.jsonl, progress.json, result.json, failure.json을 확인한다. 다른 step42 디렉터리는 완료된 원인 대조 또는 실패한 반복 수정이므로 활성 적분으로 읽지 않는다.

분류: Counterexample candidate. 저장된 실패 상태에서 이류 donor 전환과 고정 접선의 불일치를 분리했다. 고정 donor 대조의 통과를 native 해의 통과로 사용하지 않는다. gr_upwind_iteration_repair.py의 두 donor 교체 시도는 실패했고 생산 실행에 적용하지 않는다.

운영 기록: subdivision의 두 단계가 실제 통과한 경우에만 verification/gr_subdivided_full_duration.py prepare를 사용할 수 있다. 이 준비 코드는 미실행 상태이며, 99/198/396의 미래 시간 격자와 43/64/128의 수락 이력을 검증하도록 작성했다. 성공 여부를 먼저 확인하고 계획·소스 고정과 복원 검사를 거쳐 chain을 실행한다. 기존 실행이 있으면 중복하지 않는다. subdivision도 실패하면 해당 상태와 원인을 보존하고 실패 구간부터 다음 복구를 이어 간다. 과학 수락 문턱을 완화하지 않는다.

운영 기록: 같은 작업의 heartbeat gr-156-2를 1시간 점검으로 활성화했다. 실행 식별자와 진행률은 현재 파일/프로세스에서 읽고, 진단 기록만 추가하지 말고 실제 결합 적분을 이어 가는 것을 다음 목표로 삼는다.

## 2026-09-15 구간 분할 종료 후 인계

운영 기록: gr_step42_subdivision.py는 두 번째 절반 step43에서 종료됐다. step42-subdivision/failure.json이 정규화 잔차 17.560264824407337을 기록한다. 첫 절반 step-0042.npz만 native·보존·특성속도 문턱을 통과한 수락 상태이며 이력은 99b3da6c에 보존되어 있다. gr_subdivided_full_duration.py prepare/chain은 실행하지 않는다. 두 절반 성공을 요구하는 진입 조건이 실패했기 때문이다.

분류: Counterexample candidate. step42-complementarity와 step42-consistent-donor도 실패했다. 단순 donor 교체 및 그 후보 연산자는 생산 진화에 적용하지 않는다.

운영 기록: 현재 시험은 verification/gr_consistent_donor_offset.py, 출력은 step42-consistent-offset이다. 해당 실제 명령과 result.json/failure.json을 먼저 확인한다. 후보 통과 시에도 원 native 복원·수지·특성속도와 전체 기간 재실행을 구분한다. 1시간 heartbeat는 활성 상태다.
