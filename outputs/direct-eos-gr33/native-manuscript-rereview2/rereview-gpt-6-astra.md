**1. 종합 권고: 주요 수정**

**3판의 차단 지적이었던 관측 감도 배제 결론은 해소됐습니다.** 본문과 §6 추가문은 이제 척도 비교와 관측 배제를 구분합니다. 요청한 새 수치도 산술상 일치합니다.

그러나 **정오표가 두 완화 가정 아래의 작은 진폭을 다시 no-go와 붕괴 회피의 필요조건으로 바꾸고 있습니다.** 이는 본문의 신중한 결론과 충돌합니다. 초록의 가정 누락과 새 판독 정의의 적용 범위도 정리해야 하므로, 현재 상태의 통합은 권고하지 않습니다. 아래 수정에는 새로운 장기 계산이 필요하지 않습니다.

**Imported from prior work.** 선언된 결합 모형의 음의 끝점 판독값은 유지됩니다. 이번 판단은 그 수치 결과를 뒤집는 것이 아니라, 그 결과에서 도출할 수 있는 결론의 범위를 바로잡는 것입니다.

검토 HEAD는 `82ed34311`입니다. 제 초심·3판 재심 원문만 열었으며, 파일 수정·커밋·2분 이상의 계산은 하지 않았습니다.

**2. 이전 지적별 판정표**

행 번호는 [4판 초안](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/white-dwarf-free-fall-charge-section.md)과 [REQUEST284 응답](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/notes/REQUEST284_REVISION4_RESPONSE_KO.md)을 기준으로 합니다.

| 이전 지적 | 판정 | 이유 |
|---|---|---|
| ① §5 계수 구간과 백색왜성 관측 감도의 혼동 | **해소** | 초안 116·136·158행에서 템플릿 불일치, 미계산 타이밍 응답, 배제 불가를 명시했습니다. 기존 차단을 해제합니다. |
| ② 단일 완화 상한에서 A4 유지로의 도약 | **부분 해소** | 본문 114·138·144행은 두 가정과 상태 존재 여부의 미판정을 적절히 명시합니다. 그러나 응답 101–102행에 잘못된 no-go·필요조건이 남았습니다. |
| ③ 조석 감쇠식의 적용 조건 | **해소** | 본문 111행은 약감쇠 근사의 적용 실패와 공명 미배제를 유지합니다. 이전 노트의 해당 주장도 정오표에서 명시적으로 철회해야 합니다. |
| ④ GR·곡률·비선형 추정을 전하 오차 상계로 승격 | **해소** | 본문 90–93행의 제한된 수학적 진술과 물리적 추정을 구분한 수정이 유지됐습니다. |
| ⑤ 대기 보정의 방식·수렴 범위 | **해소** | 고정 핵, 개별 매개변수 변화, 제한된 광학깊이, 가중·최대 오차를 구분합니다. 비회색 1σ 범위도 수정됐습니다. |
| ⑥ 표지와 초록·논의의 결론 강도 | **부분 해소** | 조석 항목 등의 표지는 개선됐습니다. 초록 155행에는 완화 **강도** 가정이 여전히 빠져 있습니다. |
| ⑦ 원천·판독 시각과 정적 창 부호 | **해소** | 본문과 정오표가 시각, 반대 부호, 저장 표본과 최대값 격자를 구분합니다. |
| ⑧ 단위·기호·정규화 | **부분 해소** | 기호 충돌, 무질량 가정, \(T_0\) 등은 해결됐습니다. 새 \(\mathcal Q\)와 뒤의 compact 판독식을 구분해야 합니다. |
| ⑨ 재현 안내 | **부분 해소** | 결과별 manifest 대응은 개선됐지만 필수 입력의 확보·재생성 경로와 실행 연결은 아직 없습니다. |
| 3판 신규: 양의 밀도에 관한 부호 일반화 | **해소** | 보편 명제를 철회하고 보수적인 99.87% 충분조건으로 바꿨습니다. 고정 핵이라는 조건을 해당 문장에도 붙이면 됩니다. |
| 3판 신규: 3.5·4.0 ms 값의 모형과 Born 비교 노름 | **해소** | 전체별 모형으로 이동했고, Born 비교를 서로 다른 노름의 참고 비교로 제한했습니다. |
| 3판 미확인: 자기결합 힘 비 \(10^{-7}\)의 근거 | **해소 — 산술 범위** | 응답 62행의 유도는 재현됩니다. 본문의 `Conjectural` 분류도 적절합니다. |

**3. 새 지적 및 남은 문제**

**[주요] 정오표의 no-go와 붕괴 회피 필요조건은 여전히 성립하지 않습니다.**

위치: [응답 99–102행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/notes/REQUEST284_REVISION4_RESPONSE_KO.md:99), 새로 추가된 [실패 원장 2600행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/failure-ledger-dynamic-chi.md:2600).

**Proven.** \(s>0\)를 계산한 구조 감수율 규모라 두고

\[
H_{\rm struct}(\omega)=\frac{s}{2}
+\frac{s/2}{1+i\omega\tau},
\qquad \tau=\omega_{\rm orb}^{-1}
\]

를 생각하면, \(H_{\rm struct}(0)=s\), 완화 강도는 \(s/2\le s\), 완화시간은 하나입니다. 그런데 궤도 주파수에서

\[
|\operatorname{Im}H_{\rm struct}|=s/4>0
\]

입니다. **두 가정을 모두 만족하면서도 궤도 시간척도의 지연이 존재합니다.** 이는 실제 별의 상태를 발견했다는 주장이 아니라, 제시한 필요조건에 대한 대수적 반례입니다.

따라서 응답 101행의 no-go 표현과 102행의 “강도가 더 크거나 단일 완화가 아니어야 한다”는 필요조건을 삭제해야 합니다. 유지할 수 있는 분류는 **선택한 응답족의 진폭 상한에 관한 조건부 theorem progress**입니다. 관측 배제에는 별도의 타이밍 대응이 필요합니다.

정오표의 충분성도 보완해야 합니다. REQUEST272의 부호 정정과 REQUEST278의 수치 정정은 정확합니다. 반면 REQUEST274 48·60행의 공명 부재·45분 경계, REQUEST276 39행의 중성자별 기원 필요성, REQUEST277·280의 곡률 상한·포괄적 폐쇄 진술은 현 본문보다 강합니다. 역사적 원문은 보존하되, 정오표에서 해당 문장과 현재의 대체 해석을 명시해야 합니다.

본문 138행의 **“neither establishes ... nor shows that none exists”**는 적절합니다. 정오표도 이 결론과 일치시켜야 합니다.

**[경미] 완화 상한의 두 가정과 기준 진폭을 요약에도 정확히 남겨야 합니다.**

위치: [초안 114행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/white-dwarf-free-fall-charge-section.md:114), 126·138·155행.

**Proven.** 단일 Debye 형식만으로 작은 진폭은 보장되지 않습니다. 필요한 명제는

\[
|\Delta\mathcal S|\le|\mathcal S_{\rm struct}|,\qquad
|\operatorname{Im}H_{\rm rel}|
\le \frac{|\Delta\mathcal S|}{2}
\le\frac{|\mathcal S_{\rm struct}|}{2}
\]

입니다.

초록에는 강도 가정을 추가하십시오. 또한 126행의 절반이라는 표현은 **실제 주파수별 완화 응답 진폭의 절반**이 아니라 **정적 pole 강도에 대응하는 진폭의 절반**임을 명시해야 합니다. 전자의 비율은 \(\omega\tau=1\)에서 이미 \(1/\sqrt2\)입니다.

114행에서 강한 과감쇠 응답을 일괄적으로 상한 밖에 놓는 표현도 좁히십시오. 큰 \(\tau\) 자체는 Debye 상한을 깨지 않습니다. 상한 밖인 것은 선언한 강도·단일 pole 조건을 만족하지 않는 응답입니다.

**[경미] 새 질량 정규화 판독과 껍질 적분식이 같은 기호를 공유합니다.**

위치: [초안 30행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/white-dwarf-free-fall-charge-section.md:30), 63·69–82행.

**Imported from prior work.** 새 정의에는 질량 항이 있지만, [전체별 코드의 판독](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-closure-transit/scripts/phase278-transit.py:99)은 껍질 적분만 계산합니다. 질량 정규화된 결합 끝점과의 대조는 별도로 수행합니다.

**Proven.** 다음 구분이 필요합니다.

\[
\mathcal Q_c=-\Psi_{\rm out}/m_{\rm wd},\qquad
\mathcal Q=\mathcal Q_c-\alpha_0\,\delta m_{\rm wd}/m_{\rm wd}.
\]

식 63을 \(\mathcal Q_c\)로 표기하거나, 그 축약 모형에서는 \(\delta m_{\rm wd}=0\)을 취한다는 근사를 명시하십시오. 질량 보정이 작은 결합 끝점 결과를 전체별 이력의 모든 시각에 자동 적용해서는 안 됩니다.

또한 정확한 비율 미분의 질량 항 계수는 배경 천체 전하 \(a_i^{(0)}\)입니다. \(\alpha_0\) 사용은 약한 장 근사임을 밝혀야 합니다. 저장값에서 이 차이가 끝점 판독에 미치는 상대 크기는 약 \(2.0\times10^{-9}\)이므로, 제시한 두 자리 수치나 부호를 바꾸는 문제는 아닙니다.

**[경미] 부호 충분조건의 고정 핵 조건을 바로 붙이십시오.**

위치: 초안 99행, 응답 70–74행.

**Proven.** 99.87% 충분조건은 면별 응답계수와 전파기를 고정한 선형 합에 대한 것입니다. 밀도 변화에 따른 계량·배경 스칼라·전파기 변화까지 포함하는 보장은 아닙니다.

본문에 고정 핵 조건과 \(\|\delta\rho/\rho\|_\infty\)를 명시하십시오. REQUEST272 8행의 \(\delta\ln\rho\) 표기도 유한 상대변화와 구별해 정정하는 편이 정확합니다.

**[경미] 결과 목록을 실제 재현 경로까지 연결해야 합니다.**

위치: [초안 161–172행](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/docs/white-dwarf-free-fall-charge-section.md:161).

목록은 결과를 찾는 데 도움이 되지만, 실제 manifest 파일명에는 `-manifest.json`이 붙습니다. 또한 `phase278-transit.py`는 `readout268-quad64-work/gr/field-source-64.npz`를 요구하고, 실행 셸은 외부 WSL 작업 경로를 전제합니다.

주요 결과별로 **필수 입력 → 실행 스크립트·명령 → 결과 JSON**, 그리고 미공개 입력의 확보·재생성 방법 또는 독립 재실행이 불가능한 범위를 짧게 연결하십시오.

**4. 직접 확인한 새 수치**

아래 산술의 분류는 **Proven — 저장된 입력과 선언된 식에 대한 대입 확인**입니다. 실제 별의 물리적 상계나 검출 감도를 인증한다는 뜻은 아닙니다.

주요 입력은 [tides.json](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-tidal-photosphere/phase274/tides.json), [kernel.json](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-structure-eft-boundary/phase272/kernel.json), [bounds.json](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-closure-transit/phase277/bounds.json)입니다.

| 항목 | 직접 확인값·판정 |
|---|---|
| 순간 응답 | 저장 구동계수에 4를 곱하면 \(1.230431383\times10^{-9}\), \(8.124252845\times10^{-10}\). 본문의 \(1.24\times10^{-9}\), \(8.2\times10^{-10}\)은 선도차수 식에 대한 보수적 상향 반올림입니다. |
| 부호 여유 | \(m=0.998790765535947\), 즉 **99.87907655%**. 99.87% 충분조건은 유효합니다. |
| 양의 밀도 반례 | \(C_+=1.414874740\times10^{-54}\), \(C_-=-2.338701591\times10^{-51}\). 정오표의 조합은 \(+4.896330145\times10^{-55}\)입니다. |
| 구조 계수 | \(8.836576726\times10^{-9}\); \(|\beta_s|\) 대비 \(2.209144181\times10^{-9}\). 새 상한 8.84와 2.21은 적절합니다. |
| Cassini 재척도 | \(|\alpha_0|=0.003535556003\), \(|\varphi_\infty|=0.000883889001\), 제곱 배율 \(0.7812597657\), 재척도 구조 계수 \(6.903661863\times10^{-9}\). |
| 두 쌍 채널 | **\(2.128578030\times10^{-18}\)**, **\(7.525706834\times10^{-21}\)**. 2.13·7.53 표기와 일치합니다. |
| 2차 조석 | \(|dq/d\epsilon|\,\epsilon_T^2(6e_{\rm in})=3.713011755\times10^{-19}\). 제시된 균일 유효중력 대용 모형의 산술과 일치합니다. |
| 자기결합 힘 비 | \(9.915351634\times10^{-8}\); 배경 기울기 항까지 합하면 약 **\(1.33\times10^{-7}\)**입니다. |
| 질량 정규화 끝점 | compact 값 \(-2.334965720\times10^{-51}\)에 \(1+2.241582006\times10^{-5}\)를 곱해 **\(-2.335018060\times10^{-51}\)**입니다. |

**Imported from prior work.** Cassini의 채택 측정값은 원 논문 초록과도 일치합니다. 2σ 하단을 사용하는 것은 원고가 선언한 변환 절차입니다. [Bertotti·Iess·Tortora 원문](https://www.nature.com/articles/nature01997)

**Proven.** [transit.json](E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-closure-transit/phase278/transit.json)의 저장 표본을 선형 보간하면 다음과 같아 REQUEST278 정오표와 일치합니다.

| 판독 시각 | 보간값 |
|---|---:|
| 0.02 s | \(-7.037126\times10^{-47}\) |
| 0.10 s | \(-2.971978\times10^{-43}\) |
| 0.23 s | \(-8.037490\times10^{-41}\) |
| 0.2305 s | \(-8.185481\times10^{-41}\) |
| 0.35 s | \(-2.653797\times10^{-39}\) |

**Imported from prior work.** 첫 저장 부호 변화는 0.4262–0.4282 s입니다. 저장된 최대값은 0.4626 s의 \(3.489397592\times10^{-36}\)이며, 생산 코드가 0.2 ms 판독 격자에서 최대값을 추출하는 것도 확인했습니다.

순간 응답을 별도로 넣은 것은 적절합니다. §4.3의 feedback-stiffness 조건을 적용 대상으로 지목한 것도 타당하지만, **그 조건을 충족했다고 검증한 것은 아닙니다.**

현행 초안·응답·REQUEST281 및 관련 manifest 10개의 revision 해시가 일치했습니다. 별도로 선택한 산출물·스크립트 118개의 해시도 일치했습니다. 파일을 쓰지 않는 Pandoc 변환에서 목록 10개, 항목 46개, 번호식 4개와 식 참조 연결을 확인했습니다.

**5. 확인하지 못한 것**

- 결합 진화·전체별 이력·대기 계산의 독립 재실행과 전체 시간구간의 수렴.
- 실제 별의 비단열 안정성, 중성 모드 부재, 열 완화 강도·시간 및 공명.
- §4.3의 \(|\sum_j C_j/r_{pj}^{2}|\ll\kappa\) 조건의 수치적 충족.
- 두 비공통 쌍 채널의 타이밍 응답·nuisance 투영·검출 감도.
- 같은 위치·시각에서의 Born 성분과 물질 응답 비교, 전체별 이력의 질량 정규화 보정.
- 참고문헌 복원 후의 통합 TeX/PDF 및 공개 입력을 이용한 독립 재현성. 현재 참고문헌 세 항목은 여전히 통합 시 복원이 필요합니다.

