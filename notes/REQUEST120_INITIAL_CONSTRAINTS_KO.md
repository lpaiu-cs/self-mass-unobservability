# 단계120 — 실제 셀 재고를 보존하는 초기 GR·스칼라 배경

분류: Counterexample candidate. 단계119의 실제 19개 내부 셀과 512개 대기 셀에 대응하는 초기 제약 해를 만들었다. 새 질량·lapse·스칼라장·밀도와 비영 초기 광자 유속에 필요한 외재곡률을 저장했다. **이 초기 상태를 실제 광자·물질 진화기에 설치한 새 궤적은 아직 없다.** 이전 전하의 조건부 양수 하한을 새 배경에 승계하지 않는다.

분류: Imported from prior work. 방정식 기준은 [Salgado의 구면 물질·복사 RGPS 식233–235](https://arxiv.org/abs/gr-qc/0201064)와 [Novak의 Einstein-frame scalar–tensor 식2.23–2.26](https://arxiv.org/abs/gr-qc/9707041)다. lapse는 slicing 조건이며 물리 관측량 자체가 아니다. 순유속이 있는 초기 slice의 외재곡률을 0으로 놓고 정확한 정적 계량이라고 부를 수 없다.

분류: Proven. `Pi_phi=0`, `b=1−2m/r`, `A=exp(−2phi²)`, `alpha=−4phi`인 선언 모형에서 사용한 식은 다음과 같다. 에너지·압력·유속은 국소 Jordan-frame 성분이며 `m`의 단위는 cm다. 운동량 식의 두 표현과 국소 등엔트로피 구성 관계의 제1법칙은 독립 기호 검사로 확인했다.

```text
m' + r Phi² m = 4 pi r² A⁴ G E / c⁴ + r² Phi² / 2
nu' = m/(r² b) + 4 pi r A⁴ G Pr/(c⁴ b) + r Phi²/2
K^r_r = 4 pi r A⁴ G F/(c⁵ sqrt(b))
dm/dt = −c N r b K^r_r
Q' = 4 pi G N r² A⁴ alpha (E−Pr−2Pt)/(c⁴ sqrt(b))
Phi = Q/(N sqrt(b) r²),  Q(0)=0,  phi(infinity)=0.001
```

분류: Counterexample candidate. 먼저 이전 스칼라 값·미분의 C1 Hermite 재구성을 고정하고 질량·lapse 제약을 풀었다. 해당 Cauchy 상태는 스칼라 가속 잔차가 최대 `4.30059e−8 s^-2`였다. 이는 보간된 초기장의 잔차이며 물리적인 방사 신호의 발견이 아니다. 그 가속도를 0으로 가정한 채 기존 미소 전하 읽기를 재사용하지 않았다.

분류: Counterexample candidate. 다음으로 중심 정칙성·원래 무한원 스칼라 값을 고정하고 질량·lapse·순간 스칼라 균형을 동시에 풀었다. **초기 스칼라 균형은 새로 명시한 Cauchy 조건**이다. 원 고정 스칼라 Cauchy 해를 보존했고, 새 동적 전하나 부호를 계산하기 전에 조건을 선언했다. 관측 전하에 맞춰 초기 펄스를 조정하지 않았다. 광자 순유속과 외재곡률은 남으므로 정적 복사 항성이나 정수압·열평형 해가 아니다.

분류: Counterexample candidate. 첫 매끈한 원천 보간은 자체의 연속 재고를 보존했지만 실제 진화기의 셀 질량과 최대 **15.8312%** 달랐다. 내부 총질량도 실제 `1.1891854923e22 g`에 비해 `1.1928609079e22 g`였다. 이 해는 생산 입력으로 기각했다. 실제 셀의 저장 질량과 광자 에너지·반경 압력을 직접 보존하는 양의 유한체적 원천으로 바꾸어 문제를 해결했다. 원 거부 해와 검사 결과는 별도 보존한다.

분류: Counterexample candidate. 수정 모형은 실제 셀 평균의 구간별 재구성이다. 초기 proper-volume의 중심점 근사와 적분값 사이 차이만 보정하고, 계량·스칼라 보정 중에는 `A³ rho/sqrt(b)`를 보존한다. 광자는 국소 초기 광자 에너지를 고정한 Cauchy 분포로 선언하여 수·에너지·압력의 체적 배율을 함께 적용한다. 가스의 국소 등엔트로피 연장은 실제 EOS의 압력·에너지·체적 탄성률에 맞춘 Gamma 구성 관계이며, 모든 연속 상태에서의 native EOS 인증은 아니다. 알려진 진화 영역보다 깊은 내부의 복사 운동량은 인증하지 않았다. 깊은 내부의 저장 원천 보간도 그대로 전제하므로 전체 원 5,735셀 재고의 재인증과 구분한다.

분류: Counterexample candidate. 아래는 수정한 실제 유한체적 초기 입력의 독립 검사 결과다. 12/20차 구적은 같은 원천의 적분 오차 대조이며 물리적 반경 격자 수렴 대조가 아니다.

| 검사 | 최대 상대 차이 |
|---|---:|
| 원 진화기와 셀별 질량 일치 | `2.9420e−16` |
| 보정 계량에서 적분한 셀별 질량 | `3.4369e−17` |
| 보정 계량에서 적분한 광자 에너지 | `3.4478e−17` |
| 보정 계량에서 적분한 광자 반경 압력 | `3.4694e−17` |
| 질량·lapse·스칼라 등의 12/20 구적 대조 | `5.0073e−16` |
| native 압력·내부에너지 대조 | `2.9024e−7` |

분류: Counterexample candidate. native 엔트로피 역산의 첫 검사는 실패했다. 고정 이온 재고 경로에 평형 열미분 `raw[10]`을 사용한 것이 잘못이었다. 이미 저장된 같은 재고 경로의 `du/dlnT`를 사용하도록 고쳤고, 원 4회 반복과 `2e−12` 엔트로피/Cv 문턱을 유지한 88 native 호출 검사가 통과했다. 원 실패, 20호출 한도에 도달한 진단, 수정 전 코드와 추가 검사 예산을 보존했다. native 점 대조는 전역 미분 오차 상계가 아니다.

분류: Counterexample candidate. 수정 초기값의 밀도 보정은 최대 `|delta ln rho|=8.5114e−8`, 스칼라 보정은 `5.5968e−14`, lapse 중력 기울기의 상대 보정은 `3.2629e−6`다. 선언한 전반경 원천 재구성의 ADM 질량은 `29168.7540761874 cm`, `K=−116.6854204443 cm`다. 이 정적 기준값의 변화는 관측되는 동적 전하 차이가 아니다. 특히 작은 계량 보정이라고 해서 원 `~1e−27` 동적 결과에 대한 영향도 작다고 추정하지 않는다.

분류: Counterexample candidate. 성공한 원 고정 스칼라 제약 해·첫 순간 균형 해·수정 유한체적 균형 해는 각각 약 4.83/4.59/5.27초였다. 수정 native 감사는 약 5.15초였다. 실패 비용도 별도 기록했으며 유체 장기 적분은 0회다. 이 단계는 실제 원천과 초기 제약 사이의 연결 오류를 고친 loophole progress다.

분류: Conjectural. 다음 실행은 `finite-volume/balanced-initial-state.npz`를 새 초기 기준으로 삼아 계량·체적·EOS 기준·Killing 주파수의 광자 표현·중력 힘을 함께 설치하고, 같은 최소 구간의 실제 결합 응답을 계산하는 것이다. 저장된 native 온도·압력·내부에너지 감사점과 기존 EOS 표를 재사용한다. 어느 한 중력 항만 바꾸고 완전한 초기값 이전으로 세지 않는다. 그 뒤 새 실제 원천에서 전하를 다시 판정해야 한다. 반경/주파수·내부 각도, 완전한 미시 반응, 동적 계량과 관측 식별성도 남는다.

실행 코드: `verification/def_native_initial_constraints.py`, `verification/verify_native_initial_constraints.py`.
수락 입력: `outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume/balanced-initial-state.npz`.
수락 근거: 같은 디렉터리의 `balanced-result.json`, `audit.json`, `symbolic.json`.
