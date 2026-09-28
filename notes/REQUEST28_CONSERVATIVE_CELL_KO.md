# Request 28 — 반응 에너지 보존과 질량 정규화 scalar 읽기

분류: Proven. **고정 GR 보간 모형의 지정 scalar 미분에 구간 보증을 부여하고, 계량·질량 정규화에 필요한 변분식을 도출했다.** 이 보증은 실제 항성의 EOS·구조 불확실성을 포함하는 보증과 구분한다.

분류: Counterexample candidate. **26종 전체 조성과 온도, 중성미자 에너지를 함께 진행하는 고정 부피 대조를 수행했다.** 기본 EOS의 실패와 별도 HELM 근사의 결과를 각각 보존한다. 일 단위 핵 가열 상태를 바로 자유낙하 질량 변화로 읽는 해석은 사용하지 않는다.

## 1. 먼저 닫은 에너지 장부

분류: Imported from prior work. 조성이 바뀔 때 내부에너지의 조성 미분을 열수지에 포함해야 한다는 구조는 [MESA 에너지식](https://docs.mesastar.org/en/latest/reference/controls.html)에 명시되어 있다. 해당 문서를 오래된 실행 파일의 빌드 인증으로 사용하지 않는다. 이번 수치 경계는 기존 실행 파일과 저장한 입력·직접 반환값이다.

분류: Counterexample candidate. 이전 최대 가열 구역 2592에서 전체 26종 변화율을 사용하면, Q 질량 기준의 정지에너지 방출은 **732516.638203 erg/g/s**, 반응 중성미자 손실은 **48379.7343255 erg/g/s**다. 차이는 native 가열 **684136.903878 erg/g/s**와 약 **4.02×10⁻⁸ erg/g/s** 이내에서 일치한다. 이 일치를 모든 상태의 약반응 Q 보정에 대한 인증으로 확대하지 않는다. [에너지 장부](../outputs/conservative-cell28/energy-ledger.json).

분류: Counterexample candidate. 네 PP 중간 핵종과 He4 보상만 남긴 벡터의 정지에너지 방출은 약 **−4.91×10⁻¹³ erg/g/s**다. 전체의 비영 가열을 공급하는 것은 그 밖의 벌크 연료 변화다. 따라서 이전 네 상태의 준정상 반응 부분계는 닫힌 연료 계가 아니다. 그 가열 전달함수를 그대로 독립적인 총질량 원천으로 넣을 수 없다.

분류: Proven. 이번 셀 대조는 바리온 상수항을 제외한 에너지를 명시적으로 다음과 같이 정의한다.

```text
rQ_i = Qconv · qex_i / A_i
RQ(X) = Σ rQ_i X_i
E = u_native(ρ,T,X) + RQ(X)
dE/dτ = −(εν,reaction + εν,thermal)          [고정 부피, 외부 일·열 없음]
```

분류: Proven. 이 정의에서 `u_Q=u_native`로 두고, 중성 원자 W 기준으로 바꿀 때 `g=RW−RQ`, `u_W=u_Q−g`를 함께 바꾸면 총에너지는 동일하다. 이는 기준 변경의 대수적 일관성이다. native EOS의 실제 조성별 에너지 영점이 물리적으로 교정되었다는 뜻은 아니다. 정지에너지 감소와 내부 가열을 총질량 변화에 중복해서 더하지 않는다.

## 2. 기본 EOS 실패를 그대로 남긴 온도 역산

분류: Counterexample candidate. 같은 초기 밀도·온도·26종 조성에서 고유시간 **1.6294일**의 1·2·4분할 진행을 등록했다. 엄격 에너지 잔차는 초기 방출 에너지 척도의 `10⁻⁶`, 별도 해상도 진단 한계는 `0.05`다. 조성 시간 대조는 성분마다 `10⁻¹⁶+10⁻³|전체 변화량|`, 온도는 `|ΔlnT|<2×10⁻⁶`를 사용했다. [사전 계획](../outputs/conservative-cell28/plan.json).

분류: Counterexample candidate. 기본 EOS 결과는 다음과 같다. 4분할 행은 첫 단계의 중단 결과이며 최종 시각의 해가 아니다.

| 분할 | 에너지 잔차/초기 전체 방출 척도 | 엄격 에너지 | 조성 시간 대조 |
|---|---:|---|---|
| 1 | 0.0238737 | 실패 | 비교 대상 없음 |
| 2 | 0.0287696 | 실패 | 오차/허용량 42.7940, 실패 |
| 4 | 0.0841597, 첫 단계 | 실패·진단 한계 초과 중단 | 완료 해 없음 |

분류: Counterexample candidate. 실제 온도 역산 중 내부에너지가 거의 같은 값에 머물다가 점프했다. 이전 EOS 경계 측정과 일관되지만, 이 유한 반복만으로 모든 온도에 대해 근이 없다는 정리는 주장하지 않는다. 소스의 `eosdt_eval.f90`은 보간 좌표 `logRho`, `logT`, `logQ`를 기본 REAL로 저장한다. profile의 에너지는 `s%energy`를 직접 출력한다. 기록된 실패를 저장 형식만의 문제나 단순 Newton 허용치 문제로 바꾸어 설명하지 않는다. [원래 진행 결과](../outputs/conservative-cell28/continuation.json).

## 3. HELM 근사의 보존형 결합 대조

분류: Counterexample candidate. 기존 `use_eosDT_HELMEOS` 옵션으로 **별도 완전 이온화 근사**를 정의했다. 같은 초기 ρ·T·X를 사용하고 반응망·EOS 출력을 다시 계산했다. 다른 구역은 native 평가를 위한 입력 틀이며 물리적 결과는 고온의 한 구역에만 한정한다. 차가운 외피의 물리 EOS나 기존 GR 평형의 대체 해로 취급하지 않는다. [별도 계획](../outputs/conservative-cell28/HELM-plan.json).

분류: Counterexample candidate. 온도를 단계 동안 고정한 첫 HELM 대조는 에너지 역산을 통과했지만 조성 시간 기준은 실패했다. 1→2 및 2→4분할의 조성 오차/허용량은 **37.8952, 18.3059**다. 에너지 잔차가 작다는 사실만으로 시간 적분을 완료로 판정하지 않았다. [분리 적분의 실패](../outputs/conservative-cell28/HELM-continuation.json).

분류: Counterexample candidate. 다음 대조는 26종 조성·온도·누적 중성미자 에너지를 하나의 행렬 지수 단계로 결합한다. 빠른 선형 반응 과도를 손실 적분에도 포함하고, 단계 끝에서는 실제 HELM 내부에너지로 온도를 역산한다. 기존 에너지·조성·온도 기준을 유지하고 중성미자 손실의 시간 정밀화 `10⁻³` 기준도 추가했다. [결합 계획](../outputs/conservative-cell28/coupled-plan.json).

분류: Counterexample candidate. 초기 HELM 조성 에너지 미분은 H1/He4 및 C12/He4의 독립 방향으로 평가했다. 간격을 절반으로 줄였을 때 최대 정규화 차이는 **1.24×10⁻⁸ 이하**다. 완전 이온화 HELM의 `abar,zbar` 의존성을 이용한 국소 미분이며 전체 항성 EOS의 조성 미분 보증은 아니다. [미분 대조](../outputs/conservative-cell28/HELM-composition-energy.json).

분류: Proven. 알려진 비가역 연료 변환의 해석해로 조성 변화·온도 상승·중성미자 적분과 총에너지 불변량을 함께 대조했다. [양성 대조](../outputs/conservative-cell28/coupled-positive-control.json).

분류: Counterexample candidate. 결합 행렬에는 근사 Jacobian을 사용한다. 열용량의 미분, 열 중성미자의 미분, 일부 상태 의존 약반응 Q 보정 및 초기 이후 조성 에너지 미분의 변화를 생략했다. 원천 벡터와 단계 끝 EOS는 native 재평가한다. 이러한 생략을 엄밀한 전체 Jacobian이나 전역 ODE 오차 보증으로 감추지 않는다. 최종 수치 및 시간 기준은 [결합 결과](../outputs/conservative-cell28/HELM-coupled-continuation.json)에 기록한다.

분류: Counterexample candidate. **결합 계산은 모든 등록 기준을 통과했다.** 1→2 및 2→4분할의 조성 오차/허용량은 **0.122514, 0.0285460**, 온도 ln 차이는 **4.71×10⁻¹¹, 2.35×10⁻¹¹**, 중성미자 손실 상대 차이는 **2.53×10⁻⁵, 1.25×10⁻⁵**다. 전체 단계의 최대 에너지 잔차/초기 방출 척도는 **6.963×10⁻⁷ 미만**으로 `10⁻⁶` 기준 안에 있다. 이 유한 정밀화 대조는 연속 비선형 해의 전역 오차 증명과 구분한다.

분류: Counterexample candidate. 최종 4분할 상태의 `ΔlnT`는 **2.1298809×10⁻⁵**, 최대 조성 변화는 **2.3786809×10⁻⁸**이다. Q 정지에너지는 약 **1.09057008×10¹¹ erg/g 감소**, 내부에너지는 **1.00195332×10¹¹ erg/g 증가**, 중성미자 손실은 **8.86167552×10⁹ erg/g**다. 세 항의 잔차는 **29.662 erg/g**다. 가열만 총질량 증가로 읽을 때와 부호·크기가 다르지만, 이 셀 에너지 차이를 항성 ADM 질량의 실제 변화로 부르지 않는다.

분류: Counterexample candidate. 진행 중 native 가열과 `−rQ·f−νreaction`의 차이는 가열 대비 최대 **4.472×10⁻⁶ 미만**이었다. 초기의 거의 정확한 일치를 전 구간의 정확한 원천 항등식으로 승격하지 않는다. 열 중성미자 원천은 이 고온 셀에서 가열 대비 **6.662×10⁻¹¹ 미만**이었다. 실제 원천의 작은 보정과 근사 Jacobian의 생략은 기록된 모델 경계로 남는다. [최종 장부](../outputs/conservative-cell28/closure-summary.json).

## 4. scalar 읽기에 필요한 정확한 변분식

분류: Imported from prior work. 정적 구대칭 scalar 방정식과 원거리 질량 정규화는 [Damour–Esposito-Farèse의 이론 틀](https://arxiv.org/abs/gr-qc/9602056)을 따른다. 아래 변분식과 수치 대조는 이 저장소에서 별도로 도출했다.

분류: Proven. `φ(∞)=1`인 정적 해에 대해

```text
(p φ′)′ = v φ
p = r² N √f
v = 4πG β r² N (e−3P)/(c⁴√f)
φ = 1 + χ/r + O(r⁻²)

δχ = −∫[φ² δv + (φ′)² δp]dr
```

분류: Proven. 고정 면적 반지름 좌표에서 중심의 정칙성, 원거리 정규화 및 경계항 소멸을 가정한 식이다. 움직이는 불연속 경계가 있으면 해당 분포·경계항을 포함해야 한다. 원거리에서 `δp=O(r)`, `φ′=O(r⁻²)`이면 추가 경계항은 소멸한다. 따라서 고정 계량에서 `δv`만 계산하는 것은 실제 별의 전체 미분이 아니다.

분류: Counterexample candidate. 저장한 GR 보간 배경에서 지정 모양 `s(x)=16x²(1−x)²`를 적용했다. 직접 선형 변분 방정식과 적분식의 대조는 세 방향에서 약 **3.1×10⁻¹⁵ 이하**의 상대 차이를 보였다. `δv=v s` 항은 약 **251.586108 m**, 별도의 `δp=p s` 항은 약 **−0.06314893 m**다. 시험 계량 변형이 Einstein 제약식을 만족하는 실제 별의 변형이라는 주장은 하지 않는다. [계량 변분 대조](../outputs/conservative-cell28/metric-variation.json).

분류: Proven. 물리적 원거리 전개를 `φphysical=φ∞−α_A Mgeom/r+…`로 정의하면

```text
α_A = −φ∞ χ/Mgeom
δα_A = −φ∞[δχ/Mgeom − χ δMgeom/Mgeom²]      [φ∞ 고정]
```

분류: Proven. 결합에너지를 무시하는 보편 결합의 선도 한계에서 전하가 `Q=α M`이면, 내부 에너지 재분배에 대한 `δ(Q/M)`은 0이다. 유한 자기중력·압력·계량 변화는 이 한계에서 제외되어 있으므로 잔여 효과도 자동으로 0이라고 하지 않는다. 닫힌 전체 계의 내부 핵 전환은 독립적인 총에너지 원천이 아니다. 항성 부분계의 에너지 변화에는 경계를 통과하는 복사와 외부 일을 포함해야 한다. 전체 ADM 에너지와 복사 손실을 뺀 부분계/Bondi 질량을 혼동하지 않는다.

## 5. 실제로 확보한 조건부 미분 구간 보증

분류: Proven. 이번 보증의 대상은 저장한 binary64 PCHIP 계수와 매듭을 정확한 수로 해석한 **명시적 수학 모형**이다. 매개변수는 고정 계량에서 `v(ε)=v₀[1+εs(x)]`로 정의한다. 원래 연속 EOS/TOV 해의 불확실성은 이 모형의 정의에 포함되어 있지 않다.

분류: Proven. 정적 Green 적분 방정식의 연산자를 K라 하면, 절댓값 커널의 최대 적분은 중심에서 상계된다. 각 보간 구간을 구간 Horner 연산으로 감싸고, `r²` 원천과 Green 적분의 단항식 부분은 해석적으로 적분했다. Schwarzschild 외부의 Green 항도 포함했다. `||K||≤k<1/2`이면 Neumann 급수로 `|φ−1|≤k/(1−k)`를 얻고, 이를 `−∫v₁φ²dr`에 넣어 미분을 감쌌다. 수치 ODE의 유한 차분 일치를 보증의 전제로 사용하지 않는다.

분류: Proven. [mpmath 구간 산술](https://mpmath.org/doc/current/contexts.html#arbitrary-precision-interval-arithmetic-iv)의 지원되는 기본 연산·제곱근·로그·상수만으로 실제 계산을 수행했다. 40자리 구간 연산의 양 끝은 정확한 유리수로 저장한다. 일정 밀도 구의 해석해를 별도 대조했다. 보간 구간당 1개와 4개 하위 구간에서 얻은 결과는 다음과 같다.

| 하위 구간 수 | 총 구간 수 | 연산자 노름 상계 | 미분 구간, m |
|---:|---:|---:|---:|
| 1 | 5,836 | 0.000167146775 미만 | [251.119201, 251.944154] |
| 4 | 23,344 | 0.000167068151 미만 | **[251.36919, 251.70105]** |

분류: Proven. 표의 구간 끝은 바깥쪽으로 반올림한 표시다. 정확한 유리수 끝점, 계수 SHA와 재현 경계는 [구간 보증](../outputs/conservative-cell28/interval-derivative.json)에 있다. 이 구간은 고정 보간 모형의 미분 오차 보증이며 물리적 항성·timing 미분의 인증은 아니다.

분류: Proven. 별도로 `N,f≥0.99`, 물질 반지름 `≤7×10⁷ m`, `0≤e−3P≤2×10⁸c²`, `β=−4`, `|s|≤1`, `|ε|≤0.01` 영역의 1·2·3차 매개변수 미분 상계도 유리수로 도출했다. 중간 구간 차분의 절단 오차는 간격 `10⁻³`에서 **0.000324423 m 이하**, `5×10⁻⁴`에서 **0.0000811057 m 이하**다. 이것만으로 일반 부동소수점 ODE 오차가 포함되지는 않는다. [해석적 상계](../outputs/conservative-cell28/derivative-bounds.json).

분류: Proven. 두 가까운 단극 값을 빼는 대신 `[(χ(+h)−χ(−h))/(2h)]=−∫φ(+h)φ(−h)v₁dr`를 쓰는 정확한 Wronskian 항등식도 얻었다. [수치 대조](../outputs/conservative-cell28/finite-variation.json)는 이 계산의 안정성을 확인하며, 앞의 구간 보증과 구분한다.

## 6. 완료 범위와 다음 연결

분류: Proven. 이번 theorem progress는 **고정 GR 보간 모형의 scalar 미분 구간 보증**, 계량을 포함한 감수율 변분식, 질량 정규화의 상쇄 경계다. 구간 보증의 좁은 대상은 이제 “전혀 미해결”로 기록하지 않는다.

분류: Counterexample candidate. 이번 loophole progress는 **열욕을 제거한 26종 반응·열·에너지 적분의 실제 대조**다. 기본 EOS, 온도를 분리한 HELM, 온도를 결합한 HELM의 판정을 구분한다. 에너지 장부가 닫혀도 고정 부피 근사가 항성의 유체·복사 수송을 대체하지 않는다.

분류: Conjectural. 물리적 미완료 연결은 생성 핵종과 부분 이온화를 포함하는 공통 EOS의 인증 및 에너지 영점 교정, 보존된 26종 항성 초기화와 열·유체·계량 진화, 비영 scalar 배경의 실제 구동·질량 정규화 전하 읽기, 완전한 관측 전방 모형이다. 이전의 영 scalar 가지 선형 분리와 유한 carrier 보간 경계는 유지한다.

분류: Proven. 다섯 유지 문서와 역사 SHA 연결을 갱신한다. 이전 실패 판정과 Request12 원고 PDF/ZIP은 보존한다. 이번 연구 노트의 검증을 최종 투고 패키지 갱신이나 제출로 부르지 않는다.

분류: Proven. 전체 저장 대조는 **직접 추출 profile 55개, 일반 profile 1개, native 구역 호출 315425개**를 포함한다. 직접 반환 가열·중성미자와 profile의 최대 차이는 **0**이었다. 최종 상태의 별도 비추적 실행도 일치했다. 새 임시 폴더에서 선택 상태 **4개를 실제 재실행**한 최대 native 차이도 **0**이었다. 전체 결합 단계의 저장 출력 재생과 기호·변분·구간 계산 재실행을 통과했다. 모든 55개 상태를 별도 독립 물리 실험으로 세거나 모두 새로 재실행했다고 주장하지 않는다.

```bash
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 \
  python3 verification/conservative_cell.py verify
PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps OPENBLAS_NUM_THREADS=1 \
  python3 verification/conservative_cell.py recheck
```
