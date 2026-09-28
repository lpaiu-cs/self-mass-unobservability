# 단계105 — 고정 화학 EOS의 저온 실패 해결

분류: Counterexample candidate. 단계104를 중단시킨 약786K 상태의 native 계산 실패를 해결했다. 원 자유에너지·분배함수·이온 재고·밀도 구간을 유지하고, 부동소수점 연산 순서만 수정했다. 저장된34개 상태를 재사용해 원래 지정한49개 단열 상태를244.919K까지 완결했다. 실제 유한 반응·광자 열교환·유출을 함께 푼 결과나 최종 전하의 판정은 아직 아니다.

## 실제 실패 원인과 수정

분류: Counterexample candidate. 원 실패 입력을 별도 진단 빌드에서 재현했다. 첫16×16 Jacobian의 열 하나가 무한대였고, 그 열은 여기항의 전하를 가진 교란입자 합에 대한 미분이었다. native 구현이 사용하는 exp(600) 스케일 아래에서 분자 항의 `mu*qstarx/qstar`가 곱셈 단계에서 overflow했다. 실제 피연산자는 각각1.69625e243,−4.55013e65,1.22442e88이었다. 최종 값은 유한하지만 중간 곱은 배정밀도 범위를 넘는다.

분류: Proven. Q≠0일 때 `mu*x/Q = mu*(x/Q)`이며 `R*(D−A*B/Q)/Q = R*(D/Q−(A/Q)*(B/Q))`다. 이 항등식을 symbolic 검사했다. 이는 동일 수식의 평가 순서 변경이며, 모든 입력에 대한 부동소수점 오차 상계는 아니다.

분류: Counterexample candidate. 첫 수정만으로 밀도 역산의 native info는0이 되었지만 열역학 미분에는 NaN이 남았다. 공통 여기항의 혼합 이차미분70곳에 같은 중간 곱셈 문제가 있어 동일한 대수 변환을 적용했다. 이후 정확히 같은786.379K 입력의21개 출력이 모두 유한해졌다. 실제 실패 피연산자의80자리 Decimal 대조에서 수정식의 상대오차는7.57e−17 미만이었다. 원 Jacobian, 첫 불완전 수정, 실패 출력과 각각의 패치는 보존했다.

분류: Counterexample candidate. 진단 출력에 남는 H2/H2+ 저온 외삽 경고는 기존 Irwin 근사 계산에서 발생한다. 현재 native 소스는 그 뒤 결과를 명시적 분자 스펙트럼 합으로 교체한다. 따라서 그 경고만으로 이번 결과를 구 근사식의 외삽값이라고 해석하지 않았다. 이 사실이 분자 모형 전체의 물리적 인증을 뜻하지는 않는다.

## 원 구간의 완결과 검증

분류: Counterexample candidate. 초기 무친화도 상태와 기존 고정 재고 경로의5개 대표 상태에서 수정 전후21개 출력이 비트 단위로 같았다. 고정 재고 조건은 초기 원소별1e−18 이상 점유한 이온 좌표에 적용하며, 나머지 극미량 좌표와 내부 준위는 조건부 평형이다. 원래의 전체 재고 오차 검사는 유지했다.

| 분류 | 항목 | 결과 | 기존 기준 |
|---|---|---:|---:|
| Counterexample candidate | 완결 밀도점 | 49/49 | ln(rho/rho0):0부터−6 |
| Counterexample candidate | 최소 온도 | 244.919K | native100K 이상 |
| Counterexample candidate | 최대 엔트로피 근 잔차 | 1.812e−10 | 2e−10 |
| Counterexample candidate | 최대 재고 오차 | 8.171e−13 | 1e−12 |
| Counterexample candidate | 거친/미세 rapidity 상대 차이 | 2.162e−6 | 0.002 |
| Counterexample candidate | 새 저온3점의 제1법칙 최대 잔차 | 4.884e−7 | 1e−4 |
| Counterexample candidate | 두 차분 간격의 미분 차이 | 6.983e−9 | 1e−4 |

분류: Counterexample candidate. 제1법칙은 ln(rho/rho0)=−4.5,−5,−6에서 별도로 고정 재고 근을 풀고 두 차분 간격으로 확인했다. 열용량은 양수였으며 고정 화학 Gamma1은1.66670–1.66673이었다. EOS가 외부 친화도를 고정한 채 반환하는 Gamma를 그대로 고정 화학 Gamma로 쓰지 않았다. 이 표본 검사는 전 영역의 엄밀한 미분 오차 인증이 아니다.

| 분류 | ln(rho/rho0) | 고정 화학 T | LTE T | 고정 화학/LTE 압력 |
|---|---:|---:|---:|---:|
| Counterexample candidate | −2 | 3524.073K | 7229.015K | 0.51049 |
| Counterexample candidate | −4 | 928.992K | 6292.784K | 0.16604 |
| Counterexample candidate | −6 | 244.919K | 5695.252K | 0.05167 |

분류: Counterexample candidate. 같은 초기 엔트로피에서도 화학 재고의 조건에 따라 압력·온도가 크게 달라진다. 따라서 기존 LTE 유출의 전하 성분을 실제 화학을 닫지 않은 채 최종 물리적 전하로 해석할 수 없다는 단계104의 경계는 유지한다. 이 대조만으로 실제 기체가 완전히 동결됐거나 전하의 부호·크기가 바뀐다고 결론내리지 않는다.

## 계산 재사용과 남은 결정적 작업

분류: Counterexample candidate. 새15개 상태의 계산은149회 EOS 평가·2.400초, 정상 상태 대조는6회·0.681초, 독립 저온 감사는88회·1.883초였다. EOS 평가는 통상 별도의 결정론적 seed 호출을 포함하므로 이 수치를 내부 Fortran 호출 횟수와 혼동하지 않는다. 새 장기 항성 적분은 실행하지 않았다.

분류: Counterexample candidate. 49개 상태 저장 뒤 원 함수의 `abs(list)` 보고 오류가 발생했다. 계산 산출물은 온전했으므로 이를 재적분하지 않고 저장된 배열에서 특성속도와 모든 판정량을 다시 구성했다. 원 오류와 실행 소스는 보존했고 최종 코드에서는 배열 절댓값으로 수정했다.

분류: Conjectural. 다음 작업은 이 구성 관계를 실제 보존 유출의 종별 재고·에너지 갱신에 적용하는 것이다. 고정 화학 경로는 유한 반응의 비교 기준이며, 자발 복사 재결합만의 예산을 모든 중성화 경로의 상계로 간주하지 않는다. 원자 모형과 일치하는 반응·광자 에너지 교환, 안쪽 물질 응답·희박 꼬리·질량 정규화가 닫히기 전 최종 전하 완료로 표시하지 않는다.

분류: Counterexample candidate. 이번 단계는 실제 저온 구성 관계 병목을 해소한 loophole progress다. 분류: Proven. overflow를 피하는 재배열의 대수적 동등성은 theorem progress다.

근거: [단열 결과](../outputs/direct-eos-gr33/def-native-cold-population/adiabat.json), [독립 감사](../outputs/direct-eos-gr33/def-native-cold-population/audit.json), [수정 코드](../verification/def_native_cold_population.py), [감사 코드](../verification/verify_native_cold_population.py).
