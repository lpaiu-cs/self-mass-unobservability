# Request 33: EOS 값·미분 오차의 조건부 전달

분류: Proven. 이번 결과는 내부 평형을 제거한 자유에너지의 오차를 압력·엔트로피·열용량의 오차로 전달하는 조건부 정리다. 실제 FreeEOS의 물리적 오차 상계나 전체 GR 진화의 오차 인증을 이미 얻었다는 뜻은 아니다. 기호식과 정확한 유리수 대조는 `verification/eos_error_transfer.py`, 결과는 `outputs/direct-eos-gr33/eos-error-transfer/result.json`에 있다.

## 값 오차만으로는 부족하다

분류: Proven. `F0(theta,z)=mu*z^2/2`와 `F1=F0+epsilon*sin(omega*theta)`를 비교한다. 두 모형의 내부 평형은 모두 `z=0`이고 내부 Hessian은 동일한 양수 `mu`다. 최소 자유에너지의 차이는 모든 theta에서 epsilon 이하이나, 그 일차·이차 미분 차이의 최대값은 각각 `epsilon*omega`, `epsilon*omega^2`다. 따라서 기존 동위원소 병진 항의 값 상계만으로 압력이나 열역학 미분을 보증할 수 없다.

## 내부 평형에서 줄어든 자유에너지까지

분류: Proven. 보존 제약을 제거한 고정 내부 좌표 z를 사용한다. 같은 볼록 허용 집합과 공통 이웃에서 참 자유에너지 `F=F0+deltaF`가 충분히 매끄럽고, 내부 평형 `z*`가 존재하며, `F_zz >= mu I`, `mu>0`라고 가정한다. 수치점 `zhat`의 모형 잔차를 `||F0_z||<=r`, 물리적 모델 차이를 `|deltaF|<=epsilon0`, `||deltaF_z||<=epsilonz`로 제한한다. 경계 활성 집합이 바뀌는 경우는 이 정리의 적용 대상이 아니다.

분류: Proven. 강한 단조성과 강한 볼록성으로 다음을 얻는다.

```text
d = ||z* - zhat|| <= (r + epsilonz)/mu
|min_z F - F0(theta,zhat)| <= epsilon0 + (r + epsilonz)^2/(2 mu)
```

분류: Proven. 추가로 `||deltaF_theta||<=epsilon_theta`, 참 모형의 `||F_theta_z||<=L_theta_z`이면, 포락선 정리와 평균값 정리에 의해

```text
||gradient_theta(min_z F) - F0_theta(theta,zhat)||
    <= epsilon_theta + L_theta_z*d.
```

## 이차 미분

분류: Proven. 참 평형의 블록을 `A=F_theta_theta`, `B=F_theta_z`, `C=F_zz`라 하면 암시적 미분으로 `H=A-B*C^{-1}*B^T`를 얻는다. 수치점에서 구한 모형 블록 `A0,B0,C0`와의 차이를 각각 `deltaA,deltaB,deltaC`로 제한하고, `||B0||<=b0`, `C0>=mu0 I`, `mu0>0`라고 하자. 각 블록 차이에는 물리적 이차 미분 오차와 위치 오차 d에 해당하는 블록의 Lipschitz 상계를 모두 포함해야 한다.

분류: Proven. `H0=A0-B0*C0^{-1}*B0^T`에 대해

```text
||H-H0|| <= deltaA
           + deltaB*(2*b0+deltaB)/mu
           + b0^2*deltaC/(mu*mu0).
```

분류: Proven. 증명은 `B0*C^{-1}*B0^T`를 더하고 빼며, 역행렬 차이식 `C^{-1}-C0^{-1}=C^{-1}(C0-C)C0^{-1}`에 연산자 노름을 적용한다. `H0`는 수치점에서 평가한 Schur 식이며, 미수렴한 수치 최적점 곡선을 미분한 결과라고 간주하지 않는다. 정확한 이차 모형의 독립 최소화 대조에서 실제 Hessian 차이가 이 상계 안에 들어감을 유리수로 검산했다.

## 열역학량으로 변환

분류: Proven. 조성을 고정하고 `r=ln rho_B`, `t=ln T`, 바리온 1그램당 Helmholtz 자유에너지 f를 쓰면 `P=rho_B*f_r`, `s=-f_t/T`, `u=f-f_t`다. f의 값 오차 상계를 E0, 일차 성분 오차를 E1, 이차 성분 오차를 E2라 할 때 다음 상계가 성립한다. 여기의 r은 앞 절의 잔차 기호와 구분한 열역학 좌표다.

| 양 | 절대오차 상계 |
|---|---|
| P | `rho_B*E1_r` |
| `dP/dlnrho_B` | `rho_B*(E1_r+E2_rr)` |
| `dP/dlnT` | `rho_B*E2_rt` |
| s | `E1_t/T` |
| `cv*T` | `E1_t+E2_tt` |
| u | `E0+E1_t` |

분류: Conjectural. 실제 EOS에 적용하려면 전체 평형 분율·제약 잔차, 투영 Hessian의 전구간 양의 하한, 물리적 모델 차이의 일차·이차 미분, 블록 Lipschitz 상계, 부동소수점 평가 오차 및 경계·상전이 구간을 확보해야 한다. 현재 공개 래퍼의 21개 출력과 유한 두 간격 차분은 이 조건들을 모두 공급하지 않는다. 이후 EOS 인증에서 확인할 조건을 명시적으로 좁힌 정리 진전이다.
