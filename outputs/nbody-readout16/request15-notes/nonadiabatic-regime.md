# Nonadiabatic regime — exact scope

Status: Imported from prior work. Request 11.2 validates approximately 95% interval coverage within a prespecified full-nuisance, estimated-Fourier-covariance linear model (minimum 94.59% over tested cells). A diagonal covariance can under-cover severely even when the pole survives projection. See [the result and limits](../notes/REQUEST11_2_COVERAGE_RESULT.md).

Status: Proven. A passive first-order pole has no resonance peak. Its tau=Gamma/kappa need not be an inverse particle mass. Close to a scalar-charge instability, increasing susceptibility alone does not ensure negligible inertia, companion feedback or nonlinear response. The explicit physical validity inequalities are in [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md).

Status: Proven. For nonzero real beta and positive tau_chi, the settled response `c_Y+beta/(1+i omega tau_chi)` has a nonzero pole residue. No finite polynomial equals it on an open frequency interval.

Status: Proven. Finite sampling is different. A real degree-N polynomial with unrestricted shared coefficients matches K distinct positive carriers exactly if and only if `N>=2K-1`. Necessity follows from the 2K conjugate roots of `(1+tau z)(P-c_Y)-beta`; sufficiency follows from the explicit polynomial division in Theorem 3 of [the manuscript](../paper/manuscript.md).

Status: Proven. Three distinct positive carriers exclude real degree at most four; degree five suffices to match them. For unrestricted complex coefficients on positive-frequency samples alone, degree `K-1` suffices. Repeated or zero-amplitude carriers do not add independent information.

Status: Proven. Exact noninterpolation does not give a positive lower bound on the measurable residual. Nearby frequencies, small lag, and nuisance projection can make that residual small or zero. For whitened data, the required rank increment is `T_tilde^T (I-P_J) T_tilde>0`.

Status: Proven. An arbitrary independent complex projection at each carrier absorbs the response pointwise. A finite shared projection can also absorb it if its tangent directions span the signal. Rank five instead of six for a particular projection ansatz is not sufficient by itself to identify beta or tau_chi.

Status: Proven. Linear two-frequency forcing and linear readout produce only the input frequencies. Quadratic drive/readout terms can generate sum and difference sidebands, but local nonlinear static responses also generate sidebands. A sideband argument must specify and exclude that comparator.

Status: Counterexample candidate. The observable target is a shared pole relation against a bounded comparator with a justified drive and projection. The three-carrier timing benchmark is a conditional application, not a detection or a claim excluding every finite-order EFT.

Status: Imported from prior work. The stored J0337 finite-grid statistics report no detection; nuisance choices materially change conditional beta intervals. This revision reuses those artifacts without new runtime work.

Status: Proven. For known phases, the coefficient-space residual against a shared derivative comparator is invariant under block phase rotations. Multiplying it by the squared smallest singular value of the full-nuisance whitened carrier map gives a continuous all-phase information lower bound. Allowing the comparator to fit arbitrary independent quadratures is different and restores exact absorption.

Status: Imported from prior work. The item-4 map condition numbers are 98.19–124.82 under the specified covariance scenarios. Nonzero rank can coexist with a roughly 589-fold standard-error increase when fourth-order derivatives are allowed. The mathematical distinction alone therefore does not establish useful empirical precision.


## Request 12 follow-through

Status: Imported from prior work. Dedicated transient response columns pass registered amplitude-convergence gates and add at most 0.523 percent to beta standard errors in the nine tested cases. Status: Proven. Two arbitrarily close positive poles approach one pole quadratically, preventing uniform finite-precision state-count separation without additional separation/weight assumptions.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Counterexample candidate. Request 13 solves the outgoing scalar frequency response and one complex pole for a specified SLy star; changing the exterior contour and radius changes that pole by about 1.5e-8 relative. The low-frequency elastic scattering normalization and a constant-coefficient oscillator fitted to the pole are different reductions. No exact one-pole model over all frequencies or internal state count is inferred from this computation.

## Request 14 validated flow

Status: Proven. Direct interval variational integration removes finite-difference truncation from the frozen GR initial-state Jacobian calculation on its certified time domain. This does not certify dynamic-SEP response columns, physical state count or a complete nonadiabatic observable. The entire timing map still needs its own error propagation. See [certificate boundary](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).


## Request 15 후속 검증

분류: Proven. 다중 정밀도 전파, 초기화 수정·재매핑, 질량·Kepler·대수 관측식의 미분 검증을 수행했다. 이번 계산은 동적 SEP가 없는 GR 기준선에 조건부이며, 전체 기간 타이밍 인증과 비단열 검출 주장을 완성하지 않는다.

세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).
