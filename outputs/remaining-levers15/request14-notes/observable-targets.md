# Observable targets — unified paper

Current nuisance audit: [Request 11.1 result](../notes/REQUEST11_1_NUISANCE_AUDIT_RESULT.md). The full stored 90-direction space is the primary baseline for subsequent coverage work; the rank-71 result remains a sensitivity comparison. No physical prior justifying hard removal of the other directions has been established.

Status: Imported from prior work. [Request 11.2](../notes/REQUEST11_2_COVERAGE_RESULT.md) tests full nuisance and estimated Fourier covariance, obtaining minimum K=1 U coverage 94.59% within its specified family. This does not establish a universal astrophysical interval.

Status: Proven. [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md) derives a leading coplanar scalar-potential drive with unequal carrier amplitudes and phase closure 3.11837 radians for the stored parameters. The historical auxiliary physical-drive family has closure zero; no common time shift fixes it. Its beta_phys numbers are withdrawn as constraints on this realization. The unit-drive beta benchmark remains defined, but does not cover this physical drive family.

Status: Counterexample candidate. The primary target is a common transfer relation across known nonzero drives, `G(i omega)=c_Y+beta/(1+i omega tau_chi)`, tested after a specified nuisance projection. The instantaneous c_Y direction is co-fitted.

| Status | Proposed signature | Required comparison or limitation |
| --- | --- | --- |
| Proven | One-frequency quadrature | Exactly reproducible by freely fitted F and dot F coefficients. |
| Proven | Multi-frequency pole relation | K carriers exclude only real polynomial degree N below 2K-1, with known deprojection. |
| Proven | Sidebands | Absent in the linear MVP; nonlinear static comparators must also be excluded. |
| Proven | Nuisance-projected residual | Identifiability requires positive residual information, not merely a finite parameter count. |
| Counterexample candidate | Prescribed pairwise SEP readout | Requires physical matching before identification with a field-dependent mass or a universal SEP parameter. |

Status: Proven. The beta coefficient is tied to the normalized drive. Its relaxation-only carrier amplitude is `|beta d_k|/sqrt(1+omega_k^2 tau_chi^2)`. A peak is bounded by the sum of these amplitudes. A bound on total Delta also needs Delta0, c_Y and their covariance.

Status: Imported from prior work. The stored finite-origin-grid Gaussian construction gives, at tau_chi=2 days and assumed width inflation K_dyn=10, `U_beta=1.6795275e-9` for truncated nuisance and `3.5336060e-9` for the full stored nuisance construction. At the five tabulated lags the full/truncated ratio grows to about 17.35. No coverage calibration or drive-independent exclusion follows.

Status: Proven. A maximum over a registered finite origin grid is neither Bayesian phase marginalization nor a supremum over every possible relative phase. The outer orbital period is not a common period of the incommensurate carriers.

Status: Proven. A tidal `I2` correction scales as r^-7 in acceleration for a point source; a constant Nordtvedt correction scales as r^-2. Their coefficients cannot be equated without a theory-specific matching relation.

Status: Proven. The reciprocal fast-response comparator obeys S(omega_l)<=R_Lambda*S(omega_h), where S=-Im(H)/omega and R_Lambda=[1+(omega_h/Lambda)^2]/[1+(omega_l/Lambda)^2]. A positive slower pole violates it. The projected witness and a continuous known-phase information bound are derived in [Request 11.4](../notes/REQUEST11_4_COMPARATOR_RESULT.md). Free independent quadratures instead absorb all six periodic columns.

Status: Imported from prior work. After full nuisance and the declared covariance, allowing real derivatives through order four retains as little as 2.87821e-6 of the instantaneous-only information in the tested cases. Expanded independent-phase local refinement raises the unit-drive U by at most 3.72% in the tested cells; an analytic all-phase upper envelope is also reported, without treating local optimization as a certified supremum. See [Request 11.5](../notes/REQUEST11_5_PHASE_STATE_RESULT.md).


## Request 12 follow-through

Status: Proven. Inverting a joint confidence region in six carrier coefficients supports continuous phase/lag inference within its declared mean model. Status: Imported from prior work. Independent validation gives 95.17--95.43 percent inclusion in four specified covariance conditions. An omnibus carrier excess is not a relaxation detection; each tested physical lag section includes beta=0.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Imported from prior work. The specified SLy candidate has computed inner/outer carrier phase lags about 6.19e-9 and 3.08e-11 radians in the stated scattering convention. These are model response calculations, not measured timing lags. Request 13 also performs full live 28-parameter constrained timing/noise fits for alternative pulse assignments. Their stationary-point, coverage, matched-force and rigorous numerical-error gates are distinct and remain open; no relaxation detection is promoted.

## Request 14 validated flow

Status: Proven. The interval flow supplies explicit initial-state Jacobian error bounds and an outward-rounded width bound for the linear geometric-delay readout at recorded epochs. These do not yet bound a full timing residual, the 28 timing-parameter derivatives or the normalized nuisance projector. No empirical rank, exclusion or lag claim is promoted. See [validated variational certificate](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).
