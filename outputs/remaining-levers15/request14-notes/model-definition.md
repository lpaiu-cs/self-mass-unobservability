# Model definition — unified revision

The current scientific statement is [the unified manuscript](../paper/manuscript.md), Sections 3–4. Earlier runtime proposals remain in the dated REQUEST10 notes; they are not instructions to reopen runtime work.

Status: Counterexample candidate. The dimensionless prescribed-drive benchmark is

```math
\tau_\chi\dot\chi+\chi=\alpha F(t),\qquad
q(t)=c_YF(t)+c_\chi\chi(t),\qquad \beta=\alpha c_\chi,\quad\tau_\chi>0.
```

Status: Proven. With dimensionless F and chi, alpha, c_Y, c_chi and beta are dimensionless; tau_chi has units of time. For real F0 and `F=F0 cos(omega t)`, the settled response is

```math
\chi_{\rm ss}=\frac{\alpha F_0(\cos\omega t+\omega\tau_\chi\sin\omega t)}{1+\omega^2\tau_\chi^2},\qquad
G(i\omega)=c_Y+\frac{\beta}{1+i\omega\tau_\chi}.
```

Status: Proven. The full solution includes `chi_h exp(-(t-t0)/tau_chi)`. All periodic carrier results assume that term is absent or decayed. The stored templates do not estimate its independent amplitude.

Status: Counterexample candidate. The original mass readout `m_A/m0=1+q` and the timing pair potential `V_pj=-G m_p m_j(1+Delta0+q)/r_pj`, with fixed inertial masses, are separate realizations. The latter specifies the pairwise force modification for a prescribed external drive.

Status: Proven. The prescribed pair potential gives equal-and-opposite pair forces and allows energy exchange through its explicit time dependence. A varied position-dependent drive or inertial mass requires additional gradient/backreaction terms. The ODE alone does not derive these terms or a matching between the two realizations.

Status: Counterexample candidate. [Request 11.3](../notes/REQUEST11_3_MATCHING_RESULT.md) supplies a conditional force-level realization with an independent scalar charge Q_p, potential V(Q_p), inertia I and Rayleigh damping Gamma. Reciprocal pair forces follow from -(m_A m_B+Q_A Q_B)/r_AB; total orbital-plus-state energy loss is -Gamma Qdot_p^2. This is a leading Newtonian reduction, not a complete timing theory.

Status: Proven. On a stable branch kappa=V''(Q0)>0, equal fixed white-dwarf charge/mass ratios a_w give tau=Gamma/kappa and deltaDelta=B deltaU/(1+tau d/dt), B=a_w^2/(kappa m_p). For F=deltaU/Ustar, beta=B Ustar. A small inertial-error bound is needed for the one-pole approximation. Unequal companion ratios obstruct a common pair modulation.

Status: Conjectural. Actual EOS-to-body matching must still determine the coefficients and validate inertia, nonlinear response, feedback, radiation and transients for the real system. The conditional action does not supply those numerical inputs.

Status: Proven. A real local derivative comparator is `P_N(d/dt)F`, with coefficients shared across carriers. A finite-dimensional nuisance model still needs a rank test; finite dimension alone never guarantees pole observability.

Status: Imported from prior work. Explicit dynamical internal modes already appear in compact-body EFT; see Chakrabarti et al. (2013), Steinhoff et al. (2016), and Khalil et al. (2022), cited in the manuscript. The present claim is the specified comparator boundary, not priority for internal states.

Status: Counterexample candidate. Request 11.4 specifies a reciprocal overdamped comparator Gamma qdot+Kq=bF, readout b^Tq+c0F, with symmetric positive-definite matrices and all rates >=Lambda. It yields positive relaxation weights with tau_j<=1/Lambda. This physical comparator premise is not established EOS matching; nonconjugate readout, negative residues or slow/oscillatory modes need a different class.

Status: Proven. Request 11.5 shows that the six periodic response columns plus a static response cannot determine an initial-state exponential response: D*product(D^2+omega_k^2) annihilates the known inputs but not the exponential. A dedicated validated response is needed; an exponential coupling is not automatically an exponential TOA residual.


## Request 12 follow-through

Status: Proven. A positive reciprocal relaxation measure admits a two-frequency moment-variance equality identifying one observable relaxation time, with exact calibrated response; it does not count hidden physical states. Status: Imported from prior work. Request 12 constructs the corrected unequal-amplitude leading drive and dedicated exponential timing responses at 2, 52 and 500 days. Status: Conjectural. EOS coefficients and the fast-rate gap are not matched to J0337.

Details: [remaining-lever report](remaining-levers-2026-09-09.md).

## Request 13 remediation

Status: Counterexample candidate. Request 13 now fixes and numerically solves a concrete SLy/massless-DEF model (beta=-4, zero background scalar, gravitational mass 1.4378144085 solar masses). Its regular stellar scalar response and outgoing pole are calculated; it is not a choice of freely adjustable damping. The response is linear about the zero-scalar branch, where fluid and metric perturbations decouple. The published beta=-5 prototype is a separate positive control. See [stellar matching](../notes/REQUEST13_STELLAR_DERIVATION.md).

## Request 14 validated flow

Status: Proven. A separate interval implementation encloses the frozen internal GR 4-body 1PN IVP and its initial-state variational flow at recorded epochs, with four independent fractional-mass tangent columns in a local augmented run. Its input rejects nonzero dynamic-SEP and non-GR coefficients. This is a conditional numerical certificate for the GR comparator, not a completed matched scalar-tensor signal model. The physical timing-parameter-to-IVP map remains outside the certificate. See [certificate boundary](../notes/REQUEST14_VALIDATED_VARIATIONAL.md).
