# Request 10.7b/c Joint-Fit Real-Data Result: Decision Summary

Status: Imported from prior work. Contracts: the pre-registered experiment
`../notes/REQUEST10_EXTERNAL_JOINT_FIT_UPPER_LIMIT.md` (10.7b, commit
`3d9337c`, committed before the first real-data fit) and its amendment
`../notes/REQUEST10_EXTERNAL_JOINT_FIT_AMENDMENT_10_7C.md` (10.7c, commit
`92667cf`, committed before the truncated real-data fit and before either WSL
gate). Labels used exactly as required: `runtime-motivated`, `conditional`,
`collapse`.

## Verdict

Status: Counterexample candidate. **No detection; the program's first quoted
real-data upper limit, via the pre-registered conservative branch.** On the
real `J0337+1715` baseline residuals (12474 TOAs, Nancay 2013-2021), fitting
the shared-relaxation transfer `G(z) = c_Y + beta/(1 + tau_chi z)` jointly
with all 28 standard timing parameters plus an offset:

```text
No detection:  Z = 0.278 (full-rank) / 0.254 (truncated), global p = 0.93 / 0.82,
               off-carrier controls quiet in both pipelines.
QUOTED (conservative, full-rank estimator, valid-but-inefficient):
  u95(10 d) = 11.7 us   u95(26 d) = 28.6 us   u95(52 d) = 53.8 us   u95(104 d) = 102.2 us
CANDIDATE (truncated estimator; absorption-validated, held back only by Gate L):
  u95(10 d) = 1.78 us   u95(26 d) = 4.61 us   u95(52 d) = 9.18 us   u95(104 d) = 18.4 us
```

`beta`, `u95` in us of common pre-transfer drive amplitude at unit drive
(`Lambda_k F_k = 1`); physical drives rescale per 10.7a caveat (c).

## What happened (audit trail)

Status: Proven. (1) The 10.7b full-rank pipeline ran as pre-registered: no
detection (Z = 0.245, p = 0.95), all statistical checks passed
(chi2_red = 1.246, noise scale s = 1.116 < 1.75; controls quiet;
injection-recovery coverage <= 0.08 sigma at all anchors; comparator per C2
reported, not excluded, as the 10.7a counting boundary anticipated).

Status: Proven. (2) The 10.7b linearization stop rule then FIRED: at the
full-rank fitted nuisance displacement the real Nutimo response deviated from
the linear prediction by 5.9x in W-norm (`jointfit_linearization_check.json`).
Per pre-registration, no 10.7b limit was quoted.

Status: Proven. (3) Diagnosis (`jointfit_lindiag.json`): the 21 non-planet
directions are linear to 0.3% (delta_i to 5e-5); the failure is entirely the
planet block, whose 10.7a Jacobian steps (0.3 sigma_MCMC) were nonperturbative
(tasc_extra1 step = 2.8 planet periods, oman = 1 rad, acosi = 148% of value):
secants, not derivatives, with amplitude-oscillatory response.

Status: Proven. (4) The 7 planet columns were re-derived at perturbative steps
with per-column half-step linearity 0.2-2.5% (`jac_v2/rederive_report.json`,
`finite_jacobian_v2.npy`). With the corrected Jacobian: the 10.7a Stage-2
gates are UNCHANGED (principal cosines match v1 to 5e-4; rank still <= 5/6),
and the full-rank pipeline numbers shift by < 1%.

Status: Proven. (5) The full-rank beta marginal is dominated by three
near-null nuisance directions (rel SV <= 5e-4) realizable only through
period-wrapping planet excursions: removing them (SVD cut 1e-3, rank 26/29)
preserves the fitted residual structure to 1.7% but shrinks `sigma_beta(52 d)`
from 27.4 to 4.59 us. The 10.7c amendment pre-registered a truncated pipeline
plus two live-model gates before any truncated real-data statistic was
computed.

Status: Proven. (6) Gate outcomes (`jointfit_gateLA.json`):

- **Gate A (absorption reality) PASS, exactly**: a `4 sigma` chi injection
  (18.37 us) attacked with live Gauss-Newton refitting is NOT absorbed by the
  real model — surviving `z_true = 4.0001` vs truncated-linear prediction 4.0
  (estimator error 0.003%) vs full-rank prediction 0.67. The full-rank
  "absorption" is a linear fiction; the truncated estimator's calibration is
  validated end-to-end at the relevant amplitude.
- **Gate L (displacement linearity) FAIL, marginal**: relative deviation
  0.100 vs the pre-registered 0.05 at the truncated noise-fitted displacement
  (vs 5.9 at the full-rank one).

Status: Counterexample candidate. (7) Quote rule applied (fixed in advance in
the amendment): Gate L FAIL -> quote the full-rank curve as a conservative
limit; the truncated curve is NOT quoted. Note the full-rank estimator is
unbiased with correctly-propagated noise regardless of whether the near-null
directions are physical, so its u95 is valid — merely inefficient by the
factor Gate A measured (~6x).

## Request 10.7d (flagged, small and sharp)

Status: Counterexample candidate. The only obstacle to quoting the ~6x
tighter truncated curve is Gate L's 10% vs 5% at the noise-fitted
displacement — a quantity Gate A suggests is a proxy looser than needed,
since the end-to-end estimator calibration already passed at 0.003%. 10.7d
should pre-register: estimator-calibration gates at the anchor taus
(injections at u95_trunc through live-GN, plus a null run), replacing
displacement linearity as the promotion criterion, with thresholds fixed in
advance. If passed, the truncated curve is quoted.

## Files

Status: Note. Full-rank v1: `joint_fit_upper_limit.json`, `beta_limit_curve.tsv`.
Corrected Jacobian: `finite_jacobian_v2.npy` (+`_meta.json`), `jac_v2/`,
`carrier_projection_rank_v2.json`. Full-rank v2 (QUOTED):
`joint_fit_upper_limit_v2.json`, `beta_limit_curve_v2.tsv`. Truncated
(CANDIDATE): `joint_fit_upper_limit_v2trunc.json`, `beta_limit_curve_v2trunc.tsv`.
Validation: `jointfit_linearization_check.json`, `jointfit_lindiag.json`,
`jointfit_gateLA.json`, `gateLA_params.json`, `jointfit_dtheta52*.npy`,
`xb52.npy`, `xc52.npy`. Reproduction: `scripts/joint_fit_upper_limit.py`
(deterministic, seed 20260710; argv: jacobian, suffix, SV cut),
`scripts/build_jac_v2.py`, `scripts/rank_gate_v2.py`, plus the WSL-side gate
scripts described in the notes.

## Caveats

Status: Note. (a) Quoted limits use the linearized nuisance response; the
quoted branch is the conservative full-rank estimator, and Gate A bounds the
inefficiency (~6x). (b) Noise is white with release EFAC 1.1225 plus global
scale s = 1.116; chi2_red = 1.25 indicates only mild excess. (c) Unit-drive
convention; physical mapping needs the Request 10.2/10.6 source factors.
(d) The quoted-window edge tau = 10 d carries the minimum u95; the C1 window
itself derives from the sequential 10.7a pilot, and the joint fit shows
better sensitivity below 10 d (diagnostic rows of the curve files) — a
joint-fit re-derivation of the testable window belongs to 10.7d alongside the
truncated-curve promotion. (e) tau is not localizable at limit amplitudes
(expected at ~2 sigma); the 5x-amplitude diagnostic confirms joint tau
recovery works where sequential fails (0.12 d).
