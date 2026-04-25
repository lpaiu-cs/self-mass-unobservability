# Sideband Observable Targets

Status: Counterexample candidate. This note defines what counts as a sideband
observable in the dynamic-visibility loophole track.

## Does Not Count By Itself

Status: Proven. A generated frequency by itself is not a unique dynamic
signature once static nonlinear comparators are admitted.

Examples:

- Status: Proven. `F^2` can create `2\Omega`.
- Status: Proven. A two-tone input can create `\Omega_1+\Omega_2` and
  `|\Omega_1-\Omega_2|`.
- Status: Proven. An eccentric-orbit input with `n` and `2n` can create `3n`.

Status: Proven. These facts only refute a linear static or linear projection
comparator. They do not refute a static nonlinear local comparator.

## Counts As A Candidate

Status: Counterexample candidate. A sideband becomes a dynamic candidate when
the generated line carries a common relaxation factor,

```math
H_\chi(\nu)=\frac{1}{1+i\nu\tau_\chi},
```

and the same `tau_chi` also fits the linear-frequency transfer data.

Status: Counterexample candidate. The practical targets are:

- a generated sideband phase law `B_\nu/(\nu A_\nu)=\tau_\chi`;
- a ratio such as
  `\hat O(\Omega_1+\Omega_2)/[\hat O(\Omega_1)\hat O(\Omega_2)]`;
- an orbital ratio such as `\hat O(3n)/[\hat O(n)\hat O(2n)]`;
- a multi-sideband sweep whose residual cannot be killed by a finite shared
  static derivative polynomial.

Status: Proven. These targets count only when the sample design breaks the
finite comparator budget in [`sample-budget-theorem.md`](sample-budget-theorem.md).
Below that budget they are static-mimic candidates, not new observable claims.

## Projection Rule

Status: Proven. A linear time-invariant projection can multiply an existing
line by `\Lambda(\nu)`, but it cannot create a line at a frequency that is
absent from its input.

Status: Conjectural. A calibrated or finite-dimensional projection nuisance can
therefore preserve sideband-ratio tests after deprojection.

Status: Proven. If the projection is instead allowed to be an arbitrary complex
nuisance independently at every generated frequency, the ratio test collapses.

## Collapse Boundaries

Status: Proven. The sideband channel collapses as a dynamic discriminator when
any of the following holds:

- `C_dyn=0`, such as `c_chi beta_F2=0` or a vanishing forcing product;
- `tau_chi=0`, so the pole is removed;
- only one generated sideband is sampled and a free static nonlinear
  coefficient is allowed;
- static nonlinear nuisance coefficients are independent at every generated
  frequency;
- the projection nuisance is arbitrary and frequency-local.

Status: Counterexample candidate. The sideband channel remains a theorem-level
loophole candidate only when a finite shared nuisance model is enforced across
multiple lines.

Status: Proven. The current orbital and two-tone dictionaries are
underbudget for realistic first-derivative linear comparators; see
[`orbital-budget-case.md`](orbital-budget-case.md) and
[`two-tone-budget-case.md`](two-tone-budget-case.md).

## Second-Order Nonlinear Upgrade

Status: Counterexample candidate. In the second-order branch, the stronger
sideband target is not just a common relaxation factor. It is a common
quadratic denominator

```math
D_2(z)=\mu_\chi z^2+\gamma_\chi z+\omega_\chi^2
```

shared by the linear response and nonlinear generated response.

Status: Proven. The generated line remains a static-mimic candidate below its
finite nonlinear comparator budget. For a degree-`M` static sideband
polynomial and `K` shared projection nuisances, the generated-sideband sweep
needs at least `M+2+K` complex samples before it can claim a residual.

Status: Counterexample candidate. The nonlinear second-order target and its
current collapse boundaries are recorded in
[`nonlinear-second-order-mode.md`](nonlinear-second-order-mode.md).
