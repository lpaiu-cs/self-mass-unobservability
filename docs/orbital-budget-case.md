# Orbital Budget Case

Status: Counterexample candidate. This note specializes the sample-budget
theorem to the current orbital forcing dictionary.

## Current Dictionary

Status: Imported from prior work. The eccentric-orbit forcing template
`F(Y(t))=(a/r(t))^p` supplies linear input harmonics `n` and `2n` at the
current order.

Status: Imported from prior work. The minimal nonlinear drive can generate a
`3n` line through the product of the `n` and `2n` harmonics.

Status: Proven. The current orbital dictionary therefore has

```math
N_L=2,
\qquad
N_S=1 .
```

The single sideband-pair sample is `(n,2n) -> 3n`.

## Budget Classification

| Comparator class | Required `(N_L,N_S)` | Available `(N_L,N_S)` | Status | Verdict |
| --- | ---: | ---: | --- | --- |
| `N=0, M=0, K_Lambda=0` | `(2,2)` | `(2,1)` | Proven | underbudget no-go |
| `N=1, M=0, K_Lambda=0` | `(3,2)` | `(2,1)` | Proven | underbudget no-go |
| `N=1, M=1, K_Lambda=0` | `(3,4)` | `(2,1)` | Proven | underbudget no-go |
| `N=1, M=2, K_Lambda=0` | `(3,7)` | `(2,1)` | Proven | underbudget no-go |
| `N=1, M=1, K_Lambda=2` | `(3,5)` | `(2,1)` | Proven | underbudget no-go |

## Verdict

Status: Proven. The current `n,2n,3n` orbital case does not beat even the
constant static sideband comparator, because it has only one generated
sideband-pair sample and the `M=0` sideband budget already requires two.

Status: Proven. It also does not beat a realistic first-derivative linear
comparator, because `N=1` requires three linear samples and the current
dictionary has only `n` and `2n`.

Status: Counterexample candidate. The orbital case remains useful as a
forcing dictionary and ratio prototype, but not as a completed
budget-breaking observable design.

Status: Counterexample candidate. A minimum orbital-style upgrade for the
`N=1, M=1, K_Lambda=0` comparator would need at least three linear harmonics
and four generated sideband-pair samples.

Status: Counterexample candidate. The bounded higher-harmonic version of this
upgrade is worked out in [`richer-orbital-forcing-budget-audit.md`](richer-orbital-forcing-budget-audit.md).
