# Two-Tone Budget Case

Status: Counterexample candidate. This note specializes the sample-budget
theorem to a two-tone forcing dictionary.

## Current Dictionary

Status: Proven. With two input tones `Omega_1` and `Omega_2`, the linear
sample count is

```math
N_L=2 .
```

Status: Proven. Counting sum-sideband pairs gives

```math
(\Omega_1,\Omega_1),\quad
(\Omega_1,\Omega_2),\quad
(\Omega_2,\Omega_2),
```

so the conservative sideband-pair count is `N_S=3`.

Status: Proven. If the difference sideband is treated as one additional
independent pair `(Omega_1,-Omega_2)`, then `N_S=4`. This is an optimistic
count because real-signal conjugacy and projection calibration can reduce the
independent information.

## Budget Classification

| Comparator class | Required `(N_L,N_S)` | Sum-only available | Sum+difference available | Status | Verdict |
| --- | ---: | ---: | ---: | --- | --- |
| `N=0, M=0, K_Lambda=0` | `(2,2)` | `(2,3)` | `(2,4)` | Proven | beats only a trivial comparator |
| `N=1, M=0, K_Lambda=0` | `(3,2)` | `(2,3)` | `(2,4)` | Proven | underbudget no-go |
| `N=1, M=1, K_Lambda=0` | `(3,4)` | `(2,3)` | `(2,4)` | Proven | underbudget no-go |
| `N=1, M=2, K_Lambda=0` | `(3,7)` | `(2,3)` | `(2,4)` | Proven | underbudget no-go |
| `N=1, M=1, K_Lambda=2` | `(3,5)` | `(2,3)` | `(2,4)` | Proven | underbudget no-go |

## Verdict

Status: Proven. The current two-tone dictionary can beat only the trivial
`N=0, M=0, K_Lambda=0` comparator.

Status: Proven. Once a first-derivative linear comparator is admitted, the
two-tone dictionary is underbudget because it has only two linear samples and
`N=1` requires three.

Status: Proven. For the realistic `N=1, M=1, K_Lambda=0` class, the minimum
design is `N_L=3` and `N_S=4`; the current sum+difference count can reach the
sideband requirement but still fails the linear-sample requirement.

Status: Counterexample candidate. Two-tone forcing is therefore a useful
local diagnostic, but not a complete budget-breaking design unless extended
to at least one additional distinct linear tone or an externally calibrated
linear transfer constraint.

