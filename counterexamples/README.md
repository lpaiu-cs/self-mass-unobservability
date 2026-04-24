# Counterexample Classes

The repository tracks loopholes only when they are explicit and tied to a theorem assumption.

| Class | Status | Violated assumption | Minimal idea |
| --- | --- | --- | --- |
| `chi-state` | Counterexample candidate | A4: no orbital-timescale internal state variable | Add a light internal coordinate that couples to the external invariant and carries initial-condition dependence. |
| `nonanalytic-activation` | Counterexample candidate | A5: analytic coupling | Allow a threshold or cusp response that has no regular Taylor jet at the operating point. |
| `hereditary` | Counterexample candidate | A3: locality | Replace the local worldline action by a retarded memory kernel. |
| `nonmetric-clock` | Counterexample candidate | free-fall-only MVP scope | Keep free fall metric, but let the measured clock rate depend on body species. |

## Explicit Smallest Candidate

- Status: Counterexample candidate. [`chi-state/README.md`](chi-state/README.md) is the current smallest explicit loophole model because it changes only one assumption, keeps the action local, and immediately obstructs instantaneous sensitivity collapse.
