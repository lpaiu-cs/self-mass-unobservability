# Mixed-Symmetry Risk Register

## Purpose

- Status: Proven. This register records whether any mixed-symmetry or otherwise non-STF primitive-family risk remains inside the current theorem domain.

| Sector | Status | Resolution | Residual risk |
| --- | --- | --- | --- |
| Trace descendants of higher-rank ordinary tensors | Proven | Absorbed into audited scalar/vector/STF classes by [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md) | None inside the current theorem domain |
| Mixed-symmetry sectors with odd epsilon-dual parity | Proven | Excluded by the parity-even nonspinning theorem domain through `A2` | Reopens only if parity-odd or spinning sectors are admitted |
| Mixed-symmetry sectors with even epsilon-dual parity | Proven | Reduced to ordinary lower-rank tensor sectors and then absorbed into audited scalar/vector/STF classes by [`../lemmas/50-cartesian-irrep-reduction.md`](../lemmas/50-cartesian-irrep-reduction.md) and [`../lemmas/51-trace-descendant-absorption.md`](../lemmas/51-trace-descendant-absorption.md) | None inside the current theorem domain |
| Explicit unresolved mixed-symmetry primitive family | Proven | None survive after the current irrep and trace-reduction audit | None |

## Verdict

- Status: Proven. No unresolved mixed-symmetry family gate remains inside the current parity-even nonspinning local MVP theorem domain.
