> RAW BOOKKEEPING NOTE
>
> Status: Proven. This file is a supporting family-census table, not an authoritative theorem statement.

# Family-Envelope Table

## Purpose

- Status: Proven. This table records the irreducible family-envelope census under the current MVP free-fall assumptions at `\Delta \le 4`.
- Status: Proven. The table is theorem-domain specific, not universal.

| Class ID | Status | Envelope state | Family group | Irreducible character | Resolution mechanism | Current theorem role |
| --- | --- | --- | --- | --- | --- | --- |
| `Scalar` | Proven | Already audited | Scalar | Rank-0 ordinary tensor | Audited directly as `R0a` and `R0b` | Audited irreducible family class |
| `Vector` | Proven | Already audited | Vector | Rank-1 ordinary tensor | Audited directly as `R1` | Audited irreducible family class |
| `STF2` | Proven | Already audited | STF rank-2 | Rank-2 STF | Audited directly as the rank-2 STF special case | Audited irreducible family class |
| `STFge3` | Proven | Already audited | STF rank-`L \ge 3` | Rank-`L` STF | Covered by the STF self-witness threshold theorem and audited through `L = 6` | Audited irreducible family theorem class |
| `TraceDesc` | Proven | Absorbed by trace reduction | Mixed-symmetry / non-STF | Ordinary tensor with explicit delta traces | Reduced to lower-rank scalar/vector/STF classes | Not a genuinely new primitive family |
| `PseudoOdd` | Proven | Excluded by explicit assumption | Mixed-symmetry / non-STF | Odd epsilon-dual pseudo sector | Excluded by parity-even nonspinning assumption `A2` | Outside the theorem domain |
| `MixedEvenDual` | Proven | Absorbed by irrep reduction | Mixed-symmetry / non-STF | Even-dual mixed-symmetry sector | Dualized to ordinary lower-rank tensors and then reduced to scalar/vector/STF classes | Not a genuinely new primitive family |
| `State` | Proven | Excluded by explicit assumption | Loophole sector | Orbital-timescale internal state variable | Excluded by `A4` | Tracked separately as loophole class |
| `Nonlocal` | Proven | Excluded by explicit assumption | Loophole sector | Hereditary or nonlocal kernel | Excluded by `A3` | Tracked separately as loophole class |

## Reading Rule

- Status: Proven. The only genuinely new irreducible parity-even local primitive-family classes left inside the theorem domain are the already audited scalar/vector/STF classes.
- Status: Proven. No unresolved mixed-symmetry or otherwise non-STF family gate remains.
- Status: Proven. The rank-by-rank STF march is therefore superseded inside the current theorem domain.
