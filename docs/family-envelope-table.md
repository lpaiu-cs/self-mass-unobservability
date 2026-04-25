# Family-Envelope Table

## Purpose

- Status: Proven. This table records which primitive-family classes are already covered, excluded by explicit assumptions, or still unaudited under the current MVP free-fall assumptions at `\Delta \le 4`.
- Status: Proven. The table is a census tool, not a universal theorem.

| Class ID | Status | Envelope state | Tensor rank | Parity | Derivative character | Smallest expected witness type | Current theorem role |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `R2` | Proven | Audited | `2` STF | Even | Unsuppressed local family | `X2` and `EX` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R0a` | Proven | Audited | `0` | Even | Bare source allowed | `S` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R0b` | Proven | Audited | `0` | Even | Derivative-only or shift-symmetric | `dotS2` and `DtS_E2` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R1` | Proven | Audited | `1` | Even | Unsuppressed local family | `V2` audited explicitly; first mixed witness `EVV` | Audited uniqueness obstruction class and threshold class |
| `Rodd+` | Proven | Audited | Odd `\ge 3`; current representative = STF rank `3` | Even | Unsuppressed local family | `T2` audited explicitly; first mixed witness `ETT` | Audited uniqueness obstruction class and threshold class; enlarged audited-set composition now re-closed |
| `Reven4+` | Proven | Audited | Even `\ge 4`; current representative = STF rank `4` | Even | Unsuppressed local family | `Q2` audited explicitly; manual `EQQ` plus omitted exhaustive `EEQ` at the same mixed order | Audited uniqueness obstruction class and threshold class; enlarged audited-set composition now re-closed, but exhaustive bookkeeping patch is still active |
| `Rodd5+` | Proven | Audited | Odd `\ge 5`; current representative = STF rank `5` | Even | Unsuppressed local family | `U2` audited explicitly; first mixed witness `EUU` | Audited uniqueness obstruction class and threshold class; enlarged audited-set composition now re-closed |
| `Reven6+` | Proven | Still unaudited | Even `\ge 6` | Even | Unsuppressed local family | Rank-6 self-contraction type (Counterexample candidate) | Next smallest unaudited family gate after the closed `Rodd5+` composition layer |
| `Podd` | Proven | Excluded by explicit assumption | Any | Odd | Any local derivative character | Not in the MVP domain | Excluded by A2 |
| `Spin` | Proven | Excluded by explicit assumption | Any | Any | Spin- or orientation-carrying | Not in the MVP domain | Excluded by A2 |
| `State` | Proven | Excluded by explicit assumption | Not a local tensor family | Any | Orbital-timescale internal state variable | `chi_A`-type loophole | Excluded by A4; recorded as loophole class instead |
| `Nonlocal` | Proven | Excluded by explicit assumption | Not a local primitive family | Any | Hereditary or nonlocal kernel | Retarded-kernel loophole | Excluded by A3; recorded as loophole class instead |

## Reading Rule

- Status: Proven. The rows `R2`, `R0a`, `R0b`, `R1`, `Rodd+`, and `Reven4+` are now the explicitly audited primitive-family classes.
- Status: Proven. The excluded rows are outside the current theorem domain because the required assumptions already remove them.
- Status: Proven. The remaining unaudited rows are the reason family-envelope completeness still fails under the current MVP assumptions.
- Status: Proven. Among the remaining unaudited rows, `Reven6+` is now the smallest next envelope gate.
- Status: Proven. The family-envelope census still points to `Reven6+` as the next unaudited gate, but the live theorem bottleneck is the failed STF tower abstraction at rank `L = 4`.
