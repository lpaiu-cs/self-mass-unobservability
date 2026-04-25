# Family-Envelope Table

## Purpose

- Status: Proven. This table records which primitive-family classes are already covered, excluded by explicit assumptions, or still unaudited under the current MVP free-fall assumptions at `\Delta \le 4`.
- Status: Proven. The table is a census tool, not a universal theorem.

| Class ID | Status | Envelope state | Tensor rank | Parity | Derivative character | Smallest expected witness type | Current theorem role |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `R2` | Proven | Audited | `2` STF | Even | Unsuppressed local family | `X2` and `EX` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R0a` | Proven | Audited | `0` | Even | Bare source allowed | `S` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R0b` | Proven | Audited | `0` | Even | Derivative-only or shift-symmetric | `dotS2` and `DtS_E2` audited explicitly | Audited uniqueness obstruction class and threshold class |
| `R1` | Proven | Still unaudited | `1` | Even | Unsuppressed local family | `V2`-type self contraction (Counterexample candidate) | Next smallest unaudited family gate; blocks envelope completeness |
| `Rodd+` | Proven | Still unaudited | Odd `\ge 3` | Even | Unsuppressed local family | `Y2`-type self contraction (Counterexample candidate) | Later unaudited family gate |
| `Reven4+` | Proven | Still unaudited | Even `\ge 4` | Even | Unsuppressed local family | `Z2`-type self contraction (Counterexample candidate) | Later unaudited family gate |
| `Podd` | Proven | Excluded by explicit assumption | Any | Odd | Any local derivative character | Not in the MVP domain | Excluded by A2 |
| `Spin` | Proven | Excluded by explicit assumption | Any | Any | Spin- or orientation-carrying | Not in the MVP domain | Excluded by A2 |
| `State` | Proven | Excluded by explicit assumption | Not a local tensor family | Any | Orbital-timescale internal state variable | `chi_A`-type loophole | Excluded by A4; recorded as loophole class instead |
| `Nonlocal` | Proven | Excluded by explicit assumption | Not a local primitive family | Any | Hereditary or nonlocal kernel | Retarded-kernel loophole | Excluded by A3; recorded as loophole class instead |

## Reading Rule

- Status: Proven. The three audited rows are exactly the family classes whose thresholds are already jointly sufficient at `\Delta \le 4`.
- Status: Proven. The excluded rows are outside the current theorem domain because the required assumptions already remove them.
- Status: Proven. The unaudited rows are the reason audited-set sufficiency does not yet upgrade to MVP-envelope sufficiency.
- Status: Proven. Among the unaudited rows, `R1` is currently the smallest next gate.
