# Lemma 24: Family-Envelope Census

## Statement

- Status: Proven. Work under the currently stated MVP assumptions: parity-even, nonspinning, nearly spherical body data; local worldline EFT; no orbital-timescale internal state variable; fixed order `\Delta \le 4`; finite admitted external family content.
- Status: Proven. The already audited family classes now cover:
  `R2` parity-even unsuppressed rank-2 STF families,
  `R0a` parity-even rank-0 families with bare source,
  `R0b` parity-even derivative-only or shift-symmetric rank-0 families,
  `R1` genuine parity-even rank-1 vector families,
  `Rodd+` via the genuine parity-even STF rank-3 representative,
  and `Reven4+` via the genuine parity-even STF rank-4 representative.
- Status: Proven. The same assumptions explicitly exclude parity-odd, spin-carrying, stateful, and nonlocal family classes.
- Status: Proven. No current assumption excludes additional local parity-even higher-rank tensor families once the rank-3 gate is resolved.
- Status: Proven. Therefore the audited family set still does not yet exhaust the current MVP primitive-family envelope.
- Status: Proven. The smallest remaining unaudited class is now the parity-even odd-rank tensor class `Rodd5+`.

## Derivation

1. Status: Proven. A2 removes parity-odd families and spin- or orientation-carrying structures from the current theorem domain.
2. Status: Proven. A3 removes hereditary or otherwise nonlocal kernels from the primitive-family census.
3. Status: Proven. A4 removes orbital-timescale internal state variables such as `chi_A` from the local primitive-family census.
4. Status: Proven. A7 and A8 fix the order and keep the admitted family set finite, but they do not impose any tensor-rank ceiling.
5. Status: Proven. Rank `0` is already split into the audited subclasses `R0a` and `R0b`.
6. Status: Proven. Rank `1` is now resolved by the genuine-vector audit of [`26-rank1-family-admission.md`](26-rank1-family-admission.md).
7. Status: Proven. Rank `2` is already covered by the audited class `R2`.
8. Status: Proven. The minimal odd-rank tensor gate is now resolved by the genuine rank-3 audit of [`30-rank3-family-admission.md`](30-rank3-family-admission.md).
9. Status: Proven. The minimal even-rank tensor gate is now resolved by the genuine rank-4 audit of [`34-rank4-family-admission.md`](34-rank4-family-admission.md).
10. Status: Proven. The next remaining parity-even local tensor ranks not yet audited are odd ranks `\ge 5`, with `Rodd5+` the smallest.

## Coverage Map

- Status: Proven. `R0a` and `R0b` exhaust only the currently audited scalar-family subclasses.
- Status: Proven. `R1` exhausts only the currently audited genuine parity-even vector-family class.
- Status: Proven. `R2` exhausts only the currently audited parity-even rank-2 STF class.
- Status: Proven. `Rodd+` is now covered only through the current genuine rank-3 STF representative.
- Status: Proven. `Reven4+` is now covered only through the current genuine rank-4 STF representative.
- Status: Proven. The current census still does not justify a claim that all parity-even local higher-rank tensor families are already covered by those audits.

## Boundary

- Status: Proven. This lemma is a census lemma, not a composition theorem for the newly enlarged audited family set.
- Status: Proven. It identifies the next smallest unaudited family class in the census, but it does not by itself move the live theorem bottleneck past the newly audited `Rodd+` obstruction.
- Status: Proven. It only proves that the current assumptions still leave higher-rank parity-even family classes inside the MVP envelope and therefore keep the theorem family-envelope incomplete.
