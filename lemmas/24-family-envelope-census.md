# Lemma 24: Family-Envelope Census

## Statement

- Status: Proven. Work under the currently stated MVP assumptions: parity-even, nonspinning, nearly spherical body data; local worldline EFT; no orbital-timescale internal state variable; fixed order `\Delta \le 4`; finite admitted external family content.
- Status: Proven. The already audited family classes cover:
  `R2` parity-even unsuppressed rank-2 STF families,
  `R0a` parity-even rank-0 families with bare source,
  `R0b` parity-even derivative-only or shift-symmetric rank-0 families,
  and `R1` genuine parity-even rank-1 vector families.
- Status: Proven. The same assumptions explicitly exclude parity-odd, spin-carrying, stateful, and nonlocal family classes.
- Status: Proven. No current assumption excludes additional local parity-even higher-rank tensor families once the vector gate is resolved.
- Status: Proven. Therefore the audited family set does not yet exhaust the current MVP primitive-family envelope.
- Status: Proven. The smallest remaining unaudited class is now a local parity-even rank-3 tensor family.

## Derivation

1. Status: Proven. A2 removes parity-odd families and spin- or orientation-carrying structures from the current theorem domain.
2. Status: Proven. A3 removes hereditary or otherwise nonlocal kernels from the primitive-family census.
3. Status: Proven. A4 removes orbital-timescale internal state variables such as `chi_A` from the local primitive-family census.
4. Status: Proven. A7 and A8 fix the order and keep the admitted family set finite, but they do not impose any tensor-rank ceiling or any restriction to scalar and rank-2 families only.
5. Status: Proven. Rank `0` is already split into the audited subclasses `R0a` and `R0b`.
6. Status: Proven. Rank `1` is now resolved by the genuine-vector audit of [`26-rank1-family-admission.md`](26-rank1-family-admission.md).
7. Status: Proven. Rank `2` is already covered by the audited class `R2`.
8. Status: Proven. The next remaining tensor rank after `{0,1,2}` is rank `3`, and no current assumption excludes that parity-even local tensor class.
9. Status: Proven. Therefore rank `3` is now the smallest remaining unaudited family gate.

## Coverage Map

- Status: Proven. `R0a` and `R0b` exhaust only the currently audited scalar-family subclasses.
- Status: Proven. `R1` exhausts only the currently audited genuine parity-even vector-family class.
- Status: Proven. `R2` exhausts only the currently audited parity-even rank-2 STF class.
- Status: Proven. The current census does not justify a claim that all parity-even local higher-rank tensor families are already covered by those audits.

## Boundary

- Status: Proven. This lemma is a census lemma, not a direct witness audit for the remaining unaudited classes.
- Status: Proven. It does not yet prove that the remaining rank-3 tensor class survives the current reduction rules.
- Status: Proven. It only proves that the current assumptions still leave a higher-rank parity-even family class inside the MVP envelope and therefore keep the theorem family-envelope incomplete.
