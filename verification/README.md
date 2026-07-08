# Independent Verification

Independent re-derivations of the load-bearing claims of the free-fall
sensitivity-collapse theorem. These scripts were written **independently** of
`symbolic/` (not a rerun of the repo's own proof scripts); the point is to
reproduce each load-bearing quantitative claim from first principles with a
fresh implementation and compare to the theorem package. Where the two agree,
that is corroboration, not a self-consistency echo.

Dependencies: `numpy`, `sympy` only. Run each with `python verification/<file>.py`.

## What each script checks

| Script | Independently verifies | Result |
| --- | --- | --- |
| `verify_indep.py` | (A) the A5 non-analytic activation counterexample `m_A(Y)=m_0+alpha e^{-1/Y^2}Theta(Y)` has a flat Taylor jet (all derivatives vanish at 0) yet is nonzero; (B) validates the invariant-counter method: a single STF rank-2 tensor has exactly 2 SO(3) invariants (`tr E^2`, `tr E^3`) and `tr E^4 = (1/2)(tr E^2)^2`. | A5 genuine; method valid (count = 2). |
| `verify_identities.py` | the three algebraic reduction identities in `docs/reduction-rules.md`: `Tr E^4 = 1/2 (Tr E^2)^2`, `Tr B^4 = 1/2 (Tr B^2)^2` (numeric), and the mixed quartic `Tr(EBEB) = (E:B)^2 + 1/2 Tr E^2 Tr B^2 - 2 Tr(E^2 B^2)` (exact symbolic). | all three hold; mixed quartic proven exactly (simplify = 0). |
| `verify_survivors.py` | Lemma 09: the seven `Delta<=4` survivors `{E2, E3, E2^2, dotE2, gradE2, divE2, mixedGradE2}` are linearly independent (value-matrix rank), with the nontrivial content being the gradient triple `{gradE2, divE2, mixedGradE2}` (rank 3 precisely because transversality `div E = 0` is not imposed). | full rank 7; E-sector 3; gradient sector 3. Matches `symbolic/survivor_rank_check.py`. |
| `verify_ce_a3a4.py` | the A4 (chi-state) adiabatic collapse coefficients `dL_eff = g^2/(2 mu w^2) Y^2 + g^2/(2 mu w^4) Ydot^2 + ...`; the A3 (hereditary) claim that an exponential kernel is rational/finite-state while a power-law kernel `K_gamma(s) ~ s^{-gamma}` has a non-rational transfer function `~ p^{gamma-1}` (branch point at p=0), so no finite Markovian embedding exists; and that admitting the magnetic family adds an independent witness `B2`. | A4 coefficients exact; A3 power-law non-rational (`p^{-1/2}` at gamma=1/2); B2 independent. |
| `verify_completeness.py` (Tier 1) | **catalog completeness** of the `Delta<=4` electric sector: an independently written contraction enumerator sweeps EVERY building-block signature `(n_E,n_DtE,n_Dt2E,n_G,n_a)` with weight <= 4, enumerates all delta-only perfect matchings, and takes the numeric rank = exact number of independent parity-even invariants per signature; compares to the repo catalog and checks no weight<=4 signature with nonzero invariant dimension is missing. Anchored by the pure-E Hilbert series `1/((1-t^2)(1-t^3))`. | all 15 catalog signatures reproduced to matching dimension; **no missed operators**; `E^4` dim 1 confirms Cayley-Hamilton is the only relation; gradient sector dim 3. |

### Tier 1 status (catalog completeness)

The **electric sector** (`{E, DtE, Dt2E, grad E, a}`) is complete: independently
verified that the repo's `Delta<=4` enumeration misses no parity-even operator
and posits no spurious relation. For this block set, delta-contractibility
(even index count) coincides exactly with parity-even (the odd blocks `grad E`,
`a` carry the parity), so the pure-delta enumerator counts precisely the
theorem's epsilon-free sector.

Remaining Tier-1 families (magnetic `B`, scalar `S`, vector `V`, and the STF
towers `T/Q/U/Z` at ranks 3-6, plus cross-terms) reuse the same enumerator but
require an explicit **parity filter**: e.g. `B` is STF rank-2 (even index count)
but parity-odd, so `E:B` is a pseudoscalar that must be excluded. That filter,
plus each family's exact weight/parity assignment from the `*-family-ordering`
docs, is the mechanical remainder of Tier 1.

## Scope and limits

Verified independently: the electric-sector `Delta<=4` basis, the three
algebraic identities, seven-survivor independence, and the physical core of all
four boundary escapes (A5 exact, A4 exact, A3 exact via Laplace transform,
plus the magnetic no-go witness). Cross-checked by running the repo's own
`symbolic/` scripts (`survivor_rank_check`, `enumerate_contractions_delta4`,
`normal_form_reduce`, the four counterexample demos, `stf_self_witness_check`,
`r3/r4_survivor_rank_check`, `family_envelope_census`) and confirming their
output matches the memos.

**Not independently re-derived here:** catalog exhaustiveness across the full
family envelope (primitive-set adequacy) and the higher-rank STF-tower /
composition-closure chain (lemmas ~22-71). Those remain the residual
verification burden the theorem package itself flags as open; the recommended
route is a Molien/Hilbert-series count of the SO(3) invariant ring graded by
weight, cross-checked with the numeric-rank enumerator extended to all
audited families.
