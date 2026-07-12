# Independent Verification

> **CORRECTION (2026-07-12, gradient-sector kinematics).** An external
> adversarial review found that every method below modeled `grad E` as a
> generic `(STF-2) x vector` object (15 comps), dropping the Schwarz total
> symmetry of `partial_k partial_i partial_j Phi` and the vacuum trace
> condition already imposed on `E` itself. With the correct STF-3 model
> (7 comps): `divE2 = 0`, `mixedGradE2 = gradE2`, and the electric survivor
> basis is **5**, not 7 (`{E2, E3, E2^2, dotE2, gradE2}`). Every sector loses
> the two spurious electric gradient survivors, and the `E/B/S` sector
> additionally loses `divEGradS` (it contains the vanishing `div E`). The
> corrected family table — cross-checked by the exact character method and
> the repo's own symbolic modules — is
> **E,B,S,V,T,Q,U,Z = 5,16,30,15,17,23,17,21**
> (previously 7,18,33,17,19,25,19,23). Admitted families' own gradient
> blocks stay generic: they are independent primitives with no assumed
> potential structure. `GradEGradQ` survives the correction (shared l=3).

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
| `verify_survivors.py` | Lemma 09 (corrected): the five `Delta<=4` survivors `{E2, E3, E2^2, dotE2, gradE2}` are linearly independent (value-matrix rank), with `grad E` modeled as an STF-3 octupole; verifies `divE2 == 0` and `mixedGradE2 == gradE2` to machine precision, and carries a diagnostic control showing the old generic model reproduces the defective gradient rank 3. | full rank 5; E-sector 3; gradient sector 1; generic-model control 3. Matches `symbolic/survivor_rank_check.py`. |
| `verify_ce_a3a4.py` | the A4 (chi-state) adiabatic collapse coefficients `dL_eff = g^2/(2 mu w^2) Y^2 + g^2/(2 mu w^4) Ydot^2 + ...`; the A3 (hereditary) claim that an exponential kernel is rational/finite-state while a power-law kernel `K_gamma(s) ~ s^{-gamma}` has a non-rational transfer function `~ p^{gamma-1}` (branch point at p=0), so no finite Markovian embedding exists; and that admitting the magnetic family adds an independent witness `B2`. | A4 coefficients exact; A3 power-law non-rational (`p^{-1/2}` at gamma=1/2); B2 independent. |
| `verify_completeness.py` (Tier 1) | **catalog completeness** of the `Delta<=4` electric sector: an independently written contraction enumerator sweeps EVERY building-block signature `(n_E,n_DtE,n_Dt2E,n_G,n_a)` with weight <= 4 (with `G` an STF-3 octupole), enumerates all delta-only perfect matchings, and takes the numeric rank = exact number of independent parity-even invariants per signature; compares to the repo catalog and checks no weight<=4 signature with nonzero invariant dimension is missing. Anchored by the pure-E Hilbert series `1/((1-t^2)(1-t^3))`. | all 14 catalog signatures reproduced to matching dimension; **no missed operators**; `E^4` dim 1 confirms Cayley-Hamilton is the only relation; gradient sector dim 1; the `(G,a)` signature vanishes identically. |

### Tier 1 status (catalog completeness)

The **electric sector** (`{E, DtE, Dt2E, grad E, a}`) is complete: independently
verified that the repo's `Delta<=4` enumeration misses no parity-even operator
and posits no spurious relation. For this block set, delta-contractibility
(even index count) coincides exactly with parity-even (the odd blocks `grad E`,
`a` carry the parity), so the pure-delta enumerator counts precisely the
theorem's epsilon-free sector.

### Tier 1 — all audited families (`stf.py`, `verify_family_survivors.py`)

`stf.py` builds a general symmetric-trace-free rank-`L` tensor generator
(nullspace of the trace map; self-tested to dimension `2L+1` for `L=0..6`).
`verify_family_survivors.py` computes each sector's **survivor dimension**
independently — using the theorem's exact truncated block set
`{E,DtE,Dt2E,GradE} + {X,DtX,Dt2X,GradX}` and the correct total-derivative
(`D_tau`) quotient with the lower-order EOM (`a=0`) — and compares to the
repo's own survivor rank. Parity is handled exactly per family (`E,S,V,T,Q,U,Z`
parity-even; the physical magnetic `B` parity-odd, so `E:B` is excluded but
`EB2` survives).

| sector | family rank | independent survivor dim | their rank | verdict |
| --- | --- | --- | --- | --- |
| electric `E` | – | 5 | 5 | MATCH |
| `E/B` (magnetic) | 2 (odd) | 16 | 16 | MATCH |
| `E/B/S` (scalar) | 0 | 30 | 30 | MATCH |
| `E/V` (vector) | 1 | 15 | 15 | MATCH |
| `E/T` | 3 | 17 | 17 | MATCH |

So the `Delta<=4` enumeration misses no operator and posits no spurious survivor
for every family of rank `0-3`, including the physical magnetic sector and the
3-family `E/B/S` composition.

**High-rank `Q/U/Z` (ranks 4/5/6) — exact character method.**
`tier1_survivor_exact.py` computes each sector's survivor dimension with **exact
integer `O(3)` character integrals only** (no matchings, no numeric rank), via
the identity `survivor(w) = dim inv_trunc(w) - dim inv_prom(w-1)` (validated on
the electric sector = 5). The `sig_dim` gate is the load-bearing subtlety: a
delta-only scalar needs an **even total Cartesian index count**, else it would
require an epsilon (pseudoscalar, excluded) — without this gate the character
count spuriously includes gradient-block pseudoscalars. With the gate and the
corrected STF-3 `GradE` content (single `l=3` irrep), the method reproduces
**all 8 audited sectors**: `E,B,S,V,T,Q,U,Z = 5,16,30,15,17,23,17,21`,
resolving the numeric method's spin-4+ precision limit.
(`tier1_highrank_molien.py` separately confirms the exact rank-4/5/6 quartic
relations `EEQQ=3, QQQQ=2, ...` behind their nullities.)

**Finding + re-derivation — rank-4 (`Q`), two-step correction 19 -> 25 -> 23.**
Step 1: the exact character method gave `survivor = 25` for `E/Q` versus the
repo's audited `19`; the numeric brute-force independently also gave `25`.
`rederive_rank4.py` closes the gap legitimately (no new rule — only the repo's
own total-derivative / EOM / STF identities): it constructs the omitted
survivors and proves each is a nonzero, rotation-invariant, pure-primitive
(non-total-derivative) weight-4 scalar — `Q_abcd (E^2)_ab E_cd` (degree 3 in
`E`), `E_ab Q...Q...Q...` (degree 3 in `Q`), etc. The repo's `r4` mixed
candidates stopped at degree 2 in `E` (`E2Q2`), so its enumeration never
generated these higher-degree mixed operators — exactly the sector the repo's
own docs flagged as the "isolated rank-4 exception" (`family-class-table.md`,
`lemma 43`). Step 2: the 2026-07-12 gradient-kinematics correction then
removed the two spurious electric gradient survivors from the sector
(`divE2`, `mixedGradE2`), landing on the corrected value **23** = 5 (electric)
+ 18 (new); `GradEGradQ` itself survives both steps.

**Fixed in the repo scripts.** The 6 omitted survivors (`EEQ, QQQ, E3Q, EQ3,
EDtEQ, GradEGradQ`) are in `symbolic/r4_survivor_rank_check.py`
(`_r4_new_specs`, evaluator extended to `DtE`/`GradE`) and
`symbolic/r4_sector_delta4.py`; with the STF-3 `GradE` sampler the module now
reports `Total survivor rank: 23 out of 26`, cross-checked in-module by the
exact character integral (23, matches). The new-sector rank remains 18
(lifted from 12 by the six restored operators).

### Tier 2 — family-envelope closure (`tier2_irrep_census.py`, `tier2_tower.py`)

- **Irrep census / trace absorption** (`tier2_irrep_census.py`): confirms the
  STF-`L` are the complete irreducible parity-even set (dim `2L+1`), a symmetric
  rank-`r` tensor is `STF_r + STF_{r-2} + ...` (trace descendants are strictly
  lower STF), and every Cartesian rank-`r` tensor decomposes into `STF_L` pieces
  only (`sum = 3^r`) — so `TraceDesc` and `MixedEvenDual` genuinely reduce to
  STF + traces and no exotic irreducible family exists. The envelope classes are
  exhaustive up to trace/symmetry reduction.
- **Uniform STF-tower closure** (`tier2_tower.py`): the infinite tower `L>=3` is
  closed by a single uniform self-witness structure, verified via exact Molien
  for `L=2..12` (past the audited `L<=6`): for every `L>=3` the smallest self
  witness is `X2` at weight 2, `E:X` vanishes (no weight-2 mixed), and the
  smallest mixed witness is `EXX` at weight 3 — no new lower witness ever
  appears, so the tower closes with one theorem. `L=2` is the sole special case
  (`E:X` at weight 2 = the `R2` mixed witness).
- **Composition closure**: the `E/B/S` = 33 match above is a genuine 3-family
  admission-level composition with no unexpected cross-family survivor; the
  repo's `composition_attack_delta4.py` (all combos "sufficient / none") was
  cross-run and is consistent.

## Scope and limits

Verified independently: the electric-sector `Delta<=4` basis (corrected,
five survivors), the three algebraic identities, the STF-3 gradient
kinematic identities (`divE2 = 0`, `mixedGradE2 = gradE2`), five-survivor
independence, all four boundary escapes (A5/A4/A3 exact, A8), the magnetic
no-go witness, **catalog completeness for every audited family of rank 0-3
(exact survivor-dimension agreement) and the rank 4-6 relation structure
(exact Molien nullities)**, and the family-envelope closure (irrep census +
uniform tower). Cross-checked against the repo's own `symbolic/` scripts
throughout.

**Not independently re-derived here:** a fully exact (non-numeric) survivor
dimension for ranks 4-6 — the numeric quotient is precision-limited there, so
those ranks rest on exact per-signature Molien dims + the repo's exact symbolic
`survivor_rank_check`. The full 71-lemma composition-closure chain
(lemmas ~22-71) is corroborated but not re-derived line by line.
