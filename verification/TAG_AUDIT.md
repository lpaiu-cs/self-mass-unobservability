# Status-Tag Audit

Every claim in `docs/` and `lemmas/` carries a `Status: Proven` / `Status:
Conjectural` tag. A scan finds **2076 `Proven`** vs **125 `Conjectural`**. The
tag is applied so broadly — to posited assumptions, to bookkeeping notes, and
to claims that are actually open — that `Proven` no longer distinguishes
machine/hand-verified mathematics from framing prose. This will draw reviewer
skepticism and should be tidied before publication. Concrete issues and
recommended fixes:

## 1. Assumptions tagged inconsistently (`Proven` vs `Conjectural`)

- `docs/assumptions-ledger.md` marks each of `A1`-`A8` **`Conjectural`** (correct
  in spirit: they are posited hypotheses, not theorems).
- `docs/theorem-package.md:11` lists the same `A1`-`A8` under **`Status:
  Proven`** ("The domain assumptions are the active A1-A8 ledger assumptions
  ...").

An assumption is neither proven nor conjectural — it is an **assumption**.
**Fixed:** `A1`-`A8` are now tagged `Status: Assumption` in
`assumptions-ledger.md` (with a header note explaining the label), and the
`theorem-package.md:11` line that lists them now reads `Status: Assumption`.
`Proven`/`Conjectural` are reserved for derived claims. The remaining broad use
of `Proven` for bookkeeping/history lines (§3) is not yet swept.

## 2. Envelope-closure claim tagged `Proven` but is INCOMPLETE at rank 4

- `docs/theorem-package.md:12` — **`Status: Proven`**: "the irreducible
  primitive-family envelope closes on the audited scalar, vector, rank-2 STF,
  and genuine rank-`L >= 3` STF classes."

Independent exact character-integral verification
(`verification/tier1_survivor_exact.py`) reproduces the survivor dimension for
7 of 8 audited sectors exactly (`E,B,S,V,T,U,Z = 7,18,33,17,19,19,23`) but finds
the **rank-4 (`Q`) sector incomplete**: exact survivor dimension `25` vs the
repo's audited `19`. A confirmed missing survivor is `Q_abcd (E^2)_ab E_cd`
(degree 3 in `E`, 1 in `Q`): nonzero, rotation-invariant, pure-primitive (not a
total derivative), but the repo's `r4` mixed candidates cap at degree 2 in `E`
(`E2Q2`). This coincides with the repo's own flag of an "isolated rank-4
exception" needing a "manual exhaustive patch" (`docs/family-class-table.md`,
`lemmas/43-*`).

The rank-4 sector was **re-derived** (`verification/rederive_rank4.py`) and the
repo's `r4` scripts have now been **fixed**: under the theorem's own rules the
correct survivor dimension is **25, not 19**. The 6 omitted higher-degree mixed
survivors — `EEQ`, `QQQ`, `E3Q`, `EQ3`, `EDtEQ`, `GradEGradQ` — were added to
`symbolic/r4_survivor_rank_check.py::_r4_new_specs()` and
`symbolic/r4_sector_delta4.py::_r4_family_classes()` (with the evaluator
extended to the `DtE`/`GradE` blocks). Both scripts now report `Total survivor
rank: 25`, and an exact O(3) character integral cross-checks it. Ranks
`0-3,5,6` were already exact-confirmed (`7,18,33,17,19,23`); the smoke test and
the other sector scripts still pass.

**Remaining tag fix:** with the `r4` scripts now at `25`, `theorem-package.md:12`
(envelope closure) may keep `Proven` for the *dimension*, but note that the
`family-class-table.md` rank-4 row and any prose quoting "19" should be updated
to `25`. The `A1`-`A8` assumption/`Proven` inconsistency in §1 above still
stands.

## 3. `Proven` used for bookkeeping / historical notes

Many `Status: Proven` lines are not theorems but records (e.g. "This note
records the pre-M4 internal reduction attempt", `lemmas/06`, which then contains
`Conjectural` steps). Mixing "the proof step is proven" with "it is proven that
we recorded this" dilutes the tag.

**Fixed:** a `Status: Note` label is now defined in
[`../docs/status-labels.md`](../docs/status-labels.md) and applied by a
conservative sweep (353 lines): the whole `archive/` tree (superseded history)
and every self-describing meta line (`This note/file/document/table/... `) in
active docs/lemmas are now `Note`; genuine mathematical `Proven` claims are
untouched. The borderline `This lemma ...` / `The repository does not assume ...` lines
were then **individually reviewed** and swept (110 more lines): scope caveats,
disclaimers ("does not ...", "not a universal ..."), task/role descriptions,
and bookkeeping -> `Note`; genuine positive results ("proves/derives ...",
"model is/preserves ...", "minimally irreducible ... class", "self-only ...",
"necessary conditions", "EEQ class exists ...") kept `Proven`. Finally, the `assumptions-ledger.md` scope/rule sections ("Non-Assumption For
M*", "Explicit Rule For *", "Explicit Representative Choice For *", the boundary
stress-test rules) were swept as well: all 57 `Status: Proven` lines there are
scope declarations, rule records, representative choices, or methodology
bookkeeping -> `Note`. The ledger now carries only `Assumption` (A1-A8),
`Imported from prior work` (E1-E3), and `Note`; the mathematical proofs of these
facts live in the `lemmas/` files, which keep `Proven`.

## What the tags *should* read, per this verification

Confirmed `Proven` (independently reproduced in `verification/`):
- COM decoupling (Request 1), internal-structure no-go (Request 2).
- The three algebraic reductions (`Tr E^4`, `Tr B^4`, and the exact `Tr(EBEB)`
  identity); the seven-scalar `Delta<=4` electric basis and its linear
  independence; electric-sector catalog completeness.
- All four boundary escapes (`A5`, `A3`, `A4`, `A8`).
- Family completeness for ranks `0-3,5,6` (exact survivor dimensions match).
- Irrep/trace-absorption census; the uniform STF-tower self-witness structure
  (`L=2..12`).

Should NOT read `Proven` yet:
- Rank-4 (`Q`) family completeness / the unqualified envelope-closure claim
  (incomplete — see §2).
- `A1`-`A8` (assumptions, not theorems — see §1).

Everything above is reproducible from the scripts in this directory.
