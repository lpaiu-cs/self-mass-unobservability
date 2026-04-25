# collapse-theorem

Status: Conjectural. The repository now targets a theorem-first follow-up: prove or refute that self-mass-unobservability hypotheses collapse to a finite-dimensional sensitivity-manifold EFT in the free-fall sector under explicit assumptions.

## Success Modes

- Status: Conjectural. Success mode 1 is a proof that the leading body-dependent free-fall deviation reduces to sensitivity coordinates `(s_A,I, s_A,IJ, ...)` plus higher-multipole Wilson coefficients inside a local worldline EFT.
- Status: Counterexample candidate. Success mode 2 is an explicit smallest loophole model together with the exact theorem assumption it violates.

## MVP Scope

The MVP is free-fall only. The clock sector is documented only conditionally in `lemmas/04-clock-conditional-note.md`, and scalar `s_A` is treated as a corollary rather than as a starting axiom.

Out of scope for this phase:

- LLR / MLRS / PEP / Nutimo / TOA follow-up work
- runtime, build, or environment triage
- empirical claims without a derivation or a cited artifact already in this repository

## Imported Starting Point

- Status: Imported from prior work. [`request2/REQUEST2_INTERNAL_STRUCTURE.md`](request2/REQUEST2_INTERNAL_STRUCTURE.md) shows that literal internal self-unobservability fails as a self-bound equilibrium principle.
- Status: Imported from prior work. [`request1/REQUEST1_COM_DECOUPLING.md`](request1/REQUEST1_COM_DECOUPLING.md) shows that center-of-mass self-subtraction does not generate a new monopole force.
- Status: Imported from prior work. [`request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md`](request7/REQUEST7_JOINT_CONSISTENCY_SCAFFOLD.md) is only a provisional consistency scaffold, not the theorem target of this repository.

## Current Layout

- [`docs/`](docs): theorem candidates, electric and `E/B` conditional collapse notes, power counting, the primitive catalog, magnetic-ordering and primitive-set adequacy notes, dynamic sideband, shared-tau, and sample-budget comparator notes, explicit reduction rules, ledgers, and the roadmap.
- [`lemmas/`](lemmas): free-fall lemmas and reduction notes.
- [`symbolic/`](symbolic): small SymPy scripts for worldline expansion, sensitivity jets, basis enumeration, contraction enumeration, electric and `E/B` survivor-rank checks, primitive-family attacks, full `E/B`-sector audits, normal-form reduction, and dynamic-comparator checks.
- [`counterexamples/`](counterexamples): explicit loophole classes and minimal models.
- [`notes/scratch/`](notes/scratch): disposable local notes that should not become theorem claims without promotion.

Legacy request-stage material remains in `request1/` through `request7/`, `data/`, `paper/`, and `section6/` as provenance only. The root workflow should now start from the theorem documents above, not from the old pipeline entry points.

## Quick Start

Use Python `3.10+`.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Run the current symbolic scaffold:

```powershell
make worldline-expand
make enumerate-basis
make enumerate-contractions
make survivor-rank
make primitive-attack
make eb-sector
make eb-rank
make normal-form-reduce
make nonlinear-comparator
make shared-tau-ratio
make sample-budget
make orbital-harmonic-budget
make symbolic-check
```

## Working Rule

Every nontrivial scientific claim in the new theorem-first files should be labeled as one of:

- `Proven`
- `Imported from prior work`
- `Conjectural`
- `Counterexample candidate`

If a proof attempt breaks, record the exact failed step and the smallest missing assumption in [`docs/failure-ledger.md`](docs/failure-ledger.md) rather than smoothing it over in prose.
