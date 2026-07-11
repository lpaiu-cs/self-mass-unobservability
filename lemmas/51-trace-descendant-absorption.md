# Lemma 51: Trace-Descendant Absorption

## Scope

- Status: Note. This lemma applies after the Cartesian irrep reduction of Lemma 50.
- Status: Proven. It concerns only ordinary parity-even tensor sectors that still carry explicit `\delta_{ij}` traces.

## Claim

- Status: Proven. Trace descendants do not define genuinely new primitive families.
- Status: Proven. They reduce to already audited lower-rank scalar, vector, or STF classes.

## Proof

- Status: Proven. A single trace contraction lowers the ordinary tensor rank by `2`.
- Status: Proven. Repeated trace stripping therefore terminates after finitely many steps at rank `0` or rank `1`, or at a trace-free ordinary tensor.
- Status: Proven. Rank `0` and rank `1` outputs lie in the already audited scalar and vector classes.
- Status: Proven. A trace-free ordinary tensor lies in an STF class by Lemma 50.
- Status: Proven. Rank `2` trace-free outputs lie in the already audited rank-2 STF special class.
- Status: Proven. Rank `L \ge 3` trace-free outputs lie in the audited higher-rank STF threshold class.

## Required Assumptions

- Status: Proven. Primitive-family distinctness is tested after linear irrep decomposition rather than before it.
- Status: Proven. Trace descendants are not double-counted as independent primitive-family admissions once their lower-rank outputs are already present in the catalog.

## Verdict

- Status: Proven. Trace descendants of higher-rank parity-even tensors are absorbed by the already audited scalar/vector/STF classes and therefore do not obstruct irreducible-envelope closure.
