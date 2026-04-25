# High-Rank Audit Methodology

## Purpose

- Status: Proven. The rank-3, rank-4, and rank-5 family notes can no longer rely on manual survivor bookkeeping alone.
- Status: Proven. Before the family-envelope march moves upward, the current higher-rank survivor notes must be checked for exhaustiveness under the same primitive-family definitions that were used to declare the obstruction classes.

## Two Layers

- Status: Proven. `manual survivor bookkeeping` means the currently listed low-order witnesses and auxiliary survivors written directly into [`../symbolic/r3_sector_delta4.py`](../symbolic/r3_sector_delta4.py), [`../symbolic/r4_sector_delta4.py`](../symbolic/r4_sector_delta4.py), and [`../symbolic/r5_sector_delta4.py`](../symbolic/r5_sector_delta4.py).
- Status: Proven. `exhaustive contraction generation` means generating every parity-even scalar contraction class in the fixed-order candidate-survivor layer from the chosen primitive STF family `Y_L` together with the current baseline electric sector.
- Status: Proven. The present methodology patch is about exhaustiveness, not about re-proving the already recorded obstruction-class verdicts.

## Candidate-Survivor Layer

- Status: Proven. The exhaustive generator works in the explicit normal-form candidate layer rather than the unreduced raw signature layer.
- Status: Proven. In this layer, derivative-generated descendants are not treated as new primitive-family admissions.
- Status: Proven. The pure acceleration sector is not part of the current higher-rank family audit, because those terms are already assigned to lower-order EOM reduction rather than to new primitive-family survivors.
- Status: Proven. The current high-rank candidate layer keeps the following families of signatures:
  `Y_L^2`, `E Y_L^2`, `\dot Y_L^2`, `E Y_L \dot Y_L`, `E^2 Y_L^2`, `(\nabla Y_L)^2`, and `Y_L^4`, together with any additional parity-even scalar signatures at the same order that the exhaustive generator actually finds.
- Status: Proven. Total-derivative normal-form choices remove signatures such as `Y_L \dot Y_L` and `Y_L \ddot Y_L` from the independent candidate layer.

## Representative Generation

- Status: Proven. For signatures built only from fully symmetric trace-free blocks, scalar contraction classes are generated as colored multigraphs with fixed vertex degrees rather than by brute-force slot permutations.
- Status: Proven. This colored-multigraph representation is exact for the present higher-rank STF candidate layer because each STF block is symmetric in all tensor slots and forbids trace loops.
- Status: Proven. For the two-block gradient signatures `(\nabla Y_L)^2` and, when allowed, `\nabla E \nabla Y_L`, the exhaustive generator uses explicit derivative-versus-STF edge-count descriptors rather than a manually curated representative list.

## Output Rule

- Status: Proven. [`../symbolic/high_rank_family_enumerator.py`](../symbolic/high_rank_family_enumerator.py) is the new exhaustive generator for `L = 3, 4, 5`.
- Status: Proven. [`../symbolic/high_rank_diff_report.py`](../symbolic/high_rank_diff_report.py) is the machine-readable comparison layer between the current manual notes and the exhaustive generated candidate lists.
- Status: Proven. If the exhaustive generator finds an omitted contraction class, that omission must be recorded explicitly even if the obstruction-class verdict or suppression threshold does not change.
- Status: Proven. Therefore high-rank family progress now has three distinct bookkeeping layers:
  obstruction-class verdict,
  threshold formula,
  and corrected-basis or raw-survivor bookkeeping.
