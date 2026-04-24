# Scalar-Family Ordering

## Primitive Family Under Test

- Status: Proven. The scalar-like adequacy attack adjoins one parity-even external scalar family `S` to the corrected `E/B` free-fall sector.
- Status: Proven. At the same fixed-order counting level used elsewhere in the repo, the primitive blocks are
  `S`, `D_\tau S`, `\nabla_i S`, and `D_\tau^2 S`.
- Status: Proven. This is the minimal scalar-like extension that still tests whether a new rank-0 external source can enter the local worldline EFT at `\Delta \le 4`.

## Exact Physical Question

- Status: Proven. The question is not whether the current electric-only or corrected `E/B` exact-current-set theorems are internally consistent.
- Status: Proven. The question is whether a parity-even scalar external family belongs at the same fixed-order level as the currently audited tidal families once one asks for a physically justified minimal free-fall sector.
- Status: Proven. If such a scalar family is admitted unsuppressed, then the smallest new survivor appears immediately at weight `1`.

## Ordering Principles That Could Suppress The Scalar Family

- Status: Proven. Option 1: no suppression. Then `S` belongs at the same fixed-order level as the already audited families, and the first explicit new survivor is `S`.
- Status: Conjectural. Option 2: derivative-only or shift-symmetric coupling. Then bare `S`, `S2`, `S3`, and `S4` are removed, but derivative invariants such as `dotS2`, `gradS2`, `DtS_E2`, and `DtS_B2` can still survive.
- Status: Conjectural. Option 3: background restriction. If the allowed backgrounds satisfy `S = const.` or `S = 0` on the worldline, the entire scalar-family attack is removed by explicit environmental restriction rather than by operator reduction.
- Status: Conjectural. Option 4: matching-level coefficient suppression. A UV completion or matching calculation could set the scalar-family Wilson data to zero, but that must be written as an explicit assumption rather than inferred from the current symmetry list.
- Status: Conjectural. Option 5: separate ordering parameter. One may demote `S` by an additional small parameter not used for `E_ij` and `B_ij`, but that creates a stronger theorem with extra bookkeeping rather than rescuing the present MVP statement for free.

## Current Verdict

- Status: Proven. No scalar-family ordering assumption is currently active in this repository.
- Status: Proven. No shift symmetry, derivative-only coupling rule, or scalar-background restriction is being used tacitly.
- Status: Proven. [`../symbolic/es_sector_delta4.py`](../symbolic/es_sector_delta4.py) shows that the `E/B+scalar` `\Delta \le 4` extension still closes on a finite survivor list.
- Status: Proven. [`../symbolic/es_survivor_rank_check.py`](../symbolic/es_survivor_rank_check.py) shows that the corrected `E/B+scalar` survivor list is linearly independent.
- Status: Proven. Therefore the scalar-like family currently obstructs adequacy claims for the corrected `E/B` exact-current-set theorem, not finite fixed-order closure by itself.
