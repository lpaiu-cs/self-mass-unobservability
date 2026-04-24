# Lemma 08 Draft: Mixed Time-Derivative Audit

## Question

- Status: Conjectural. Do `Delta<=4` mixed `E-E-D_tau E` contractions survive, or do they reduce to total derivatives?

## Enumerated Mixed Time-Derivative Classes

- Status: Proven. The exhaustive contraction generator [`../symbolic/enumerate_contractions_delta4.py`](../symbolic/enumerate_contractions_delta4.py) finds exactly one cubic mixed class built from `E`, `E`, and `DtE`.

| Operator | Classification | Status | Meaning |
| --- | --- | --- | --- |
| `TrE2DtE` | Proven reducible | Proven | `Tr(E^2 D_tau E)` |

## Reduction

- Status: Proven. By cyclicity of the trace,

```math
D_\tau Tr(E^3) = 3 Tr(E^2 D_\tau E).
```

- Status: Proven. Therefore

```math
Tr(E^2 D_\tau E) \simeq 0
```

modulo worldline total derivatives.

## Related Lower-Weight Time-Derivative Terms

| Operator | Classification | Status | Reduction |
| --- | --- | --- | --- |
| `E_DtE` | Proven reducible | Proven | `D_\tau(E2) = 2 E_DtE` |
| `E_Dt2E` | Proven reducible | Proven | `D_\tau(E_DtE) = dotE2 + E_Dt2E` |
| `dotE2` | Surviving candidate | Conjectural | remains in normal form |

## Audit Result

- Status: Proven. No `Delta<=4` mixed cubic time-derivative obstruction survives from the exact current primitive blocks.
- Status: Proven. The mixed `E-E-D_tau E` sector does not block the `Delta<=4` corrected normal-form path.
