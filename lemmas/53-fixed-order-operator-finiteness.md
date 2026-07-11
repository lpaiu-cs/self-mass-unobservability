# Lemma 53: Fixed-Order Operator Finiteness

## Statement

- Status: Proven. Fix the current theorem domain and let `{\cal F} = \{X_a\}_{a=1}^{N_F}` be any explicit finite admitted primitive-family catalog after irreducible family-envelope closure.
- Status: Proven. Assume the fixed-order counting rule of [`../docs/power-counting.md`](../docs/power-counting.md): each primitive species has positive intrinsic weight `w(X_a) \ge 1`, while `w(D_\tau) = w(\nabla) = 1`, and only operators with total weight `\Delta \le 4` are retained.
- Status: Proven. Then the candidate local parity-even free-fall scalar operator space is finite before any reduction is applied.

## Proof

1. Status: Proven. Because the admitted primitive-family catalog is finite, there are only finitely many primitive species `X_a`.
2. Status: Proven. For each species `X_a`, the decorated descendants `D_\tau^p \nabla^q X_a` that can contribute at `\Delta \le 4` satisfy `p + q \le 4 - w(X_a)`.
3. Status: Proven. Since `w(X_a) \ge 1`, the set of allowed pairs `(p, q)` is finite for each `X_a`, so the total decorated block set `{\cal B}_{\le 4}` is finite.
4. Status: Proven. Every retained scalar operator is a contraction of blocks from `{\cal B}_{\le 4}` whose total weight is at most `4`.
5. Status: Proven. Because every block weight is positive, such an operator contains at most `4` factors.
6. Status: Proven. For each finite multiset of factors, the total tensor rank is finite because each admitted primitive family has finite spatial rank and only finitely many derivatives are allowed.
7. Status: Proven. A finite-rank tensor multiset admits only finitely many complete contraction patterns.
8. Status: Proven. Therefore only finitely many parity-even scalar contraction classes can appear at `\Delta \le 4`.

## Role Of Irreducible-Envelope Closure

- Status: Proven. Irreducible-envelope closure is not the source of finiteness by itself.
- Status: Proven. Its role is to ensure that trace descendants and mixed-symmetry sectors are not double-counted as extra primitive families before the finite-family count is done.
- Status: Proven. Once that closure is in place, the only genuinely new primitive-family types in the current theorem domain are the already audited scalar, vector, and STF classes.

## What Remains After This Lemma

- Status: Proven. This lemma proves finiteness before reduction.
- Status: Proven. The next step is to quotient this finite candidate space by the explicit total-derivative, lower-order-EOM, algebraic, and linear-dependence relations already fixed in the repo.
- Status: Proven. That quotient step is the content of [`54-normal-form-basis-closure.md`](54-normal-form-basis-closure.md).
