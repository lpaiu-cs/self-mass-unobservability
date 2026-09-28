# Lemma 55: finite-order response expansion

Corrected 2026-09-09; the earlier assertion that analyticity is necessary for a finite jet is withdrawn.

Status: Proven. Let f be `C^(D+1)` in a neighborhood of zero in finitely many arguments x_i. Suppose `x_i=epsilon^(w_i) xbar_i` with bounded xbar and positive integer weights w_i. Then

```math
f(x)=\sum_{\sum_i n_iw_i\le D}\frac{\partial^{\mathbf n}f(0)}{\mathbf n!}x^{\mathbf n}
+O(\epsilon^{D+1}).
```

Status: Proven. Taylor expansion through total degree D has remainder `O(norm(x)^(D+1))`. Since norm(x)=O(epsilon), this is the stated order. Each discarded weighted monomial has integer weight at least D+1, so it also belongs to the remainder. There are finitely many retained multi-indices. At D=4, C5 suffices. The constant term appears exactly once in this sum.

Status: Proven. If the arguments are redundant invariant representatives, such as I2 and I2^2, their algebraic relations must be reduced before counting independent response coefficients. Calling the coefficients sensitivities does not make them different from the corresponding EFT action coefficients without an explicit matching convention.

Status: Proven. Analyticity A5 suffices but is unnecessary for the finite-order statement. A smooth-flat response has a valid zero jet with a beyond-all-orders remainder. Analyticity instead addresses exact reconstruction by the convergent Taylor series. See [Lemma 58](58-nonanalytic-jet-failure.md).

Status: Proven. The threshold `sqrt(Y) Theta(Y)` lacks the required regularity at zero and is not approximable by an ordinary polynomial with `O(Y^5)` error there. This is an expansion failure at the threshold, not a change in the specified polynomial catalog's size.
