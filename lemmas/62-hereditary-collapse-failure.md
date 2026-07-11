# Lemma 62: Hereditary Collapse Failure

## Statement

- Status: Proven. Admit the one-coordinate power-law hereditary model of [`60-hereditary-kernel-construction.md`](60-hereditary-kernel-construction.md) and keep the primitive-family envelope, the fixed-order local operator catalog, and the local scalar normal-form basis `Y^I` otherwise unchanged.
- Status: Proven. Then the first broken layer of the positive finite-family collapse theorem is not fixed-order operator-space finiteness.
- Status: Proven. The first broken layer is the reduction of the monopole response to a local function `m_A(Y(\tau))` of finitely many instantaneous normal-form coordinates.
- Status: Proven. Consequently Lemma 55 cannot even be started in its local form, and Lemma 56 fails secondarily because the monopole bookkeeping is no longer a finite local sensitivity jet.

## Proof

1. Status: Proven. The candidate local operator space remains the one already counted in [`53-fixed-order-operator-finiteness.md`](53-fixed-order-operator-finiteness.md). The kernel does not add a new primitive family or a new local scalar monomial.
2. Status: Proven. The reduced local normal-form quotient also remains intact as a local operator statement, because the hereditary kernel is introduced only after that local quotient has already been fixed.
3. Status: Proven. Let `\tau_0` be the readout time and let `Y_1` and `Y_2` be two histories with the same instantaneous value `Y_1(\tau_0)=Y_2(\tau_0)` but different past profiles. The explicit demo in [`../symbolic/hereditary_kernel_demo.py`](../symbolic/hereditary_kernel_demo.py) constructs such a pair.
4. Status: Proven. For the power-law kernel one finds

```math
\int_{-\infty}^{\tau_0} d\tau' \,
K_{\gamma}(\tau_0-\tau')\,Y_1(\tau')
\neq
\int_{-\infty}^{\tau_0} d\tau' \,
K_{\gamma}(\tau_0-\tau')\,Y_2(\tau').
```

5. Status: Proven. Therefore the monopole response cannot be written as a single-valued function of the instantaneous coordinate `Y(\tau_0)` alone, or of any finite Taylor jet built only from instantaneous local normal-form coordinates.
6. Status: Proven. The exact failing step relative to [`55-monopole-jet-collapse.md`](55-monopole-jet-collapse.md) is even earlier than analyticity: the prerequisite local form `m_A(Y)` is absent once `A3` is dropped.
7. Status: Proven. Since [`56-sensitivity-vs-wilson-split.md`](56-sensitivity-vs-wilson-split.md) assumes a local worldline action with monopole function `m_A(Y)` plus residual local operators, its sensitivity side also ceases to apply in the genuinely hereditary branch.

## Boundary Classification

- Status: Proven. The exponential kernel control does not give the same sharp failure, because it can be rewritten with a finite auxiliary local state and therefore moves the stress point into an `A4`-type extension instead.
- Status: Proven. The power-law kernel is therefore the first explicit pure hereditary failure of the local-collapse branch.
