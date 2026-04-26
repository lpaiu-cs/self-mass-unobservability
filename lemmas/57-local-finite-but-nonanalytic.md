# Lemma 57: Local And Finite Without Analyticity

## Statement

- Status: Proven. Keep `A2`, `A3`, `A4`, `A6`, `A7`, and `A8`, and keep the irreducible family-envelope closure already proved in [`../docs/irreducible-envelope-theorem.md`](../docs/irreducible-envelope-theorem.md).
- Status: Proven. Drop only analyticity `A5`.
- Status: Proven. Then the following statements remain intact:
  irreducible family-envelope closure,
  finite admitted primitive-family content,
  fixed-order candidate operator finiteness,
  and finite-dimensional scalar normal-form closure.
- Status: Proven. Therefore the nonanalytic counterexample program is aimed only at the monopole-jet layer and its bookkeeping descendants.

## Proof

1. Status: Proven. The irreducible family-envelope theorem uses only the parity-even nonspinning local tensor decomposition, trace-descendant absorption, and explicit theorem-domain exclusions. It does not use analyticity `A5`.
2. Status: Proven. Finite admitted primitive-family content is the explicit assumption `A8`. Dropping `A5` does not enlarge the primitive-family catalog.
3. Status: Proven. [`53-fixed-order-operator-finiteness.md`](53-fixed-order-operator-finiteness.md) uses local weight-spectrum finiteness below the cutoff, positive operator weights, and fixed `\Delta \le 4`. It does not use analyticity `A5`.
4. Status: Proven. [`54-normal-form-basis-closure.md`](54-normal-form-basis-closure.md) identifies the reduced scalar operator space as a finite-dimensional quotient of that already finite candidate space. It also does not use analyticity `A5`.
5. Status: Proven. Therefore locality and finite operator closure survive the isolated `A5` drop.

## Boundary

- Status: Proven. This lemma does not preserve the finite Taylor sensitivity jet.
- Status: Proven. It preserves only the local finite operator scaffold on which the nonanalytic monopole counterexample is built.
