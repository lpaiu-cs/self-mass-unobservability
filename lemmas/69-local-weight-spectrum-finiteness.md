# Lemma 69: Local Weight-Spectrum Finiteness

## Statement

- Status: Proven. Fix the current theorem domain, a cutoff `\Delta_{\max} = 4`, and the positive intrinsic-weight counting rule of [`../docs/power-counting.md`](../docs/power-counting.md).
- Status: Proven. Replace the old finite-total-family hypothesis by the weaker condition:
  only finitely many primitive-family species have intrinsic weight `w \le \Delta_{\max}`.
- Status: Proven. Then the candidate local parity-even free-fall scalar operator space at `\Delta \le \Delta_{\max}` is finite before reduction.

## Proof

1. Status: Proven. Any primitive-family species with intrinsic weight `w > \Delta_{\max}` cannot appear in a scalar operator of total weight `\Delta \le \Delta_{\max}`, because all derivative decorations have nonnegative weight.
2. Status: Proven. Therefore only the finite subcatalog

```math
{\cal F}_{\le \Delta_{\max}}
:=
\{X_a \mid w(X_a) \le \Delta_{\max}\}
```

can contribute at the chosen order.
3. Status: Proven. For each such species `X_a`, the decorated descendants `D_\tau^p \nabla^q X_a` that survive at `\Delta \le \Delta_{\max}` satisfy `p+q \le \Delta_{\max} - w(X_a)`, so each contributing family supplies only finitely many decorated blocks.
4. Status: Proven. Because `{\cal F}_{\le \Delta_{\max}}` is finite, the entire contributing decorated-block set is finite.
5. Status: Proven. Every retained scalar operator is a parity-even complete contraction of finitely many blocks of positive weight whose total weight is at most `\Delta_{\max}`.
6. Status: Proven. Positive block weights bound the number of factors in any such operator, and finite-rank tensor multisets admit only finitely many complete contraction patterns.
7. Status: Proven. Therefore only finitely many scalar contraction classes can appear at the chosen cutoff.

## Bridge Consequences

- Status: Proven. Lemma 53 survives with this weaker hypothesis after replacing “finite total admitted family content” by “local weight-spectrum finiteness below `\Delta_{\max}`.”
- Status: Proven. Lemma 54 still goes through, because the reduced scalar operator space is still a quotient of a finite candidate space.
- Status: Proven. Lemma 55 still goes through, because it only needs a finite scalar normal-form basis together with analyticity `A5`.
- Status: Proven. Lemma 56 still goes through, because it only needs the finite reduced scalar basis, locality `A3`, and the separate finite residual higher-multipole operator sector.

## Verdict

- Status: Proven. Local weight-spectrum finiteness below the chosen cutoff is sufficient for the fixed-order operator-finiteness layer and for the rest of the positive collapse bridge inside the current theorem domain.
