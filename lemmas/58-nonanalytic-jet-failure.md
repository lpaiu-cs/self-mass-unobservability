# Lemma 58: nonanalyticity and exact reconstruction

This title's historical filename is retained. The finite-jet failure previously attributed to a smooth-flat function is corrected here.

Status: Proven. Define phi(Y)=exp(-1/Y^2) for Y>0 and zero otherwise. Every positive-side derivative is a polynomial in 1/Y times the exponential. Each derivative tends to zero at the origin, so phi is C-infinity and its Taylor jet at every finite order is zero.

Status: Proven. For every finite n, `phi(Y)/|Y|^n -> 0` at zero. This follows from `x^(n/2) exp(-x) -> 0` as x=1/Y^2 tends to infinity. Thus every finite-order approximation with a polynomial-order remainder remains valid. In particular, phi=0+O(Y^5).

Status: Proven. Since phi is positive for Y>0, its all-zero Taylor series does not recover the exact germ. This is a counterexample to exact analytic reconstruction, not to [Lemma 55](55-monopole-jet-collapse.md).

Status: Proven. For the threshold function sqrt(Y) Theta(Y), a polynomial matching the zero value is either zero or begins at O(Y). The error therefore retains a leading Y^(1/2) term and cannot be O(Y^5). This distinct example fails finite-order differentiability.

Status: Proven. Neither example enlarges the polynomial operator catalog. Its physical use would require a specified domain and coupling, and does not by itself establish a measurable loophole.
