# Request 8: Binary Theorem Closure

## Status

This memo closes the current binary leading-order clock branch as a conditional
classification result.

It does not claim a final all-binary or all-sector theorem. The statement below
holds only under the assumptions listed here.

## Assumptions

- free-fall sector = GR
- propagation sector = GR
- carrier orbit = Newtonian Kepler binary
- clock sector = decoupled EFT only
- clock order = leading `O(c^-2)`
- eccentricity treatment = exact in `e` for the Newtonian Kepler projector
- timing-basis comparison = DD-style Einstein-delay basis plus standard spin
  phase/frequency nuisance structure

No triple systems, source chasing, external estimators, or legacy timing
runtimes enter this closure.

## Theorem Candidate

Under GR free-fall and GR propagation, a decoupled clock-sector EFT at leading
`O(c^-2)` on a Newtonian Kepler binary produces no distinguishable periodic
timing basis beyond the DD Einstein-delay `sin E` basis and secular
spin-redshift renormalization.

Equivalently, the decoupled clock terms remain symbolically present, but after
projection into timing residuals they only renormalize:

- spin phase / spin frequency through secular terms,
- Einstein-delay gamma through `sin E` periodic terms.

## Derivation Closure

The active clock EFT is

```text
d tau_A / dt
= 1 - v_A^2/(2 c^2)
  - (1 + zeta_1 s_A + zeta_2 s_A^2) U_ext / c^2 .
```

For the Newtonian binary carrier,

```text
r = a (1 - e cos E)
U_ext = G m_B / r
v_A^2 = (m_B/M)^2 G M (2/r - 1/a)
dt/dE = (1 - e cos E) / n .
```

The exact-in-`e` projector rewrites each contribution as

```text
d tau_i/dE = (d tau_i/dt) (dt/dE),
Delta_i(E) = integral (d tau_i/dE) dE .
```

The machine table in
`outputs/tsv/binary_pk_dictionary_exact_e.tsv` gives the complete
classification at this order.

## DD Absorption Map

The timing model already fits a spin phase and spin frequency. Therefore any
residual piece of the form

```text
A (t - T0) + B
```

is not an independent binary timing observable at this order. It is absorbed by
spin phase and spin-frequency redefinition. In the exact projector,

```text
t - T0 = (E - e sin E) / n,
```

so secular absorption must be performed against coordinate time, not against
`E` alone.

The DD Einstein delay has the periodic basis

```text
Delta_E = gamma sin E .
```

Any residual piece of the form

```text
C sin E
```

is therefore not a new timing basis. It is absorbed as

```text
gamma -> gamma + C
```

up to the timing-residual sign convention used by the final estimator.

## Exact-in-e Periodic Coefficients

After removing the secular spin-redshift trend, the exact-in-`e` periodic
pieces are:

```text
baseline kinetic:   -G e m_B^2 / (M a c^2 n) sin E
baseline potential: -G e m_B / (a c^2 n) sin E
zeta_1 term:        -G e m_B s_A zeta_1 / (a c^2 n) sin E
zeta_2 term:        -G e m_B s_A^2 zeta_2 / (a c^2 n) sin E
```

The baseline kinetic and potential pieces are GR carrier terms. The decoupled
clock-sector gamma shift candidate is therefore

```text
delta gamma_clock
= -G e m_B (zeta_1 s_A + zeta_2 s_A^2) / (a c^2 n)
```

again modulo the residual sign convention.

No independent periodic basis remains after this projection.

## Counterexample Search Status

Within the assumptions of this memo, the binary search did not find a
counterexample to the gamma-only classification. The exact-in-`e` projection
also removes the main concern that the earlier `O(e^1)` result was only a
truncation artifact.

The current binary branch therefore closes as:

```text
observable_status = absorbed or degenerate
independent leading-order binary clock observable = none found
```

## Limits

This closure does not cover:

- post-Newtonian corrections to the carrier orbit,
- non-GR free-fall or propagation sectors,
- higher clock orders beyond leading `O(c^-2)`,
- spin-orbit, Shapiro, aberration, or propagation-sector modifications,
- hierarchical triples or outer-potential clock terms,
- source-specific estimator covariance.

Those are separate questions. They should not be folded back into the binary
leading-order branch.

## Next Boundary

The next scientifically clean question is not more binary source searching.
It is:

```text
Does a hierarchical triple outer potential open a clock-sector timing basis
that is absent in the closed binary leading-order projector?
```

That question should begin only after this binary closure memo and the exact
machine-readable dictionary are treated as the binary branch baseline.
