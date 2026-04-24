# Request 8: Triple Outer-Potential Clock Dictionary

## Scope

This memo starts the hierarchical-triple clock dictionary after the binary
leading-order branch has closed as an absorbed/degenerate classification result.

Assumptions:

- free-fall sector = GR,
- propagation sector = GR,
- clock sector = decoupled EFT only,
- carrier = Newtonian nested-Kepler hierarchical triple,
- clock order = leading `O(c^-2)`,
- hierarchy order = through leading nontrivial
  `epsilon = a_in / a_out`.

This is not a full three-body timing code, not a Nutimo/Tempo2 runtime task, and
not a source-driven analysis.

## Notation

Inner binary:

- `A`: clock-carrying pulsar.
- `B`: inner companion.
- `m_A`, `m_B`: inner masses.
- `M_in = m_A + m_B`.
- `r_in = a_in (1 - e_in cos E_in)`.
- `n_in`: inner mean motion.
- `rho_A = beta_A r_in`, with `beta_A = m_B / M_in`, is the pulsar offset
  from the inner-binary barycenter along the inner separation vector.

Outer orbit:

- `C`: outer body.
- `m_C`: outer-body mass.
- `R_out = a_out (1 - e_out cos E_out)`.
- `n_out`: outer mean motion.
- `epsilon = a_in / a_out`.
- `Phi`: instantaneous angle between the inner separation direction and the
  outer separation direction, so `cos Phi = rhat_in dot Rhat_out`.

Clock-sector coefficient:

```text
alpha_A = 1 + zeta_1 s_A + zeta_2 s_A^2 .
```

## Clock Rate

The pulsar clock rate is written as

```text
d tau_A/dt
= 1 - v_A^2/(2 c^2)
  - alpha_A (U_inner + U_outer,A) / c^2 .
```

The current triple outer-potential dictionary focuses on the scalar potentials:

```text
U_inner = G m_B / r_in .
```

For the outer body, expand about the inner-binary barycenter:

```text
U_outer,A
= G m_C / |R_out - rho_A|
= G m_C/R_out
 + G m_C (rho_A dot Rhat_out)/R_out^2
 + O(rho_A^2/R_out^3).
```

Using `rho_A = beta_A r_in rhat_in`,

```text
U_outer,A
= G m_C/R_out
 + G m_C beta_A r_in cos Phi/R_out^2
 + O(epsilon^2).
```

The first term is the outer monopole at the inner-binary barycenter. The second
term is the leading nontrivial `O(epsilon)` outer-gradient/dipole clock term.

## Projection

The inner monopole projects exactly as the closed binary result:

```text
dt/dE_in = (1 - e_in cos E_in)/n_in .
```

After secular spin-frequency absorption, the nonzero periodic inner piece is
proportional to `sin E_in`; it is an inner Einstein-delay gamma shape.

The outer monopole similarly uses

```text
dt/dE_out = (1 - e_out cos E_out)/n_out .
```

After secular spin-frequency absorption, the nonzero periodic outer-monopole
piece is proportional to `sin E_out`; it is an outer Einstein/redshift delay
shape.

The leading outer-gradient/dipole term has rate kernel

```text
delta dot_tau_dipole
= -G m_C beta_A epsilon factor
   [(1 - e_in cos E_in)/(1 - e_out cos E_out)^2]
   cos Phi
   /(a_out c^2),
```

where `factor` is one of

```text
1, zeta_1 s_A, zeta_2 s_A^2 .
```

Its timing projection is the formal residual

```text
Delta_dipole = integral dt delta dot_tau_dipole .
```

This is a combination structure because it depends on both the inner phase and
the outer phase through `r_in`, `R_out`, and `Phi`.

## Classification

The machine table uses the required harmonic labels:

- `secular`: coordinate-time trend absorbed by spin phase/frequency,
- `inner-only`: inner `sin E_in` Einstein-delay basis,
- `outer-only`: outer `sin E_out` Einstein/redshift-delay basis,
- `combination`: inner-outer phase combination from the leading outer gradient.

Classification at this bounded step:

- Inner monopole secular pieces are `absorbed`.
- Inner monopole periodic pieces are `degenerate` with inner Einstein delay
  gamma.
- Outer monopole secular pieces are `absorbed`.
- Outer monopole periodic pieces are `degenerate` with an outer
  Einstein/redshift-delay parameter.
- The baseline outer dipole is a GR carrier combination term and is marked
  `degenerate` with the GR triple clock carrier or fitted triple parameters.
- The `zeta_1 s_A` and `zeta_2 s_A^2` outer-dipole pieces are marked
  `observable` at the dictionary level because they are not pure spin,
  inner-only gamma, or outer-only gamma bases.

The last classification is not a final detectability claim. It flags the first
clock-sector structure that survives the binary gamma-only collapse and should
later be tested against full triple parameter covariance.

## Machine Outputs

The generated outputs are:

```text
outputs/tsv/triple_clock_dictionary_leading.tsv
outputs/json/triple_clock_dictionary_leading.json
```

They contain the notation, projection rules, and one classified row per
leading-order contribution.

## Stop-Rule Status

The leading triple dictionary does not collapse completely into spin,
inner-only gamma, and outer-only gamma nuisance structure. The outer-gradient
combination rows survive as clock-sector observable candidates at the dictionary
level.

Therefore this is not a no-go closure yet, and it is also not a final theorem.
The next bounded task should be a covariance-aware or basis-rank check of the
combination kernel, not source chasing or runtime provisioning.
