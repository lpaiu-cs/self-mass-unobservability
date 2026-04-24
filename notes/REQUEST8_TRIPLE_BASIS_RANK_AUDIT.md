# Request 8: Triple Basis-Rank Audit

## Goal

Determine whether the leading `O(epsilon)` outer-gradient clock-sector
combination kernel in a hierarchical triple lies outside the span of standard
triple timing nuisance structure.

This memo is a basis-rank gate, not a covariance fit and not a numerical timing
runtime.

## Scope

Assumptions:

- free-fall sector = GR,
- propagation sector = GR,
- clock sector = decoupled EFT only,
- carrier = Newtonian nested-Kepler hierarchical triple,
- clock order = leading `O(c^-2)`,
- hierarchy order = leading `O(epsilon = a_in/a_out)`.

This pass does not use Nutimo, Tempo2, external estimators, source chasing, or
full three-body numerical timing.

## Kernel Under Audit

The previous triple dictionary found that the inner and outer monopoles collapse
into spin, inner-only Einstein-delay, or outer-only redshift-delay structure.
The only remaining candidate was the outer-gradient/dipole term

```text
delta dot_tau_dipole
= -G m_C beta_A epsilon factor
   [(1 - e_in cos E_in)/(1 - e_out cos E_out)^2]
   cos Phi
   /(a_out c^2),
```

where

```text
factor in {1, zeta_1 s_A, zeta_2 s_A^2}.
```

The baseline `factor = 1` term is the GR outer-dipole clock carrier. The
decoupled clock terms are scalar multiples of the same carrier shape.

## Harmonic Expansion

For the first rank audit, use the coplanar phase model

```text
cos Phi ~= cos(lambda_in - lambda_out)
```

and expand the eccentricity modulation through first order in `e_in` and
`e_out`:

```text
[(1 - e_in cos lambda_in)/(1 - e_out cos lambda_out)^2]
cos(lambda_in - lambda_out)
```

which gives

```text
cos(lambda_in - lambda_out)
- (e_in/2) cos(lambda_out)
- (e_in/2) cos(2 lambda_in - lambda_out)
+ e_out cos(lambda_in)
+ e_out cos(lambda_in - 2 lambda_out).
```

The true combination frequencies are therefore:

```text
lambda_in - lambda_out
2 lambda_in - lambda_out
lambda_in - 2 lambda_out
```

The `O(e)` expansion also carries inner-only and outer-only sidebands, but those
sit in the already-known nuisance span.

## Nuisance Spans

The minimal nuisance span includes:

- secular spin phase/frequency,
- inner-only Einstein-delay basis,
- outer-only redshift-delay basis.

Against this minimal span, the outer-dipole kernel is a genuine combination
vector.

The expanded standard triple span includes:

- the minimal span,
- GR outer-gradient clock carrier amplitude,
- GR outer-gradient clock carrier phase quadratures,
- inner/outer orbital phase and amplitude redefinitions.

Against this expanded span, the decoupled outer-dipole kernels do not add rank:

```text
K_zeta1 = (zeta_1 s_A) K_GR_dipole
K_zeta2 = (zeta_2 s_A^2) K_GR_dipole .
```

Thus both EFT pieces are basis-degenerate with an effective rescaling of the GR
outer-dipole clock carrier.

## Classification Update

The previous triple dictionary correctly marked the `zeta_1` and `zeta_2`
outer-dipole terms as `observable` at dictionary level: they are not pure spin,
inner-only gamma, or outer-only gamma terms.

This rank audit uses a larger nuisance span. Under that stricter span, the same
rows are reclassified as:

```text
observable_status = degenerate
absorption_target = GR outer-gradient clock amplitude / geometry nuisance
```

The two EFT coefficients are also mutually degenerate in a single-clock timing
basis audit:

```text
delta alpha_A = zeta_1 s_A + zeta_2 s_A^2 .
```

At this level, timing sees one effective outer-dipole clock amplitude, not two
separately identifiable EFT coefficients.

## Machine Outputs

The audit table is written to:

```text
outputs/tsv/triple_basis_rank_audit.tsv
outputs/json/triple_basis_rank_audit.json
```

The table records:

- the explicit harmonic vector,
- the minimal-span result,
- the expanded-span result,
- the updated observable classification.

## Result

The leading triple outer-gradient clock kernel is an algebraic escape hatch from
the binary gamma-only structure only if the nuisance span excludes the GR
outer-dipole carrier amplitude.

Once the standard GR triple carrier/geometry nuisance span is included, the
decoupled clock-sector kernel is a scalar rescaling of an existing carrier
shape. Therefore this gate produces a leading triple no-go candidate:

```text
new combination harmonic exists,
but it is rank-degenerate with the GR outer-dipole carrier at this order.
```

This is not a final theorem. It is a stricter basis-rank result that should
prevent premature runtime or source-driven work. A later closure would still
need to state whether external priors on masses/geometry can break the carrier
amplitude degeneracy.
