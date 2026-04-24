# Request 8.2: Triple Prior-Breaking / Identifiability Audit

## Goal

Determine whether the leading `O(epsilon)` outer-gradient clock-sector
combination kernel in a hierarchical triple remains degenerate with the GR
outer-dipole carrier after external priors are introduced.

This is a theory-side identifiability gate. It is not a Nutimo/Tempo2 runtime
task, not source chasing, and not a full three-body numerical timing fit.

## Starting Point

The triple basis-rank audit found:

```text
K_zeta1 = (zeta_1 s_A) K_GR,dipole
K_zeta2 = (zeta_2 s_A^2) K_GR,dipole .
```

Thus the timing basis sees one effective clock amplitude,

```text
delta alpha_A = zeta_1 s_A + zeta_2 s_A^2,
```

multiplying the same outer-dipole carrier shape as the GR clock term. Without
external information, this is an exact carrier-amplitude degeneracy.

## Normalized Amplitude Model

Use the two-parameter amplitude basis

```text
theta = [A_GR_dipole, delta alpha_A] .
```

At leading order the data columns are identical:

```text
dR/dA_GR_dipole      = K_GR,dipole
dR/d(delta alpha_A)  = K_GR,dipole .
```

After normalizing the carrier vector to unit information, the timing-only Fisher
block is

```text
F_data = [[1, 1],
          [1, 1]] .
```

External priors are represented as a prior precision `p` on the GR carrier
amplitude:

```text
F = [[1 + p, 1],
     [1,     1]] .
```

The Schur-complement information for the clock amplitude is

```text
I_clock = p / (1 + p).
```

This makes the logic explicit:

- `p = 0`: exact degeneracy.
- small `p`: weak conditional identifiability.
- large `p`: practical identifiability in this normalized model.

## Prior Ladder

The machine table evaluates five bounded prior scenarios:

```text
none
mass_only
mass_plus_geometry_weak
outer_amplitude_moderate
idealized_strong_prior
```

Interpretation:

- `none`: no external prior; rank remains 1.
- `mass_only`: masses alone do not fix the projected carrier amplitude if
  geometry remains free; rank remains 1 in this model.
- `mass_plus_geometry_weak`: finite but weak carrier-amplitude prior; rank is
  formally restored but clock information is small.
- `outer_amplitude_moderate`: carrier-amplitude prior comparable to timing
  information; conditional identifiability appears.
- `idealized_strong_prior`: strong external carrier prior; motivates full
  runtime only if such a prior is realistic.

## Machine Outputs

The audit outputs are:

```text
outputs/tsv/triple_prior_identifiability.tsv
outputs/json/triple_prior_identifiability.json
```

The key columns are:

```text
prior_set
carrier_amplitude_degeneracy
clock_amplitude_identifiability
carrier_prior_precision
fisher_rank
marginal_clock_information
condition_number
verdict
notes
```

## Result

The current leading triple result is:

```text
without an external carrier-amplitude prior: no-go
with weak external priors: conditional but weak
with strong carrier priors: runtime work becomes scientifically motivated
```

This does not prove a final triple theorem. It shows that the next real
scientific gate is prior realism, not runtime availability.

## Consequence

Full triple timing runtime work should remain blocked unless a realistic prior
set can move the audit from `no-go` or `conditional-no-go` into the
`conditional` or `motivates-runtime` regime.

If realistic priors cannot constrain the GR outer-dipole carrier amplitude,
then the leading triple clock branch closes as carrier renormalization rather
than a new observable.
