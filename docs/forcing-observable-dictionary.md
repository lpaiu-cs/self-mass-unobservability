# Forcing And Observable Dictionary

Status: Counterexample candidate. M3 lowers the abstract transfer obstruction
to explicit forcing and readout templates without making empirical claims.

## Dictionary Scope

Status: Proven. The dynamic result to project is

```text
q_A(t) = delta m_A(t) / m_A^(0)
G(z) = q_hat(z) / F_hat(z) = c_Y + beta/(1 + tau_chi z),
beta = alpha c_chi,  z = i Omega.
```

Status: Conjectural. A concrete observable dictionary is a map from a physical
drive `F(Y(t))` to a measured linear channel

```text
O_hat(Omega_k) = Lambda(Omega_k) G(i Omega_k) F_hat(Omega_k).
```

Status: Counterexample candidate. The M3 dictionary is successful only when
`Lambda(Omega_k)` is known, calibrated, or restricted to its own finite
nuisance basis. If `Lambda(Omega_k)` is an arbitrary complex nuisance at every
frequency, the pole relation is absorbed and the dictionary fails.

## Concrete Forcing Class 1: Orbital Harmonics

Status: Counterexample candidate. Let the scalar drive be a normalized
free-fall-sector invariant with a power-law radial dependence along a bound
orbit:

```text
F(Y(t)) = Y(r(t)) / Y_*,
Y(r) / Y_* = (a/r)^p.
```

Status: Proven. In a small-eccentricity Kepler template, with mean motion `n`
and mean anomaly `M = n t`, the expansion through `O(e^2)` is

```text
(a/r)^p =
  1
  + p e cos(M)
  + [p(p-1)/4] e^2
  + [p(p+3)/4] e^2 cos(2M)
  + O(e^3).
```

Status: Counterexample candidate. Therefore a single eccentric system provides
at least two positive drive frequencies, `n` and `2n`, whenever

```text
e != 0,  p != 0,  p != -3.
```

Status: Proven. Those two frequencies are already enough to distinguish the
relaxation pole from a real-coefficient degree-`1` or degree-`2` derivative EFT,
provided the same projection factor does not introduce an arbitrary
frequency-by-frequency nuisance.

## Observable Projection

Status: Counterexample candidate. A minimal readout template is

```text
O(t) = O_static + Lambda q_A(t),
```

where `Lambda` is a known or separately modeled linear projection from the
body-dependent mass response into the chosen free-fall observable channel.

Status: Proven. For a harmonic drive `F_k cos(Omega_k t)`, the observable
contains

```text
O_k(t) =
  Lambda F_k [c_Y + beta/(1 + Omega_k^2 tau_chi^2)] cos(Omega_k t)
  + Lambda F_k [beta Omega_k tau_chi/(1 + Omega_k^2 tau_chi^2)] sin(Omega_k t).
```

Status: Proven. If `Lambda` is nonzero and frequency-independent over the
tested band, the calibrated ratio `O_hat_k/(Lambda F_hat_k)` recovers the same
transfer `G(i Omega_k)`.

Status: Counterexample candidate. If `Lambda(Omega)` is a known finite-order
projection polynomial, the distinguishability test must increase the derivative
order by the projection nuisance order and then apply the same frequency-sweep
counting.

## Concrete Forcing Class 2: Multi-System Sweep

Status: Conjectural. A set of systems with different orbital frequencies can
provide a frequency sweep only if the same `tau_chi` is shared or its
body-dependence is modeled by an independently specified relation.

Status: Proven. Without that shared-parameter assumption, each system can be
fit by its own effective derivative coefficients, so the pole relation is not
tested.

## Concrete Forcing Class 3: Independent Two-Tone Forcing

Status: Counterexample candidate. If one physical system has two independently
known drive frequencies `Omega_1` and `Omega_2`, the linear one-state model
still creates no sidebands, but it does provide two transfer samples.

Status: Proven. This is enough to beat a real-coefficient degree-`1` or
degree-`2` derivative EFT, but not enough to beat arbitrary higher-order
finite derivative EFTs.

## Failure Conditions

Status: Proven. The orbital-harmonic dictionary fails as a frequency-sweep
test if the orbit is circular in the chosen invariant, `e = 0`.

Status: Proven. The `O(e^2)` two-harmonic dictionary fails if the selected
invariant has no first harmonic, `p = 0`, or no second harmonic at this order,
`p = -3`.

Status: Proven. The observable projection fails if `Lambda = 0`.

Status: Proven. The pole test fails if `Lambda(Omega_k)` is allowed to be an
unconstrained complex nuisance at every sampled frequency.

Status: Counterexample candidate. The next required step after this dictionary
is to choose a specific free-fall observable channel and decide whether its
projection factor is known, finite-dimensional, or an arbitrary nuisance.
