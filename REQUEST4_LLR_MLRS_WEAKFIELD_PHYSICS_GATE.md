# Request 4: MLRS Weak-Field Physics Gate

## Scope

This is the first gate that no longer uses an arbitrary synthetic displacement
class.

Instead, it takes the weak-field analytic structure already fixed in
[REQUEST3_LLR_MOCK.md](/Users/lpaiu/vs/lab/self-mass-unobservability/REQUEST3_LLR_MOCK.md:1)
and asks whether an externally computed, weak-field-inspired `delta_SEP` state
correction can still be injected through the local `jjreadnp` seam.

It is still **not** a full dynamical rederivation of Nordtvedt physics inside
`MLRS`. It is a bounded observation-operator gate.

## External Weak-Field Model

Using the Request 3 weak-field benchmark,

```math
\Delta_{\rm SEP}
=
\sigma_1 (s_E-s_M) + \sigma_2 (s_E^2-s_M^2),
```

with

```math
\delta r_{\cos D} \approx A_N \Delta_{\rm SEP},
\qquad
A_N = \frac{13.1\ \mathrm{m}}{s_E-s_M},
```

this gate constructs an externally parameterized relative state correction

```math
\delta \mathbf r_{\rm rel}(t) = A_N \Delta_{\rm SEP}\,\hat{\mathbf u}_{ES}(t),
\qquad
\delta \mathbf v_{\rm rel}(t) = A_N \Delta_{\rm SEP}\,\frac{d\hat{\mathbf u}_{ES}}{dt},
```

where `u_ES(t)` is the instantaneous Earth-to-Sun unit vector from the JPL
state already present at the `jjreadnp` bridge.

The correction is then split between Earth and Moon to preserve the
Earth-Moon barycenter.

## Implementation

The gate runner is
[request4_llr_mlrs_weakfield_physics_gate.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_weakfield_physics_gate.py:1).

The edit scope stays local:

- [jjreadnp.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f:1)
- [smu_probe_adjust_states.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/smu_probe_adjust_states.f:1)

No edits were made to:

- `PLEPH`,
- `jeulpkg.f`,
- the batch scripts,
- or the ephemeris files.

## Probe Configuration

Two bundled public lunar samples were used:

- `s25y11d016t0231#103`
- `s25y11d042t0234#103`

The gate used `sigma_2 = 0` and four symmetric `sigma_1` amplitudes:

- `-1.0e-3`
- `-5.0e-4`
- `+5.0e-4`
- `+1.0e-3`

These map to relative displacement amplitudes

- `-1.31 cm`
- `-6.55 mm`
- `+6.55 mm`
- `+1.31 cm`

through the Request 3 analytic coefficient.

## Outputs

- summary:
  [request4_llr_mlrs_weakfield_physics_gate_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_weakfield_physics_gate_summary.json:1)
- figure:
  [request4_llr_mlrs_weakfield_physics_gate_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_weakfield_physics_gate_summary.svg)
- runner:
  [request4_llr_mlrs_weakfield_physics_gate.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_weakfield_physics_gate.py:1)

## Result

The local seam survives the weak-field-inspired correction.

For `s25y11d016t0231#103`:

- all four runs keep `1456 / 1761` reacting `93` residual lines,
- `frfin` stays at `12` differing lines,
- `frd` remains exact match,
- `npt` stays a one-line drift,
- and the response is nearly perfectly linear and odd-symmetric:
  `Δ93 / sigma1 ≈ 58.1 ns`,
  `ΔPmM / sigma1 ≈ 5.81e4 ps`.

For `s25y11d042t0234#103`:

- all four runs keep `256 / 714` reacting `93` residual lines,
- `frfin` stays at `9` differing lines,
- `frd` remains exact match,
- `npt` still stays a one-line drift,
- and the response is weaker but still approximately linear, with a sign flip
  relative to the first sample:
  `Δ93 / sigma1 ≈ -3.8` to `-4.1 ns`.

So the same weak-field correction class propagates through both bundled
sessions without widening the code-change surface.

## Interpretation

This is the strongest Request 4 result so far.

What is now justified:

- `MLRS` is no longer just a live software hand-off candidate,
- it is a plausible **weak-field observation-operator** candidate,
- because an externally computed weak-field-style `delta_SEP` state correction
  can be injected through the local seam and survive the downstream
  residual/normal-point pipeline.

What is still **not** justified:

- that `MLRS` internally generates the correct SEP dynamics,
- or that this local state-remapping approximation is already a complete
  physical model of real LLR Nordtvedt physics.

So the remaining question is now narrow and explicit:

- is this local observation-operator approximation physically adequate enough
  for a Request 4 weak-field remapping study,
- or does a faithful implementation still require ephemeris-level dynamics that
  should be handed off to a more mature external estimator.
