# Request 4: MLRS Sun-Driven State Gate

## Scope

This is the last bounded gate before declaring that the remaining `MLRS`
question is purely physical rather than architectural.

It asks:

- can a more physical, Sun-driven perturbation class still fit through the
  local `jjreadnp` state bridge,
- without widening into `PLEPH` edits, ephemeris regeneration, or broad batch
  stack surgery.

It still does **not** implement the full weak-field `delta_SEP` dynamics.

## Implementation

The gate is implemented in
[request4_llr_mlrs_sundriven_gate.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_sundriven_gate.py:1).

The edit scope remains local:

- [jjreadnp.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f:1)
- [smu_probe_adjust_states.f](/Users/lpaiu/vs/lab/self-mass-unobservability/data/request4_llr/mlrs_handshake_lab/src/llr_npt/smu_probe_adjust_states.f:1)

The only new architectural move is:

- one extra `PLEPH` fetch for the Sun state inside `jjreadnp`,
- followed by a local Moon displacement along the instantaneous
  Earth-to-Sun direction.

No edits were made to:

- `PLEPH`,
- `jeulpkg.f`,
- the batch scripts,
- or the ephemeris files themselves.

## Probe Configuration

The probe used the first bundled public sample:

- `s25y11d016t0231#103`

and four symmetric amplitudes:

- `-0.010 m`
- `-0.005 m`
- `+0.005 m`
- `+0.010 m`

The synthetic state class is:

- Moon displacement along the instantaneous Earth-to-Sun line.

This is more physical than the earlier Earth-to-Moon radial probe because it is
explicitly Sun-driven, but it is still only a bounded surrogate for the real
`delta_SEP` class.

## Outputs

- summary:
  [request4_llr_mlrs_sundriven_gate_summary.json](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_sundriven_gate_summary.json:1)
- figure:
  [request4_llr_mlrs_sundriven_gate_summary.svg](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_sundriven_gate_summary.svg)
- runner:
  [request4_llr_mlrs_sundriven_gate.py](/Users/lpaiu/vs/lab/self-mass-unobservability/request4_llr_mlrs_sundriven_gate.py:1)

## Result

The Sun-driven local state class still stays bounded.

For all four amplitudes:

- the pipeline runs end-to-end,
- the same `1456 / 1761` `93` residual lines respond,
- `frfin` changes only `12` lines,
- `frd` remains an exact match,
- and `npt` changes only one `11` record.

The response is also close to linear and odd-symmetric:

- `Δ93 / amplitude` is `4.429` to `4.440 ns / m`,
- `ΔPmM / amplitude` is exactly `4440 ps / m`,
- the odd-symmetry sums stay small:
  `|frrec sum| <= 1.03e-5 ns`,
  `npt TOF sum = 0`,
  `npt PmM sum = 0.1 ps`.

Numerically:

- `-0.010 m` gives mean `Δ93 = -4.43908e-2 ns`, `ΔPmM = -44.3 ps`
- `-0.005 m` gives mean `Δ93 = -2.21470e-2 ns`, `ΔPmM = -22.1 ps`
- `+0.005 m` gives mean `Δ93 = +2.21463e-2 ns`, `ΔPmM = +22.2 ps`
- `+0.010 m` gives mean `Δ93 = +4.44011e-2 ns`, `ΔPmM = +44.4 ps`

## Interpretation

This is the key gate result:

- the `MLRS` branch has now survived not just a recalc-layer tweak,
- and not just an arbitrary local state perturbation,
- but also a bounded **Sun-driven** local state class.

So the software-architecture stop rule stays cleared.

What remains open is not software scope but physics adequacy:

- does the real weak-field `delta_SEP` effect reduce cleanly to a local
  Sun-driven state remapping at this bridge,
- or does a faithful implementation require ephemeris-level dynamics that
  should be handed off to a more mature estimator.

That is now the remaining decision.
