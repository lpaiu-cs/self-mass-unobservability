# Request 4: MLRS Interface Probe

## Scope

This stage is the first follow-on after the `MLRS` handshake passed sample
replay.

It does **not** implement the real `delta_SEP` physics yet.
It asks two narrower questions:

1. is the residual mismatch against the bundled reference only a formatting
   issue, or does it already appear inside the recalc state/residual layer,
2. can a strictly local correction be injected into the `recalc` layer and
   propagated through the existing batch stack without touching the downstream
   filtering executables.

## Implementation

The probe lives in
[request4_llr_mlrs_interface_probe.py](request4_llr_mlrs_interface_probe.py#L1).

It performs two runs on the first bundled lunar sample:

- `baseline`: no synthetic correction,
- `offset`: a local `SMU_OMC_OFFSET_NS=0.005` correction injected inside
  the recalc/normalpoint executable.

The local hook is intentionally minimal and lives in the lab copy only:

- [np_crd.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/np_crd.f#L441)
  now calls a bounded helper immediately after computing `OMC`,
- [smu_probe_adjust_omc.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/smu_probe_adjust_omc.f#L1)
  reads `SMU_OMC_OFFSET_NS` and, by default, is a no-op.

This is a software-architecture probe, not the final EFT.

## Outputs

- probe summary:
  [request4_llr_mlrs_interface_probe_summary.json](request4_llr_mlrs_interface_probe_summary.json#L1)
- probe figure:
  [request4_llr_mlrs_interface_probe_summary.svg](request4_llr_mlrs_interface_probe_summary.svg)
- probe script:
  [request4_llr_mlrs_interface_probe.py](request4_llr_mlrs_interface_probe.py#L1)

## Result 1: The Existing Drift Is Not Pure Formatting

The baseline replay shows that the current `MLRS` lab still differs slightly
from the bundled reference, but the location of the difference is now clearer.

For the first sample:

- `frrec` differs from reference in `79` lines,
- every one of those differences is in `93` residual records,
- the maximum observed numeric drift is only about `0.001 ns`,
- `frfin` then inherits only `1` residual-line difference,
- `frd` remains an exact match,
- and the final `npt` difference is just one line in the higher-moment /
  return-rate fields.

So the current mismatch is **not** "just text formatting at the very end".
It appears first in the recalc residual stream, but only at the
picosecond-scale residual level and largely collapses before the final
stripped full-rate output.

## Result 2: A Bounded Local Interface Exists

The synthetic probe then adds `+0.005 ns` through the local `OMC` hook.

That probe shows:

- the full batch chain still runs end-to-end,
- `frrec` changes in `1456 / 1761` residual (`93`) lines, with observed
  quantized deltas clustered around `+0.005 ns`,
- `frfin` changes in only `12` residual lines,
- `frd` stays **exactly unchanged**,
- `npt` changes by one line, with the normal-point time-of-flight shifting by
  `5e-12 s` and `PmM` shifting by `+5.0 ps`.

That is the critical architectural result.

It means a correction inserted at the recalc prediction/residual layer can
propagate through normalpoint construction **without** requiring changes to:

- `Poisson_crd`,
- `lun_auto_crd`,
- `lun_sdqi_crd`,
- or the downstream stripped full-rate layout.

## Interface Reading

After this probe, the bounded insertion map is narrower:

1. nearest local correction seam:
   [np_crd.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/np_crd.f#L441)
   at the `OMC` computation / `write_93` path
2. deeper state-level seam:
   [jjreadnp.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/jjreadnp.f#L30)
   through the `PLEPH` state bridge
3. geometry/prediction consumer:
   [jeulpkg.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/jeulpkg.f#L341)
   and
   [jeulpkg.f](data/request4_llr/mlrs_handshake_lab/src/llr_npt/jeulpkg.f#L1216)

So the next real question is no longer "can MLRS be made to run?".
It is:

- can the actual `delta_SEP` correction be written as a bounded prediction /
  residual perturbation at one of these seams,
- or would that require broad ephemeris/batch-stack surgery.

## Decision

This probe upgrades the `MLRS` status again.

- It remains **not** a weak-field posterior.
- But it is now more than a sample-replay candidate.
- It exposes a **bounded recalc-layer correction interface**.

That is enough to justify one more short step on the MLRS branch.

The stop rule also becomes sharper:

- if the real `delta_SEP` remapping can be expressed through this bounded
  interface, MLRS can carry the weak-field hand-off,
- if it instead expands into broad `PLEPH` / ephemeris / batch-stack surgery,
  this branch should stop here and Request 4 should fall through to an
  external mature estimator.
