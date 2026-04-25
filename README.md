# dynamic-chi-observable

Status: Imported from prior work. This worktree treats the previous collapse-theorem repository as background and attacks assumption A4 directly: no orbital-timescale internal state variable in the free-fall sector.

Status: Conjectural. The central question is whether a one-state internal-visibility variable `chi_A` can generate an observable that cannot be absorbed into a static finite-dimensional sensitivity-manifold EFT.

## Success Modes

Status: Counterexample candidate. Success mode 1 is an explicit dynamic loophole: `chi_A` produces phase lag, quadrature response, sidebands, or frequency-dependent transfer that is not reducible to a finite static sensitivity/Wilson basis.

Status: Conjectural. Success mode 2 is a sharp no-go boundary for the minimal model: the notes must state the exact limit or degeneracy in which the one-state `chi_A` model collapses back to the static description.

## MVP Model

Status: Counterexample candidate. The first model is the one-state relaxation system

```text
tau_chi * d chi_A / dt + chi_A = alpha * F(Y)
m_A(Y, chi_A) = m_A^(0) * [1 + c_Y F(Y) + c_chi chi_A]
```

Status: Conjectural. The default drive is the simplest deterministic scalar drive, `F(t) = F0 cos(Omega t)`, standing in for one selected free-fall-sector invariant or normalized invariant.

Status: Proven. The periodic response has transfer function

```text
H(Omega) = chi_hat / F_hat = alpha / (1 + i Omega tau_chi)
```

with an in-phase component and a quadrature component.

## Current MVP Classification

Status: Proven. In the strict instantaneous limit `Omega tau_chi = 0`, or when `alpha c_chi = 0`, the dynamic state produces no new observable and collapses to a static sensitivity redefinition.

Status: Proven. In the adiabatic regime `Omega tau_chi << 1`, the model is locally expandable as an infinite derivative series and is absorbable order-by-order by a derivative EFT to any fixed truncation order.

Status: Counterexample candidate. At finite `Omega tau_chi`, the full rational transfer function and its associated phase lag cannot be represented exactly by any finite static instantaneous sensitivity basis. If the admissible static EFT also allows arbitrary time-derivative Wilson coefficients, the one-frequency quadrature alone is not sufficient; novelty then requires frequency dependence across more than one drive frequency or the full pole structure.

Status: Counterexample candidate. The M2 frequency-sweep target is sharper: a real-coefficient degree-`N` finite derivative EFT can absorb at most `floor((N+1)/2)` positive frequency samples before the next sample exposes the relaxation pole. Even a complex-coefficient polynomial comparator fails at `N + 2` distinct samples if `alpha c_chi tau_chi != 0`.

Status: Counterexample candidate. The M3 forcing dictionary lowers that target to a concrete template: a normalized power-law invariant `(a/r)^p` on a small-eccentricity orbit produces `n` and `2n` harmonics through `O(e^2)`, and a calibrated linear observable projection preserves the same pole-bearing transfer relation.

Status: Counterexample candidate. The M4 projection audit fixes two concrete free-fall-style readouts: an acceleration-like channel with constant projection and a range-like channel with `Lambda_R(z)=Gamma/(kappa^2+z^2)`. The relaxation pole survives both unless the channel is blind, singular, or allowed an arbitrary frequency-local nuisance.

Status: Counterexample candidate. The triple shared-tau bridge gives a cleaner
carrier inventory: existing hierarchical-triple inner and outer GR carriers
provide two distinct samples, enough to test real shared-coefficient degree-1
and degree-2 derivative comparators when both carrier amplitudes are nonzero.

Status: Counterexample candidate. The triple GR carrier inventory adds the
outer-dipole combination carrier `|Omega_in-Omega_out|` as a conditional third
sample. In a nonresonant hierarchy with finite projection nuisance, this
pressures real degree-`N<=4` and complex degree-`N<=1` comparators.

Status: Counterexample candidate. The triple projection-nuisance realism gate
keeps the three-carrier bridge alive for calibrated or finite shared projection
models, but downgrades the result to conditional whenever the projection
parameters consume rank.

Status: Proven. If the timing projection is granted arbitrary complex nuisance
freedom independently at each carrier, the shared-`tau_chi` bridge collapses
point by point and is not runtime-motivated.

Status: Counterexample candidate. The physical projection-manifold gate shows
that a phase-locked outer-dipole combination projection has generic rank `5`
inside the six-real-dimensional three-carrier complex-amplitude vector. This
keeps dynamic-chi runtime-worthy if the timing model enforces the phase link.

Status: Proven. If the combination carrier amplitude and phase are floated
independently from the inner and outer carrier phases, the projection rank
becomes `6` and the bridge collapses as arbitrary per-carrier complex nuisance.

Status: Counterexample candidate. The named timing-model projection audit
places the standard public `Nutimo` triple timing core on the finite
state/geometry side of the gate, so dynamic-chi is conditionally
runtime-worthy for a configuration that preserves those links.

Status: Proven. The same source inspection finds an explicit harmonic
special case (`RN_PL`) with per-harmonic amplitudes and phases; enabling that
kind of nuisance on the target carriers collapses the shared-`tau_chi` bridge.

Status: Counterexample candidate. Request 10.7 promotes this to a narrow
external runtime-worthiness pilot contract: configuration closure first,
finite Jacobian rank gate second, and a minimal synthetic shared-`tau_chi`
injection only if the first two gates pass.

Status: Proven. A build-only or dependency-repair exercise does not count as
Request 10.7 progress. Runtime is scientifically motivated only if it tests
the named projection class and excludes carrier-local harmonic soak-up.

Status: Counterexample candidate. The external handoff packet fixes what a
pre-provisioned `Nutimo` environment must return: configuration manifest,
finite fitted-parameter Jacobian, carrier projection rank, dynamic-chi test
column, synthetic injection recovery, and a decision summary.

Status: Counterexample candidate. The M5 nonlinear sideband test adds the
minimal nonlinear drive/readout terms. These generate frequencies absent from
the linear input, including `Omega1+Omega2`, `|Omega1-Omega2|`, and the orbital
`3n` sideband from M3's `n,2n` forcing.

Status: Proven. For the linear readout above and a two-frequency linear drive, no sidebands are generated. Sidebands require a nonlinear drive/readout, a hereditary kernel with mixing, a second-order internal mode coupled nonlinearly, or a threshold/hysteretic extension.

## Layout

- `docs/model-definition.md`: one-state and deferred second-order model definitions.
- `docs/observable-targets.md`: criteria for new versus absorbed observables.
- `docs/adiabatic-limit.md`: low-frequency collapse derivation.
- `docs/nonadiabatic-regime.md`: finite-frequency response and observable classification.
- `docs/frequency-sweep-distinguishability.md`: finite-derivative EFT distinguishability theorem candidate.
- `docs/forcing-observable-dictionary.md`: concrete forcing templates and observable projection boundaries.
- `docs/projection-channel-audit.md`: concrete projection-channel audit for acceleration/range readouts.
- `docs/nonlinear-sideband-test.md`: minimal nonlinear drive/readout sideband test.
- `docs/failure-ledger.md`: exact failure modes and next escalation paths.
- `symbolic/chi_relaxation_response.py`: monochromatic response formulas.
- `symbolic/chi_two_frequency_response.py`: two-frequency superposition and sideband check.
- `symbolic/frequency_sweep_distinguishability.py`: polynomial interpolation and Taylor-residual checks.
- `symbolic/forcing_observable_dictionary.py`: harmonic forcing and projection checks.
- `symbolic/projection_channel_audit.py`: projection deconvolution and pole-survival checks.
- `symbolic/triple_shared_tau_bridge.py`: inner/outer carrier bridge to the shared-`tau_chi` frequency-sweep test.
- `symbolic/triple_gr_carrier_inventory.py`: three-carrier GR triple inventory and comparator-count audit.
- `symbolic/triple_projection_nuisance_gate.py`: projection-nuisance realism gate for the three-carrier bridge.
- `symbolic/triple_projection_manifold_gate.py`: physical projection-manifold rank and runtime-worthiness gate.
- `symbolic/named_timing_model_projection_audit.py`: named `Nutimo` source-inspection projection-class audit.
- `symbolic/nutimo_runtime_worthiness_pilot.py`: external `Nutimo` runtime-worthiness pilot contract.
- `symbolic/nutimo_external_handoff_packet.py`: external runtime handoff artifact checklist.
- `symbolic/nonlinear_sideband_test.py`: minimal nonlinear sideband generation checks.
- `symbolic/checks/test_symbolic.py`: deterministic symbolic checks.
- `counterexamples/README.md`: loophole framing and escalation list.

Status: Imported from prior work. Other inherited theorem files remain as background provenance only; the current worktree entry points are the dynamic-chi files listed above.

## Quick Start

Use Python `3.10+`.

```powershell
python -m pip install -r requirements.txt
python symbolic/checks/test_symbolic.py
```

If GNU Make is available:

```powershell
make symbolic-check
```
