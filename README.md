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

Status: Proven. For the linear readout above and a two-frequency linear drive, no sidebands are generated. Sidebands require a nonlinear drive/readout, a hereditary kernel with mixing, a second-order internal mode coupled nonlinearly, or a threshold/hysteretic extension.

## Layout

- `docs/model-definition.md`: one-state and deferred second-order model definitions.
- `docs/observable-targets.md`: criteria for new versus absorbed observables.
- `docs/adiabatic-limit.md`: low-frequency collapse derivation.
- `docs/nonadiabatic-regime.md`: finite-frequency response and observable classification.
- `docs/frequency-sweep-distinguishability.md`: finite-derivative EFT distinguishability theorem candidate.
- `docs/forcing-observable-dictionary.md`: concrete forcing templates and observable projection boundaries.
- `docs/projection-channel-audit.md`: concrete projection-channel audit for acceleration/range readouts.
- `docs/failure-ledger.md`: exact failure modes and next escalation paths.
- `symbolic/chi_relaxation_response.py`: monochromatic response formulas.
- `symbolic/chi_two_frequency_response.py`: two-frequency superposition and sideband check.
- `symbolic/frequency_sweep_distinguishability.py`: polynomial interpolation and Taylor-residual checks.
- `symbolic/forcing_observable_dictionary.py`: harmonic forcing and projection checks.
- `symbolic/projection_channel_audit.py`: projection deconvolution and pole-survival checks.
- `symbolic/triple_shared_tau_bridge.py`: inner/outer carrier bridge to the shared-`tau_chi` frequency-sweep test.
- `symbolic/triple_gr_carrier_inventory.py`: three-carrier GR triple inventory and comparator-count audit.
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
