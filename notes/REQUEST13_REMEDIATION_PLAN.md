# Request 13: complete the three unfinished research objectives

Baseline: feede76; pre-task checkpoint: 9098434. User instruction: plan an appropriate order and execute all three remediations. This explicitly authorizes the numerical runtime work needed here despite the repository's older default restriction. Preserve REQUEST10 and REQUEST12 artifacts and verdicts.

## Order and dependencies

1. **Numerical reliability.** Trace the full parameter-to-residual path. Separate integration, interpolation, floating-point cancellation, and finite-difference truncation. Repair identified implementation defects in an isolated build. Compare independent tolerances/meshes and derivative constructions. A rigorous certificate requires an actual upper bound on the residual/derivative error, including all stages; step agreement alone never passes this gate.
2. **Physical matching.** Fix a published scalar-tensor theory/EOS and its units, stable branch, mass definition and domain before evaluating it. Reproduce published matched coefficients as a positive control, then calculate a stellar structure and response for a specified tabulated EOS. Do not insert a gravitational mass into a baryon-mass fit. Compute dynamical response/damping from the selected theory or published mode data, not a freely chosen relaxation time. A benchmark prototype is not an empirical J0337 constraint.
3. **Nonlinear observation inference.** Use physical-domain constraints and the complete live timing residual evaluator. Profile or marginalize explicitly stated linear nuisance terms; jointly optimize nonlinear timing parameters and noise hyperparameters with their determinant terms. Include the zero pulse assignment and competing assignments, multiple physical starts, convergence diagnostics and injection checks. Finite local displacement tests do not pass this gate. Tie any physical signal conclusion to the forward-model matching and the numerical error budget.
4. **Integration.** Re-run relevant symbolic and numerical checks, update the five maintained dynamic-chi notes, record outstanding gates without relabeling them completed, and update the manuscript only for validated findings. Checkpoint before and after major tasks.

The final physical signal fit depends on both steps 1 and 2. Read-only physical-source work and preliminary constrained timing fits may proceed while numerical calculations run, but do not receive a certified scientific verdict early.

## Completion gates (all initially open)

| Gate | Required artifact | Status |
| --- | --- | --- |
| D1 | Independent numerical settings, all derivative columns, precision/mesh provenance | Open |
| D2 | Instantiated rigorous error enclosure for the continuous timing model and propagated projector/likelihood error | Open |
| E1 | Published EOS/theory matching reproduced with dimensions and branch checks | Met for the published beta=-5 prototype; distinct from the SLy candidate |
| E2 | Numerical stellar equilibrium and scalar response at the target mass with convergence controls | Met numerically for the specified SLy/beta=-4 zero-scalar branch; no rigorous continuum certificate |
| E3 | Matched dynamics, companion response and complete physical force/readout mapping | Open |
| N1 | Physical-domain live nonlinear joint timing/noise fit; no clamped evaluations accepted | Open |
| N2 | Competing pulse assignments, multiple starts, convergence and injection/coverage checks | Open |
| N3 | Physical signal inference using matched forward model and numerical error budget | Open |

Status: Conjectural. These are work objectives, not promised positive detections or guarantees that the chosen candidate produces orbital-scale relaxation. A failed gate requires a corrective attempt; writing a no-go note does not complete an empirical gate.

## Initial concrete choices

Status: Imported from prior work. The existing engine uses long-double integration, integration tolerance 1e-16 and 250 interpolation samples per inner orbit. The archived fit has 28 nonlinear timing parameters and 12,474 TOAs. The previous weak-gap compensation exits the eccentricity domain.

Status: Imported from prior work. The initial matching positive control is Khalil et al., Phys. Rev. D 106, 104016, https://arxiv.org/html/2206.13233v2, equations 28--31 and their beta=-5 two-piece-polytrope model. It is a prototype, not an observationally accepted J0337 parameter choice. Its baryon mass is explicitly distinct from the observed gravitational mass.

Further numerical settings and candidate-domain choices are recorded before their production runs. Prior failure cases remain in the audit trail.

## Execution status after the first remediation cycle

Status: Imported from prior work. D1 has four independent setting pilots, a second residual arithmetic implementation, and a fresh full 28-column Jacobian at the nonlinear fit checkpoint. D2 remains open: neither tighter settings nor the algebraically stable residual implementation provides a rigorous enclosure of the integrated continuous model. E1/E2 now contain actual computations, including a matched complex stellar pole and an independent EOS/TOV crosscheck. N1 has physical-domain live joint fits and a fresh-Jacobian continuation; stationarity and complete inference are not certified. E3, N2 and N3 remain open.

The next required implementation is a differentiable/validated forward calculation with propagated integration, interpolation and floating-point bounds on an explicitly bounded physical parameter region. It must be connected to the full force/readout model before the final physical likelihood can pass. Continuing ordinary finite-difference step scans alone is not this implementation. The original all-remediation objective remains unfinished; this plan has not redefined it as completion of the pilots.
