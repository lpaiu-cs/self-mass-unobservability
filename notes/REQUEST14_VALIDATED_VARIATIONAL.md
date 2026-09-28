# Request 14: validated variational integration

Status: Conjectural. The acceptance target is an interval enclosure of the actual 4-body 1PN continuum flow and its initial-state Jacobian, followed by a documented chain to timing-parameter derivatives and residual errors. Finite-difference agreement does not satisfy this target. Request 10, 12 and 13 evidence remains frozen.

Status: Imported from prior work. Request 13 left D2 open because it supplied no continuous ODE error enclosure. Its initial-state map, interpolation and observation delays also lack certified derivative bounds. Those omissions cannot be covered by an ODE-only result.

Status: Conjectural. Execute in order: export the live engine's internal initial-value problem and native RHS samples; bind a CAPD automatic-differentiation vector field to those equations; verify exact-solvable positive controls; integrate interval state and variational equations without reinitializing uncertainty; record the achieved time domain and any enclosure failure; then propagate only certified quantities to observables. A failure or a partial observation map must leave the full D2 gate false.

Status: Imported from prior work. CAPD provides interval ODE solvers and C1 doubleton sets for validated flow and variational integration. The dependency is isolated under the existing WSL research runtime, built from a recorded commit with rounding-aware compiler flags. Sources: https://github.com/CAPDGroup/CAPD and https://ww2.ii.uj.edu.pl/~wilczak/papers/capd-review.pdf .

## Exact certificate boundary

Status: Proven. The declared IVP treats the exported hexadecimal long-double state and coefficients as exact real numbers. Conversion brackets every such value between adjacent binary64 numbers. The validated flow therefore includes that frozen internal IVP; it does not include an unspecified error in converting a physical timing parameter file into the internal IVP. The exported input reader rejects any nonzero dynamic coupling, non-GR pair coefficient, gamma or beta coefficient. This certificate is for the existing GR 1PN force truncation, not an enclosure of omitted higher-PN physics.

Status: Proven. For each successfully recorded epoch, the CAPD C1 set encloses the continuum state and its 24 by 24 initial-state Jacobian on the declared initial box, with coefficient intervals carried throughout. The method integrates variational equations using automatic differentiation and validates Taylor remainder inclusion; it does not infer an error from sampled finite differences. Identity is the initial Jacobian. No intermediate set is replaced by its midpoint or restarted to narrow its uncertainty.

Status: Proven. An additional augmented run sets M_i=M_i,reference*(1+delta_i), delta_i'=0, initially delta_i=0. Its 28 by 28 variational system includes the four independent fractional-mass derivatives as well as the 24 state derivatives. Both the linear and quadratic mass factors in the 1PN force are differentiated. These 28 IVP coordinates are not the 28 fitted timing parameters. Units and reference Jacobi coefficients are held fixed; the physical parameter initialization, barycenter convention and parameter-dependent unit conversion still require a separate chain rule.

Status: Proven. For a returned interval matrix J with a binary64 midpoint matrix M, let e_ij=max(|M_ij-J_ij.lower|,|J_ij.upper-M_ij|). Then ||J_true-M||_2 <= sqrt(sum e_ij^2). The audit computes each radius as an exact rational number and rounds the final square root upward, checking its square against that rational sum. This gives an actual error upper bound for this Jacobian, rather than a step-difference diagnostic.

Status: Proven. For the geometric delay d=n dot r_p/c with any fixed unit line-of-sight n, a position-box diameter vector w gives a delay diameter at most ||w||_2/c. The native length/time units are restored with outward interval arithmetic. The stored geometric-delay width is the width contribution of this linear readout alone; it is not the error of a computed timing residual, nor a bound on Einstein, Shapiro, aberration, observer-motion, interpolation or emission-time inversion errors.

## Conditioning and exact transformations

Status: Proven. The dyadic change s=128 t, w=v/128 represents the same continuous IVP. If x=S y, the physical initial-state Jacobian is S J_y S^-1. These pullbacks are computed as interval operations before reporting errors in the original coordinates.

Status: Proven. Jacobi relative positions and a weighted center can use fixed dyadic mass-ratio coefficients as a coordinate choice. The forward and inverse affine maps are exactly inverse over the reals even when those coefficients only approximate physical mass fractions. The audit checks this identity symbolically for arbitrary coefficients. The transformed initial interval box may be wider than the exact image; this conservatively enlarges the certificate domain.

Status: Proven. The native RHS uses only pairwise position differences, so the common position cancels identically. Removing it algebraically before automatic differentiation preserves the entire RHS and makes its common-position derivative exactly zero. Common velocity must remain because the 1PN RHS is not Galilean invariant. The transformed RHS is checked against the original specialization both at the exported state and with an added common velocity.

Status: Imported from prior work. Generic Jacobi reconstruction initially failed the CAPD remainder-inclusion test at the first step, both with adaptive steps and fixed-step refinement. Algebraically eliminating common position allowed the same validated inclusion test to pass. These failed attempts and their exact source versions are retained; the solver's acceptance test was not weakened.

## Reproduction and interpretation

Status: Imported from prior work. CAPD 6.1.0 was built from commit `731079217a9254ea2948d742df2b170895effe7f` in the isolated WSL runtime. The producer stores the compile command, `-frounding-math` flags, source and library hashes. No `-ffast-math` is used. The exact-rational oscillator control and nonlinear x'=x^2 interval-box control validate both state and variational outputs. The GR term specialization and coordinate identities are checked by SymPy.

Status: Conjectural. The full D2 gate remains false until the entire observation span, timing-parameter-to-IVP derivatives including the mapping to masses, all delay terms and emission-time inversion are enclosed and propagated into the normalized nuisance projector and likelihood. Successful short-time state and independent-mass derivative enclosures cannot substitute for any of these links. Failure of the width ceiling establishes a limit of this numerical representation, not instability or chaos of the physical system.

Classification: theorem progress through conditional continuum-flow and variational enclosures. No new empirical SEP exclusion or submission-ready full-inference claim follows from this lever.

## Recorded results

Status: Proven. At the binary64 internal epoch 0.07 (approximately 1.68542119 days), the Jacobi run encloses the 24-state flow with midpoint state error norm at most 4.1463103E-11 and initial-state Jacobian operator error at most 0.0000055316127, in the original dimensionless coordinates. The geometric-delay box width is at most 0.000044770472 microseconds. These displayed upper bounds are rounded upward. Raw IVP Jacobian norms cannot be compared directly to singular values of the normalized timing nuisance matrix.

Status: Proven. The Cartesian, scaled and Jacobi state/Jacobian intervals at this common epoch intersect component by component. The augmented mass calculation also agrees with the scaled state block. The four independent fractional-mass column error norms, in pulsar/inner/outer/extra order, are bounded by 1.8732067E-7, 3.6273306E-8, 7.6514719E-11, 1.0275596E-19.

Status: Imported from prior work. The full observation interval extends to approximately 2989.511 days after the reference epoch. Every long-run representation below stopped at the declared Jacobian-width ceiling, without restarting its uncertainty. The listed last epochs are not useful-accuracy guarantees over the whole preceding interval.

| Representation | Last returned epoch (days) | Full horizon completed |
|---|---:|---|
| full-forward | 6.09714964 | False |
| ho-forward | 6.12213234 | False |
| scaled-forward | 31.94608159 | False |
| jacobi-translation-forward | 45.78006494 | False |

Status: Conjectural. The next dependency is a long-span variational representation that avoids this interval overestimation, with the same inclusion test retained. In parallel mathematical preparation, the timing-parameter-to-IVP chain and full delay/inverse-time map must be defined. Extra numerical precision or Kepler-based coordinates are candidate remedies, not established solutions. D2 and full physical inference remain open.
