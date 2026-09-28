# Request 11.4 — comparator identifiability registration

Input revision b1d57e7; before-task checkpoint 992d48b. This continues the user's original item 4: justify the comparator and test information after nuisance fitting. Preserve all REQUEST10 and completed Request 11.1–3 outputs.

Use the full stored rank-90 nuisance span and the existing extra-Fourier covariance family. Compare fixed covariance amplitudes a=0,0.08316109496333533,1 at unit white scale; the middle value is an already observed conditional REML value, not a new globally fitted noise model. Work on sufficient coordinates and verify direct projection agreement.

For lags 2,5,18,52,200,500 days, compare real local derivative orders N=0,...,5 and an even-only {0,2,4} conservative comparator. Use three declared phase vectors: zero, (0.7,1.1,2.2), and the leading physical phases from Request 11.3, with unit carrier amplitudes in every case. The latter is a phase stress, not a matched physical-drive amplitude dictionary. Compute the projected beta information relative to co-fitting only the instantaneous column. Degree-five collapse must be checked as an algebraic identity; a floating-point residual is not new information.

Justify a stronger comparator condition without choosing an arbitrary polynomial order: a reciprocal overdamped gradient system Gamma qdot+Kq=bF with symmetric positive-definite Gamma,K, output b^Tq, and all rates at least Lambda. Derive its positive relaxation weights and a two-frequency quadrature inequality. For a predeclared illustrative gap Lambda=10*max(omega), compute the slow-pole violation and noise-weighted witness after full nuisance. The gap is an explicit hypothetical premise, not established neutron-star microphysics. Also state failure when reciprocity, positivity or the gap is absent.

Derive a continuous known-phase lower bound using the smallest singular value of the whitened six-column map and the phase-invariant coefficient-space residual. Check it against every evaluated phase. This does not make unknown drive phases known. Symbolic controls cover odd conservative quadratic terms, the monotone quadrature ratio, and exact fifth-degree interpolation; numerical controls compare QR/SVD and spectral/direct gradient-system transfer.

Outcome is restricted-comparator progress or an exact absorption boundary. Do not declare a detection, a physical coupling limit, or an EOS-derived gap from these calculations.
