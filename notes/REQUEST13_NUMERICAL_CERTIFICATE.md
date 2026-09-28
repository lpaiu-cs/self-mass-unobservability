# Numerical certificate: supplied and missing parts

Status: Proven. Let r(theta) be the continuous-model residual vector, with computed values r_hat and a proved weighted function-error upper bound epsilon at each endpoint theta+-h e_j. If ||W partial_j^3 r||<=M3_j throughout that interval, the central-difference column has weighted error at most

    E_j <= [epsilon_plus + epsilon_minus]/(2h) + h^2 M3_j/6 + rounding_j.

Different endpoint bounds and column normalization must be carried explicitly. Frobenius aggregation bounds the matrix operator error by sqrt(sum_j E_j^2) after the stated normalization. For a full-column-rank computed matrix with smallest singular value sigma_hat, a sufficient projector perturbation bound is E/(sigma_hat-E), provided E<sigma_hat. The cubic truncation coefficient and the stable residual polynomial identity are checked symbolically by verification/remediation_audit.py.

Status: Proven. A finite set of sampled function values cannot supply the required third-derivative or true-error upper bound without an additional regularity/enclosure argument: a smooth function vanishing at all sampled points can have an arbitrarily large derivative at the central point. Agreement of two finite-difference estimates therefore cannot close this certificate.

Status: Imported from prior work. Request 13 corrects index-test ordering in forward/backward dense integration and exposes the inverse-time tolerance. Control, tighter integration (1e-18 versus 1e-16), denser interpolation (500 versus 250) and tighter inverse-time termination are evaluated separately. The control reproduces the frozen baseline exactly. Tighter inverse-time termination leaves the tested residuals and derivatives identical. Integration and interpolation changes remain visible, and do not consistently reduce step dependence.

Status: Proven. For fixed turns and profiled offset, the spin polynomial can be rearranged to retain the small emission delays separately, cancel the integer turns with a fused multiply-add, and avoid forming a large emission epoch before subtracting its small delay. The new expression equals the old centered model over the real numbers. The isolated stable build falls back to the original code for fractional-phase calculations, unprofiled offsets and the two additional timing special cases.

Status: Imported from prior work. The stable arithmetic baseline differs from the previous numerical implementation by at most 4.8459e-5 microseconds and recovers exactly after zero shifts. Step dependence in the three pilot columns remains approximately 0.104%, 0.0662% and 0.665%. This implementation change is not sufficient to certify the physical derivatives.

Status: Conjectural. Missing inputs to an actual certificate are the integrated continuous-model error enclosure, interpolation remainder bounds, parameter-to-initial-state derivative bounds, higher derivatives or validated variational equations, and accumulated floating-point error bounds. These must hold on a declared region excluding orbit singularities and collisions. The output audit deliberately leaves these inputs null and the certificate gate false. No result in this note marks D2 complete.

Classification: numerical repair and theorem progress, with the empirical certificate unresolved. New producers, patches, raw residuals and provenance are under outputs/research-remediation; REQUEST10 and REQUEST12 remain unchanged.
