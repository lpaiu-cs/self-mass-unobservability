# Request 11.4 — physically restricted comparators and projected information

Design fb269b9; command `python verification/comparator_audit.py` (four BLAS threads). Evidence: `outputs/research-completion/comparator-audit.json`. The run passes eight algebra/metric control groups. All input arrays and the earlier outcomes are preserved.

## Which restriction has a physical reason?

Status: Proven. For a real scalar drive with vanishing boundary variations, a local conservative quadratic action (1/2) integral F P(d/dt)F contributes the self-adjoint kernel [P(D)+P(-D)]/2 to its conjugate response. Odd derivatives cancel as boundary terms. Thus a real even polynomial is justified for this particular time-reversal-invariant conservative response. Dissipation, nonconjugate drive/readout, a nonstationary background or unknown phase calibration invalidates that restriction. It is not a restriction on every static EFT or timing nuisance parameter.

Status: Counterexample candidate. A broader physically specified comparator is a reciprocal overdamped gradient system

    Gamma qdot + K q = b F,      output = b^T q + c0 F,

with real symmetric positive-definite Gamma and K, and a lower bound Lambda on every generalized rate. This follows from a positive quadratic state energy with linear conjugate forcing and positive quadratic Rayleigh dissipation. The lower rate bound is an explicit scale-separation assumption, not measured J0337 microphysics.

Status: Proven. Diagonalizing Gamma^(-1/2) K Gamma^(-1/2) yields

    H(s)=c0+sum_j a_j/(1+s*tau_j),   a_j>=0,   0<tau_j<=1/Lambda.

For omega_h>omega_l>0, define S(omega)=-Im H(i omega)/omega. Nonnegative weights give

    S(omega_l) <= R_Lambda S(omega_h),
    R_Lambda=[1+(omega_h/Lambda)^2]/[1+(omega_l/Lambda)^2].

This follows because (1+omega_h^2*tau^2)/(1+omega_l^2*tau^2) increases with tau. A positive single pole with tau>1/Lambda violates the inequality. The arbitrary real instantaneous c0 cancels from S. This test avoids choosing an arbitrary Taylor-polynomial order, but needs reciprocal nonnegative relaxation weights and the gap. A general dissipative response, nonconjugate readout, oscillatory modes or negative residues need not obey it. The gap and readout assumptions must be independently matched before an empirical claim.

## Does the distinction survive the stored nuisance?

Status: Proven. Let A be the covariance-whitened six-carrier map after removing the full nuisance, W_phi the coefficient columns of the declared derivative comparator, and w_phi the pole coefficient. If A has smallest singular value s_min>0, then

    min_c ||A(w_phi-W_phi c)||^2 >= s_min^2 min_c ||w_phi-W_phi c||^2.

Changing known carrier phases rotates each cosine/sine pair orthogonally. Therefore the coefficient-space distance on the right is phase invariant. This gives a lower bound over the continuous three-phase domain, not just the sampled origins. Unknown independently fitted phases are a different comparator: freely fitting both quadratures spans all six stored columns and absorbs the periodic signal exactly.

Status: Imported from prior work. On the frozen arrays, the six-column map remains rank six after all 90 nuisance directions. Its singular-value condition number is 98.19 for diagonal covariance, 111.69 at the earlier conditional REML amplitude a=0.0831611, and 124.82 under a=1 extra-Fourier stress. Direct/sufficient-coordinate whitening and QR/SVD comparator residuals agree within their registered controls. These are numerical properties of stored derivatives, not a derivative-accuracy certificate.

| Status | Allowed derivative powers | Minimum information / instantaneous-only information | Maximum |
| --- | --- | ---: | ---: |
| Imported from prior work | 0,1 | 0.082888 | 0.999988 |
| Imported from prior work | 0,1,2 | 8.35514e-6 | 0.191680 |
| Imported from prior work | 0,1,2,3 | 3.18699e-6 | 0.185620 |
| Imported from prior work | 0,1,2,3,4 | 2.87821e-6 | 0.134298 |
| Imported from prior work | 0,2,4 only | 0.00303486 | 0.643531 |
| Proven | 0,1,2,3,4,5 | Exact absorption | Exact absorption |

Status: Imported from prior work. Minima/maxima range over the six registered lags, three known-phase stress vectors and three fixed covariance amplitudes. A sampled fourth-order case retains only 2.87821e-6 of the instantaneous-only information, so its unit-noise standard error grows by about 589 times. Mathematical noninterpolation does not mean useful experimental precision. Fifth-order numerical relative residual information lies between 9.29e-30 and 1.38e-24 and is classified as floating-point error around the algebraic zero, not measurable information.

Status: Proven. For the fast-spectrum inequality, construct a coefficient-space vector v extracting S_l-R_Lambda*S_h. If v^T theta<=0 for the comparator and v^T w>0 for a positive unit pole, the covariance-weighted distance of beta*w from the entire fast-spectrum cone is at least beta*(v^T w)/sqrt(v^T(A^T A)^(-1)v). A freely fitted real instantaneous term is annihilated. Thus the witness tests a restricted physical class after nuisance fitting, rather than simply counting frequencies.

Status: Imported from prior work. For the explicitly illustrative Lambda=10*omega_in, every registered slow-pole lag violates that inequality, and the run records its positive noise-weighted witness for all phases/covariances. No data significance, rate-gap measurement or EOS limit is inferred. Unit carrier amplitudes are retained; the phase-only stress using Request 11.3 values is not a restoration of its withdrawn physical amplitude map.

Decision: include the reciprocal fast-spectrum criterion and the severe projected-information loss in the paper. Keep the exact fifth-degree no-go. Neither arbitrary polynomial freedom nor freely fitted independent phases can be excluded by this dataset alone. A conservative or gapped reciprocal comparator is an explicit conditional hypothesis, not a post-result prior chosen for a detection.

Classification: theorem progress (physical-comparator inequality and continuous-phase bound) and loophole progress (measured conditional information and exact absorption cases).
