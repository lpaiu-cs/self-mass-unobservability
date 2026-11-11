# Finite-frequency identifiability of relaxation responses in free-fall timing

**Author:** Juneyoung Kim\thanks{Corresponding author: \href{mailto:lpaiu.cs@gmail.com}{lpaiu.cs@gmail.com}}\\{\small Independent researcher, Seoul, Republic of Korea}
**Date:** 30 September 2026

## Abstract

**Proven.** A relaxation pole can be distinct from an instantaneous response yet indistinguishable in a timing measurement. We formulate a procedure that separates three questions: whether a declared comparator can interpolate the response at the measured carriers, whether the difference survives nuisance fitting, and whether the numerical and noise models support inference. For a known drive at \(K\) distinct positive frequencies, a real derivative comparator exactly matches a nonzero single pole if and only if its degree is at least \(2K-1\). Below that boundary, a projected distance determines local amplitude identifiability; its singular-value bounds expose the loss between coefficient-space distinction and measurement precision. Exact controls exhibit surviving, erased and arbitrarily attenuated differences. Positivity and an independently known relaxation-rate gap give a separate comparator test. These results combine classical interpolation and linear estimation in a free-fall response setting. **Imported from prior work.** A stored PSR J0337+1715 analysis illustrates the procedure's limits: fourth-degree comparators retain as little as \(4.63\times10^{-6}\) of the reference information in the tested cases, a physical-drive phase condition invalidates our earlier auxiliary mapping, and recomputed timing derivatives change the nuisance geometry substantially. **Conjectural.** The archived application supports a conditional methodological case study, not a detected relaxation or an empirical strong-equivalence-principle constraint. The accompanying code checks the analytic controls and replays the stored inference without rerunning the timing engine.

---

**Keywords:** dynamical identifiability; relaxation response; nuisance projection; pulsar timing; tests of gravity

## 1. Introduction

**Imported from prior work.** Internal structure enters compact-body effective descriptions through response coefficients and, when needed, explicit dynamical variables \cite{goldberger2006eft,porto2016eft}. Internal modes produce response poles \cite{chakrabarti2013response,steinhoff2016dynamical}; dynamical scalarization supplies a monopolar example \cite{khalil2022scalarization}. Free-fall tests constrain specified couplings within a timing model \cite{will2014confrontation}. PSR J0337+1715 provides a particularly precise setting for static strong-equivalence-principle (SEP) tests \cite{archibald2018universality,voisin2020sep}, with a later timing analysis and public release available for methodological work \cite{voisin2025planet,voisin2025release}.

**Proven.** A finite timing record raises a narrower question than whether a response contains a pole. Only certain carrier coefficients are sampled, and fitted timing parameters can absorb their effects. Thus a transfer function can fail to lie in a specified instantaneous family while its measured response lies in that family's nuisance-projected image. Even when the distinction survives, it need not be large enough to measure.

**Imported from prior work.** The mathematical ingredients are established. The real interpolation count is a special case of polynomial and rational interpolation \cite{mayo2007realization} and agrees with persistent-excitation counting \cite{ljung1999system}. Nuisance removal and information after co-fitting are linear-estimation constructions, also used with correlated noise in pulsar timing \cite{vanhaasteren2013noise}. Positive-spectrum tests connect to the Stieltjes moment problem \cite{schmudgen2017moment}; distinguishing close relaxation times is an ill-conditioned inverse problem \cite{honerkamp1989relaxation}.

**Proven.** The contribution here is an explicit sequence of compatibility and separation tests for a declared free-fall response measurement. It connects the finite-frequency comparator boundary to the actual response matrix after nuisance removal, provides exact controls for the failure of that connection, and identifies when a physical-drive or numerical-validity check must stop the inference. The projection identity and interpolation theorem are not claimed as new general mathematics. Sections 2--4 give the complete method and analytic controls; Section 5 applies the sequence to stored results, including its unresolved validity checks.

**Imported from prior work.** The J0337 application audits our own unpublished repository analysis. It does not correct the published SEP results cited above. Supplementary Information (SI) preserves the full derivations, static operator catalog, numerical audits and conditional white-dwarf calculation, with its own numbering. Supplementary Code provides the new exact controls and the records needed for the compact replay. Neither a source hash nor a successful replay certifies the physical adequacy of the model that generated those records.

Scientific statements carry one of four labels. **Proven** denotes mathematics under stated premises; **Imported from prior work** denotes literature or archived numerical results; **Counterexample candidate** denotes a proposed physical realization; **Conjectural** denotes an interpretation whose physical or numerical premises remain unestablished. The labels do not replace the premises or their validation.

## 2. A response model and its measurement

### 2.1 Declared drive, response and comparator

**Counterexample candidate.** Let a known dimensionless drive F excite one internal variable and a scalar readout q:

\begin{equation}
\tau_\chi\dot\chi+\chi=\alpha F(t),\qquad
q=c_YF+c_\chi\chi,\qquad \beta=\alpha c_\chi,\quad \tau_\chi>0.
\label{eq:model}
\end{equation}

The readout corrects a specified observable coupling. Assigning q to an inertial mass or to a pair coupling are different modeling choices and need not give the same forces. The settled transfer is

\begin{equation}
G(i\omega)=c_Y+\frac{\beta}{1+i\omega\tau_\chi}.
\label{eq:response}
\end{equation}

The periodic model excludes, or assumes the decay of, the homogeneous term \(\chi_h e^{-(t-t_0)/\tau_\chi}\).

**Proven.** Write the retained drive as \(F(t)=\sum_k d_k\cos(\omega_k t+\phi_k)\), with distinct positive \(\omega_k\), known \(d_k>0\) and known phases. Use the carrier convention \(\cos\omega_k t,\cos(\omega_k t+\pi/2)=-\sin\omega_k t\). Stacking real and imaginary parts gives the single-pole vector w and real derivative-comparator matrix \(V_N\):

\begin{equation}
\begin{split}
w_{2k-1}+i w_{2k}&=\frac{d_ke^{i\phi_k}}{1+i\omega_k\tau_\chi},\\
(V_N)_{2k-1,n+1}+i(V_N)_{2k,n+1}
 &=d_ke^{i\phi_k}(i\omega_k)^n,\quad 0\le n\le N.
\end{split}
\label{eq:carrier-vectors}
\end{equation}

Column n represents a shared real coefficient multiplying \(F^{(n)}\); column zero includes the co-fitted instantaneous term. Rescaling a nonzero column for numerical conditioning leaves its span unchanged. Freely chosen quadratures at each carrier are a larger comparator, which can absorb every vector in this carrier space.

### 2.2 Projected separation and its assumptions

**Proven.** Consider a fixed linear measurement model

\begin{equation}
y=J\eta+C\theta+\varepsilon,\qquad
\operatorname{Cov}(\varepsilon)=\Sigma\succ0.
\label{eq:measurement}
\end{equation}

Here J contains nuisance derivatives, C contains the 2K carrier-response columns, and \(\theta=V_N a+\beta w\). For the following statements, J, C, the drive and \(\Sigma\) are specified, independent of the fitted amplitudes. Choose a whitening matrix L with \(L^TL=\Sigma^{-1}\); let \(P_J\) project orthogonally onto \(\operatorname{col}(LJ)\), and define

\begin{equation}
A=(I-P_J)LC,\qquad P_{AV}=\operatorname{proj}_{\operatorname{col}(AV_N)}.
\label{eq:projected-map}
\end{equation}

The projector uses the actual column space, including when columns are dependent. No statistical prior deleting weak nuisance directions is implicit in this definition.

**Proven. Proposition 1 (projected separation).** At fixed \(\tau_\chi\), the amplitude \(\beta\) is locally identifiable against this comparator precisely when

\begin{equation}
\mathcal I_{\beta,N}
=\|(I-P_{AV})Aw\|^2
=\min_a\|A(w-V_Na)\|^2>0.
\label{eq:information}
\end{equation}

With a fixed Gaussian covariance, this is the Fisher information after co-fitting the comparator and nuisance coefficients, and the unconstrained linear estimator has variance \(\mathcal I_{\beta,N}^{-1}\). With a non-Gaussian error distribution the rank statement still holds, but this alone gives no confidence level. The claim concerns one amplitude at a specified lag, not global identification of an unknown lag, unknown drive or nonlinear physical model.

**Proven.** Whiten Equation (\ref{eq:measurement}) and minimize first over \(\eta\), then over a. These two least-squares eliminations leave the response \((I-P_{AV})Aw\) for \(\beta\). It adds a direction if and only if it is nonzero, proving Equation (\ref{eq:information}) and the Gaussian information statement. Equivalently, the rank of \([AV_N,Aw]\) must exceed that of \(AV_N\).

**Proven.** Let \(\delta_N=\min_a\|w-V_Na\|\). If A has full column rank, its extreme singular values give

\begin{equation}
\sigma_{\min}(A)^2\delta_N^2
\le \mathcal I_{\beta,N}
\le \sigma_{\max}(A)^2\delta_N^2.
\label{eq:map-bound}
\end{equation}

For the lower bound apply \(\|Au\|\ge\sigma_{\min}\|u\|\) before minimization; for the upper bound evaluate at a coefficient-space minimizer. A rank-deficient map has no positive lower singular bound on the full carrier space. The exact rank test (\ref{eq:information}) remains valid, but noninterpolation before projection no longer suffices. The bounds are deterministic statements for the supplied A; they do not bound errors in J, C or the covariance.

## 3. Comparator boundaries

### 3.1 Exact finite-frequency boundary

**Proven. Theorem 1.** For real \(c_Y,\beta\), \(\beta\ne0\), \(\tau_\chi>0\), and K distinct positive frequencies, a real polynomial of degree at most N matches \(G(i\omega_k)\) at every carrier if and only if \(N\ge2K-1\).

**Proven.** If \(P_N\) matches, the real polynomial
\(R(z)=(1+\tau_\chi z)[P_N(z)-c_Y]-\beta\)
has the 2K distinct roots \(\pm i\omega_k\). When \(N<2K-1\), its degree is smaller than 2K, so R would vanish identically, contradicting \(R(-1/\tau_\chi)=-\beta\). Conversely, put

\begin{equation}
Q_K(z)=\prod_{k=1}^K(z^2+\omega_k^2),\qquad
P_{2K-1}(z)=c_Y+
\frac{\beta[1-Q_K(z)/Q_K(-1/\tau_\chi)]}{1+\tau_\chi z}.
\label{eq:interpolant}
\end{equation}

The numerator vanishes at \(-1/\tau_\chi\), so division gives a real polynomial of degree 2K-1. It agrees at every carrier.

**Proven.** Since all drive multipliers in Equation (\ref{eq:carrier-vectors}) are known and nonzero, this is also the boundary for \(w\in\operatorname{col}(V_N)\). Thus N at least 2K-1 implies \(\mathcal I_{\beta,N}=0\) for every A. Below the boundary, a full-column-rank A guarantees strictly positive information, but no useful lower magnitude without its conditioning and \(\delta_N\). Three carriers exclude degrees up to four before projection; degree five interpolates them exactly.

**Proven.** The theorem admits unrestricted shared real coefficients. It does not justify a derivative cutoff or establish natural coefficient sizes. In the adiabatic domain \(|\omega\tau_\chi|\le\rho<1\), truncating the geometric series in Equation (\ref{eq:response}) gives error at most \(|\beta|\rho^{N+1}\). A mathematically distinct pole can therefore be arbitrarily close to a low-order comparator even before measurement loss.

### 3.2 A physically restricted alternative

**Proven.** Even derivative powers can be justified for a time-reversal-invariant conjugate quadratic response: integration by parts replaces \(P(D)\) by \([P(D)+P(-D)]/2\), eliminating odd powers when boundary terms vanish. This restriction need not hold for a dissipative or nonconjugate readout. SI Section 3.4 states the underlying action assumptions.

**Counterexample candidate.** An alternative to selecting a polynomial degree is a reciprocal dissipative response

\begin{equation}
\Gamma\dot q+Kq=bF,\qquad q_{\rm out}=b^Tq+c_0F,
\label{eq:gradient}
\end{equation}

with real symmetric positive-definite K and \(\Gamma\), positive quadratic energy and Rayleigh dissipation. Suppose, independently, that all eigenvalues of \(\Gamma^{-1/2}K\Gamma^{-1/2}\) are at least \(\Lambda>0\).

**Proven.** Diagonalization gives \(H(z)=c_0+\sum_j a_j/(1+z\tau_j)\), with \(a_j\ge0\) and \(0<\tau_j\le\Lambda^{-1}\). For \(\omega_h>\omega_l>0\) and \(S(\omega)=-\operatorname{Im}H(i\omega)/\omega\),

\begin{equation}
S(\omega_l)\le R_\Lambda S(\omega_h),\qquad
R_\Lambda=\frac{1+(\omega_h/\Lambda)^2}{1+(\omega_l/\Lambda)^2}.
\label{eq:gap}
\end{equation}

The ratio \((1+\omega_h^2\tau^2)/(1+\omega_l^2\tau^2)\) increases with \(\tau>0\); summing its upper bound with nonnegative weights proves the result, also for a convergent positive spectrum. A positive pole slower than the gap violates it; the real instantaneous term cancels.

**Proven.** When A has full column rank, let u extract \(S_l-R_\Lambda S_h\) from the drive-calibrated carrier coefficients. Every comparator satisfies \(u^T\theta\le0\). If \(u^Tw>0\) and \(\beta>0\), the distance in measurement space from \(\beta w\) to this entire comparator class is at least

\begin{equation}
\frac{\beta u^Tw}{\sqrt{u^T(A^TA)^{-1}u}}.
\label{eq:witness}
\end{equation}

Indeed \(u^T(\beta w-\theta)\ge\beta u^Tw\), and Cauchy--Schwarz in the metric \(A^TA\) supplies the bound. This is a deterministic separation, not a detection significance. Negative residues, nonreciprocal readout, oscillatory modes or an unestablished rate gap invalidate this comparator. The gap is not inferred from an equilibrium susceptibility.

## 4. A reproducible decision procedure

### 4.1 Exact controls for measurement loss

**Proven.** A dimensionless example isolates the distinction without any timing-engine calculation. Choose two carriers \((1,2)\), \(\tau_\chi=1\), unit drive amplitudes and zero phases, and compare powers zero and one. Then

\begin{equation}
w=(1/2,-1/2,1/5,-2/5)^T,\qquad
V_1=\begin{pmatrix}1&0\\0&1\\1&0\\0&2\end{pmatrix}.
\label{eq:example}
\end{equation}

Let \(P=V_1(V_1^TV_1)^{-1}V_1^T\), \(Q=I-P\). Direct rational arithmetic gives \(\delta_1^2=w^TQw=117/1000>0\). Table 1 applies Equation (\ref{eq:information}) to three maps.

\begin{table}[htbp]
\centering
\caption{Proven. Exact controls with the same response and comparator. The parameter satisfies \(0<\epsilon\le1\). These are constructed dimensionless maps, not fitted astronomical measurements.}
\label{tab:controls}
\begin{tabular}{lll}
\toprule
Measurement map A & Information & Interpretation\\
\midrule
\(I\) & \(117/1000\) & Difference survives\\
\(P\) & \(0\) & Difference erased\\
\(P+\epsilon Q\) & \(117\epsilon^2/1000\) & Full rank, arbitrarily weak\\
\bottomrule
\end{tabular}
\end{table}

**Proven.** In the last row, AV equals V and the map multiplies its orthogonal complement by \(\epsilon\). Its smallest singular value is \(\epsilon\), so it saturates the lower bound in Equation (\ref{eq:map-bound}). At degree three the interpolation theorem instead makes the information zero for every map. Supplementary Code checks these statements with exact symbolic arithmetic. This example verifies the algebra, not an astrophysical likelihood or numerical derivative accuracy.

### 4.2 Application protocol and stopping conditions

**Proven.** For a fixed linear model, the following ordering keeps transfer-function and measurement claims separate:

1. Declare the physical readout, drive amplitudes and phases, carrier set, transient treatment and comparator class. Failure of a model's phase or coupling relations excludes that mapping before a coefficient interval is interpreted.
2. Form the whitened nuisance-projected map A and the comparator image. Evaluate Equation (\ref{eq:information}), its rank condition and conditioning; use Equation (\ref{eq:witness}) only with its additional physical assumptions. Zero information stops amplitude identification within the declared model.
3. Treat finite nonzero information as a design diagnostic. Confidence regions additionally require a justified error model and reliable response and nuisance derivatives. A covariance estimated from the same data needs calibration of that procedure; the fixed-covariance variance formula alone is insufficient.
4. For lag or phase selection, propagate a simultaneous carrier-coefficient region through the physical map rather than interpreting a collection of pointwise intervals as simultaneous coverage.

**Proven.** For the last step, if E covers the true carrier vector \(\theta\) with probability at least \(1-\alpha\), then
\(\{(b,\tau,\phi):W(\tau,\phi)b\in E\}\)
covers any true physical parameter tuple whose specified map gives that vector with at least the same probability. This follows directly from set inclusion. It does not repair an incorrect W, an unreliable carrier fit, or omitted mean terms; an empty physical section can indicate model inadequacy rather than a useful upper bound on its amplitude.

**Proven.** Two further limits remain even for perfect periodic measurements. A positive reciprocal spectrum can be tested for one observable relaxation time using a two-sample moment-variance equality (SI Section 3.5), but arbitrarily close positive poles approach a single pole quadratically, preventing uniform finite-precision discrimination without separation assumptions. Periodic carriers also do not determine a transient response: the stable proper filter
\(s\prod_k(s^2+\omega_k^2)/(s+\lambda)^{2K+2}\), \(\lambda>0\), vanishes at zero and every retained carrier yet changes the response to an exponential input. The full proofs and state-counting qualifications are retained in SI Sections 3.5 and 5.8.

### 4.3 Computational scope and AI assistance

The revision adds only the exact analytic controls above. The numerical application below reuses frozen records; it does not rerun stellar evolution, the timing engine, noise simulations or the finite-difference derivatives. Supplementary Code contains \path{verification/projected_identifiability.py}, the stored comparator results and the compact conditional-inference replay, with commands and hashes. Symbolic control identities are evaluated with SymPy; replay and table arithmetic use NumPy and the Python standard library.

AI agents wrote and ran most of the original analysis and verification code, drafted scientific text, and assisted internal review under the author's direction. They also assisted the present reorganization and the analytic control script. The disclosure in the Acknowledgements identifies the tools and the validation performed. Internal AI reviews and reproducibility checks are not independent experimental confirmation.

## 5. An archived timing analysis as a conditional case study

### 5.1 What the measurement map retains

**Imported from prior work.** The stored application uses 12,474 public Nançay pulse times over approximately 2987.9 days in 2013--2021, a published planet-model baseline and released timing software \cite{voisin2025planet,voisin2025release}. Six carrier-response columns are finite differences of a modified pair-coupling model at \(\omega_{\rm in}=3.856137\) rad/day, \(\omega_{\rm out}=0.019200\) rad/day and their difference. The high-frequency carriers are close; distinct frequencies are not equally informative measurements. SI Sections 4.2 and 5.1 specify the modification and the archived inputs.

**Imported from prior work.** The full nuisance construction retains 90 directions, including 28 timing derivatives, an offset, 30 low-frequency sine/cosine pairs and a static-SEP guard. A truncated construction retains 71. With the stored derivatives, interval widths depend on this choice by factors up to 27.5 in the audited cases (SI Sections 5.3 and 5.5). No physical prior excluding the discarded subspace was established, so the full span is the reference.

**Imported from prior work.** A stored comparator audit fixes six lags, three known phase sets and covariance amplitudes \(a=0,0.0831611,1\). It uses unit carrier amplitudes, including in the phase-only stress drawn from the physical model; it is not a scan over all physical drives. The projected maps have rank six and condition numbers 98.19--124.82. Table 2 reports minima of the information ratio relative to co-fitting only the instantaneous term. They reproduce archived results, not new fits or noise-independent lower bounds.

\begin{table}[htbp]
\centering
\caption{Imported from prior work, except the exact zero marked Proven. Minimum information retained in the stored comparator audit. The ratio is \(\mathcal I_{\beta,N}/\mathcal I_{\beta,0}\) for the same lag, phase and covariance. SI Section 5.6 specifies the grid and continuous-phase bound.}
\label{tab:retained}
\begin{tabular}{lc}
\toprule
Comparator powers & Minimum retained information\\
\midrule
\(\{0,1\}\) & \(0.005389\)\\
\(\{0,1,2\}\) & \(6.29\times10^{-5}\)\\
\(\{0,1,2,3\}\) & \(1.17\times10^{-5}\)\\
\(\{0,1,2,3,4\}\) & \(4.63\times10^{-6}\)\\
\(\{0,1,2,3,4,5\}\) & \(0\), exact (Proven)\\
\(\{0,2,4\}\) & \(0.0311\)\\
\bottomrule
\end{tabular}
\end{table}

**Imported from prior work.** The fourth-degree comparator widens the unit-noise standard error by about 465 times in the worst tested case. **Proven.** At degree five, reported floating-point residuals are rounding error around the exact zero, not a surviving measurement channel. **Imported from prior work.** The gap witness is positive at every tested lag for the illustrative \(\Lambda=10\omega_{\rm in}\); this rate gap was not physically matched or measured.

### 5.2 Drive compatibility precedes physical interpretation

**Counterexample candidate.** A restricted scalar-charge realization has a stable pulsar charge \(Q_p\) with inertia I, damping \(\Gamma\) and stiffness \(\kappa>0\). With fixed companion charges and constant leading orbital masses, reciprocal forces follow from the interaction \(\sum_{A<B}Q_AQ_B/r_{AB}\), while
\(I\ddot Q_p+\Gamma\dot Q_p+V'(Q_p)=\sum_jQ_j/r_{pj}\).
This uses the dynamical-charge framework of \cite{khalil2022scalarization}, not a new gravity theory. SI Sections 4.3--4.5 retain the action and full reduction conditions.

**Proven.** If both companion charge/mass ratios equal \(a_w\), the small-response overdamped limit gives the common pair modulation

\begin{equation}
\delta\Delta(\omega)=\frac{B\,\delta U(\omega)}{1+i\omega\tau_\chi},\qquad
B=\frac{a_w^2}{\kappa m_p},\quad \tau_\chi=\Gamma/\kappa,
\label{eq:matching}
\end{equation}

where \(\delta U=\sum_jm_j[1/r_{pj}-\langle1/r_{pj}\rangle]\) and \(\beta=BU_*\ge0\) for \(F=\delta U/U_*\). The reduction requires small nonlinear force terms, decayed homogeneous modes, \(I\kappa/\Gamma^2\ll1\), and \(I\omega^2/|\kappa+i\Gamma\omega|\ll1\) at each carrier. Companion susceptibility must not appreciably shift \(\kappa\). Unequal companion charge/mass ratios do not in general produce the common modulation. Fixed equilibrium stiffness alone leaves \(\Gamma\), hence the relaxation time, undetermined.

**Proven.** In an aligned coplanar expansion, linear separately in eccentricities and the inner-to-outer separation ratio, let \(\lambda_p,\lambda_b\) be the inner and outer mean longitudes and \(M_A=\lambda_A-\varpi_A\). The leading potential drive has positive amplitudes and signs

\begin{equation}
\delta U=A_{\rm in}\cos M_p+A_{\rm out}\cos M_b
-A_{\rm dif}\cos(\lambda_p-\lambda_b)+\cdots .
\label{eq:physical-drive}
\end{equation}

Its phases therefore obey

\begin{equation}
\phi_{\rm dif}-\phi_{\rm in}+\phi_{\rm out}
=\pi+\varpi_p-\varpi_b\pmod{2\pi}.
\label{eq:closure}
\end{equation}

This combination is invariant under a common origin shift because \(\omega_{\rm dif}=\omega_{\rm in}-\omega_{\rm out}\). It tests compatibility with a common physical drive without fitting a relaxation coefficient.

**Imported from prior work.** The released orbital parameters follow the ELL1 convention \(\eta=e\sin\varpi\), \(\kappa=e\cos\varpi\) \cite{lange2001ell1,voisin2025release}. With those parameters the closure is 3.16481295 radians. Our earlier auxiliary drive dictionary instead had zero closure, from a pericenter convention and difference-term sign error. **Proven.** No common origin shift repairs this discrepancy. **Imported from prior work.** The corresponding coupling interpretation of our earlier intervals was withdrawn. This is a correction of our unpublished analysis, not of published SEP measurements. SI Section 4.4 gives amplitudes, phases and the provenance of the correction.

### 5.3 The statistical and numerical stopping conditions

**Imported from prior work.** A subsequent stored analysis uses the corrected leading drive and a simultaneous six-coefficient region. Its quadratic threshold, 12.8241766, was calibrated with 8192 draws at each of four specified covariance amplitudes before independent validation; inclusion was 95.17--95.43 percent within those four generating conditions. The archived data statistic is 16.3525. At six evaluated lags from 2 to 500 days, the instantaneous-plus-relaxation physical sections are empty relative to that frozen threshold. This is conditional on the archived derivatives and the specified mean and covariance models. The selected covariance family is not an independently validated description of astrophysical noise.

**Imported from prior work.** The 2-day minimum, about 12.84, exceeds the threshold only marginally; individual calibration quantiles have Monte Carlo standard errors near 0.13. For the nonnegative-\(\beta\) realization the smallest section minimum is approximately 14.59 in the same stored fit, but this greater conditional separation does not bound derivative or model errors. No relaxation is detected, and the excess is not attributed. An empty section is not an EOS-matched upper limit. Supplementary Code reconstructs the drive columns and solves all six sections from the exported coefficients and precision matrix.

**Imported from prior work.** Recomputing the 28 timing derivatives at half the archived steps preserves rank 90 but changes the nuisance span with maximum principal sine 0.999913. Physical-drive standard errors change by factors 0.7858--0.8722 and fitted coefficients shift. A rigorous derivative-error budget remains unavailable. The leading-drive approximation also omits input terms with RMS 3.37427 percent of its retained RMS; that is not a bound on omitted timing error. The stored nonlinear pulse-count follow-through does not exclude all alternative nonlinear solutions. These limitations are reported in SI Sections 5.9--5.10.

**Conjectural.** Consequently, step 3 of the protocol has not passed for empirical interpretation of this case. Timing-model mismatch, derivative error, omitted forces or drive terms, and covariance misspecification remain possible explanations for the excess. The case demonstrates how the declared method diagnoses the limits of an analysis; it does not establish a new astrophysical exclusion.

## 6. Discussion and conclusion

**Proven.** The useful boundary is a sequence, not a frequency count alone. Interpolation determines whether the declared comparator can absorb the response before measurement. Projected distance determines whether any difference survives the actual linear measurement map. Its conditioning determines how little information may remain. A positive result at these stages is still conditional on the validity of the drive, derivatives and error model. The exact examples separate these logical steps without treating a generic linear-algebra identity as a new physical law.

**Imported from prior work.** The archived application supplies a concrete warning: three carriers can retain very little information against a fourth-degree comparator; a phase mismatch can invalidate the physical interpretation of an interval; and unchanged nuisance rank does not imply stable nuisance geometry. These are reproducible properties of the stored analysis, with the derivative uncertainty left visible rather than folded into an empirical bound.

**Conjectural.** A detected internal-state response would additionally need a physically established coupling, timescale and comparator restriction. The conditional white-dwarf calculation retained in SI Section 4.6 does not close that gap: it concerns a different pair-response pattern and a transient signal, its thermal-lag sums omit a deeper-layer contribution for which no bound was established, and no corresponding timing response was computed. It is background on the physical matching problem, not evidence validating Section 5 or the method's acceptance conditions.

**Conjectural.** The present result is therefore a methods contribution with exact analytic controls and a conditional negative application. An empirical SEP statement requires a reliable derivative space, adequate nonlinear and noise inference, and completed physical matching. No amount of coefficient-space separation by itself supplies those missing premises.

## Acknowledgements {-}

The author received no external funding for this research and declares no competing interests. J.K. conceived and directed the study, set its questions and acceptance criteria, and takes responsibility for the work.

### AI disclosure {-}

This work used AI agents substantively in research and manuscript preparation. Claude Opus 5.5 (Anthropic, through Claude Code) and GPT-6-Astra (OpenAI, through Codex) wrote and ran most of the analysis, stellar-structure and verification code, carried out the recorded computations, and drafted and revised scientific claims and explanations. Those models and Claude Fable 5.1 (Anthropic, as a Claude Code subagent) assisted internal review of earlier versions. OpenAI Codex also assisted the present methods-focused revision and exact-control implementation. The author set the questions, acceptance criteria and claim labels, and checked earlier outputs against numerical records, symbolic identities and verification scripts. This revision adds exact rational controls and a replay of frozen numerical results; their successful execution does not certify the remaining derivative, noise or physical-matching assumptions. Section 4.3 describes the computational scope. The author is responsible for the interpretation and conclusions; AI tools are not authors or independent empirical validators.

## Data and code availability {-}

The public timing release is cited as a dataset \cite{voisin2025release}. Supplementary Code contains the exact controls, a compact copy of the stored comparator results, and the scripts and records for the conditional six-coefficient replay. Its README gives the environment, commands, expected outputs and the distinction between analytic controls and empirical validation. Supplementary Information preserves the detailed derivations and earlier audits, with the historical section numbering used above.

The earlier research sources and audit history are available at \url{https://github.com/lpaiu-cs/self-mass-unobservability}, frozen tag \path{grg-submission-2026-09-28}. That tag identifies the earlier version, not this revision. Its manifests identify roughly 50 GB of runtime arrays that are not deposited because of hosting limits and are available from the author on request. The compact replay does not require those arrays, but cannot revalidate their derivatives or the full inference pipeline. Historical commit identifiers in SI refer to the full local history; the public repository is a snapshot. The present source files and supplemental records accompany this submission.
