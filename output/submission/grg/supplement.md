# Online Resource 1: Supplementary Information for ``Identifying a relaxing internal state in free-fall timing: finite-frequency boundaries and an application to PSR J0337+1715''

**Author:** Juneyoung Kim\thanks{Contact author: \href{mailto:lpaiu.cs@gmail.com}{lpaiu.cs@gmail.com}}\\{\small Independent researcher, Seoul, Republic of Korea}
**Date:** 28 September 2026

## Abstract

This Supplementary Information (Online Resource 1, abbreviated SM) is the complete technical account of the paper named in the title. It is the full earlier version of the work. Apart from this title and abstract, its text is unchanged except for corrections made in this version: the periastron convention of the J0337 drive (Sections 4.4, 5.6, 5.9 and 5.10, Figures 3 and 4, and the scale list of Section 4.6), wording aligned with the main text in Sections 3.3, 3.5 and 5.3 and in the Use of AI tools section, the thermal-scan cutoff and uncomputed tail in Section 4.6, and the reproduction notes of the data-availability section. Main-text Sections 2 and 3.1--3.3 correspond to Sections 3.1--3.5 here, main-text Section 3.4 to Section 5.8, main-text Section 4 to Sections 4.3--4.5, main-text Section 5 to Sections 5.1--5.10, and the white-dwarf paragraph of the main-text discussion to Section 4.6. It contains the finite static operator catalog (Section 2 and Appendices A--C), the nuisance-projected rank and the pair-coupling benchmark (Sections 4.1--4.2), the complete scalar-charge matching (Sections 4.3--4.5), the worked calculation for the inner white dwarf of PSR J0337+1715 (Section 4.6), and all audits of the stored timing analysis with their provenance and reproduction commands (Sections 5.1--5.10 and Appendices D--E). Its section, equation, table and figure numbers are its own; Theorem 3 here is Theorem 1 of the main text. The status labels are defined in Section 1.

---

## 1. Introduction

**Imported from prior work.** Compact-body structure is represented in worldline effective field theory (EFT) by body-dependent coefficients and, when appropriate, explicit internal degrees of freedom. Tidal operators, sensitivities to additional fields, and violations of the strong equivalence principle (SEP) are related topics, but their parameters are not interchangeable without a specified theory and matching calculation \cite{goldberger2006eft,damour1992tensor,porto2016eft,will2014confrontation}.

**Imported from prior work.** Dynamical response is already established in this framework. Chakrabarti, Delsate and Steinhoff represent compact-object multipoles through response functions whose poles describe internal modes \cite{chakrabarti2013response}. Steinhoff et al. construct a relativistic worldline action with dynamical quadrupoles \cite{steinhoff2016dynamical}. Khalil et al. study nonadiabatic dynamical scalarization, including a monopolar internal mode \cite{khalil2022scalarization}. Introducing an internal state or identifying a response pole is therefore not, by itself, a new physical principle.

**Proven.** Three questions must be separated. Does a specified local action have finitely many independent operators at the chosen order? Can a dynamical response be reproduced on the available frequencies by a chosen local comparator? Does any difference survive the nuisance parameters of the measurement? A finite-dimensional operator space answers only the first question. It implies neither an observational null nor complete absorption into an orbital fit.

**Counterexample candidate.** We use a scalar relaxation state as a minimal benchmark for the latter two questions. The timing benchmark is a prescribed time-dependent pairwise SEP coupling. We construct a conditional scalar-charge force realization and identify its approximation domain, without claiming numerical neutron-star microphysical matching or a unique extension of static response. A separate audit tests nuisance and covariance assumptions on the stored linear arrays.

**Proven.** Our analytic contribution is an explicit separation of finite-order representation, exact functional reconstruction, finite-frequency interpolation, and nuisance-projected identifiability. An electric-tidal quotient gives a concrete static reference space. A root-counting argument gives a sharp sampling boundary for a nonzero single pole. These results do not establish a new theory of gravity or a general unobservability theorem.

**Imported from prior work.** The pulsar triple PSR J0337+1715 has provided sensitive tests of a static SEP parameter \cite{archibald2018universality,voisin2020sep}. The later planet/noise analysis and public data release provide the setting for the stored dynamic-template calculation used here \cite{voisin2025planet,voisin2025release}. We report that calculation as a conditional application and expose its nuisance-space dependence. A measured degeneracy of one static-SEP column is not a theorem about every static sensitivity.

Each scientific paragraph or result carries a status label. **Proven** denotes conditional mathematics, including numerical evaluation of explicit algebraic identities; **Imported from prior work** denotes literature or numerical experiments recorded in the repository, including the new audits reported here; **Counterexample candidate** denotes a proposed physical/model realization; **Conjectural** denotes an interpretation requiring additional work. A theorem's premises are hypotheses of its conditional statement, not independently established facts.

## 2. A finite-order static reference space

### 2.1 Scope and counting convention

**Proven.** Suppose a local, nonspinning, parity-even response is built from primitives with positive integer weights and only finitely many species of weight at most D. Assign each spatial or worldline derivative weight one. At fixed D there are finitely many decorated blocks, products, and tensor contractions. A quotient by linear reduction relations remains finite. No measured template rank enters this argument.

**Proven.** For the explicit D=4 example, impose the leading-Newtonian representation and block set

\[
E_{ij}=\partial_i\partial_j\Phi_{\rm ext},\qquad
\nabla^2\Phi_{\rm ext}=0,\qquad
\mathcal D_E=\{E,D_tE,D_t^2E,\nabla E,a\}.
\]

E has weight one, its first derivatives weight two, its second worldline derivative weight three, and the acceleration insertion a weight one. We quotient by total worldline derivatives and the leading free-fall equation of motion for a. Spatial integration by parts is not a worldline reduction. This is an explicit block-set calculation, not an unrestricted relativistic Weyl-tensor EFT.

**Proven.** This premise makes \(\partial_kE_{ij}=\partial_k\partial_i\partial_j\Phi_{\rm ext}\) fully symmetric and trace-free on every pair. It is an STF rank-three tensor, so

\[
\partial_iE_{ij}=0,\qquad
(\partial_kE_{ij})(\partial_iE_{kj})=(\partial_kE_{ij})(\partial_kE_{ij}).
\]

**Imported from prior work.** A general relativistic electric-Weyl gradient has richer gravitoelectric and gravitomagnetic constraints \cite{danehkar2022gem}. Our Newtonian benchmark and the relativistic timing application are not claimed to share an identical tidal operator basis.

### 2.2 The electric scalar quotient

**Proven. Proposition 1 (explicit static quotient).** For the block set and reductions of Section 2.1, the scalar action quotient through D=4 has five representatives:

\begin{equation}
\mathcal B_E=\{I_2,I_3,I_2^2,I_t,I_g\},\qquad
\begin{aligned}
I_2&=E_{ij}E_{ij}, & I_3&=E_{ij}E_{jk}E_{ki},\\
I_t&=(D_tE_{ij})(D_tE_{ij}), & I_g&=(\partial_kE_{ij})(\partial_kE_{ij}).
\end{aligned}
\label{eq:basis}
\end{equation}

**Proven.** Cayley-Hamilton gives \(\operatorname{tr}E^4=I_2^2/2\). Worldline differentiation gives \(E:D_tE=D_tI_2/2\), \(\operatorname{tr}(E^2D_tE)=D_tI_3/3\), and \(E:D_t^2E\equiv-I_t\). The STF rank-three gradient has one quadratic norm and vanishing trace contractions. Acceleration-bearing terms vanish in the declared equation-of-motion quotient. Appendix A completes the degree-by-degree argument and independence check.

**Proven.** These are five linearly independent operator representatives, not five algebraically independent scalar coordinates: the third is the square of the first. A nonredundant action at this order is

\begin{equation}
\delta L_A=C_{A,2}I_2+C_{A,3}I_3+C_{A,4}I_2^2
+C_{A,t}I_t+C_{A,g}I_g+O(\epsilon^5).
\label{eq:static-action}
\end{equation}

A weight-w block scales as O(\(\epsilon^w\)). Expanding a function of all five representatives creates no additional independent coefficients after reduction. Coefficients may depend on body structure, but calling them sensitivities does not remove their role as EFT operator coefficients.

### 2.3 Smooth finite-order expansion and exact reconstruction

**Proven. Proposition 2 (finite-order response).** Let \(x^1,\ldots,x^m\) be finitely many local arguments with integer weights \(w_i\ge1\), and \(x^i=\epsilon^{w_i}\bar x^i\) with bounded \(\bar x\). If f is \(C^{D+1}\) near zero, then

\begin{equation}
f(x)=\sum_{\sum_i n_iw_i\le D}
\frac{\partial^{\mathbf n}f(0)}{\mathbf n!}x^{\mathbf n}
+O(\epsilon^{D+1}).
\label{eq:jet}
\end{equation}

The sum is finite. Algebraically dependent arguments require a subsequent reduction rather than interpretation of all displayed coefficients as independent data.

**Proven.** Multivariate Taylor expansion through total degree D has remainder O(\(\|x\|^{D+1}\)). Because \(\|x\|=O(\epsilon)\), it obeys the stated bound. Any Taylor monomial discarded by the integer weighted cutoff has weight at least D+1 and also obeys that bound. Analyticity is sufficient, but unnecessary.

**Proven.** The function

\[
f_{\rm flat}(Y)=\begin{cases}e^{-1/Y^2},&Y>0,\\0,&Y\le0\end{cases}
\]

has all Taylor coefficients zero and satisfies \(f_{\rm flat}(Y)=o(|Y|^n)\) for every finite n. It defeats exact recovery from the full Taylor series but not Equation (\ref{eq:jet}) at any finite order. A threshold \(\sqrt{Y}\Theta(Y)\) instead lacks the required differentiability at the threshold. Neither example changes the specified polynomial catalog's count. Appendix B separates these failure layers.

### 2.4 What the static calculation does not identify

**Proven.** An additional independent primitive can add an operator while preserving fixed-order finiteness. A freely specifiable STF tensor X of weight one has a weight-two norm \(X:X\). On constant fields this cannot be a worldline total derivative or an operator of E alone. This limited catalog nonuniqueness result does not select the fields a physical theory contains. The optional family census in Appendix A is an algebraic extension, not part of a purely electric physical background.

**Proven.** The \(I_2\) coefficient is not generally the Nordtvedt SEP parameter. For an external point mass M,

\[
I_2=6(GM)^2/r^6,\qquad |\delta a_{I_2}|\propto r^{-7}.
\]

The ordinary constant Nordtvedt mass-ratio correction multiplies the external monopole acceleration, proportional to \(r^{-2}\). A specific matching calculation may relate coefficients, but their scalar character supplies no such relation. We therefore do not translate lunar-ranging bounds into constraints on Equation (\ref{eq:static-action}).

**Imported from prior work.** The standard SEP mass-ratio parameterization is reviewed by Will \cite{will2014confrontation}. Finite-size worldline actions and relativistic tidal extensions provide the comparator for tidal coefficients \cite{goldberger2006eft,bini2012tidal}. These are distinct applications of effective descriptions.

## 3. A one-state dynamical benchmark

### 3.1 Model and solution

**Counterexample candidate.** Let a known dimensionless drive excite a dimensionless state and scalar readout q:

\begin{equation}
\tau_\chi\dot\chi+\chi=\alpha F(t),\qquad
q(t)=c_YF(t)+c_\chi\chi(t),\qquad
\beta=\alpha c_\chi,\quad\tau_\chi>0.
\label{eq:model}
\end{equation}

The readout may correct a specified observable coupling. Assigning it to a mass or a pairwise gravitational parameter are different modeling choices. Relaxation presupposes an effective dissipative setting; it is not derived from an isolated conservative one-variable action.

**Proven.** With \(F(t)=F_0\cos\omega t\) and real \(F_0\), the settled periodic solution is

\begin{equation}
\chi(t)=\frac{\alpha F_0[\cos\omega t+\omega\tau_\chi\sin\omega t]}{1+\omega^2\tau_\chi^2},
\qquad G(i\omega)=c_Y+\frac{\beta}{1+i\omega\tau_\chi}.
\label{eq:response}
\end{equation}

The general solution also contains \(\chi_h e^{-(t-t_0)/\tau_\chi}\). The periodic analysis assumes this transient has decayed or is absent by initial conditions. An appreciable transient needs an additional amplitude outside the stored templates. For positive real pole strength the relaxation phase is \(-\arctan(\omega\tau_\chi)\); a negative strength adds a phase of \(\pi\).

### 3.2 Adiabatic and finite-frequency collapse

**Proven.** For z=\(i\omega\), the degree-N derivative expansion obeys

\begin{equation}
G(z)-\left[c_Y+\beta\sum_{n=0}^{N}(-\tau_\chi z)^n\right]
=\frac{\beta(-\tau_\chi z)^{N+1}}{1+\tau_\chi z}.
\label{eq:remainder}
\end{equation}

On \(|\omega\tau_\chi|\le\rho<1\), the error is at most \(|\beta|\rho^{N+1}\). Zero lag gives the settled readout \((c_Y+\beta)F\). Zero frequency gives a settled constant shift, but does not eliminate a transient. Zero pole strength eliminates the driven-state contribution to the settled q; an independently initialized transient can remain when \(c_\chi\) is nonzero.

**Proven.** At one known frequency, \(a_0F+a_1\dot F\) matches Equation (\ref{eq:response}) exactly with \(a_0=c_Y+\beta/(1+\omega^2\tau_\chi^2)\) and \(a_1=-\beta\tau_\chi/(1+\omega^2\tau_\chi^2)\). A quadrature at one frequency is insufficient evidence for an internal state against this derivative comparator.

**Proven.** Linear superposition of two drives through Equation (\ref{eq:model}) produces only the input frequencies. Sum/difference sidebands require nonlinearity in the drive, readout, or observable dynamics. Local nonlinear static comparators can also generate sidebands; their existence alone does not establish a state variable.

### 3.3 Exact sampling theorem

**Proven. Theorem 3 (real finite-frequency boundary).** Suppose \(c_Y,\beta\) are real, \(\beta\ne0\), \(\tau_\chi>0\), and \(\omega_1,\ldots,\omega_K\) are distinct positive frequencies at which the drive is known and nonzero, so that the values \(G(i\omega_k)\) are available. A polynomial \(P_N(z)\) of degree at most N with freely chosen shared real coefficients matches \(G(i\omega_k)\) at all K frequencies if and only if \(N\ge2K-1\).

**Proven.** Define \(R(z)=(1+\tau_\chi z)[P_N(z)-c_Y]-\beta\). Real coefficients imply that a root at \(i\omega_k\) is accompanied by \(-i\omega_k\). Matching thus gives 2K distinct roots, but R has degree at most N+1. If \(N<2K-1\), R would be identically zero, contradicting \(R(-1/\tau_\chi)=-\beta\ne0\).

**Proven.** Conversely, set

\begin{equation}
Q_K(z)=\prod_{k=1}^K(z^2+\omega_k^2),\qquad
P_{2K-1}(z)=c_Y+
\frac{\beta[1-Q_K(z)/Q_K(-1/\tau_\chi)]}{1+\tau_\chi z}.
\label{eq:interpolant}
\end{equation}

The numerator vanishes at \(-1/\tau_\chi\), giving a real polynomial of degree 2K-1 upon division. At every carrier it equals G. The construction allows unrestricted coefficients; power-counting bounds, passivity or external calibration would define a smaller comparator class.

**Proven.** The first exact obstruction to a real degree-N comparator is \(K=\lfloor(N+1)/2\rfloor+1\). Three positive carriers exclude N up to four, but a fifth-degree real comparator interpolates them. With complex coefficients and only positive-frequency samples, interpolation permits K up to N+1 and the first obstruction occurs at N+2. No finite polynomial equals a nonzero single-pole response on an open frequency interval.

**Proven.** Exact noninterpolation does not guarantee detectability. Nearly coincident carriers or small \(\omega\tau_\chi\) can make the difference arbitrarily small. The theorem concerns equality on specified samples; precision and nuisance geometry are additional requirements.

### 3.4 Physically specified comparator restrictions

**Proven.** A real local conservative quadratic action \(\frac12\int F P(D)F\,dt\), with \(D=d/dt\) and vanishing boundary variations, contributes the self-adjoint response \([P(D)+P(-D)]/2\). Odd derivatives cancel. An even polynomial is therefore justified for this particular time-reversal-invariant conjugate response. Dissipation, nonconjugate readout, time-dependent background or unknown drive phase can invalidate this restriction. It is not a condition on every static EFT or timing nuisance parameter.

**Counterexample candidate.** A different comparator class follows from positive quadratic state energy, reciprocal forcing/readout and positive Rayleigh dissipation:

\begin{equation}
\boldsymbol\Gamma\dot{\boldsymbol q}+\boldsymbol K\boldsymbol q=\boldsymbol b F,
\qquad q_{\rm out}=\boldsymbol b^T\boldsymbol q+c_0F.
\label{eq:fast-gradient}
\end{equation}

Both matrices are real symmetric positive definite. Impose a lower bound \(\Lambda\) on every eigenvalue of \(\boldsymbol\Gamma^{-1/2}\boldsymbol K\boldsymbol\Gamma^{-1/2}\). This describes modes that relax faster than a specified scale; the gap is a hypothesis requiring independent physical matching.

**Proven.** Orthogonal diagonalization gives

\begin{equation}
H(z)=c_0+\sum_j\frac{a_j}{1+z\tau_j},\qquad
a_j\ge0,\quad 0<\tau_j\le\Lambda^{-1}.
\label{eq:positive-spectrum}
\end{equation}

For \(\omega_h>\omega_l>0\), put \(S(\omega)=-\operatorname{Im}H(i\omega)/\omega\). Then

\begin{equation}
S(\omega_l)\le R_\Lambda S(\omega_h),\qquad
R_\Lambda=\frac{1+(\omega_h/\Lambda)^2}{1+(\omega_l/\Lambda)^2}.
\label{eq:gap-witness}
\end{equation}

Indeed \((1+\omega_h^2\tau^2)/(1+\omega_l^2\tau^2)\) increases with positive \(\tau\); summing its bound with the nonnegative weights \(a_j\tau_j/(1+\omega_h^2\tau_j^2)\) proves the inequality. A positive single pole with \(\tau_\chi>\Lambda^{-1}\) violates it. The real instantaneous coefficient cancels from S. The argument also applies to a convergent nonnegative distribution of relaxation times with the same support.

**Proven.** This supplies a comparator restriction tied to energy, reciprocity and a rate gap rather than an arbitrary Taylor cutoff. Negative residues, a different drive/readout pair, oscillatory modes or removal of the gap can evade it. Rejecting this class would not uniquely identify one internal state: multiple slow states or a slow memory distribution can also violate the inequality. Neither the rate gap nor reciprocal readout has been empirically established for J0337 here.

### 3.5 Observable-pole uniqueness and its finite-precision boundary

**Proven.** Suppose \(H(i\omega)=c_0+\int a(d\tau)/(1+i\omega\tau)\), with real \(c_0\) and a nonnegative finite measure on \(\tau>0\). At two calibrated frequencies \(0<l<h\), let \(S_\omega=-\operatorname{Im}H(i\omega)/\omega\) and define

\begin{equation}
R=\frac{h^2S_h-l^2S_l}{h^2-l^2},\qquad
D=\frac{\operatorname{Re}H(il)-\operatorname{Re}H(ih)}{h^2-l^2},\qquad
Q=\frac{S_l-S_h}{h^2-l^2}.
\label{eq:positive-moments}
\end{equation}

These are the zeroth, first and second moments of \(d\mu=\tau a(d\tau)/[(1+l^2\tau^2)(1+h^2\tau^2)]\). Therefore

\begin{equation}
RQ-D^2=R^2\operatorname{Var}_{\mu/R}(\tau)\ge0.
\label{eq:positive-single-pole}
\end{equation}

For \(R>0\), equality holds exactly for one observable relaxation time \(\tau_0=D/R\). Its amplitude is \(a_0=R(1+l^2\tau_0^2)(1+h^2\tau_0^2)/\tau_0\), and either real part fixes \(c_0\). A relaxation-time atom at \(\tau=0\) would be indistinguishable from \(c_0\). This assumes exact common drive/readout calibration; it is not an observed equality in J0337.

**Proven.** One observable pole does not count hidden, unexcited or degenerate internal variables. Two positive poles at \(\tau_0\pm\epsilon\), each of weight \(a_0/2\), differ from one pole by

\begin{equation}
\frac{a_0}{2}\left[\frac{1}{1+s(\tau_0-\epsilon)}+\frac{1}{1+s(\tau_0+\epsilon)}\right]-\frac{a_0}{1+s\tau_0}
=\frac{a_0s^2\epsilon^2}{(1+s\tau_0)[(1+s\tau_0)^2-s^2\epsilon^2]}.
\label{eq:close-poles}
\end{equation}

For \(0<\epsilon<\tau_0\) both modes are stable, and the difference vanishes quadratically. A finite-noise test cannot uniformly distinguish exactly one pole from arbitrarily close two-pole alternatives without a separation or weight condition. Without positivity, the stable strictly proper addition \(\eta\prod_k(s^2+\omega_k^2)/(s+\lambda)^{2K+1}\), \(\lambda>0\), vanishes at every stored carrier but changes the response elsewhere. It need not be a positive reciprocal spectrum.

## 4. From response functions to an observable

### 4.1 Nuisance-projected rank

**Proven.** Linearize a measurement as \(y=J\theta+T\beta+\varepsilon\), with positive-definite covariance C. Whiten by \(C^{-1/2}\) and let P project onto \(\operatorname{col}(\widetilde J)\). The information for an unconstrained signal amplitude is

\begin{equation}
\mathcal I_\beta=\widetilde T^T(1-P)\widetilde T.
\label{eq:projection}
\end{equation}

Local structural identifiability requires \(\mathcal I_\beta>0\): T must add rank beyond J. A co-fitted instantaneous direction must also enter the nuisance projection. Finite dimension alone provides no lower bound on this information.

**Proven.** A finite shared projection can erase a signal when its tangent directions span it. Arbitrary independent complex projections \(\Lambda_k\) give a stronger exact degeneracy: for nonzero \(G_kF_k\), the choice \(\Lambda_k=O_k/(G_kF_k)\) fits any carrier observations. Survival of a common pole requires more than counting state parameters.

### 4.2 A transparent SEP benchmark

**Counterexample candidate.** Define a phenomenological pair potential on pulsar-companion pairs,

\begin{equation}
V_{pj}(t)=-\frac{Gm_pm_j}{r_{pj}}[1+\Delta_0+q(t)],
\label{eq:pair}
\end{equation}

with fixed inertial masses and the prescribed drive in Equation (\ref{eq:model}). Equal-and-opposite pair forces acquire the same coupling factor in the respective accelerations. This defines the modified pair term used to interpret the stored dynamic columns, not a replacement for the additional dynamics of the full timing code.

**Proven.** For externally prescribed q(t), the potential is explicitly time-dependent and mechanical energy can be exchanged with the drive. If F is instead a field or position-dependent function varied in an action, additional derivatives and backreaction must be included. A field-dependent inertial mass also changes the equations beyond rescaling a pair potential. Equation (\ref{eq:pair}) is therefore a specified phenomenological readout, not a general derivation from \(m_A(Y,\chi)\).

### 4.3 A conditional scalar-charge realization

**Imported from prior work.** Dynamical scalar-charge worldline models contain a charge potential and mode inertia; compact-body matching and monopole radiation damping have been developed by Khalil et al. \cite{khalil2022scalarization}. We use that framework as prior art for the following restricted reduction.

**Counterexample candidate.** In units G=c=1, hold the companion charges fixed, retain leading constant orbital masses, and let only the pulsar charge \(Q_p\) evolve. Consider

\begin{equation}
L=\sum_A\frac{m_Av_A^2}{2}+\frac{I\dot Q_p^2}{2}-V(Q_p)
+\sum_{A<B}\frac{m_Am_B+Q_AQ_B}{r_{AB}},\qquad
\mathcal R=\frac{\Gamma\dot Q_p^2}{2},
\label{eq:charge-action}
\end{equation}

with positive I and \(\Gamma\). A constant scalar background is included in V. The Rayleigh function \(\mathcal R\) specifies dissipation. This leading orbital/state model omits higher post-Newtonian and radiation terms of a complete timing theory.

**Proven.** Variation at independent positions and charge gives reciprocal pair forces with \(\Delta_{pj}=Q_pQ_j/(m_pm_j)\), and

\begin{equation}
I\ddot Q_p+\Gamma\dot Q_p+V'(Q_p)=\sum_j\frac{Q_j}{r_{pj}},
\qquad \dot E_{\rm orbit+state}=-\Gamma\dot Q_p^2.
\label{eq:charge-eom}
\end{equation}

No derivative of an already substituted Q(r) is added: Q is an independent coordinate during variation. Linearizing a stable branch with \(\kappa=V''(Q_0)>0\) gives transfer \(H=1/(\kappa-I\omega^2+i\Gamma\omega)\). For \(H_0=1/(\kappa+i\Gamma\omega)\),

\begin{equation}
\epsilon_I=\frac{I\omega^2}{|\kappa+i\Gamma\omega|}<1
\quad\Longrightarrow\quad
\left|\frac{H-H_0}{H_0}\right|\le\frac{\epsilon_I}{1-\epsilon_I}.
\label{eq:inertial-bound}
\end{equation}

The error bound must be small at every retained carrier. Separated overdamped roots require \(I\kappa/\Gamma^2\ll1\); neglected homogeneous modes must also have decayed. Dissipation alone supplies neither condition.

**Proven.** If both companion charge/mass ratios equal \(a_w\), define \(\delta U=\sum_jm_j[1/r_{pj}-\langle1/r_{pj}\rangle]\). Then \(\delta\varphi=a_w\delta U\), and the common pulsar-pair response is

\begin{equation}
\delta\Delta(\omega)=\frac{B\,\delta U(\omega)}{1+i\omega\tau_\chi},\qquad
B=\frac{a_w^2}{\kappa m_p},\qquad \tau_\chi=\frac{\Gamma}{\kappa}.
\label{eq:matching}
\end{equation}

For \(F=\delta U/U_*\), the benchmark coefficient is \(\beta=BU_*\). This stable equal-charge branch has \(B\ge0\), whereas the generic signed-beta fit admits a larger model. The minimal slow-mode realization has \(c_Y=0\); co-fitting a free instantaneous coefficient allows a fast response or a broader comparator.

**Proven.** If \(a_i\ne a_o\), the pair responses differ by \((a_i-a_o)\delta Q_p/m_p\), obstructing a common nonzero modulation. A responsive companion with susceptibility \(C_j\) produces a leading feedback stiffness shift \(-\sum_j C_j/r_{pj}^2\). Neglect requires this small relative to \(\kappa\). Linearization also requires small nonlinear force ratios, for example \(|V^{(3)}(Q_0)\delta Q|/(2\kappa)\ll1\) and \(|V^{(4)}(Q_0)\delta Q^2|/(6\kappa)\ll1\). Large susceptibility near an instability does not guarantee any of these inequalities. A first-order pole has no resonance peak, and \(\Gamma/\kappa\) is not a model-independent inverse particle mass.

### 4.4 Physical-drive phases provide an additional gate

**Proven.** Let \(\boldsymbol r=\boldsymbol x_p-\boldsymbol x_i\), \(\boldsymbol R=\boldsymbol x_b-\boldsymbol x_o\), b be the inner center of mass, and \(f=m_i/(m_p+m_i)\). Then \(r_{po}=|\boldsymbol R+f\boldsymbol r|\). In the aligned coplanar limit, retain terms linear separately in eccentricity and \(f a_{\rm in}/a_{\rm out}\), dropping their products. With mean longitudes \(\lambda_p=n_{\rm in}(t-t_{{\rm asc},p})\), \(\lambda_b=n_{\rm out}(t-t_{{\rm asc},b})\), and \(M_A=\lambda_A-\varpi_A\),

\begin{equation}
\begin{split}
\delta U={}& A_{\rm in}\cos M_p+A_{\rm out}\cos M_b
-A_{\rm dif}\cos(\lambda_p-\lambda_b)+\cdots,\\
A_{\rm in}={}&\frac{m_i e_{\rm in}}{a_{\rm in}},\quad
A_{\rm out}=\frac{m_o e_{\rm out}}{a_{\rm out}},\quad
A_{\rm dif}=\frac{m_o f a_{\rm in}}{a_{\rm out}^2}.
\end{split}
\label{eq:physical-drive}
\end{equation}

The minus sign follows from the first-order expansion of \(1/|\boldsymbol R+f\boldsymbol r|\). Higher eccentricity, mixed and inclination terms require separate control.

**Imported from prior work.** The published timing convention has \(t_{\rm asc}=t_{\rm pericenter}-P\varpi/(2\pi)\) \cite{voisin2025planet}. In the released timing code the parameters eta and kappa are \(e\sin\varpi\) and \(e\cos\varpi\), the ELL1 convention \cite{lange2001ell1}; the frozen parameters then place the pericenters close to those of Ref. \cite{ransom2014triple}. An earlier version of this analysis read eta as \(e\cos\varpi\); the values below and in Sections 5.6, 5.9 and 5.10 are corrected in this version.

**Proven.** For positive carrier amplitudes the phase closure

\begin{equation}
\mathcal C=\phi_{\rm dif}-\phi_{\rm in}+\phi_{\rm out}
=\pi+\varpi_p-\varpi_b\pmod{2\pi}
\label{eq:closure}
\end{equation}

is invariant under a common origin shift because \(n_{\rm dif}=n_{\rm in}-n_{\rm out}\). The frozen parameters give \(\varpi_p=1.69406454\), \(\varpi_b=1.67084424\) and \(\mathcal C=3.16481295\) radians. The archived auxiliary potential-drive dictionary instead has \(\mathcal C=0\). Its inner, outer and difference phases differ from Equation (\ref{eq:physical-drive}) by 0.12326821, 0.10004792 and \(\pi\) radians. No common time shift repairs the mismatch. The unit-drive family also has zero closure; its envelope does not contain this leading physical family.

**Counterexample candidate.** Equations (\ref{eq:charge-action})--(\ref{eq:matching}) supply conditional force-level matching. They do not establish its realization in J0337. We withdraw the archived auxiliary physical-beta numbers as constraints on this realization. A corrected physical analysis must prescribe the unequal amplitudes, phases, neglected harmonics and additional forces, then validate its inference. Rescaling the unit-drive interval cannot do this. Numerical EOS-to-body matching and a complete astrophysical likelihood remain outside the present result.

### 4.5 What equilibrium matching cannot determine

**Proven.** In the admitted damped-charge EFT class, fixing the potential, equilibrium charge and drive coupling fixes the static susceptibility \(1/\kappa\), but the family \(\Gamma=\kappa\tau_0\), \(\tau_0>0\), leaves that equilibrium unchanged and permits arbitrary relaxation time. Positive inertia can vary independently of equilibrium as well. Equilibrium information alone therefore does not establish a fast-rate gap. A specified microscopic theory may relate these parameters; this statement does not make damping free within every fixed theory.

**Imported from prior work.** Khalil et al. obtain potential coefficients from stellar equilibrium sequences and kinetic coefficients using scalar-led modes \cite{khalil2022scalarization}. Their leading coupled monopole damping force is proportional to the negative sum of charge velocities. **Proven.** Its two-charge all-ones damping matrix has a null direction, so the positive-definite gradient comparator of Section 3.4 is not automatic. Additional inertia, dissipation or a justified reduction must be considered. **Conjectural.** Numerical J0337 matching still requires a gravity coupling, EOS/branch, companion charges and dynamical response data not fixed by the stored orbital parameters.

### 4.6 A worked white-dwarf calculation: a transient retarded scattering signal and orbital-timescale bounds

**Imported from prior work.** Spectroscopy gives the inner white dwarf of PSR J0337+1715 \(T_{\rm eff}=15{,}800\pm100\) K, \(\log g=5.82\pm0.05\), \(R=0.091\pm0.005\,R_\odot\) and a hydrogen (DA) atmosphere \cite{kaplan2014j0337}. Masses and semimajor axes below are those of the Section 5.9 parameter set (inner white dwarf \(0.19754\,M_\odot\)). Orbital eccentricities are from \cite{ransom2014triple}.

**Counterexample candidate.** We examine one internal state of this star: the scalar-driven displacement of its own matter. From Section 4.3 we reuse the companion charge-to-mass ratios \(a_i\) (inner white dwarf) and \(a_o\) (outer white dwarf), and we write \(a_p=Q_p/m_p\) for the pulsar. All other symbols are local to this section.

- **Theory.** The declared theory has one massless scalar with \(A(\varphi)=e^{\beta_s\varphi^2/2}\) and \(\beta_s=-4\), so the matter coupling is \(\alpha_s(\varphi)=\beta_s\varphi\). The asymptotic value \(\varphi_\infty=10^{-3}\) is a linearization point, and \(\alpha_0=\alpha_s(\varphi_\infty)\).
- **Background.** The static background metric is \(ds^2=-\mathcal N^2c^2dt^2+\mathcal L^2dr^2+r^2d\Omega^2\). The optical coordinate is \(x=-\int_r^{R_d}(\mathcal L/\mathcal N)\,dr\), with \(x'=dx/dr\). The background scalar is \(\varphi_0(r)\), the density \(\rho_0(r)\), and \(T_0\) is the trace of the background matter stress-energy tensor. The field variable is \(\Psi=r\,\delta\varphi\).
- **Masses.** \(M_{\rm wd}\) is the white dwarf's mass and \(m_{\rm wd}=GM_{\rm wd}/c^2\) its geometric mass.

**Counterexample candidate.** The background is a static weak-field model with a helium core and a hydrogen-rich envelope (\(X\simeq0.89\)). It has \(M_{\rm wd}=0.198\,M_\odot\), \(R=6.91\times10^9\) cm, \(T_{\rm eff}=15{,}898\) K and \(\log g=5.740\), with a gray Eddington envelope built with Rosseland tables. An ingoing pulse drives the coupled evolution of baryons, gas, photon transport and one returned GR metric until \(t_e=2t_p\):
\[
\Psi_{\rm in}=\eta R_d\,\mathcal P\big((t+x/c)/t_p\big),\qquad \mathcal P(s)=(4s(1-s))^4\ (0<s<1),\ 0\ \text{otherwise},
\]
with \(t_p=1.717\) ms, \(\eta=10^{-30}\) and \(R_d=6.91\times10^9\) cm. The compact readout at \(R_d\) is \(\mathcal Q_c=-\Psi_{\rm out}/m_{\rm wd}\), and the readout is
\[
\mathcal Q=\mathcal Q_c-a_i^{(0)}\,\delta m_{\rm wd}/m_{\rm wd},
\]
where \(\delta m_{\rm wd}\) is the change of the exterior mass and \(a_i^{(0)}\) the background charge per unit mass, for which we use its weak-field value \(\alpha_0\). \(\mathcal Q\) is the matter-response part of the change of \(a_i\). It excludes the instantaneous polarization \(\beta_sT_0\delta\varphi\), which belongs to the wave operator; in the long-wavelength limit that polarization gives the instantaneous response \(\delta a_i=\beta_s\delta\varphi\) discussed below. The pulse length \(ct_p\simeq515\) km is far below \(R\), so this is a finite-size response outside the point-particle regime \(\lambda\gg R\).

**Imported from prior work.** For orientation only, since the two numbers use different norms: the stored maximum ratio of the first-return Born field to the incident field is \(1.7\times10^{-10}\), while the endpoint readout below corresponds to an outgoing amplitude of about \(10^{-26}\) of the incident field at \(R_d\).

**Imported from prior work.** At the protocol endpoint the readout is negative. We extrapolate \(\mathcal Q_c\) over interior refinements by 1, 2 and 4 with the Richardson method (observed order 2.00; change from 2\(\times\) to 4\(\times\) of \(-1.2\%\)). The mass term, measured on the unrefined solution, multiplies it by \(1+2.2\times10^{-5}\). The result is

\begin{equation}
\mathcal Q(t_e)=-2.3\times10^{-51}.
\label{eq:wd-readout}
\end{equation}

Its composition is as follows.

- The direct part, the instantaneous product of the pulse and the background geometry, contributes \(10^{-8}\) of it.
- The rest-mass redistribution of baryons contributes 0.999996, and 99.9% comes from layers the pulse has already left.
- With the background fixed, disabling the Planck–Larkin occupation probabilities in the equation of state and opacities changes it by at most a relative \(2.2\times10^{-10}\).

**Conjectural.** With the \(\varphi_\infty^2\) scaling discussed below, Equation (\ref{eq:wd-readout}) corresponds to \(\mathcal Q(t_e)/(\eta\varphi_\infty^2)\simeq-2.3\times10^{-15}\).

**Proven.** For pressureless matter at linear order, the radial displacement

\begin{equation}
\xi=4c^2\frac{\mathcal N^2}{\mathcal L^2}\eta R_dt_p^2\Big[\Big(\frac{\varphi_0'}{r}-\frac{\varphi_0}{r^2}\Big)\mathcal P_2(s)+\frac{\varphi_0x'}{rct_p}\mathcal P_1(s)\Big],\qquad s=\frac{t+x/c}{t_p},
\label{eq:free-fall}
\end{equation}

with \(\mathcal P_1(s)=\int_0^s\mathcal P\) and \(\mathcal P_2(s)=\int_0^s\mathcal P_1\), solves \(\ddot\xi=-c^2(\mathcal N^2/\mathcal L^2)\partial_r[\alpha_s(\varphi_0)\delta\varphi]\), with \(\xi=\dot\xi=0\) before the pulse arrives. **Imported from prior work.** The resulting face mass fluxes \(4\pi r^2\mathcal L\rho_0\xi\), passed through the retarded propagator of the coupled calculation, reproduce its continuum limit to \(5.2\times10^{-5}\).

**Proven.** In flat space the scattered outgoing monopole of spherical shells depends only on the retarded time \(u=t-r/c\). A shell at radius \(r_e\) carrying the mass perturbation \(\delta M_e(t)\) contributes

\begin{equation}
\mathcal Q_c(u)=\frac{1}{M_{\rm wd}}\sum_e w_e\,\frac{c}{2r_e}\int_{u-r_e/c}^{u+r_e/c}\delta M_e(t')\,dt',\qquad w=\alpha_s(\varphi_0)\,\mathcal N,
\label{eq:wd-window}
\end{equation}

whose static limit is \(\sum_ew_e\delta M_e/M_{\rm wd}\). A mass-conserving redistribution therefore contributes only in two ways: through the radial variation of \(w\), and through light-crossing windows that are not yet complete.

**Imported from prior work.** We also evolved a whole-star linear adiabatic radial model with Newtonian self-gravity. It is driven by the scalar force behind Equation (\ref{eq:free-fall}), including the regular pulse reflected through the center, and it evaluates \(\mathcal Q_c\) through Equation (\ref{eq:wd-window}); the mass term is not applied to its history.

- It reproduces the endpoint readout to \(-0.81\%\) and the free-fall displacement in the charge layers to \(1.2\times10^{-6}\).
- In this model the static-window value at \(t_e\) is \(+3.6\times10^{-54}\), of opposite sign and 1/640 of the magnitude.
- The readout \(\mathcal Q_c\) rises steeply after \(t_e\): \(-2.8\times10^{-51}\) at 3.5 ms and \(-8.7\times10^{-51}\) at 4.0 ms.
- In source time, the pulse reaches the center at 0.2305 s and leaves the star at 0.4627 s.

**Imported from prior work.** In readout time \(t_R\) at \(R_d\), the model gives the following history. Values are interpolated from the stored samples (0.25 ms apart up to 20 ms, 2 ms apart after); the peak comes from the run's 0.2 ms grid.

- The readout stays negative until \(t_R\simeq0.43\) s. It grows from \(-2.3\times10^{-51}\) at \(t_e\) to \(-3.0\times10^{-43}\) at 0.1 s, \(-8.0\times10^{-41}\) at 0.23 s and \(-2.7\times10^{-39}\) at 0.35 s.
- The signals from the center arrival and from the outgoing pass reach the readout together near \(2R_d/c\simeq0.46\) s. There \(|\mathcal Q_c|\) peaks at \(3.5\times10^{-36}\), an outgoing amplitude of \(1.5\times10^{-11}\) of the incident one.
- After the pulse has left, the star rings in radial p modes, with dominant periods of 44–59 s and an rms readout of \(8.1\times10^{-45}\).

**Imported from prior work.** Only the endpoint is cross-checked against the coupled calculation. The model is linear, adiabatic and Newtonian. The photosphere, whose thermal times are at most 0.4 s (about 1 s by an estimate that includes ionization energy), is not adiabatic during the ringing.

**Conjectural.** The endpoint readout is therefore the retarded monopole of a transient redistribution, Equation (\ref{eq:wd-window}), and its value is tied to the protocol. It is not a static charge of the displaced configuration.

**Proven.** Consider a discrete linear system with positive-definite mass and stiffness matrices, driven by a force of finite duration. After the force ends, the response is a finite sum of modes of nonzero frequency. Its time average vanishes, and it tends to zero if every mode is positively damped. **Imported from prior work.** The radial model's lowest eigenvalue is \(\omega_0^2=1.25\times10^{-3}\) s\(^{-2}\), a period of 178 s. **Conjectural.** Suppose the star has no neutral modes and all its modes are damped, including the thermal (entropy) modes that the adiabatic model omits; we made no non-adiabatic stability calculation. Then the pulse-induced change of \(a_i\) vanishes at first order in \(\eta\) as \(t\to\infty\). Thermal modes can decay as slowly as the cooling time, so a first-order residual may persist over any observing span. A change of order \(\eta^2\) from absorbed energy remains.

The endpoint value is subject to the following corrections.

- **Proven.** In the Newtonian radial channel, self-gravity feedback changes the displacement over \(t_e\) by at most a fraction \(2\pi G\rho_{\max}t_e^2=4.5\times10^{-18}\), where \(\rho_{\max}\) is the largest \(\rho_0\) in the evolved layers. **Conjectural.** With post-Newtonian terms neglected, the once-returned GR solution differs from the self-consistent one at this order.
- **Conjectural.** Curvature backscatter changes the flat-space readout at order \(m_{\rm wd}/R_d=4.2\times10^{-6}\). This is an order estimate, not a bound.
- **Proven.** Because \(\alpha_s\) is linear in \(\varphi\), the second-order force and source terms are at most \(|\delta\varphi|/\varphi_0=1.0\times10^{-27}\) of the first-order ones. This compares those terms only; it does not bound the full readout error.
- **Proven.** Invariance under \(\varphi\to-\varphi\) makes \(\mathcal Q/\eta\), at first order in \(\eta\), even in \(\varphi_\infty\). **Conjectural.** Since the force and the readout weight are both proportional to \(\alpha_s(\varphi_0)\), the dynamic part scales as \(\varphi_\infty^2\) at leading order in the weak background field. The exponent of the \(10^{-8}\) direct part is not fixed.
- **Imported from prior work.** The ADM mass differs from the mass at the readout radius by a relative \(6.8\times10^{-11}\).

**Imported from prior work.** The photosphere sets the magnitude of the readout.

- The kernel centroid lies 363 km below the \(P_{\rm gas}=1\) dyn cm\(^{-2}\) cut, at gray optical depth 1.06.
- The kernel is nearly single-signed. With the face kernel held fixed (response coefficients and propagator unchanged), no relative change of the face densities with \(\|\delta\rho/\rho\|_\infty\) below 99.87% reverses the sign, and the tested smooth structure families keep it.
- A spherical gray envelope rebuilt with the same equation of state and tables reproduces the declared density with a kernel-weighted deviation of 0.25% (maximum 0.38%).
- At fixed mass, and with the pulse kernel of the declared geometry, varying \(\log g\) or \(T_{\rm eff}\) alone within Kaplan's 1\(\sigma\) ranges scales \(|\mathcal Q(t_e)|\) by 1.6–6.2 (3.2 at the central values).
- An LTE non-gray correction multiplies this by a further 1.9–2.4 over those 1\(\sigma\) cases (2.1 at the central values). It uses H and He continuum opacity whose Rosseland mean is 0.55–0.74 of the tables, is plane-parallel, and is converged to \(5.6\times10^{-5}\) in flux for \(\tau_R\le100\).

**Conjectural.** These magnitudes only show how strongly the protocol readout depends on surface gravity; no conclusion below uses them.

**Imported from prior work.** At orbital timescales, the white dwarf's radial (\(l=0\)) restoring times, below 180 s, are far shorter than the lags \(\tau_\chi=2\)–500 d scanned in Section 5.

- **Imported from prior work.** *Radial.* The lowest radial mode has a period of 178 s, and the acoustic cutoff periods of the charge layers are 28–47 s. At the inner orbit the in-phase correction is \((\omega/\omega_0)^2=1.6\times10^{-6}\). At equal amplitude, a drive longer than the star exerts a monopole force about \(3\times10^{-8}\) of the pulse's through the \(\varphi_0'\) term. **Conjectural.** Bounding instead the full change of the effective gravity under a uniform scalar shift, at most 0.048 per unit shift, gives a ratio of order \(10^{-7}\).
- **Proven.** *Tidal.* On a static spherical background, linear response couples \(a_i\) only to the \(l=0\) part of the drive, so tides do not enter it at linear order.
- **Imported from prior work.** *Tidal magnitudes.* The equilibrium tide has apsidal constant \(k_2=3.3\times10^{-4}\) and relative tidal force \(9.1\times10^{-12}\); its apsidal rate is \(1.5\times10^{-5}\) of the 1PN rate.
- **Conjectural.** *Second order and resonance.* At second order the tide has an \(l=0\) part. Treating the squared tidal parameter as a uniform change of the effective gravity gives a modulation of \(a_i\) of \(3.7\times10^{-19}\) at \(\varphi_\infty=10^{-3}\) for the equilibrium tide, scaling as \(\varphi_\infty\). A one-way radiative damping-depth estimate for the dynamical tide (\(8.9\times10^6\)) uses a weak-damping formula outside its validity, so discrete g-mode resonances, which could raise this, are not excluded.
- **Imported from prior work.** *Structural and thermal.* The adiabatic static structural susceptibility \(\mathcal S_{\rm struct}\) is the change of \(a_i\) per uniform scalar shift from hydrostatic readjustment. We find \(|\mathcal S_{\rm struct}|\le8.84\times10^{-9}\) at \(\varphi_\infty=10^{-3}\), at most \(2.21\times10^{-9}\) of \(|\beta_s|\), scaling as \(\varphi_\infty^2\). The thermal time of the overlying layers is only 4 s at the envelope base, 624 km deep (up to about 2.5 times longer by the ionization-energy estimate).
- **Imported from prior work.** *Deeper layers.* We estimated the thermal response layer by layer. Each layer relaxes on the thermal time of the layers above it, and its fully relaxed limit is bracketed by the isothermal exponent \(P_{\rm gas}/P\). Solving the static structural response again with all layers faster than a given time relaxed gives the relaxation strength of each deeper shell. The layers whose thermal time equals \(1/\omega\) lie about 2,600 km deep for the inner orbit and 6,600 km for the outer one, above about \(2\times10^{-9}\) and \(1.4\times10^{-7}\) of the mass. The scan uses 161 thermal-time cuts from \(10^{-3}\) to \(10^{13}\) s. The last cut encloses the outer 1.70 percent of the mass, whereas the central thermal time is \(2.01\times10^{15}\) s. Summing the absolute strengths of the scanned shells with the Debye weight \(\omega\tau/(1+\omega^2\tau^2)\) gives partial lag estimates of \(4.0\times10^{-9}\,|\mathcal S_{\rm struct}|\) at the inner orbital frequency and \(3.3\times10^{-7}\,|\mathcal S_{\rm struct}|\) at the outer one. **Conjectural.** This estimate is not a non-adiabatic calculation, and its relaxed limit can differ from thermal equilibrium by factors of order unity. It does not assume a single pole or the strength bound below, but it leaves the relaxation strength of layers with \(\tau>10^{13}\) s uncomputed.
- **Proven.** *Uncomputed tail.* Within this shell model, let \(T_{\rm tail}=\sum_{\tau_j>\tau_{\rm cut}}|s_j|/|\mathcal S_{\rm struct}|\), with \(\tau_{\rm cut}=10^{13}\) s and a finite sum. Since \(\omega\tau/(1+\omega^2\tau^2)\le1/(\omega\tau)\), the tail quadrature obeys \(|Q_{\rm tail}|/|\mathcal S_{\rm struct}|\le T_{\rm tail}/(\omega\tau_{\rm cut})\). **Imported from prior work.** The multiplying factors are \(2.24\times10^{-9}\) and \(4.50\times10^{-7}\) at the inner and outer frequencies. **Conjectural.** No bound on \(T_{\rm tail}\) was computed. A total lag below one percent would follow, in this model, from \(T_{\rm tail}\le2.22\times10^4\), but this is an additional condition rather than a measured property. The partial sums alone therefore do not bound the total lag.
- **Conjectural.** *Lag bound.* Write the thermally relaxing part of the structural response to a drive at angular frequency \(\omega\) as \(\Delta\mathcal S/(1+i\omega\tau)\). Two assumptions, neither computed, bound it: \(|\Delta\mathcal S|\le|\mathcal S_{\rm struct}|\), and a single relaxation time \(\tau\) (Debye form). Its quadrature (lagged) part is then at most \(|\Delta\mathcal S|/2\le|\mathcal S_{\rm struct}|/2\) per unit drive at every frequency. The bound limits the size of the lag; it does not remove it. Responses that violate either assumption, such as resonant, two-pole (inertial) or mixed-sign multi-relaxation dynamics, fall outside it.

**Proven.** A change \(\delta a_i\) modulates the pulsar–inner pair factor by \(a_p\,\delta a_i\) and the inner–outer pair factor by \(a_o\,\delta a_i\), and leaves the pulsar–outer pair unchanged, so it is not the common modulation of the Section 4.3 template. **Conjectural.** For J0337 we therefore compare scales only; the coefficient envelopes of Table 1 do not translate into it (Section 5.4). **Imported from prior work.** To leading order in the eccentricities and in \(a_{\rm in}/a_{\rm out}\), the orbital modulation of the companions' scalar field at the white dwarf obeys

\begin{equation}
\delta\varphi_{\rm mod}\le3.08\times10^{-10}|a_p|+2.04\times10^{-10}|a_o|,
\label{eq:wd-modulation}
\end{equation}

where the first term comes from the inner eccentricity and the second from the outer eccentricity and the inner-orbit motion. The charge responds in two ways.

- **Proven.** *Instantaneous.* In the weak-field limit \(a_i=\alpha_s(\varphi)\) at the white dwarf, so \(\delta a_i=\beta_s\,\delta\varphi_{\rm mod}\) follows the field without lag; by Equation (\ref{eq:wd-modulation}) it is at most \(1.24\times10^{-9}|a_p|+8.2\times10^{-10}|a_o|\). **Conjectural.** In the language of Section 3 this is a zero-lag coefficient, not an internal state. The fixed-companion reduction of Section 4.3 omits this in-phase modulation of both pair factors. For the pulsar–inner pair factor it is at most \(1.24\times10^{-9}a_p^2+8.2\times10^{-10}|a_pa_o|\), which reaches only the smallest stored Section 5 scales, for \(|a_p|\gtrsim0.5\); its timing response was not computed. The published static SEP limits for J0337, at most about \(2.6\times10^{-6}\) \cite{archibald2018universality,voisin2025planet}, bound the leading-order difference \(\Delta_{po}-\Delta_{io}=a_o(a_p-a_i)\) of the Section 4.3 pair factors. With \(|a_o|\simeq|\alpha_0|\) near the Cassini limit they give \(|a_p|\lesssim4\times10^{-3}\) and a modulation below \(4\times10^{-14}\); \(|a_p|\gtrsim0.5\) would require \(|a_o|\lesssim5\times10^{-6}\). The outer white dwarf has the analogous zero-lag response in the pulsar–outer pair factor, of comparable size by a leading-order estimate; we did not evaluate it further. The same susceptibility also shifts the pulsar-charge stiffness, which Section 4.3 requires to be small relative to \(\kappa\); we did not evaluate that condition either.
- **Conjectural.** *Structural.* The adiabatic structural part \(\delta a_i^{\rm struct}=\mathcal S_{\rm struct}\,\delta\varphi_{\rm mod}\) is in phase, up to the correction \((\omega/\omega_0)^2\). Under the two assumptions of the lag bound, the lagged part of the thermal relaxation is at most \(|\mathcal S_{\rm struct}|\,\delta\varphi_{\rm mod}/2\).

**Imported from prior work.** The Cassini measurement \(\gamma-1=(2.1\pm2.3)\times10^{-5}\) \cite{bertotti2003cassini}, with \(\gamma-1=-2\alpha_0^2/(1+\alpha_0^2)\) \cite{damour1992tensor} taken at its 2\(\sigma\) lower limit, gives \(|\alpha_0|\le3.54\times10^{-3}\), or \(\varphi_\infty\le8.84\times10^{-4}\) for \(\beta_s=-4\). The declared \(\varphi_\infty=10^{-3}\) lies outside this range, so we rescale \(\mathcal S_{\rm struct}\) as \(\varphi_\infty^2\). **Conjectural.** With \(|a_o|\simeq|\alpha_0|\) and \(|a_p|\le1\), the structural scale is \(|\mathcal S_{\rm struct}|\,\delta\varphi_{\rm mod}\le2.13\times10^{-18}\). The in-phase adiabatic part and, under the two assumptions, the thermally relaxing part therefore each change the pulsar–inner pair factor by at most \(2.13\times10^{-18}\) (\(4.3\times10^{-18}\) together) and the inner–outer pair factor by at most \(7.53\times10^{-21}\); the lagged thermal part is at most half of its bound. These lie about eight orders of magnitude below the stored Section 5 scales at \(\tau_\chi=2\) d, the shortest stored lag, where the Table 1 envelopes are smallest:

- \(2.8\times10^{-10}\) (K=1, truncated)
- \(1.7\times10^{-9}\) (K=10, truncated)
- \(3.5\times10^{-9}\) (K=10, full)
- \(1.57\times10^{-7}\) (\(K\approx934\))

**Conjectural.** We did not compute the timing response of either part in these two non-common channels, including the nuisance projection. Reading the structural comparison as a sensitivity statement would require their timing sensitivity to lie within about eight orders of magnitude of the common template's, which was not established. The outer white dwarf's own internal states were not computed.

**Conjectural.** In summary, a scalar pulse leaves in the displaced matter of the inner white dwarf a transient retarded monopole signal, with no first-order permanent charge as \(t\to\infty\) if all its modes are damped. At orbital timescales the charge follows the field instantaneously through \(\beta_s\), and the displacement state adds a computed adiabatic structural term that is small. Its lagged part is bounded under the two stated assumptions on the thermal relaxation. The layer-by-layer calculation gives partial contributions of \(4.0\times10^{-9}\) and \(3.3\times10^{-7}\) of that term at the inner and outer orbital frequencies, from layers scanned up to \(\tau_{\rm cut}=10^{13}\) s. A bound on the total also requires control of the uncomputed tail \(T_{\rm tail}\), as specified above; a non-adiabatic calculation of the deep thermal response was not made. The calculation therefore neither establishes nor excludes an orbital-timescale internal state with a coupling at the Section 5 scales, as the Section 3 benchmark posits. The remaining assumptions are:

- linear response and the weak-field charge \(a_i=\alpha_s(\varphi)\);
- \(|a_p|\le1\) in the scale comparison;
- a Newtonian radial model with a flat-exterior readout;
- no neutral modes and positive damping of all modes;
- neglect of interior Born scattering, where the core compactness is about \(10^{-4}\);
- for the lag, either a thermal relaxation strength bounded by \(\mathcal S_{\rm struct}\) with a single relaxation time, or the layer-by-layer model with an explicit bound on its uncomputed tail;
- the unfixed \(\varphi_\infty\) exponent of the \(10^{-8}\) direct part;
- a non-gray correction with simplified LTE continuum opacity and a fixed pulse kernel;
- a slowly rotating white dwarf;
- a mixed hydrogen–helium envelope, whereas the observed star has a hydrogen (DA) atmosphere \cite{kaplan2014j0337}.

## 5. Conditional application to PSR J0337+1715

### 5.1 Data and template provenance

**Imported from prior work.** The stored analysis uses 12,474 public Nançay pulse times spanning approximately 2987.9 days in 2013--2021, a published planet-model baseline, and the released Nutimo implementation \cite{voisin2025planet,voisin2025release}. Its baseline weighted residual RMS is 1.887 microseconds. Sections 5.1--5.4 retrieve the fixed artifact record in Appendix D. Sections 5.5--5.9 add registered linear-array audits and conditional inference; Section 5.10 adds isolated live transient/derivative evaluations and bounded nonlinear displacement checks. No complete nonlinear astrophysical likelihood or arbitrary pulse reconnection is performed.

**Imported from prior work.** Six stored response columns are central finite differences of the modified integrator: cosine and plus-\(\pi/2\) drives at

\[
\omega_{\rm in}=3.856137\ {\rm rad/day},\quad
\omega_{\rm out}=0.019200\ {\rm rad/day},\quad
\omega_{\rm dif}=\omega_{\rm in}-\omega_{\rm out}.
\]

They are numerical derivatives of a nonlinear model, not solutions free of numerical and linearization error. The two high-frequency carriers are close; three distinct samples are not three equally informative frequency measurements.

**Proven.** The plus-\(\pi/2\) column \(C_{k,s}\) responds to \(-\sin\omega_k t\). Put \(g_{1k}=1/(1+\omega_k^2\tau_\chi^2)\), \(g_{2k}=\omega_k\tau_\chi/(1+\omega_k^2\tau_\chi^2)\). For drive phase \(\phi_k\),

\begin{equation}
\begin{split}
T_\beta=\sum_k d_k\big[&(g_{1k}\cos\phi_k+g_{2k}\sin\phi_k)C_{k,c}\\
&+(g_{1k}\sin\phi_k-g_{2k}\cos\phi_k)C_{k,s}\big].
\end{split}
\label{eq:stencil}
\end{equation}

Zero phase gives \(g_1C_c-g_2C_s\). Reversing g2 defines the advanced comparison template. The instantaneous column is \(T_{c_Y}=\sum_kd_k(\cos\phi_k C_{k,c}+\sin\phi_k C_{k,s})\).

### 5.2 Estimator and interval definition

**Imported from prior work.** The nuisance block includes 28 timing-parameter derivative columns, an offset, 30 low-frequency sine/cosine pairs, and a static-SEP guard. Column normalization and a relative singular-value cut of \(10^{-3}\), with the guard residual retained, yield rank 71. The full stored construction retains all 90 singular vectors. The 19 discarded directions have not been established to be physically forbidden. The co-fitted instantaneous signal column is additional to this nuisance block.

**Proven.** For a chosen nuisance projector, the two-column fit X=\([T_{c_Y},T_\beta]\) has

\begin{equation}
\widehat{\boldsymbol b}=(X_w^TP_\perp X_w)^{-1}X_w^TP_\perp y_w,
\qquad
\sigma_F^2=s^2[(X_w^TP_\perp X_w)^{-1}]_{22}.
\label{eq:fit}
\end{equation}

The subscript w denotes the stored diagonal error weighting; s is the residual scale for that nuisance space. For a fixed Gaussian covariance and a flat amplitude prior, define U by

\begin{equation}
\Phi\!\left(\frac{U-\widehat\beta}{K\sigma_F}\right)
-\Phi\!\left(\frac{-U-\widehat\beta}{K\sigma_F}\right)=0.95,
\qquad U\ge0.
\label{eq:interval}
\end{equation}

This encloses 95 percent of the specified Gaussian in [-U,U]; it is not generally K times the K=1 value. Uncertain covariance requires more than an arbitrary width multiplier.

**Imported from prior work.** K=1 uses the Fisher width with its residual scale. K=10 was selected as a safety factor, not inferred as a noise parameter. K approximately 934 was motivated by a discrepancy in the static-SEP direction, which does not calibrate the dynamic direction. Its stored two-day interval is \(1.57\times10^{-7}\), retained as a sensitivity scenario rather than a guaranteed conservative bound.

**Imported from prior work.** The scan contains 65 lag values from 2 to 500 days and 4821 origins \(t_{\rm off}\in[0,P_{\rm out})\), spaced by \(P_{\rm in}/24\), with \(\phi_k=\omega_k t_{\rm off}\) and \(d_k=1\). For each lag and K the maximum U on this finite grid is recorded. We call it a *registered-grid envelope*, not Bayesian phase marginalization or a certified continuous-phase supremum.

**Proven.** A noninteger inner/outer frequency ratio means that advancing the origin by one outer period does not return every carrier phase to its initial value. This interval is an analysis domain, not an exact common period. Independent carrier phases enlarge the model further. Neither the finite-grid envelope nor its calibrated scan statistic automatically extends to these larger domains.

### 5.3 Stored results and their scope

**Imported from prior work.** The registered-grid causal maximum is \(\max|z|=2.2851\), with simulation-based global p=0.26. The advanced comparison gives 2.3398 and p=0.282. The stored detection flag is false under the registered rule. These probabilities refer to the recorded Gaussian-noise simulations, nuisance space and finite grid; they are not an independent validation of the noise model.

**Imported from prior work.** Table 1 reports envelopes for the normalized coefficient \(\beta\). The full-space scenario precedes the truncated result to expose the latter's stronger apparent sensitivity. Neither K=10 column is a calibrated astrophysical 95 percent exclusion. The K=1 column is conditional on the truncated Gaussian model, not a model-independent floor.

\begin{table}[htbp]
\centering
\caption{Imported from prior work. Registered-grid envelopes U for the unit-drive coefficient \(\beta\). Full and truncated use the stored 90- and 71-direction nuisance spaces. K is the width multiplier in Equation (\ref{eq:interval}).}
\label{tab:intervals}
\begin{tabular}{r r r r r}
\toprule
\(\tau_\chi\) (day) & Full, K=10 & Trunc., K=10 & Trunc., K=1 & Full/trunc.\\
\midrule
2 & \(3.534\times10^{-9}\) & \(1.680\times10^{-9}\) & \(2.794\times10^{-10}\) & 2.10\\
5 & \(9.673\times10^{-9}\) & \(1.852\times10^{-9}\) & \(2.884\times10^{-10}\) & 5.22\\
18 & \(3.166\times10^{-8}\) & \(1.976\times10^{-9}\) & \(2.933\times10^{-10}\) & 16.02\\
52 & \(4.488\times10^{-8}\) & \(2.600\times10^{-9}\) & \(3.669\times10^{-10}\) & 17.26\\
200 & \(1.241\times10^{-7}\) & \(7.155\times10^{-9}\) & \(9.955\times10^{-10}\) & 17.35\\
\bottomrule
\end{tabular}
\end{table}

**Imported from prior work.** Increasing the truncated model's Fourier block from 30 to 60 pairs changes the stored K=10 intervals by about 0.52--0.60 percent. This does not address the factor 2.1--17.4 full/truncated dependence over the 65 stored lags, or replace covariance inference with correlated achromatic and chromatic noise. Figure 1 displays the nuisance-space comparison explicitly.

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figures/conditional-intervals.pdf}
\caption{Imported from prior work. Stored interval envelopes and nuisance-space dependence. Left: K=1 and K=10 truncated constructions and the K=10 full-space scenario. Right: the full/truncated K=10 ratio. Connecting lines guide the eye; they do not certify unsampled lag or phase values. These are normalized-coefficient intervals, not universal SEP exclusions.}
\label{fig:intervals}
\end{figure}

**Imported from prior work.** Corrected live-integrator gates reported null offsets no larger than 0.039 Fisher sigma, detection-amplitude recovery deviations of about 0.10 sigma or less, and relative recovery deviations of about 0.13 percent or less at tested limit amplitudes. The two-day, 104.349-day-origin gate targeted the truncated-grid maximum. These local tests do not validate coverage of either K scenario or full-space maxima at every lag and origin.

**Imported from prior work.** A finite turn-number lattice tested smooth linear-plus-quadratic reassignments. Seven viable alternatives per tested amplitude survived its fit criterion with small estimator displacements. Arbitrary isolated slips in observing gaps were not tested. The causal/advanced overlap study sampled two origins per lag and found both strong and weaker overlaps. Similar maxima of the detection statistics do not establish collinearity throughout the origin domain. No lag-sign detection is claimed.

### 5.4 What amplitude is constrained?

**Proven.** The relaxation-only amplitude of carrier k is

\begin{equation}
A_{\chi,k}=\frac{|\beta d_k|}{\sqrt{1+\omega_k^2\tau_\chi^2}}.
\label{eq:carrier-amplitude}
\end{equation}

Thus \(\beta\), a carrier amplitude, a peak and an RMS are distinct. For the specified drive, \(\sup_t|q_\chi(t)|\le\sum_k A_{\chi,k}\). The co-fitted instantaneous contribution requires its coefficient and covariance too. A beta interval alone does not bound the total \(\Delta(t)=\Delta_0+c_YF+c_\chi\chi\).

**Proven.** At 200 days and the stored frequencies, the unit-drive amplitude multipliers are approximately 0.001297, 0.2520 and 0.001303 for inner, outer and difference carriers. The table therefore cannot be relabeled as an undifferentiated amplitude of SEP oscillation.

**Imported from prior work.** An auxiliary stored analysis used a potential-based drive dictionary whose phase map fails the gate in Section 4.4. We withdraw its physical-coupling interpretation; the historical unit-drive table retains only its explicitly prescribed meaning. Published static-SEP bounds also depend on planet/noise assumptions \cite{voisin2025planet}; comparing their magnitudes with Table 1 does not compare the same parameter or statistical construction.

### 5.5 Nuisance and coverage audits on the frozen linear arrays

**Imported from prior work.** The subsequent SVD audit finds numerical rank 90 with minimum relative singular value \(4.18\times10^{-7}\). Independent full-space QR/SVD constructions agree to \(1.3\times10^{-10}\) in subspace residual. The weak modes mix timing, astrometric, planet and Fourier directions; nineteen dimensions are discarded, not nineteen individually identified parameters. No physical prior excluding that subspace has been established. All historical full/truncated K=10 anchors reproduce within the registered relative tolerance \(10^{-5}\). Full rank is therefore the primary baseline below. This numerical audit does not independently certify finite-difference derivative accuracy or nonlinear nuisance curvature.

**Proven.** For a truncated beta estimator l and orthonormal omitted residual directions D, the largest bias in units of its nominal white-noise width, for omitted mean norm at most R, is \(R\|D^Tl\|/\|l\|\). Full nuisance fitting removes this specific mean bias. At fixed nominal width, no finite multiplier protects uniformly against an unbounded overlapping omitted mean. That statement does not cover every data-dependent scale rule; the following experiments test the actual residual-scale procedure separately.

**Proven.** With an unbiased Gaussian estimate and known covariance, the K=1 symmetric interval in Equation (\ref{eq:interval}) has frequentist coverage at least 0.95. For positive true beta exceeding \(1.959964\sigma\), define h by
\(\Phi((\beta-h)/\sigma)+\Phi((\beta+h)/\sigma)=1.95\).
Failure is \(|\widehat\beta|<h\), of probability
\(2\Phi((\beta+h)/\sigma)-1.95\le0.05\).
Smaller absolute beta is always included and symmetry handles negative beta. An envelope containing the true grid cell cannot reduce this coverage. Estimated or misspecified covariance needs separate validation.

**Imported from prior work.** A registered experiment uses 8192 independent realizations per condition, 18 preselected lag/origin cells and nine amplitudes (0, \(\pm2,\pm5,\pm20,\pm50\) full-space unit-noise sigma). The two signal coefficients are co-fitted. Common random numbers correlate comparisons, so rows are not independent repeated experiments. The main stresses are omitted means of norm 3, 10 or 30 and extra Fourier modes 31--60 with spectral variance proportional to \(j^{-4}\), normalized to weighted RMS 0.25 or 1. A direct/compressed projection check passes. Table 2 reports pointwise coverage minima over the tested cells and amplitudes.

\begin{table}[htbp]
\centering
\caption{Imported from prior work. Minimum pointwise coverage of the K=1 U construction across tested cases. Each condition has 8192 realizations. Minima are descriptive, not simultaneous guarantees. Estimated covariance is a separately registered follow-through; dashes denote conditions not rerun for that estimator.}
\label{tab:coverage}
\begin{tabular}{l r r r r}
\toprule
Generating condition & Trunc. diag. & Full diag. & Oracle GLS & Est. GLS\\
\midrule
White noise & 0.9468 & 0.9481 & 0.9481 & 0.9473\\
Omitted mean, norm 3 & 0.0956 & 0.9481 & 0.9481 & --\\
Omitted mean, norm 10 & 0 & 0.9481 & 0.9481 & --\\
Omitted mean, norm 30 & 0 & 0.9481 & 0.9481 & --\\
Extra Fourier RMS 0.25 & 0.9310 & 0.7469 & 0.9485 & 0.9471\\
Extra Fourier RMS 1 & 0.8132 & 0.5825 & 0.9490 & 0.9459\\
\bottomrule
\end{tabular}
\end{table}

**Imported from prior work.** The worst full diagonal K=1 cell gives 4772/8192 coverage, with Wilson 95 percent interval [0.57180,0.59316]. Its actual estimator noise is 9.21 times the nominal unit-noise sigma, while its median fitted residual scale is only 1.19. A truncated pointwise K=10 interval also has a 0/8192 stress case. This does not demonstrate failure of the historical K=10 grid envelope: separately registered complete-grid tests, with 512 realizations at lags 2 and 200 days, give K=10 coverage at least 511/512 in every tested scenario. These finite successes do not certify all nuisance means, noise or phases.

**Imported from prior work.** The oracle GLS control restores approximately nominal coverage for the specified generating covariance. A follow-through registered after those outcomes estimates
\(C=\sigma^2(I+a^2LL^T)\), using the same fixed extra-Fourier family and an a grid consisting of zero and 100 logarithmic values from 0.01 to 4. The restricted maximum-likelihood (REML) objective includes \(\log\det C_0\), \(\log\det(X^TC_0^{-1}X)\) and \(\nu\log({\rm RSS}/\nu)\), where \(C_0=I+a^2LL^T\), the full-nuisance contrasts are understood, and \(\nu=12474-90-2\). No grid is expanded after observing a result.

**Imported from prior work.** Across 486 follow-through conditions, minimum K=1 U coverage is 7749/8192=0.945923, Wilson interval [0.940813,0.950615]. No fit reaches the upper covariance-amplitude endpoint. Zero-amplitude ordinary-fit and signal-translation controls pass. This is Monte Carlo support within the prescribed covariance family, not an exact coverage theorem or an independently validated astrophysical noise model. Figure 2 compares the noise stresses.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.90\linewidth]{figures/coverage-validation.pdf}
\caption{Imported from prior work. Minimum pointwise K=1 U coverage over the preselected lag/origin/amplitude cells under three noise conditions. Dots are descriptive minima from 8192 realizations per condition, not independent global confidence statements. Oracle and estimated GLS use the generating covariance family; estimated GLS uses a separately registered simulation seed. The dashed line is 0.95.}
\label{fig:coverage}
\end{figure}

**Imported from prior work.** Applying that REML specification to stored residuals at the same preselected cells gives covariance amplitude 0.07368 or 0.08316. At lag 2 days and the full-space reference origin 174.27780 days, the local K=1 U is \(5.25094\times10^{-10}\). This is a conditional local fit, not a newly searched envelope, a detection or a limit for the physical drive of Section 4.4. Chromatic/instrumental noise, nonlinear timing errors and other covariance families remain unvalidated here.

### 5.6 Projected comparison with derivative and fast-response models

**Proven.** Let A be the six-carrier response map after full nuisance removal and covariance whitening. For pole coefficients w and comparator columns W, the residual information is

\begin{equation}
I_{\beta,W}=\min_{\boldsymbol c}\|A(w-W\boldsymbol c)\|^2
\ge s_{\min}(A)^2\min_{\boldsymbol c}\|w-W\boldsymbol c\|^2.
\label{eq:map-bound}
\end{equation}

Known phase changes rotate each cosine/sine pair orthogonally. The coefficient-space distance on the right is phase invariant and therefore bounds the continuous three-phase domain. This is different from allowing a comparator to fit unknown phases and arbitrary amplitudes independently: free quadratures span all six response columns and absorb any signal in that span.

**Imported from prior work.** A registered follow-through compares derivative powers 0 through N, for N=0,...,5, and the even-only set \(\{0,2,4\}\). It uses six lags, three specified known-phase stresses and fixed covariance amplitudes \(a=0,0.0831611,1\). The middle value is inherited from earlier local REML fits, not a new global noise estimate. Unit carrier amplitudes are used even for the phase-only stress borrowed from Section 4.4. The whitened maps have rank six and condition numbers 98.19, 111.69 and 124.82. Direct/coordinate whitening and QR/SVD controls agree.

**Imported from prior work.** Relative to co-fitting only the instantaneous column, the minimum retained information over these cases is 0.005389 at N=1, \(6.29\times10^{-5}\) at N=2, \(1.17\times10^{-5}\) at N=3, and \(4.63\times10^{-6}\) at N=4. The last case widens the unit-noise standard error by about 465 times. The even-only comparator retains at least 0.0311 in the tested cases. These minima are descriptive; Equation (\ref{eq:map-bound}) separately supplies a continuous-phase lower bound. Mathematical noninterpolation can leave very little measurable information.

**Proven.** At N=5 the information is exactly zero by Theorem 3. Numerical relative residual information between \(5.4\times10^{-30}\) and \(5.3\times10^{-25}\) is rounding error around that identity, not a weak detection channel.

**Proven.** The fast-spectrum inequality also gives a nuisance-aware distance witness. Let v extract \(S_l-R_\Lambda S_h\) from the dephased coefficients. For the comparator \(v^T\theta\le0\), whereas a positive slow pole has \(v^Tw>0\). Its covariance-weighted distance from the entire comparator cone is at least

\begin{equation}
\frac{\beta\,v^Tw}{\sqrt{v^T(A^TA)^{-1}v}},\qquad\beta>0.
\label{eq:projected-witness}
\end{equation}

This follows by Cauchy--Schwarz in the coefficient covariance metric. A freely fitted real instantaneous term is annihilated. The audit evaluates a positive witness for every tested lag with the explicitly illustrative \(\Lambda=10\omega_{\rm in}\). This is a conditional discrimination calculation, not a data significance, a measured gap or an EOS constraint.

### 5.7 Expanded phases and an analytic envelope

**Proven.** For the instantaneous-plus-pole fit, define \(d_\tau=\min_c\|w_\beta-cw_0\|\). Equation (\ref{eq:map-bound}) and \(U\le|\widehat\beta|+z_{0.975}\sigma_\beta\) imply the continuous all-phase upper envelope

\begin{equation}
\sup_{\boldsymbol\phi}U(\boldsymbol\phi)
\le\frac{\|\operatorname{Proj}_{\operatorname{col}(A)}y\|+z_{0.975}s}
{s_{\min}(A)d_\tau}.
\label{eq:phase-envelope}
\end{equation}

The estimator direction belongs to col(A), and \(d_\tau\) is unchanged by phase rotations. The residual scale s is the pre-signal-fit scale of the specified construction, computed in the fixed covariance metric. This deterministic inequality does not establish coverage for misspecified or estimated covariance. Its numerical evaluation is not an interval-arithmetic rounding certificate or necessarily a tight supremum.

**Imported from prior work.** The audit compares the historical 4821-origin grid, a \(72\times72\) grid of two orbital longitudes with zero closure, and a \(24^3\) grid of independent carrier phases. They are not nested. In particular, the coarse three-phase grid misses narrow long-lag extrema already sampled in the historical domain. A separately registered refinement starts from each grid's maximum and uses 17 successively halved angular steps. It preserves all seeds, reaches no move cap, and has final-step relative improvements below \(2.6\times10^{-9}\).

**Imported from prior work.** The local-refined U is 1.00073--1.03721 times the historical-grid value across the tested lag/covariance cells. This is not certification of the global optimum; Equation (\ref{eq:phase-envelope}) provides the all-phase upper bound. At fixed \(a=0.0831611\), local-refined U at lags 2, 200 and 500 days is \(6.3066\times10^{-10}\), \(3.3933\times10^{-8}\) and \(8.2658\times10^{-8}\), while the analytic upper bounds are \(1.4226\times10^{-8}\), \(4.6943\times10^{-8}\) and \(1.1317\times10^{-7}\). Figure 3 shows the differing tightness. These are unit-drive constructions; expanded-domain detection significance and estimated-covariance coverage are not newly calibrated here.

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figures/comparator-phase-validation.pdf}
\caption{Imported from prior work. Left: ranges of retained information over the specified lag, phase and covariance cases, relative to the instantaneous-only comparator. E4 allows only powers 0,2,4; unrestricted P5 has exact zero information. Right: local-refined three-phase U and the analytic all-phase upper bound at fixed covariance amplitude a=0.0831611. The refinement is not a global-optimum certificate, and neither curve is a physical-drive SEP limit.}
\label{fig:comparator-phase}
\end{figure}

### 5.8 Initial state, observing gaps and numerical scope

**Proven.** A homogeneous amplitude decays to fraction \(\epsilon\) only after \(\tau_\chi\log(1/\epsilon)\); without an initial-amplitude bound this gives no uniform absolute signal bound. At \(\tau_\chi=500\) days, one-per-mille settling requires 3453.88 days, longer than the stored observing span. More fundamentally, the causal differential operator

\begin{equation}
\mathcal A(D)=D\prod_k(D^2+\omega_k^2)
\label{eq:transient-obstruction}
\end{equation}

annihilates the constant and all six periodic inputs, while its action on \(e^{-t/\tau_\chi}\) is nonzero. Linear operators L and \(L+\eta\mathcal A(D)\) can agree on every stored carrier/static response and differ arbitrarily on a transient. This proves insufficiency of the finite response record even within causal linear maps; it does not assert that every such operator is a realization of the timing theory. A validated forward model or a dedicated transient response is needed. An exponential in the coupling must not be substituted directly as an exponential TOA residual.

**Imported from prior work.** On the stored sampling, 0.7984--0.9973 of an exponential input's norm lies outside the constant-plus-six-periodic input span across the tested lags. This is input independence, not a transient timing-response computation. The original periodic-state assumption is therefore retained explicitly.

**Imported from prior work.** Every one of the 565 adjacent-TOA gaps longer than one day is tested with an idealized permanent single-cycle residual step, using the frozen spin frequency. Both step signs are evaluated after fitting all 90 nuisance directions and all six free harmonic coefficients. Under diagonal covariance the weakest step leaves norm 6.4337 and minimum change in residual chi-square 27.9406. Under the a=1 correlated-noise stress these fall to 1.1149 and 0.1529, at the 223.3708-day gap between relative days 2066.3729 and 2289.7438. Thus this linear stress test cannot uniformly rule out a single-cycle error. It is not evidence of an actual slip: the stress covariance is not the measured best noise model. Arbitrary multiple slips and nonlinear pulse-number reconnection remain untested.

**Imported from prior work.** QR/SVD agreement for the fixed full nuisance matrix remains \(1.3\times10^{-10}\). The superseded coarse-step v1 and corrected v2 matrices have maximum principal sine 0.999776, illustrating why v1 is not an equally valid alternative. That difference does not estimate the remaining v2 error. Stored half-step summaries range from 0.0016004 to 0.0248096 and do not supply all projected error vectors. For full-rank B, a bound \(\|E\|\le\epsilon<s_{\min}(B)\) would control projector rotation by \(\epsilon/[s_{\min}(B)-\epsilon]\), capped at one; numerical agreement at a fixed B supplies no such derivative-error budget. The new audits therefore retain their frozen-linear-model scope.

### 5.9 Corrected leading drive and simultaneous coefficient inference

**Imported from prior work.** Applying the runtime's parameter-set convention gives Newtonian masses \((1.43781441,0.19753639,0.41010271)M_\odot\) and relative semimajor axes \((4.77619153\times10^9,1.76487508\times10^{11})\) meters. In this convention the inner inclination is obtained from the outer inclination plus \(\delta_i\) in degrees. With \(F=(\delta U/c^2)/U_\star\) and \(U_\star=\sum_k A_k=1.74999156\times10^{-10}\), the normalized positive amplitudes are \((0.24394925,0.69195682,0.06409393)\); phases follow Equation (\ref{eq:closure}). Exact coplanar Kepler-potential grids at \(128^2\) and \(256^2\) points both give omitted-input RMS 3.37427 percent of the leading-drive RMS. This is not an omitted timing-error bound or complete scalar-force matching.

**Proven.** Let the nuisance-projected data have mean \(C\theta\), with six real coefficients, and write \(C=Q_6R_6\). Fitting all six coefficients before selecting a lag/phase makes the fitted covariance and residual scale invariant under translations by any signal in \(\operatorname{col}(C)\). A confidence region \(E\) for \(\theta\) can be inverted through \(\theta=W(\phi,\tau)b\): whenever it contains the true coefficient vector, its full inverse image contains the true physical parameter point. This controls continuous phase/lag selection within the declared mean model; a few numerical sections do not equal the full projection.

**Proven.** A conservative alternative is available if the white scale is known and \(\Sigma(a)=I+a^2LL^T\) has a declared \(a\le a_{\max}\). Then \(\Sigma(a)\preceq\Sigma(a_{\max})\). GLS with the latter covariance has true estimation-error covariance no larger than its nominal value. Its six-dimensional quadratic error is stochastically bounded by \(\chi^2_6\). A 95-percent region therefore covers every amplitude in the declared interval, without a phase-grid correction. Unknown white scale, unbounded amplitude or covariance misspecification is outside this theorem.

**Imported from prior work.** A six-coefficient REML experiment calibrates a quadratic threshold from 8192 draws at each \(a=0,0.25,1,4\), then freezes it before independent validation. The threshold is 12.8241766; independent inclusion rates are 0.951660, 0.954346, 0.952759 and 0.954102. A separate shallower Fourier-spectrum stress yields 0.959229. These finite results support the specified procedure, not every astrophysical noise process or every intermediate covariance. Figure 4 shows the validation.

**Imported from prior work.** The recorded data's six-coefficient omnibus statistic is 16.3525, above that threshold (nominal chi-square p of about 0.012 for six degrees of freedom). Its null is that all six carrier coefficients vanish, not \(\beta=0\) with a physical instantaneous term fitted. The minimum statistic over the physical plane is 12.84 at 2 days and 14.39--14.59 at the other five evaluated lags, all above the threshold, so every evaluated physical lag section is empty: at these six lags, no instantaneous term plus single relaxation driven by the leading physical drive reproduces the carrier coefficients at the calibrated level. The 2-day rejection is marginal; with the sign \(\beta\ge0\) of the equal-charge realization every section minimum is at least 14.59. The region, its threshold and the omnibus statistic do not depend on the drive phases, so only the sections were recomputed. No beta interval or EOS-matched bound follows, and the excess is not attributed. With the earlier periastron error, every section contained \(\beta=0\); that result is withdrawn.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.95\linewidth]{figures/remaining-levers-validation.pdf}
\caption{Imported from prior work. Independent six-coefficient region inclusion, with Wilson 95-percent intervals from 8192 draws per condition (left); beta standard-error ratio after adding a measured transient at the fixed middle covariance (right). The latter is a local fit diagnostic, not validation of a seven-coefficient confidence region.}
\label{fig:remaining-validation}
\end{figure}

### 5.10 Live transient and derivative checks, and pulse-count limits

**Imported from prior work.** An isolated copy of the dynamic-SEP engine adds an exponential drive referenced to the integrator epoch. Both archived and rebuilt engines reproduce the stored zero-drive residual array exactly; source and binary hashes are retained. At \(\tau=2,52,500\) days, central differences at amplitudes \(10^{-8}\) and \(5\times10^{-9}\) give full-nuisance-projected derivative changes 0.389791, 0.005752 and 0.000732 percent. All pass the registered 5-percent convergence gate, with maximum perturbations below 58 microseconds and successful final zero-drive recovery. These dedicated columns supply the response missing from the earlier finite periodic record; they do not invalidate its insufficiency theorem.

**Imported from prior work.** Co-fitting the measured transient with the corrected leading instantaneous and relaxed drive changes beta standard errors by factors 0.999998--1.00742 across the nine tested lag/covariance combinations. This resolves a specific local transient question. It does not certify every lag, large initial amplitude, or coverage after adding that coefficient to the six-dimensional mean model.

**Imported from prior work.** All 28 timing derivatives were recomputed at half the archived corrected steps, and seven planet derivatives at quarter steps. The largest half-step weighted change is 0.665177 percent. Some quarter-step changes increase, so uniform convergence is not established. Both nuisance matrices have rank 90, yet their maximum principal sine is 0.999913. The normalized matrix change norm 0.006857 exceeds the archived minimum singular value \(1.12207\times10^{-6}\). Replacing the derivatives changes corrected-drive diagonal-noise standard errors by factors 0.7858--0.8722 and shifts estimates. The half-step matrix is a sensitivity diagnostic, not a certified replacement. A rigorous derivative-error budget remains unavailable.

**Imported from prior work.** All 201 assignments with at most two nonzero \(\pm1\) steps on the ten longest gaps were evaluated. The weakest nonzero assignment remains the 223.37-day single step: minimum changes in chi-square are 27.9406, 11.6183 and 0.152887 at \(a=0,0.0831611,1\). Its unconstrained linear nuisance compensation requires \(\delta\eta_{\rm extra1}=26.0400\). Fractions 0.25, 0.5 and 1 yield eccentricities exceeding one, where the runtime clamps the requested orbit; that displacement is therefore rejected as a physical nonlinear fit. This does not exclude a different nonlinear solution for the pulse assignment.

**Imported from prior work.** A separately registered local follow-through uses fractions 0.001, 0.003 and 0.01 of the same direction. Weighted nonlinear discrepancies are 0.4797, 1.5042 and 4.0799 percent, below the 5-percent gate. These probes reach only one percent of the full compensation and cannot rescue it. Arbitrary pulse reconnection and a constrained nonlinear likelihood search remain incomplete.

## 6. Discussion

**Proven.** The two parts answer one question at different levels. Positive weights and a finite cutoff make a specified local description finite; sufficient smoothness permits its weighted expansion. Neither fact removes its observable effects. A settled single pole escapes a restricted derivative comparator, but finitely many frequencies can be interpolated once enough coefficients are admitted. Nuisance freedom can erase a remaining functional distinction.

**Counterexample candidate.** The physical target is a shared transfer relation with known drive and constrained projection, tested against an explicitly bounded comparator. The one-state model is one realization. Hereditary kernels, multiple states, nonlinear thresholds and other sectors remain separate possibilities.

**Conjectural.** Section 4.6 works through one internal state of the inner white dwarf, the scalar-driven displacement of its own matter. That state leaves a transient retarded monopole signal. Under the stated assumptions, its structural modulation of the pair factors lies about eight orders of magnitude below the stored Section 5 scales. The white dwarfs' instantaneous zero-lag responses could reach the smallest of those scales only for \(|a_p|\gtrsim0.5\), which the static SEP limits allow only if \(|a_o|\lesssim5\times10^{-6}\); they are not lags. Neither timing response was computed, so this is a scale comparison, not an exclusion. Other internal states of the white dwarfs, and the neutron-star charge of Section 4.3, remain open.

**Imported from prior work.** Stored J0337 columns and local validation gates support a conditional timing-template application. Its sensitivity changes substantially with the retained nuisance space. That dependence is part of the result, not a small correction to the narrowest interval.

**Proven.** The completed audits close several conditional questions: the full stored nuisance span removes the identified mean-bias mechanism; known-Gaussian U coverage has an analytic guarantee; a scalar-charge model supplies reciprocal force matching and a bounded inertial reduction. Its leading physical phase closure also rejects the archived auxiliary interpretation. Reciprocal fast-response comparison, continuous phase bounds and the causal transient obstruction further specify what the periodic evidence can identify. These conclusions preserve unsuccessful gates rather than rescuing a numerical headline.

**Conjectural.** A physically calibrated exclusion still requires numerical EOS-to-body matching, control of omitted drive/force terms, a reliable derivative space and an adequate astrophysical likelihood. Corrected leading-drive and joint-region calculations address defined conditional questions and reject the leading physical-drive plane at every evaluated lag; the live transient responses cover three lags. The measured derivative sensitivity and unphysical full pulse compensation prevent empirical promotion. The fast-rate gap remains an unestablished physical premise. The finished result is a conditional analytic and methodological paper, not a universal empirical SEP bound.

## Use of AI tools

This work used AI agents substantively, under the author's direction and responsibility. Coding agents built on Claude Opus 5.5 (Anthropic, run through Claude Code) and GPT-6-Astra (OpenAI, run through the Codex command-line tool) wrote and ran most of the analysis, stellar-structure and verification code, carried out the computations reported here, and drafted and revised manuscript text, including scientific claims and explanations. GPT-6-Astra, Claude Opus 5.5 and Claude Fable 5.1 (Anthropic, run as a Claude Code subagent) independently reviewed drafts of Section 4.6 in several rounds and the complete manuscript before submission; their reports and the responses are archived in the repository. The author set the research questions, the acceptance criteria and the claim labels, and checked the AI output against the stored numerical records, the SHA-256-bound manifests, symbolic checks, the verification scripts listed below and these reviews.

## Statements and declarations

### Competing interests

The author declares no competing interests.

### Funding

The author received no external funding for this research.

### Author contributions

J.K. conceived and directed the study, set the research questions and acceptance criteria, checked the results against archived numerical records and verification tests, and reviewed and approved the manuscript. AI assistance with coding, computation, drafting and internal review is disclosed in the Use of AI tools section. J.K. takes responsibility for the work.

## Data and code availability

The public timing release is cited as a dataset \cite{voisin2025release}. The manuscript source, code, notes, manifests and numerical results are available in the public snapshot tagged \path{grg-submission-2026-09-28} at \url{https://github.com/lpaiu-cs/self-mass-unobservability}. The accompanying \path{paper/revision-manifest.json} identifies this revision's inputs by SHA-256, and \path{output/submission/grg/manifest.json} binds the submission files. Binary runtime arrays, about 50 GB in total, exceed the hosting limits and are not deposited; the manifests identify them by SHA-256, and they are available from the author on request. Commit identifiers cited in this paper, such as the unified input state `4897038`, refer to the author's full repository history, of which the public repository holds a snapshot; \path{PUBLIC_SNAPSHOT.md} there states what the snapshot omits. The ordered follow-through designs and results are recorded in the revision manifest.

The compact public record \path{outputs/research-completion/public-inference.json} contains the six fitted coefficients, their precision matrix, drive parameters and frozen threshold. Running \path{verification/replay_public_inference.py} with NumPy reproduces the omnibus and six physical lag sections without private arrays. Its exported quantities are conditional on the archived timing derivatives and covariance fit; this replay does not certify those derivatives, rerun the engine, or recalibrate coverage. The broader \path{verification/verify_unified_paper.py} requires the omitted arrays, including the baseline and 35 finite-difference Jacobian files.

The white-dwarf calculation of Section 4.6 was added after the unified input state. It is recorded in Korean-language notes \path{notes/REQUEST244_*} through \path{notes/REQUEST288_*} and in manifests under \path{outputs/direct-eos-gr33/}; all are bound in the revision manifest. Each manifest lists the scripts and results of its step by hash. Verdict, classification and derived-bound fields written into manifests before the independent reviews are superseded by the errata in \path{notes/REQUEST284_*} through \path{notes/REQUEST286_*}. The thermal-lag verdict of step 287 applies to the scanned shells only; \path{notes/REQUEST295_SUBMISSION_CORRECTIONS_KO.md} records the cutoff and tail qualification added here. The main results map to these manifests:

- endpoint readout and Born first-return ratio: \path{native-quad-refined-primary-manifest.json}
- equation-of-state test: \path{native-eos-sensitivity-manifest.json}
- mechanism and free-fall reproduction: \path{native-charge-mechanism-manifest.json}
- density kernel and radial boundary: \path{native-structure-eft-boundary-manifest.json}
- tides, structural susceptibility and field modulation: \path{native-tidal-photosphere-manifest.json}
- Cassini rescaling: \path{native-final-closure-manifest.json}
- corrections, mass normalization and history: \path{native-closure-transit-manifest.json}
- photosphere: \path{native-atmosphere-reconstruction-manifest.json}
- layer-by-layer thermal relaxation: \path{native-thermal-relaxation-manifest.json}
- pair-channel and instantaneous-response arithmetic: \path{notes/REQUEST284_*}

\path{notes/REQUEST285_*} maps each result to its inputs, scripts and outputs. Steps that read runtime arrays of the endpoint run, which are listed by hash in \path{native-quad-refined-primary-manifest.json} but not deposited, cannot be re-run from the deposit alone. \path{verification/verify_unified_paper.py} does not check Section 4.6.

Reproduction commands and their scope are in \path{paper/README.md}. \path{verification/verify_unified_paper.py} checks the smooth-flat correction, finite-carrier interpolation, tidal/SEP scaling, and table arithmetic without invoking Nutimo. Existing character and ODE checks remain separate. New nuisance, coverage and force-matching programs have fixed-input results under \path{outputs/research-completion/}. Historical runtime records remain provenance; the separately hashed Request 12 returns contain the new live evaluations described in Section 5.10. Outputs superseded by the periastron-convention correction are kept in \path{outputs/research-completion/withdrawn-periastron-convention/}. Rerunning the registered six-coefficient validation with the same seeds under a different linear-algebra configuration, such as another thread count, reproduces its inclusion rates only within Monte Carlo error, because the simulated noise is drawn in an eigenbasis that is not unique for repeated eigenvalues; the registered simulation rows are therefore retained.

\appendix

## Electric-sector reduction and optional family census

**Proven.** After acceleration reduction, pure-E invariants through weight four are \(I_2,I_3,\operatorname{tr}E^4,I_2^2\); trace-free E has no linear scalar. Its characteristic identity is

\[
E^3-\tfrac12\operatorname{tr}(E^2)E-\tfrac13\operatorname{tr}(E^3)1=0.
\]

Multiplying by E and tracing eliminates \(\operatorname{tr}E^4\). The possible single-time-derivative contractions are \(E:D_tE\) and \(\operatorname{tr}(E^2D_tE)\), both total derivatives. With two time derivatives, \(E:D_t^2E\equiv-I_t\).

**Proven.** A single rank-three gradient paired with weight-allowed rank-two blocks has an odd Cartesian index count and no delta-only scalar without acceleration. Two gradients at weight four give the unique norm of the STF rank-three irrep. Internal traces vanish. No spatial integration by parts is used. This exhausts the specified acceleration-free block set.

**Proven.** Independence in the action quotient can be tested by integration over periodic worldline data, for which total derivatives vanish. First take constant E and zero spatial gradient. Scaling E separates degrees two, three and four; a trace-free E with nonzero cubic trace then forces the coefficients of \(I_2,I_3,I_2^2\) to vanish. Next use a time-independent harmonic cubic potential at its stationary origin: E vanishes there but its gradient has a nonzero norm, isolating \(I_g\). Finally take a harmonic quadratic potential with \(E(t)=A\cos\nu t\) at its stationary origin. With the previous coefficients zero, the positive period integral of \(I_t\) isolates the last coefficient. These trajectories satisfy the leading free-fall condition. Thus the five representatives are independent modulo the declared relations; the exact character census corroborates completeness.

**Imported from prior work.** An optional algebraic census admits independent additional primitives. The exact-character program reproduces dimensions 5, 16, 30, 15, 17, 23, 17 and 21 for E, E/B, E/B/S, E/V and E/X at STF ranks three through six. B has its declared magnetic parity. Additional primitives' gradients are generic unless separately constrained. The census is not a claim that B is present in the purely electric reference background, or that every parity-even relativistic scalar lies in the chosen delta-only catalog. It is not required for Theorem 3.

## Boundary statements at their proper level

**Proven.** Repeated differentiation of the flat function on Y>0 gives a polynomial in 1/Y times \(e^{-1/Y^2}\). All such expressions tend to zero, establishing smoothness with a zero jet. Substituting x=\(1/Y^2\) turns the remainder ratio into \(x^{n/2}e^{-x}\), also tending to zero for every finite n. Exact analytic-germ reconstruction fails; Proposition 2 does not.

**Proven.** For \(\sqrt Y\Theta(Y)\), a polynomial matching the value at zero is either zero or starts at O(Y). Its approximation error retains leading order \(Y^{1/2}\) as Y approaches zero from above and cannot be O(\(Y^5\)). This failure is at the differentiability/expansion layer; a physical branch model must specify threshold and domain.

**Proven.** A causal exponential memory kernel is finitely realizable by a relaxation state. A finite linear time-invariant state system has rational transfer \(C(z1-A)^{-1}B+D\). A noninteger power-law transfer with a branch point cannot equal it on an open domain. This concerns exact linear time-invariant realization, not finite approximation on a measured band or arbitrary nonlinear encodings.

**Proven.** Finite explicit states preserve a finite local description but invalidate an instantaneous external-variable-only readout when state/history matters. Infinitely many independent low-weight primitives can instead make the candidate catalog infinite before reduction. These are different failure layers, and do not select one unique two-parameter escape.

## Interpolation example and checks

**Proven.** At \(c_Y=0\), \(\beta=\tau_\chi=1\), \(\omega=1,2,3\), Equation (\ref{eq:interpolant}) gives

\[
P_5(z)=-\frac{z^5}{100}+\frac{z^4}{100}-\frac{3z^3}{20}
+\frac{3z^2}{20}-\frac{16z}{25}+\frac{16}{25}.
\]

It equals \(1/(1+z)\) at \(z=\pm i,\pm2i,\pm3i\). These real coefficients define a local derivative comparator, not an internal state. This is exact finite-sample collapse away from the adiabatic limit; Theorem 3 supplies the general construction.

**Proven.** The verification checks all six equalities, the polynomial's degree and reality, the ODE solution and derivative-expansion residual. It separately tests the flat remainder and threshold obstruction. The output is mathematical verification, not an observational detection.

## Stored-analysis provenance and revision boundary

Table 1 and Figure 1 read \path{request10_external/sep_dynamic/sep_phase_marg_10_8e.json}. The Fourier comparison uses \path{sep_rn_robustness_10_8h.json}; selected gate records are \path{sep_gateG2.json} and \path{sep_gateG2wp.json}; the two-origin overlap is \path{sep_quadrature_overlap_10_8g.json}. Their checksums and those of the stored numerical inputs appear in the accompanying manifest.

**Imported from prior work.** The chain includes failed gates, a causal/advanced sign correction, and correction of a harness missing the retained static-guard direction. Corrected successes are not attributed to earlier failed implementations. Request notes preserve amendments and results. This revision is not a new blinded analysis or an independent audit of every preregistration timestamp.

**Proven.** At the stored two-day truncated-grid maximum, \(\widehat\beta=5.51335\times10^{-12}\) and \(\sigma_F=8.56900\times10^{-11}\). Equation (\ref{eq:interval}) with K=10 reproduces \(U=1.6795275\times10^{-9}\). This arithmetic identity uses stored estimates, not a new raw-data inference. A different nuisance prior, phase domain or covariance would require a separately specified analysis and could change the result.


## Follow-through registration and reproduction

The nuisance plan was committed at e790643, with the pre-execution standard-library CDF amendment at b8d931b; its result is 3506522. The main coverage design is 9412142. Estimated covariance was registered at 12d9f9a after the main coverage outcomes; both results are recorded at 6217da2. The force/drive matching design is 1e1cd54 and its result is 8edbe84. Repository notes retain the designs, failures and exact scope. This is an internally registered follow-through, not an externally blinded replication.

**Imported from prior work.** Main coverage seed 2026090902 and estimated-covariance seed 2026090903 are fixed. Full-space contrasts remove the nuisance means; low-dimensional sufficient coordinates retain the signal/correlated-noise span, with the orthogonal white-noise residual drawn as a chi-square variable of its declared dimension. Direct quadratic-form and subspace controls validate the compression. The oracle regression screen is a debugging control, not a criterion selected to prove astrophysical adequacy.

**Proven.** The force-matching check verifies reciprocal forces, total energy loss, branch stiffness, the inertial transfer residual, common-charge factorization, companion feedback and the leading inverse-radius expansion. It also checks phase-closure invariance against the frozen parameter dictionary. It does not simulate a star or certify all terms in the relativistic timing observable.


The comparator and phase/state designs were registered at fb269b9, with the comparator result at 65f4b65. Coarse phase, gap and derivative diagnostics were preserved at 5b9dba2 before the phase-refinement follow-through; the completed phase/state result is 7636d72. These deterministic follow-throughs use the fixed recorded covariance amplitudes and retain their non-detection scope. The repository programs are \path{verification/comparator_audit.py}, \path{verification/phase_state_audit.py} and \path{verification/phase_refinement.py}; their output JSON files are separately hashed.


### Remaining-lever follow-through and live provenance

The Request 12 plan was committed at abb973e before its analyses. Corrected leading-drive output and joint-region calibration were frozen at 6d7f3f3 before independent validation; their base seeds are 2026090912 and 2026090913, with per-condition offsets fixed in the program. Signal translations leave covariance selection unchanged in the controls. The known-scale covariance-envelope result is separate from empirical REML calibration.

The archived live baseline and isolated exponential build were recorded at 150fa43. All zero-drive preflights reproduce the recorded residuals exactly. The new source changes only the prescribed dynamic coupling, selecting an exponential when its positive time constant is provided; it does not add a fully matched gravity theory. The original full nonlinear displacement was rejected outside the bound-orbit domain; its local follow-through was registered at a08c36a before execution. Numerical results are recorded at 05d7af9. Neither old REQUEST10 verdicts nor files were overwritten.

The new symbolic, frozen-array and external-runtime producers and return arrays are identified in the revision manifest. Runtime reproduction uses the documented external WSL source/data release and isolated directories. The manuscript source package compiles without that engine. This internally registered follow-through is not a blinded astrophysical replication or independent verification of the original data reduction.
