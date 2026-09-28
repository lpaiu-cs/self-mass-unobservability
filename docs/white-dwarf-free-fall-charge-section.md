# Draft section (revision 7): a worked white-dwarf calculation — a transient retarded scattering signal and orbital-timescale bounds

Status (2026-09-28): revision 6 passed the confirmation review of revision 5 (fable5.1: accept; gpt-6-astra and opus5.5: accept after minor revision); revision 7 adds the layer-by-layer thermal-relaxation estimate of phase 287 after a single self-review (user decision). The review records are `notes/REQUEST282_INDEPENDENT_REVIEW_REVISION_KO.md` through `notes/REQUEST288_RELAXATION_ESTIMATE_SECTION46_KO.md`. The draft goes after Section 4.5 of `paper/manuscript.md`. The earlier integration (a322c6462) was reverted (b40864a22). The proposed companion edits are at the end.

Labels follow the manuscript convention:

- **Imported from prior work**: literature and repository computations
- **Counterexample candidate**: the proposed model realization
- **Proven**: conditional mathematics
- **Conjectural**: interpretations

---

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
- **Imported from prior work.** *Deeper layers.* We estimated the thermal response layer by layer. Each layer relaxes on the thermal time of the layers above it, and its fully relaxed limit is bracketed by the isothermal exponent \(P_{\rm gas}/P\). Solving the static structural response again with all layers faster than a given time relaxed gives the relaxation strength of each deeper shell. The layers whose thermal time equals \(1/\omega\) lie about 2,600 km deep for the inner orbit and 6,600 km for the outer one, above about \(2\times10^{-9}\) and \(1.4\times10^{-7}\) of the mass. Summing the absolute shell strengths with the Debye weight \(\omega\tau/(1+\omega^2\tau^2)\) gives a lag of \(4.0\times10^{-9}\,|\mathcal S_{\rm struct}|\) at the inner orbital frequency and \(3.3\times10^{-7}\,|\mathcal S_{\rm struct}|\) at the outer one. **Conjectural.** This estimate does not use the two assumptions below. It is not a non-adiabatic calculation, and its relaxed limit can differ from thermal equilibrium by factors of order unity.
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
- \(4.1\times10^{-10}\) (physical-drive endpoint)
- \(1.7\times10^{-9}\) (K=10, truncated)
- \(3.5\times10^{-9}\) (K=10, full)
- \(1.57\times10^{-7}\) (\(K\approx934\))

**Conjectural.** We did not compute the timing response of either part in these two non-common channels, including the nuisance projection. Reading the structural comparison as a sensitivity statement would require their timing sensitivity to lie within about eight orders of magnitude of the common template's, which was not established. The outer white dwarf's own internal states were not computed.

**Conjectural.** In summary, a scalar pulse leaves in the displaced matter of the inner white dwarf a transient retarded monopole signal, with no first-order permanent charge as \(t\to\infty\) if all its modes are damped. At orbital timescales the charge follows the field instantaneously through \(\beta_s\), and the displacement state adds a computed adiabatic structural term that is small. Its lagged part is bounded under the two stated assumptions on the thermal relaxation, and a layer-by-layer estimate puts it at \(4.0\times10^{-9}\) and \(3.3\times10^{-7}\) of that term at the inner and outer orbital frequencies; a non-adiabatic calculation of the deep thermal response was not made. The calculation therefore neither establishes nor excludes an orbital-timescale internal state with a coupling at the Section 5 scales, as the Section 3 benchmark posits. The remaining assumptions are:

- linear response and the weak-field charge \(a_i=\alpha_s(\varphi)\);
- \(|a_p|\le1\) in the scale comparison;
- a Newtonian radial model with a flat-exterior readout;
- no neutral modes and positive damping of all modes;
- neglect of interior Born scattering, where the core compactness is about \(10^{-4}\);
- for the lag, either a thermal relaxation strength bounded by \(\mathcal S_{\rm struct}\) with a single relaxation time, or the layer-by-layer relaxation estimate;
- the unfixed \(\varphi_\infty\) exponent of the \(10^{-8}\) direct part;
- a non-gray correction with simplified LTE continuum opacity and a fixed pulse kernel;
- a slowly rotating white dwarf;
- a mixed hydrogen–helium envelope, whereas the observed star has a hydrogen (DA) atmosphere \cite{kaplan2014j0337}.

---

## Proposed companion edits

**Abstract, after the damped scalar-charge sentence:**
**Conjectural.** A worked calculation for the inner white dwarf of PSR J0337+1715 finds that the scalar-driven displacement of its own matter leaves a transient retarded monopole signal, with no first-order permanent charge as \(t\to\infty\) if all modes, including slow thermal ones, are damped. At orbital timescales its computed adiabatic structural term is small. A bound on its lagged part assumes a thermal relaxation strength no larger than that term and a single relaxation time, and the corresponding timing sensitivity was not computed.

**Section 6, after the paragraph on the physical target:**
**Conjectural.** Section 4.6 works through one internal state of the inner white dwarf, the scalar-driven displacement of its own matter. That state leaves a transient retarded monopole signal. Under the stated assumptions, its structural modulation of the pair factors lies about eight orders of magnitude below the stored Section 5 scales. The white dwarfs' instantaneous zero-lag responses could reach the smallest of those scales only for \(|a_p|\gtrsim0.5\), which the static SEP limits allow only if \(|a_o|\lesssim5\times10^{-6}\); they are not lags. Neither timing response was computed, so this is a scale comparison, not an exclusion. Other internal states of the white dwarfs, and the neutron-star charge of Section 4.3, remain open.

**Data and code availability, new sentences:**
The white-dwarf calculation of Section 4.6 was added after the unified input state. It is recorded in Korean-language notes \path{notes/REQUEST244_*} through \path{notes/REQUEST288_*} and in manifests under \path{outputs/direct-eos-gr33/}; all are bound in the revision manifest. Each manifest lists the scripts and results of its step by hash. Verdict, classification and derived-bound fields written into manifests before the independent reviews are superseded by the errata in \path{notes/REQUEST284_*} through \path{notes/REQUEST286_*}. The main results map to these manifests:

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

**References:** restore \texttt{kaplan2014j0337}, \texttt{ransom2014triple} and \texttt{bertotti2003cassini} as in commit a322c6462.

**README and revision manifest:** record that the PDF and submission archive predate Section 4.6 until TeX is rebuilt.

**Introduction:** no sentence. The introduction has no section roadmap, and the abstract carries the summary.
