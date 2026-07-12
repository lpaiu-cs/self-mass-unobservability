# Finite-Family Collapse of Free-Fall Self-Energy Couplings: A Fixed-Order Worldline-EFT Theorem, its Uniqueness No-Go, and Sharp Boundary Escapes

**Status:** draft manuscript (Paper A of the two-paper split; static free-fall sector)
**Repository:** `lpaiu-cs/self-mass-unobservability` (main branch)
**Author:** Juneyoung, Kim
**Note:** Draft 2026-07-12, revision 1 (adversarial-review corrections: the gradient block \(\nabla E\) is an STF-3 octupole, collapsing the electric normal-form dimension from seven to five and shifting every family dimension — see Sections 3.2, 3.3, 7; framing and citation fixes throughout).

---

## Abstract

We ask a sharp version of a simple question — *how much of a body's gravitational
self-energy can be observable in its free fall?* — and answer it as a structural
theorem about the worldline effective field theory (EFT) of a self-gravitating
body. Two elementary no-go results first rule out the literal reading of
"self-mass unobservability": at quadrupole order the center-of-mass equation is
independent of any self-subtraction parameter, and removing self-gravity inside a
Newtonian barotrope decompactifies the star rather than producing a new
equilibrium. What survives is a restricted EFT in which the body's self-energy can
enter only through body-dependent sensitivity couplings of an already-formed
object. We then prove, in the parity-even, nonspinning, local free-fall sector at
fixed operator order \(\Delta \le 4\) and under an explicit assumption ledger
\(A1\)–\(A8\): **(i)** the admissible parity-even scalar sensitivity operator space
is finite, reduces to an explicit finite normal-form basis of computed dimension (five for the electric sector), and — under locality
and analyticity — the monopole response collapses to a finite Taylor jet in
finitely many scalar coordinates (the *positive finite-family collapse theorem*);
**(ii)** this does *not* single out a unique minimal sector — admitting any
unsuppressed primitive family (magnetic, scalar, vector, or higher-rank
symmetric-trace-free) produces a new low-order witness that survives the reduction
quotient and is linearly independent of the electric basis, a *class-limited
family-admission no-go*; and **(iii)** each of the four load-bearing assumptions
\(A3, A4, A5, A8\) is sharp: dropping it admits an explicit counterexample, with
\(A4\) additionally admitting a finite state-augmented salvage. All load-bearing
counts and identities are reproduced by independent computational methods, and the
verification process itself exposed and corrected two errors in earlier
bookkeeping — an undercount of the rank-4 survivor dimension, and an over-count of
the gradient sector in every dimension quoted previously (the tidal gradient
\(\nabla E = \partial^3\Phi\) is a totally symmetric, vacuum trace-free octupole,
not a generic tensor) — so the dimensions quoted here are the twice-corrected
values. We are explicit about the theorem's boundaries: nothing
here is claimed for parity-odd, spinning, nonlocal, clock, or orbital-timescale-state
sectors, for a unique physically privileged minimal sector, or beyond the fixed
cutoff.

---

## 1. Introduction

### 1.1 The question

Let \(A\) be a compact, self-gravitating body carrying gravitational self-energy
\(\Omega_A<0\), and write its dimensionless self-gravity fraction as
\[
s_A \;\equiv\; \frac{|\Omega_A|}{m_A c^2}\;>\;0 .
\]
For the Earth \(s_\oplus\sim 5\times10^{-10}\); for a neutron star
\(s_{\rm NS}\sim 0.1\)–\(0.2\). In general relativity the strong equivalence
principle (SEP) makes \(s_A\) invisible to free fall: passive, inertial, and active
masses coincide regardless of self-energy content. A violation is classically
parameterized in the Nordtvedt manner,
\[
\frac{m_{P,A}}{m_{I,A}} \;=\; 1 + \sigma_1\,s_A + \sigma_2\,s_A^2 + \cdots ,
\tag{1}
\]
with \(\sigma_1\) essentially the Nordtvedt parameter and \(\sigma_2\) a
strong-field extension. The motivating slogan of this program — *"matter cannot
sense its own mass"* — is the statement that the self-energy contribution to
(1) is unobservable. Our question is structural rather than observational: **within
a controlled effective description, how many independent self-energy-dependent
free-fall couplings can there be, and is any one of them singled out?**

### 1.2 Two elementary no-gos force a restricted EFT

The literal reading of the slogan fails twice, and both failures are used only as
imported inputs here (Appendix A):

* **COM decoupling.** Model self-subtraction by weighting the internal
  self-interaction by \((1-\lambda)\). Because the self-energy is a
  translation-invariant functional of internal coordinates, it never appears in the
  center-of-mass Euler–Lagrange equation; through quadrupole order,
  \(M\ddot X^i = -M\,\partial_i\Phi_{\rm ext} - \tfrac12 I^{jk}\partial_i\partial_j\partial_k\Phi_{\rm ext} + O(\ell^3\nabla^4\Phi)\),
  with \(\lambda\) absent. "Not counting" self-energy internally does nothing to
  free fall.
* **Structural inconsistency.** Weighting self-gravity by \(\alpha=1-\lambda\) inside
  Newtonian hydrostatic equilibrium and holding \((M,K)\) fixed for an \(n=1\)
  polytrope gives \(R\propto\alpha^{-1/2}\), \(\rho_c\propto\alpha^{3/2}\),
  \(P_c\propto\alpha^{3}\); as \(\lambda\to1\) the object decompactifies into an
  unbound cloud rather than reaching a new equilibrium.

Together these motivate the restricted setting adopted here: an effective
description in which self-energy enters through *sensitivity couplings of an
already-formed body*, i.e. the coefficients of (1) and their higher analogues. We
emphasize that the two no-gos rule out alternatives but do not by themselves
entail this description; the restricted EFT is *adopted as premises* — the
worldline-EFT existence, sector, and equilibrium assumptions \(A2, A3, A6\) of
Section 2 — and everything downstream is conditional on that ledger. The
remainder of this paper classifies those couplings.

### 1.3 Results

Working in the worldline EFT (Section 2) at fixed order \(\Delta\le 4\) under the
ledger \(A1\)–\(A8\):

1. **Positive finite-family collapse (Section 3).** The parity-even scalar
   sensitivity operator space is finite; after explicit reductions it has a finite
   normal-form basis; under locality \(A3\) and analyticity \(A5\) the monopole
   response is a finite Taylor jet; higher multipoles are carried by finitely many
   Wilson coefficients, separate from the monopole sensitivity data.
2. **Uniqueness no-go (Section 4).** No unique minimal sector is singled out:
   admitting any unsuppressed primitive family yields a new low-order witness.
3. **Sharp boundary escapes (Section 5).** \(A3, A4, A5, A8\) are each necessary —
   dropping any one admits an explicit counterexample — with a finite
   state-augmented salvage for \(A4\).
4. **Independent verification (Section 7).** Independent computational methods
   reproduce every load-bearing count and identity. The verification corrected two
   errors in earlier bookkeeping: a rank-4 survivor undercount (six omitted
   higher-degree mixed operators), and a gradient-sector over-count — the earlier
   enumeration modeled \(\nabla_k E_{ij}\) as a generic tensor, missing the Schwarz
   total symmetry of \(\partial_k\partial_i\partial_j\Phi\) and the vacuum trace
   condition already imposed on \(E\) itself, which had inflated the electric
   dimension from five to seven and every family dimension with it.

---

## 2. Worldline setup and assumptions

We treat the body as a point worldline dressed by finite-size operators built from
the external field evaluated on the worldline. The leading tidal object is the
electric tidal tensor \(E_{ij}=\partial_i\partial_j\Phi_{\rm ext}\), a symmetric
trace-free (STF) rank-2 tensor; its worldline time derivative \(D_\tau E_{ij}\) and
spatial gradient \(\nabla_k E_{ij}\) enter at higher order, as do higher-multipole
and additional primitive families. The gradient block carries kinematic structure
that the enumeration must respect: \(\nabla_k E_{ij}=\partial_k\partial_i\partial_j
\Phi_{\rm ext}\) is totally symmetric in \((k,i,j)\) by equality of mixed partials,
and in the external vacuum (\(\nabla^2\Phi_{\rm ext}=0\), the same condition that
makes \(E\) itself trace-free) it is trace-free on every index pair — an STF-3
octupole with seven components, not a generic \(({\rm STF\text{-}2})\otimes
{\rm vector}\) object with fifteen. The **monopole sensitivity** is the piece of the
worldline action that renormalizes how the body's mass responds to the external
field — the operator analogue of the coefficients in (1) — and is a scalar function
of finitely many local invariants \(Y^I\).

Operators are graded by a weight \(\Delta\) (Appendix A): each primitive field
insertion and each worldline/spatial derivative carries positive weight, and only
\(\Delta\le\Delta_{\max}=4\) is retained. The theorem domain is fixed by the
assumption ledger, whose entries are **premises**, not derived results:

| ID | Assumption | Sharp role |
| --- | --- | --- |
| \(A1\) | quasi-static body on the orbital timescale | removes fast internal readout variables |
| \(A2\) | nearly spherical, nonspinning, parity-even minimal sector | fixes the operator search space |
| \(A3\) | a local worldline EFT exists | supplies the operator language |
| \(A4\) | no orbital-timescale internal state variable | rules out hidden \(\chi\) coordinates |
| \(A5\) | analytic monopole response near the reference background | defines a finite Taylor jet |
| \(A6\) | self-bound equilibrium before external perturbation | separates formation from coupling |
| \(A7\) | fixed operator order \(\Delta\le4\) | prevents a silent all-orders claim |
| \(A8\) | locally finite primitive-family spectrum below the cutoff | finitely many species with \(w\le\Delta_{\max}\) |

The classification criterion we use is explicit: an assumption is *load-bearing*
if removing it changes a named conclusion, and the interesting load-bearing cases
are those whose failure modes are structurally informative. By that criterion
\(A3, A4, A5, A8\) are load-bearing and individually non-redundant — dropping any
one admits the explicit counterexample of Section 5. \(A7\) is also load-bearing
for finiteness, but trivially so: removing the cutoff reopens the infinite
all-orders tower, and there is no finite-basis claim left to test, so it is kept
as a stated scope choice rather than a Section 5 entry. \(A1, A6\) delimit the
sector, and \(A2\) is the sector restriction whose cost is quantified by Theorem 2:
each family excluded by \(A2\) would, if admitted unsuppressed, contribute new
survivors.

---

## 3. The positive finite-family collapse theorem

**Theorem 1 (finite-family collapse).** *Inside the stated domain at fixed order
\(\Delta\le4\), the admissible parity-even local free-fall scalar operator space is
finite. After the explicit reduction rules are imposed, the reduced scalar operator
space is finite-dimensional and admits a finite normal-form basis \(\{Y^I\}\). Under
locality \(A3\) and analyticity \(A5\), the monopole response collapses to a finite
Taylor jet in those finitely many scalar coordinates. The residual higher-multipole
sector is carried by finitely many Wilson coefficients that remain separate from the
monopole sensitivity coefficients.*

The proof has four layers: (a) the irreducible family envelope is fixed (Section
3.1); (b) fixed-order candidate finiteness follows from \(A7\)–\(A8\) and positive
weights; (c) the candidate space is reduced to a finite normal-form basis by
explicit rules (Section 3.2); (d) analyticity turns the monopole response into a
finite jet, and the sensitivity/Wilson split isolates it. We state the load-bearing
content and the numbers we independently verified.

**Remark (where the content is).** Finiteness of the candidate space is, by
itself, the weakest part of Theorem 1: given a fixed order \(A7\) and a locally
finite primitive spectrum \(A8\), boundedness of a graded invariant space is
expected on general invariant-theory grounds. The theorem's nontrivial content is
therefore *not* "the space is finite" but three sharper facts that finiteness alone
does not supply: **(i)** the reduced space is a *specific* normal-form quotient of
computed dimension — five for the electric sector (2), and the tabulated
dimensions for each admitted family (Section 3.3) — obtained by quotienting the raw
contractions by total worldline derivatives, the lower-order equation of motion,
the kinematic identities of the gradient block, and the STF identities (3);
**(ii)** that quotient is exact, i.e. the redundancy is
removed by a *closed, finite* set of algebraic identities rather than an open-ended
reduction; and **(iii)** the irreducible family envelope itself closes (Section
3.1). Assumptions \(A7\)–\(A8\) deliver only boundedness; layers (a), (c), (d)
deliver the actual collapse, and it is the explicit dimensions and identities — not
the bare word "finite" — that the verification in Section 7 is designed to certify.

### 3.1 The irreducible family envelope

In three spatial dimensions the only \(SO(3)\) irreducibles are the STF-\(L\)
(spin-\(L\), dimension \(2L+1\)) tensors; every Cartesian tensor decomposes into
STF-\(L\) pieces, and a symmetric rank-\(r\) tensor is \({\rm STF}_r\oplus{\rm
STF}_{r-2}\oplus\cdots\) so that its trace descendants are strictly lower-rank STF
(Appendix B). Consequently the primitive-family envelope closes on the parity-even
**scalar (rank 0), vector (rank 1), rank-2 STF, and genuine rank-\(L\ge3\) STF**
classes; trace-descendant and even-dual mixed-symmetry sectors reduce to STF +
traces, and parity-odd/pseudo sectors are excluded by \(A2\). The rank-\(L\ge3\)
tower is closed uniformly rather than rank-by-rank, by elementary rank arithmetic
valid for all \(L\): the full contraction \(X2=X\!\cdot\!X\) exists at weight 2 for
every rank (smallest self witness); a \(\delta\)-contraction scalar of two STF
tensors requires equal rank, so the mixed quadratic \(E\!:\!X\) exists only at
\(L=2\) (no weight-2 mixed witness for \(L\ne2\)); and \(X\otimes X\) contains an
STF-2 for every \(L\ge1\), so the mixed cubic \(E X X\) exists at weight 3
(smallest mixed witness). We additionally confirmed this uniform structure by
exact character integrals for \(L=2,\dots,12\), well past the audited \(L\le6\)
(Section 7).

### 3.2 The electric-sector normal-form basis

For the minimal electric sector \(\{E, D_\tau E, D_\tau^2 E, \nabla E\}\) at
\(\Delta\le4\), all parity-even scalar contractions reduce — modulo total
worldline derivatives, the lower-order equation of motion (\(a_i=0\)), the
kinematic identities of the gradient block, and the STF algebraic identities — to
the five-dimensional normal-form basis
\[
\{\,E2,\; E3,\; E2^2,\; \dot E2,\; \nabla\!E2\,\},
\tag{2}
\]
where \(E2=E_{ij}E^{ij}\), \(E3=E_i{}^jE_j{}^kE_k{}^i\), and
\(\nabla\!E2=(\nabla_kE_{ij})(\nabla^kE^{ij})\). The gradient sector is
one-dimensional for a kinematic reason: as noted in Section 2, \(\nabla_kE_{ij}\)
is an STF-3 octupole — totally symmetric by equality of mixed partials and
trace-free in the external vacuum — and a single \(SO(3)\) irreducible carries
exactly one quadratic invariant. The two further contractions admitted by a
generic \(({\rm STF\text{-}2})\otimes{\rm vector}\) model,
\({\rm div}E2=(\nabla_iE^{ij})(\nabla^kE_{kj})\) and
\({\rm mixedGrad}E2=(\nabla_kE_{ij})(\nabla^iE^{kj})\), respectively vanish
(vacuum trace) and coincide with \(\nabla\!E2\) (Schwarz symmetry). An earlier
version of this enumeration modeled the gradient block generically and quoted a
seven-dimensional electric basis; the error — declining for \(\nabla E\) the very
vacuum condition already used to make \(E\) trace-free, and overlooking the Schwarz
symmetry entirely — was caught in adversarial review and corrected throughout
(Section 7). The only genuinely algebraic
reductions are the STF identities
\[
{\rm Tr}\,E^4 = \tfrac12({\rm Tr}\,E^2)^2,\qquad
{\rm Tr}\,B^4 = \tfrac12({\rm Tr}\,B^2)^2,\qquad
{\rm Tr}(EBEB) = (E\!:\!B)^2 + \tfrac12{\rm Tr}\,E^2\,{\rm Tr}\,B^2 - 2\,{\rm Tr}(E^2B^2),
\tag{3}
\]
the third of which we verified as an exact symbolic identity. The basis (2) is both
**complete** (no weight-\(\le4\) parity-even scalar is missed) and **linearly
independent** (the coefficient matrix has rank 5); we verified completeness by two
independent methods, independence by explicit polynomial rank, and the two
kinematic identities both exactly (on the explicit 7-parameter STF-3
parametrization) and numerically to machine precision (Section 7).

### 3.3 Survivor dimensions of the audited families

Admitting a primitive family enlarges the surviving normal-form set. Under the
theorem's own reduction rules the survivor dimension of each audited sector is
finite and equals, for the electric baseline plus one family:

| Sector | family rank / parity | survivor dimension |
| --- | --- | --- |
| electric \(E\) | — | 5 |
| \(E\) + magnetic \(B\) | rank 2, parity-odd | 16 |
| \(E\) + \(B\) + scalar \(S\) | rank 0 | 30 |
| \(E\) + vector \(V\) | rank 1 | 15 |
| \(E\) + rank-3 STF \(T\) | rank 3 | 17 |
| \(E\) + rank-4 STF \(Q\) | rank 4 | **23** |
| \(E\) + rank-5 STF \(U\) | rank 5 | 17 |
| \(E\) + rank-6 STF \(Z\) | rank 6 | 21 |

The table reflects the two corrections found in verification (Section 7). First,
the original hand-built enumeration capped the mixed \(E/Q\) sector at degree 2 in
\(E\), omitting six genuine higher-degree survivors (\(EEQ, QQQ, E^3Q, EQ^3,
E\,D_\tau E\,Q, \nabla E\,\nabla Q\)); restoring them lifts the new-sector rank
from 12 to 18 (Appendix C). Second, the gradient-kinematics correction of Section
3.2 removes the two spurious electric-gradient entries from every sector, and
additionally removes the mixed invariant \((\nabla_iE^{ij})(\nabla_jS)\) from the
scalar sector (the \(\ell=3\) octupole shares no irreducible with the scalar
gradient); the cross-gradient \(\nabla E\,\nabla Q\) itself survives, because
\(\nabla Q\) does contain an \(\ell=3\) piece. Admitted families' own gradient
blocks are kept generic: they are independent primitives with no assumed potential
structure, so no analogue of the Schwarz/vacuum collapse applies to them. Every
number shown is reproduced independently by an exact \(O(3)\) character identity
(Section 7). In all cases the space is finite — the content of
Theorem 1's finiteness at the family level.

### 3.4 Monopole jet and the sensitivity/Wilson split

Given the finite scalar basis \(\{Y^I\}\), analyticity \(A5\) makes the monopole
response \(m_A(Y)\) a finite Taylor jet in the \(Y^I\); locality \(A3\) is what
allows the response to be a function of the *instantaneous* invariants at all. The
higher-multipole Wilson coefficients remain separate from the monopole sensitivity
jet and are not merged into the same coefficient family. This completes Theorem 1.

---

## 4. The uniqueness no-go

Theorem 1 does **not** identify a unique, physically privileged minimal sector.

**Theorem 2 (class-limited family-admission no-go).** *Fix the domain and
reduction rules of Theorem 1, and admit alongside \(E\) any one primitive family
from the audited classes — scalar, vector, rank-2 STF (magnetic), or genuine
rank-\(L\ge3\) STF — at its natural weight assignment and with untuned
coefficients. Then the enlarged sector contains, already at \(\Delta\le4\), at
least one new operator that survives the full reduction quotient and is linearly
independent of the electric basis (2). Hence no theorem can single
out a unique minimal sector without adding explicit suppression, ordering, or
background-restriction assumptions family by family.*

The existence of *some* new invariant upon admitting a new field is, of course,
immediate (the family supplies its own quadratic norm); the theorem's content is
that this invariant **survives the reduction quotient** — it is not a total
derivative, not EOM-reducible, and not algebraically dependent on the electric
basis — which we verified by explicit rank computation for every audited class
(Section 7). The smallest new survivors are:

| Admitted family | smallest new witness | weight |
| --- | --- | --- |
| magnetic rank-2 STF \(B_{ij}\) | \(B2\) | 2 |
| bare scalar \(S\) | \(S\) | 1 |
| derivative-only scalar | \(\dot S2\) | 4 |
| vector \(V_i\) | \(V2\) | 2 |
| rank-\(L\ge3\) STF \(X_L\) | \(X2\) | 2 |

Because nothing in the minimal-sector setup forbids these families, and each
admission adds a genuine low-order coupling, minimal-sector *uniqueness* is
obstructed while fixed-order *finiteness* (Theorem 1) survives. This is the precise
formal counterpart of the informal statement that "self-mass unobservability" does
not pick out one canonical coupling.

---

## 5. Sharp boundary escapes

Each of \(A3, A4, A5, A8\) is individually non-redundant: dropping it, with the
rest of the ledger held fixed, admits an explicit counterexample. (This
establishes sharpness of each assumption separately, not minimality or
completeness of the set; and \(A7\) is absent here because its removal breaks
finiteness trivially rather than informatively — see Section 2.) We verified all
four counterexamples (Section 7).

| Dropped | Smallest explicit counterexample | Layer broken | Replacement data |
| --- | --- | --- | --- |
| \(A5\) | \(m_A(Y)=m_0+\alpha\,e^{-1/Y^2}\Theta(Y)\) | finite analytic Taylor jet | non-Taylor monopole germ |
| \(A3\) | one-coordinate causal power-law kernel \(K_\gamma(s)\propto s^{-\gamma}\), \(0<\gamma<1\) | local reduction to instantaneous \(m_A(Y)\) | memory kernel / spectral data |
| \(A4\) | one-state local analytic \(\chi\) model, \(\ddot\chi+\omega_\chi^2\chi=(g_\chi/\mu_\chi)Y\) | \(Y\)-only monopole reduction | finite state-space data \((Y^I,\chi^a)\) |
| \(A8\) | infinite low-weight STF tower \(\{T^{(n)}\}\), each \((T^{(n)})^2\) at weight 2 | candidate operator finiteness before reduction | no finite pre-reduction catalog |

* **\(A5\).** The smooth-flat model is \(C^\infty\) with every derivative vanishing at
  \(Y=0\) yet nonzero for \(Y>0\); its Taylor jet is trivial while the function is
  not, so no finite analytic jet captures the response. Physically this escape is
  not merely a mathematical sharpener: a non-analytic monopole germ is what a sharp
  phase transition in the dense-matter equation of state would produce — near a
  strong first-order or near-critical transition the bulk response to an external
  tidal control parameter can turn non-analytic — so \(A5\) marks a candidate
  physical regime for compact stars with a stiff-to-soft transition, not only a
  contrived function.
* **\(A3\).** An exponential memory kernel is finitely Markovianizable (a single local
  auxiliary state), but a power-law kernel has a non-rational transfer function
  \(\sim p^{\gamma-1}\) with a branch point at \(p=0\); no finite local state-space
  reproduces it.
* **\(A4\).** A driven internal oscillator makes the response depend on history and
  initial data, not on the instantaneous \(Y\). This escape has a **salvage**: if
  finitely many local state variables are kept explicit, a finite *state-augmented*
  collapse theorem survives on \((Y^I,\chi^a)\); adiabatically eliminable states
  collapse back to a \(Y\)-only model and do not count as sharp counterexamples. This
  is the seam Paper B develops dynamically.
* **\(A8\).** If infinitely many primitive species survive below the cutoff, the
  candidate operator space is infinite before any reduction; local weight-spectrum
  finiteness below \(\Delta_{\max}\) is exactly what Theorem 1 needs, and it is
  weaker than a globally finite catalog.

---

## 6. Scope: what is and is not claimed

**Claimed.** A fixed-order (\(\Delta\le4\)), parity-even, nonspinning, local
free-fall theorem: finite-family collapse of the scalar sensitivity sector (Theorem
1); a class-limited no-go against minimal-sector uniqueness (Theorem 2); and the
sharpness of \(A3, A4, A5, A8\) (Section 5), with a finite state-augmented salvage
for \(A4\).

**Not claimed.** No theorem for parity-odd, spinning, nonlocal, clock, or
orbital-timescale-state sectors; no unique physically privileged minimal sector; no
all-orders closure beyond the fixed cutoff; no universal mixed-pattern theorem for
all higher-rank families; and no empirical, runtime, or observational claim — those
belong to the separate data-side program and to Paper B (the dynamic loophole).

---

## 7. Independent verification

Every load-bearing count and identity in Sections 3–5 was reproduced by methods
independent of the repository's own enumeration scripts. The verification is not
merely confirmatory: it exposed and corrected two errors, and the second of them
carries a methodological lesson we state explicitly.

* **Algebraic identities (3).** Verified numerically over random STF tensors and, for
  the mixed quartic, as an exact symbolic identity.
* **Gradient kinematic identities.** \({\rm div}E2=0\) and \({\rm mixedGrad}E2
  =\nabla\!E2\) verified exactly on the explicit 7-parameter STF-3 parametrization
  and numerically to machine precision over random octupoles; a diagnostic control
  confirms that the *generic* gradient model reproduces the defective rank-3
  gradient sector, isolating the old error to the modeling choice.
* **Electric basis (2).** *Completeness* confirmed two ways — a from-scratch
  \(\delta\)-only contraction enumerator over every weight-\(\le4\) signature
  (fourteen signatures with nonzero invariant content), and an
  \(O(3)\) character (Molien) count anchored by the pure-\(E\) Hilbert series
  \(1/((1-t^2)(1-t^3))\) — and *independence* by explicit polynomial rank 5.
* **Family completeness.** Survivor dimensions computed by an exact character
  identity, \({\rm survivor}(w)=\dim{\rm inv}_{\rm trunc}(w)-\dim{\rm inv}_{\rm
  prom}(w-1)\), reproducing all eight audited sectors
  (\(5,16,30,15,17,23,17,21\)) after the two corrections below were applied to the
  enumeration; the numeric \(D_\tau\)-quotient method independently confirms the
  brute-force-feasible sectors.
* **Correction 1 (rank-4 undercount).** The original hand-built rank-4 candidate
  list capped the mixed \(E/Q\) sector at degree 2 in \(E\); the character identity
  exposed six omitted operators, which were constructed explicitly and shown
  nonzero, rotation-invariant, and non-total-derivative, lifting the rank-4
  new-sector rank from 12 to 18.
* **Correction 2 (gradient-sector over-count).** An external adversarial review
  found that *every* method — including the nominally independent ones — had
  modeled \(\nabla E\) as a generic \(({\rm STF\text{-}2})\otimes{\rm vector}\)
  object, dropping the Schwarz symmetry of \(\partial^3\Phi\) and the vacuum trace
  condition already imposed on \(E\) itself. Agreement between two methods is no
  protection against a wrong ansatz they share; the corrected kinematics is now
  encoded in all methods, which again agree, and the diagnostic control above
  keeps the failure mode visible. This moved the electric dimension from seven to
  five and shifted every family dimension in Section 3.3.
* **Family-envelope closure.** The STF-\(L\) irrep census and the uniform tower
  self-witness structure (\(X2\) at weight 2, \(EXX\) at weight 3, no weight-2 mixed
  witness for \(L\ge3\)) proven by rank arithmetic (Section 3.1) and confirmed by
  exact character integrals for \(L=2,\dots,12\).
* **Boundary escapes.** All four counterexamples verified from their explicit models
  (flat Taylor jet for \(A5\); non-rational transfer function for \(A3\); adiabatic
  collapse coefficients and genuine state dependence for \(A4\); infinite pre-reduction
  count for \(A8\)).

The verification scripts and a written audit are in `verification/` and
`symbolic/` in the repository; a DOI-archived snapshot of the verification tree
will accompany submission so that referees need no repository access. This
manuscript's status labels follow the taxonomy in
`docs/status-labels.md`: `Proven` marks a computationally or analytically established
claim, `Assumption` marks the \(A1\)–\(A8\) premises, and bookkeeping/provenance is
tagged `Note`.

Finally, because the character/Molien count is order-agnostic — the same
Hilbert-series machinery applies at any \(\Delta_{\max}\) — extending the
electric-sector count to \(\Delta\le6\) is mechanical, and is the natural
robustness check that the \(\Delta\le4\) cutoff (\(A7\)) is a scope choice rather
than an artifact of the order. We flag it as the first planned extension.

---

## 8. Discussion

The theorem reframes "self-mass unobservability" from a slogan into a bounded,
falsifiable EFT statement: for a static, parity-even, nonspinning body the space of
self-energy-dependent free-fall couplings is *finite at each fixed order* but *not
unique*, and the finiteness rests on exactly four load-bearing assumptions, each of
which is sharp. The two most physically suggestive escapes are \(A4\) (an internal
dynamical variable on the orbital timescale) and \(A8\) (an unbounded low-weight
family tower); the former is the subject of the companion dynamic paper, where a
second-order internal mode generates orbital sidebands that evade the static
collapse. That companion analysis (Paper B, `paper/paper-B-dynamic-sep-limit.md`
on the dynamic branch) is now complete: instrumenting the salvaged \(A4\) state
data \((\beta,\tau_\chi)\) on PSR J0337+1715, it bounds the amplitude of any
lag-responding SEP oscillation at \(|\delta\Delta|<1.7\times10^{-9}\)
(a 95% interval evaluated at a systematically inflated width, anchor tier
\(K_{\rm dyn}=10\); \(\tau_\chi=2\) d, worst drive phase, over the window
\(\tau_\chi\in[2,500]\) d). The uniqueness no-go is, we think, the honest headline: the idea does not
select a canonical coupling, and any stronger "minimal-sector" claim must import an
explicit suppression principle.

**Relation to existing tidal worldline-EFT bases.** The operator language used
here is the standard worldline EFT of extended bodies (Goldberger–Rothstein), and
finite-size operator towers built from \(E_{ij}\), \(B_{ij}\), and their
derivatives have long been organized for post-Newtonian and effective-one-body
tidal dynamics — notably by Bini, Damour and Faye; Levi and Steinhoff;
Henry, Faye and Blanchet; and, within modern scattering-amplitude/worldline
programs, Kälin and Porto. Those constructions classify *tidal response*
operators (Love-number couplings and their spin/derivative refinements) for
waveform modeling. What this paper adds is orthogonal to that goal: a
*sensitivity* reading of the same operator space — which couplings can carry
dependence on the body's own self-energy in free fall — together with the exact
fixed-order reduction quotient with computed dimensions, the family-admission
no-go, and the sharp assumption-boundary map. We are not aware of a prior
classification of this sector under this question; the operator kinematics,
however, deliberately match the established tidal-EFT conventions so that the
two programs can be composed.

**The empirical anchor.** The theorem is structural, but its leading coordinate is
already measured. The linear coefficient \(\sigma_1\) in (1) is, to leading
post-Newtonian order, the Nordtvedt parameter, which lunar laser ranging constrains
at the \(|\eta|\sim10^{-4}\) level (Will 2014). In the normal-form language
\(\sigma_1\) is the coefficient of the lowest scalar coordinate — the
linear-in-\(E2\) monopole response — so existing Nordtvedt/LLR data already pins the
*first* survivor. What the collapse theorem adds is the structure above it: a finite
but non-unique tower (\(\sigma_2\) and the tabulated family survivors) that current
weak-field data does not separately resolve, and whose strong-field entries
(\(s_{\rm NS}\sim0.1\)–\(0.2\)) are where the additional couplings first become
numerically relevant. Translating the full \(\{Y^I\}\) basis into a data-facing
bound table is the natural bridge to the empirical program and is left to that
track.

---

## Appendices (support)

* **A. Worldline setup and imported exclusions.** COM decoupling and the
  internal-structure no-go (imported); power counting and the operator grading;
  worldline reduction. Sources: `../lemmas/01`–`03`, `docs/power-counting.md`.
* **B. Irreducible family-envelope closure.** Cartesian irrep reduction, trace
  absorption, mixed-symmetry gate. Sources: `docs/irreducible-envelope-theorem.md`,
  `../lemmas/50`–`52`.
* **C. Family-admission no-go and thresholds.** Audited witness classes, the uniform
  STF self-witness threshold, and the corrected rank-4 basis. Sources:
  `docs/family-admission-theorem.md`, `docs/stf-self-witness-theorem.md`,
  `docs/family-class-table.md`, `../symbolic/r4_survivor_rank_check.py`,
  `../verification/rederive_rank4.py`.
* **D. Fixed-order operator closure.** Candidate finiteness, normal-form quotient,
  monopole jet, sensitivity/Wilson split. Sources:
  `docs/finite-family-collapse-theorem.md`, `../lemmas/53`–`56`.
* **E. Boundary counterexamples.** The four assumption-drop failures and the
  state-augmented salvage. Sources: `docs/boundary-escape-map.md`, `../lemmas/57`–`71`,
  `counterexamples/`.

## Selected references (draft; convert to BibTeX before submission)

1. K. Nordtvedt, *Equivalence principle for massive bodies*, Phys. Rev. **169**, 1017 (1968).
2. C. M. Will, *The confrontation between general relativity and experiment*, Living Rev. Relativity **17**, 4 (2014).
3. W. D. Goldberger and I. Z. Rothstein, *An effective field theory of gravity for extended objects*, Phys. Rev. D **73**, 104029 (2006). arXiv:hep-th/0409156.
4. T. Damour and G. Esposito-Farèse, *Tensor-multi-scalar theories of gravitation*, Class. Quantum Grav. **9**, 2093 (1992).
5. É. É. Flanagan and T. Hinderer, *Constraining neutron-star tidal Love numbers with gravitational-wave detectors*, Phys. Rev. D **77**, 021502 (2008).
6. B. Kol and M. Smolkin, *Black hole stereotyping: induced gravito-static polarization*, JHEP **02**, 010 (2012).
7. R. A. Porto, *The effective field theorist's approach to gravitational dynamics*, Phys. Rep. **633**, 1 (2016).
8. T. Damour, *The problem of motion in Newtonian and Einsteinian gravity*, in *300 Years of Gravitation*, eds. Hawking and Israel (Cambridge, 1987).
9. D. Bini, T. Damour, and G. Faye, *Effective action approach to higher-order relativistic tidal interactions in binary systems and their effective one body description*, Phys. Rev. D **85**, 124034 (2012). arXiv:1202.3565.
10. M. Levi and J. Steinhoff, *Spinning gravitating objects in the effective field theory in the post-Newtonian scheme*, JHEP **09**, 219 (2015). arXiv:1501.04956.
11. Q. Henry, G. Faye, and L. Blanchet, *Tidal effects in the gravitational-wave phase evolution of compact binary systems to next-to-next-to-leading post-Newtonian order*, Phys. Rev. D **102**, 044033 (2020). arXiv:2005.13367.
12. G. Kälin and R. A. Porto, *Post-Minkowskian effective field theory for conservative binary dynamics*, JHEP **11**, 106 (2020). arXiv:2006.01184.
