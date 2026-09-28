# Observable-state count and physical matching

## Observable-pole uniqueness in a restricted positive class

Status: Proven. Let H(iw)=c0+integral a(dt)/(1+iwt), with real c0 and a nonnegative finite relaxation measure on t>0. At two known frequencies 0<l<h, define S_w=-Im H(iw)/w and

    R=(h^2*S_h-l^2*S_l)/(h^2-l^2),
    D=(Re H(il)-Re H(ih))/(h^2-l^2),
    Q=(S_l-S_h)/(h^2-l^2).

These are respectively the zeroth, first and second moments of the positive measure dmu=t*a(dt)/[(1+l^2*t^2)(1+h^2*t^2)]. Hence RQ-D^2=R^2 Var_mu(t)>=0. If R>0, equality holds exactly when the observable finite relaxation measure is supported at one t0=D/R. Its amplitude is a0=R*(1+l^2*t0^2)*(1+h^2*t0^2)/t0; c0 follows from either real part. Instantaneous mass at t=0 is inseparable from c0. This is stronger than merely rejecting a fast spectrum. It requires exact calibrated complex response data; arbitrary unknown per-frequency drive or readout phases invalidate the construction.

Status: Proven. This identifies one observable relaxation time, not the number of physical internal variables. Decoupled or unexcited states and degenerate modes can be added. At finite noise, two positive poles at t0+-epsilon with weights a0/2 approach the single-pole response with exact difference

    a0*s^2*epsilon^2 / [(1+s*t0)*((1+s*t0)^2-s^2*epsilon^2)].

For 0<epsilon<t0 this is a stable positive model and tends to zero quadratically. No uniformly separated test of exactly one versus arbitrarily close two poles exists at finite precision without a minimum separation/weight condition. The two-atom moment determinant is w1*w2*(t1-t2)^2.

Status: Proven. Without positive-residue restrictions, a stable strictly proper addition eta*product_k(s^2+w_k^2)/(s+lambda)^(2K+1), lambda>0, vanishes at all measured +/-iw_k but alters the response elsewhere. This is a causal stable finite-frequency indistinguishability construction; it is not generally a reciprocal positive relaxation spectrum. No unrestricted state-count uniqueness follows from the stored carriers.

## Matching what equilibrium information actually supplies

Status: Proven. Within the damped scalar-charge candidate, fixing V(Q), the drive coupling and equilibrium solution leaves the family Gamma=kappa*t0, t0>0 with identical static susceptibility 1/kappa and arbitrary relaxation time t0. Equilibrium information alone therefore does not imply a fast-rate gap. Positive inertia may also vary without changing the equilibrium solution. Once a specific radiation/matter theory is fixed, Gamma need not be freely adjustable: the counterexample addresses insufficiency of equilibrium information in the admitted EFT class.

Status: Imported from prior work. Khalil et al. match potential coefficients from stellar equilibrium sequences and the kinetic coefficient using scalar-led radial modes. Their leading coupled monopole radiation force is proportional to minus the sum of charge velocities (their Eq. 18). The corresponding two-charge damping matrix has a null direction; the positive-definite fast-gradient comparator is therefore not automatic even for this explicit scalar-field realization. Source: [Khalil et al., Phys. Rev. D 106, 104016, Sections III--IV](https://arxiv.org/html/2206.13233v2).

Status: Proven. The 2x2 all-ones damping matrix has eigenvalues 0 and 2. Additional inertia, dissipation channels or reduction of the undamped sector must be addressed before applying the SPD-gradient spectral-gap theorem to that coupled model.

Status: Conjectural. Actual J0337 matching still requires a specified gravity coupling, EOS and stellar equilibrium branch, companion charges, complex response/mode data, and radiation/feedback matching. None is selected by the recorded orbital timing parameters. We have not fabricated these inputs or declared an EOS-derived bound. The completed lever is an exact insufficiency theorem and a concrete matching route; numerical J0337 body matching remains incomplete.

Checks: symbolic/state_identifiability.py. Classification: theorem progress and loophole boundary.
