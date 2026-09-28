# Request 11.3 — scalar-charge matching and physical-drive obstruction

Design commit: 1e1cd54. `python symbolic/physical_matching.py` passes 14 exact/algebraic and frozen-parameter checks; output: `outputs/research-completion/physical-matching.json`. No timing integration or stellar-structure calculation was performed.

## Action, force and dissipation

Status: Imported from prior work. A dynamical scalar-charge worldline action and its compact-body matching are established in [Khalil et al. (2022)](https://arxiv.org/html/2206.13233v2), especially Sections II–IV. Their leading monopole radiation force couples charge derivatives; for a single varying charge and fixed companion charges it reduces to a damping term. The inertia and potential coefficients require stellar-structure and mode information. These are prior ingredients, not a new gravity theory here.

Status: Counterexample candidate. In units G=c=1, consider the conditional leading Newtonian reduction

    L = sum_A m_A^0 v_A^2/2 + I Qdot_p^2/2 - V(Q_p)
        + sum_{A<B} (m_A^0 m_B^0 + Q_A Q_B)/r_AB,
    R = Gamma Qdot_p^2/2,  I>0, Gamma>0.

Only Q_p varies; white-dwarf charges are held fixed as a stated approximation. The background scalar can be absorbed into a linear term in V. This is a leading conservative orbital plus dissipative internal-state model, not a complete relativistic timing theory.

Status: Proven. Varying positions and Q independently gives reciprocal central pair forces with fractional coupling Delta_pj=Q_p Q_j/(m_p m_j), and

    I Qddot_p + Gamma Qdot_p + V'(Q_p) = sum_j Q_j/r_pj.

The total orbital-plus-state energy obeys Edot=-Gamma Qdot_p^2 for time-independent parameters. There is no additional derivative from an already substituted Q(r): Q is an independent coordinate during variation. Substituting an assumed mass readout into a force law without this step is not equivalent. Leading constant orbital masses omit the higher-order mass/redshift effects of the relativistic worldline action.

## Controlled one-pole reduction

Status: Proven. At a stable equilibrium Q0, let kappa=V''(Q0)>0. Linearization gives I deltaQddot+Gamma deltaQdot+kappa deltaQ=delta phi. For V=c2 Q^2/2+c4 Q^4/24, kappa=c2+c4 Q0^2/2. The full harmonic transfer and its reduced form are

    H(omega) = 1/(kappa-I omega^2+i Gamma omega),
    H0(omega) = 1/(kappa+i Gamma omega).

If epsilon_I=I omega^2/|kappa+i Gamma omega|<1, then |(H-H0)/H0|<=epsilon_I/(1-epsilon_I). A chosen small bound must hold at every retained carrier. The one-pole approximation is a controlled band approximation, not a consequence of dissipation alone. Globally separated overdamped roots additionally require I*kappa/Gamma^2 << 1, and the rapid homogeneous mode must have decayed. At omega*tau~1 these criteria have comparable scale, but neither follows just from choosing tau.

Status: Proven. If a_i=Q_i/m_i=a_o=Q_o/m_o=a_w, define deltaU=sum_j m_j(1/r_pj-mean(1/r_pj)). Then delta phi=a_w deltaU, the two pulsar-pair readouts agree, and

    tau = Gamma/kappa,
    deltaDelta(omega) = B deltaU(omega)/(1+i omega tau),
    B = a_w^2/(kappa m_p),
    F=deltaU/Ustar  =>  beta=B Ustar.

B>=0 on this stable equal-charge branch. The generic signed phenomenological beta is a larger model. The minimal slow-charge realization has c_Y=0; fitting a free c_Y allows a fast response or a broader nuisance comparator and is not an independently derived coefficient of this single-mode action.

Status: Proven. Unequal companion charge/mass ratios give deltaDelta_pi-deltaDelta_po=(a_i-a_o)deltaQ_p/m_p. A common nonzero pair modulation is then impossible unless a_i=a_o. Separate pair-response columns are needed; the existing six common-coupling columns cannot be assumed to supply them. A responsive companion with susceptibility C_j introduces a leading feedback stiffness shift -sum_j C_j/r_pj^2. Neglect requires its magnitude small relative to kappa; otherwise effective stiffness, coupled modes and possibly instability change the transfer.

Status: Proven. Small nonlinear corrections require, for example, |V'''(Q0)deltaQ|/(2 kappa)<<1 and |V''''(Q0)deltaQ^2|/(6 kappa)<<1. Close to kappa=0, a large static susceptibility alone does not ensure these conditions, small inertia or negligible companion feedback. A one-pole response has monotonically decreasing squared magnitude 1/(1+omega^2 tau^2): it has no orbital quasi-resonance. Tau=Gamma/kappa is not a model-independent inverse particle mass.

## Physical drive and the archived phase error

Status: Proven. Let r=x_p-x_i, R=x_b-x_o, where b is the inner center of mass, f=m_i/(m_p+m_i), and r_po=|R+f r|. For aligned coplanar orbits, retain terms linear separately in eccentricity and f*a_in/a_out and discard their products. With lambda_p=n_in(t-tasc_p), lambda_b=n_out(t-tasc_b), M_p=lambda_p-varpi_p and M_b=lambda_b-varpi_b,

    deltaU = A_in cos M_p + A_out cos M_b
             - A_dif cos(lambda_p-lambda_b) + higher terms,
    A_in=m_i e_in/a_in, A_out=m_o e_out/a_out,
    A_dif=m_o f a_in/a_out^2.

The minus sign follows directly from expanding 1/|R+f r|. The vector orientation must be kept consistent with the published pulsar and inner-center-of-mass ascending-node definitions. The positive-amplitude difference carrier therefore has phase lambda_p(0)-lambda_b(0)+pi. Higher eccentricity terms, mixed terms and inclination corrections are not validated by the three-carrier truncation.

Status: Imported from prior work. [Voisin et al. (2025), footnote 3 and Table 4](https://arxiv.org/html/2411.10066v2) define tasc=t_pericenter-P*varpi/(2*pi) and list e*cos(varpi), e*sin(varpi). Matching the frozen parameter values identifies this release's eta as e*cos(varpi) and kappa as e*sin(varpi); names must not be interpreted using a different implementation's convention.

Status: Proven. Applying that definition to the stored parameters gives varpi_p=-0.12326821 and varpi_b=-0.10004792 radians. Relative to the archived auxiliary F_PHYS dictionary, the leading physical phases differ by 1.69406454, 1.67084424 and pi radians for inner, outer and difference carriers. The phase closure C=phi_dif-phi_in+phi_out is invariant under a common origin shift because n_dif=n_in-n_out. Its physical value is pi+varpi_p-varpi_b=3.11837236 radians modulo 2*pi, whereas the archived auxiliary and unit-drive families have C=0. No common time shift repairs this discrepancy in this model.

Decision: withdraw the archived auxiliary beta_phys values as bounds on this scalar-potential-driven realization. Preserve the raw historical output; the unit-drive beta table remains a defined phenomenological benchmark, not a physical drive translation. A future corrected physical analysis must prescribe unequal amplitudes, phases, omitted harmonics and model-error tolerances before recomputing intervals and their coverage. Replacing only phase labels or rescaling the old beta cannot perform this correction.

## Completed result and residual physical boundary

Status: Counterexample candidate. The force-level realization is now explicit, dissipative and reciprocal, with a coefficient map and controlled reduction conditions. This completes the conditional EFT matching task. It does not establish that a real J0337 neutron star and both white dwarfs realize this regime.

Status: Conjectural. A numerical EOS-to-body matching must determine the branch, susceptibility, inertia, damping and companion charges at fixed baryon number, and bound higher post-Newtonian, radiation, retardation and nonlinear corrections over the observing span. The current data do not identify those inputs uniquely. An astrophysical exclusion would additionally need the corrected drive and a validated timing likelihood. These are explicit limits on the completed conditional paper, not manufactured positive outcomes.

Classification: theorem progress (force/transfer and phase-closure boundaries) and loophole progress (conditional scalar-charge realization; historical physical interpretation fails its phase gate).
