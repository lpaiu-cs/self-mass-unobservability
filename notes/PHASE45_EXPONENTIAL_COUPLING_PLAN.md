# Phase45 working derivation — not yet an implemented coupled solver

Status: Conjectural. The Phase44 native paths failed their time-response gates.
The fixed-matter carrier reproduces their odd scalar to relative ~1e-8, while
the finest midpoint carrier differs from exact finite-grid propagation by 31%.
Do not rerun those histories on a finer grid before resolving the fast carrier.
The following is a concrete implementation route; its actual energy identity
must be derived and checked before a production run.

## Canonical scaling and the old discrete gradient

Status: Conjectural. Let `R=star.rf[-1]`, `tau=t*c/R`,
`w_i=star.scalar_volume_i/(4*pi*R^3)`, `x_i=sqrt(w_i)*phi_i`,
`v_i=sqrt(w_i)*R*Pi_i`, and use dimensionless mass `mu=m_surface/R`.
The scalar canonical skew matrix is `Q=[[0,I],[-I,0]]` in `(x,v)`.
The old midpoint field equation has `Delta x = h_tau*k*v_mid`,
`k=H*(b_old+b_new)/2`. Its gradient in v is therefore `k*v_mid`.

Status: Conjectural. With conductances
`c_face=scalar_area*kf/distance`, form the tridiagonal negative Laplacian `D`
on field cells. The scalar-gradient Hessian in x is
`-D/(4*pi*R*sqrt(w_i*w_j))`, and the matter term in its discrete gradient is

`-G_over_c4/R * H_i*V_i*beta/(a_old+a_new)`
` * (a_old*trace_old*phi_old+a_new*trace_new*phi_new)/sqrt(w_i)`.

Status: Conjectural. Holding each current metric/trace evaluation fixed makes
this affine in the midpoint y: `g=Gmat*y_mid+bvec`.
The matter diagonal is
`-2*G_over_c4/R * H*V*beta*a_new*trace_new/((a_old+a_new)*w)`;
the affine term is
`-G_over_c4/R * H*V*beta*(a_old*trace_old-a_new*trace_new)`
` * phi_old/((a_old+a_new)*sqrt(w))`.
These factors must be checked against `def_resolved_scalar_pulse.wave` rather
than accepted from this note. The old global mass identity uses the same H,
kf and reciprocal `gstar`; it should generalize if the new scalar update is
skew with respect to exactly this discrete gradient.

## Exponential discrete gradient, without inverting a singular Hamiltonian

Status: Proven. For a constant symmetric reference Hessian `M` and constant
skew `Q`, let `A=exp(h Q M)` and
`B=integral_0^h exp(t Q M) Q dt`. Compute A and B from the upper blocks of
`expm([[h Q M,h Q],[0,0]])`; no inverse of M is needed. Where `A+I` is
invertible, `S=2 solve(A+I,B)` is skew in exact arithmetic.
The exponential discrete-gradient equation can be written

`y_new-y_old = S*g`.

For a true discrete gradient, `g^T Delta y=Delta H`, this conserves H.
When `g=M*y_mid`, it gives exactly `y_new=A*y_old`.
This is an algebraic statement under the given matrix assumptions; an
implementation must monitor `A+I` conditioning and matrix skew error.

Status: Conjectural. For the affine frozen current gradient above, solve
`(I-S*Gmat/2)*y_new=(I+S*Gmat/2)*y_old+S*bvec`.
Use binary64 factorization with long-double residual refinement, as in the
existing wave solve. Enforcing `(S-S.T)/2` may retain the energy algebra but
does not certify propagation accuracy. The nonlinear outer scalar/metric
iteration and the full native matter residual remain active.

## Autonomous forcing reservoir for the prescribed boundary

Status: Conjectural. A plain replacement of midpoint fields with an exact
carrier breaks the existing boundary-work formula. Instead, augment the
canonical variables with nine clock coordinates and their nine momenta:
one constant coordinate, and `(cos(2*pi*j*tau),sin(2*pi*j*tau))`, j=1..4.
Their quadratic clock Hamiltonian is
`sum_j 2*pi*j*(q_cos*p_sin-q_sin*p_cos)`.
Clock coordinates then follow the prescribed oscillations, independently of
their reservoir momenta; field backreaction changes only those momenta.

Status: Conjectural. Write the boundary as `b^T q_clock`, with coefficients
`[phi_boundary0+amplitude*35/128, amplitude*(-56)/128,0,`
`amplitude*28/128,0, amplitude*(-8)/128,0, amplitude/128,0]`.
At tau=0 all cos coordinates and the constant are one, all sin coordinates
and momenta are zero. Include the full outer gradient energy
`c_boundary*(b^T q_clock-phi_last)^2/(8*pi*R)` in M, not just its
field-clock cross term. The constant clock makes M singular; the block
exponential construction above handles this.

Status: Conjectural. During a current nonlinear evaluation the clock physical
gradient comes from the same outer-face energy using current kf. Add the
clock Hamiltonian gradient at the midpoint. The nonlinear remainder has no
clock momentum dependence; check that the new clock coordinates reproduce
the analytic pulse even while physical field/metric coefficients change.

Status: Conjectural. The expected integrated boundary work is now
`-R/G_over_c4 * Delta H_clock`. It should tend to the physical surface
scalar work as h tends to zero and exactly reproduce the known forced
quadratic carrier. The combined desired identity is
`Delta(m_surface/G_over_c4)-work_clock-sum(H*V*heat0*matter_residual_energy)=0`.
Prove this with the actual discrete mass chain before interpreting it as
physical conservation. Do not assign a remaining scalar defect to matter heat.

## Minimal implementation and budget sequence

Status: Conjectural. Reuse Phase44 initial Cauchy data and unprojected momentum,
the full native EOS/transport/species residual and nonlinear solver. A subclass
can bind a replacement `wave` into `PulseStar.evaluate`, retain the old frozen
tangent only as an iteration preconditioner, and save clock states separately.
All saved Phase40–44 sources and verdicts must stay immutable. Check the
initial `at` placeholder issue and refuse overwriting finished/failed outputs.

Status: Conjectural. First test a fixed-background forced carrier against the
existing exact propagator, then a manufactured nonzero metric with the energy
identity, then at most two native steps under 90 seconds. Only after those
tests and a measured budget should a finite matched comparison be planned.
No larger grid or longer physical duration is authorized by this design note.
Physical free-surface/radiative-exterior coupling remains a separate unmet
requirement even if the new finite-cavity integrator succeeds.
