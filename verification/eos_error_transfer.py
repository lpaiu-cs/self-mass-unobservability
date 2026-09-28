"""Proven conditional EOS error transfer; no measured plasma certificate.

The constraints are eliminated in fixed coordinates before applying the
strong-convexity argument. Boundary active-set changes need separate analysis.
"""
import json
from fractions import Fraction as Q
import sympy as s
import direct_eos_gr as g

OUT=g.OUT/'eos-error-transfer'


def main():
    assert not OUT.exists();OUT.mkdir()
    theta,z,mu,eps,omega=s.symbols('theta z mu eps omega',positive=True)
    f=mu*z*z/2;perturbed=f+eps*s.sin(omega*theta)
    assert s.diff(perturbed,z,2)==mu
    assert s.diff(perturbed,theta).subs(theta,0)==eps*omega
    assert s.diff(perturbed,theta,2).subs(theta,s.pi/(2*omega))==-eps*omega**2
    # Exact independent quadratic minimization with one internal coordinate.
    f0=Q(3)*theta**2/2+theta*z+Q(4)*z*z/2
    f1=f0+Q(1,20)*theta**2+Q(1,10)*theta*z+Q(1,10)*z*z
    z0=-theta/4;z1=-Q(11,42)*theta
    assert s.diff(f0,z).subs(z,z0)==0 and s.diff(f1,z).subs(z,z1)==0
    h0=s.diff(f0.subs(z,z0),theta,2);h1=s.diff(f1.subs(z,z1),theta,2)
    # Constant A,B,C blocks: d=state shift, no spatial block Lipschitz term.
    da,db,dc,b0,m,m0=Q(1,10),Q(1,10),Q(1,5),Q(1),Q(4),Q(4)
    bound=da+db*(2*b0+db)/m+b0*b0*dc/(m*m0)
    assert abs(h1-h0)<=bound
    gradient_error=abs(s.diff(f1.subs(z,z1),theta).subs(theta,1)-
        s.diff(f0.subs(z,z0),theta).subs(theta,1))
    # The full interval |z|<=1 covers both minimizers and their segment.
    ez=Q(3,10);shift_bound=ez/m
    eg=Q(1,5)+Q(11,10)*shift_bound
    assert gradient_error<=eg and abs(z1.subs(theta,1)-z0.subs(theta,1))<=shift_bound
    r,t=s.symbols('r t',real=True);F=s.Function('F')(r,t)
    P=s.exp(r)*s.diff(F,r);entropy=-s.exp(-t)*s.diff(F,t);u=F-s.diff(F,t)
    assert s.simplify(s.diff(P,r)-s.exp(r)*(s.diff(F,r)+s.diff(F,r,2)))==0
    assert s.simplify(s.diff(P,t)-s.exp(r)*s.diff(F,r,t))==0
    assert s.simplify(s.diff(u,t)-(s.diff(F,t)-s.diff(F,t,2)))==0
    assert s.simplify(s.diff(entropy,t)+s.exp(-t)*(s.diff(F,t,2)-s.diff(F,t)))==0
    report=dict(classification='Proven',passed=True,
        assumptions=['A fixed convex feasible set in independent internal coordinates z; all required derivatives exist on a common convex neighbourhood.',
            'The true constrained free energy has F_zz >= mu I with certified mu>0 throughout that neighbourhood, and an interior minimizer.',
            'The saved numerical point is feasible and has approximate-model stationarity residual norm at most r.',
            'Physical-model perturbation derivative bounds and block Lipschitz constants are independently supplied, not estimated from two finite steps.'],
        stationary_state_bound='d=||z_true-z_hat|| <= (r+epsilon_z)/mu',
        free_energy_value_bound='|f_true-F0(theta,z_hat)| <= epsilon_0+(r+epsilon_z)^2/(2 mu)',
        reduced_gradient_bound='||gradient f_true-F0_theta(theta,z_hat)|| <= epsilon_theta+L_theta_z*d',
        blocks='A=F_theta_theta; B=F_theta_z; C=F_zz; H=A-B*C^{-1}*B^T. Bound deltaA,deltaB,deltaC by perturbation second derivatives plus the corresponding spatial block Lipschitz constant times d.',
        reduced_hessian_bound='||H-H0|| <= deltaA+deltaB*(2*b0+deltaB)/mu+b0^2*deltaC/(mu*mu0), with ||B0||<=b0 and C0>=mu0 I.',
        proof='Strong monotonicity of the z-gradient gives the stationary-state bound. Strong convexity gives the value suboptimality bound. The envelope theorem and mean-value theorem give the gradient bound. Implicit differentiation gives the Schur complement; add and subtract B0*C^{-1}*B0^T and use C^{-1}-C0^{-1}=C^{-1}(C0-C)C0^{-1}.',
        thermodynamic_coordinates='F is Helmholtz energy per baryon gram, r=ln rho_B and t=ln T, at fixed nuclear composition. P=exp(r)*F_r, s=-exp(-t)*F_t, u=F-F_t.',
        component_error_transfer=dict(P='rho*E1_r',dP_dlnrho='rho*(E1_r+E2_rr)',
            dP_dlnT='rho*E2_rt',entropy='E1_t/T',cvT='E1_t+E2_tt',energy='E0+E1_t'),
        counterexample='F0=mu*z^2/2 and F1=F0+epsilon*sin(omega*theta) have the same strictly convex internal minimizer and uniformly bounded value difference epsilon, while the first and second parameter derivative differences grow as epsilon*omega and epsilon*omega^2. Value-only model bounds cannot certify EOS derivatives.',
        exact_quadratic_control=dict(H0=str(h0),H1=str(h1),Hessian_error=str(abs(h1-h0)),
            Hessian_error_bound=str(bound),gradient_error=str(gradient_error),gradient_error_bound=str(eg),
            minimizer_shift_bound=str(shift_bound)),
        numerical_roundoff='Certified floating evaluation errors must be added to the residual, derivative blocks and reported outputs; this algebra does not provide them.',
        unverified_for_actual_EOS=['Uniform projected strong convexity','Complete equilibrium fractions and stationarity residuals',
            'Physical perturbation first and second derivatives','Block Lipschitz bounds','Active-set/phase-boundary coverage'],
        actual_physical_or_continuous_EOS_certificate=False)
    (OUT/'result.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('PASS conditional EOS value/gradient/Hessian transfer and exact counterexample',flush=True)


if __name__=='__main__': main()
