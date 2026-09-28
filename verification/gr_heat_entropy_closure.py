"""Heat-only quadratic entropy closure, retaining the full transport term.

Two spatially constant trial proper heat times are constructed from the
initial cone bounds. They are not measurements or microscopic calibrations.
"""
import json,sys,urllib.request
from fractions import Fraction as F
import numpy as np
import sympy as sp
import gr_heat_initial_tangent as old

g=old.g;OUT=g.OUT/'gr-heat-entropy-closure'
URL='https://arxiv.org/pdf/astro-ph/9609119'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();old.verify()
    with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'maartens-9609119.pdf').write_bytes(response.read())
    assert (OUT/'maartens-9609119.pdf').read_bytes().startswith(b'%PDF')
    paths=[g.ROOT/'verification/gr_heat_entropy_closure.py',old.OUT/'manifest.json',old.char.OUT/'coefficients.npz',
        g.OUT/'gr-heat-initial-constraints/initial-data.npz',g.OUT/'gr-heat-fluid-frame/rate-coefficients.npz',
        g.OUT/'gr-heat-relaxation-boundary/thresholds.npz',g.OUT/'gr-opacity/evaluation.npz',OUT/'maartens-9609119.pdf']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='7c9cc12',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source=dict(classification='Imported from prior work',url=URL,author='Roy Maartens',year=1996,
            equations=['2.17','2.20','2.22','2.24'],scope='Heat-only specialization, no viscosity/heat cross coupling or vorticity. The lecture explicitly cautions that neglected terms require checking.'),
        time_rule='Find max_i h_i*lambda_star_i using the previously defined local zero-flux cone rule. Use 2 and 3 times this single frozen maximum as spatially and temporally constant proper times in separate declared models.',
        root_bracket_halfwidth_powers=[45,40,35,30,25,20],
        initial_target='Q_F=Q from the existing discrete face-to-cell target, as in the previous tangent. No independent continuous gradient certificate.',
        physical_heat_time_calibrated=False,full_GR_evolution=False))
    symbolic()


def symbolic():
    q,K,T,tau,dq,force,divbeta=sp.symbols('q K T tau dq force divbeta',real=True)
    # Natural units for this identity. beta=tau/(K*T^2).
    entropy=-q*force/T**2-tau*q*dq/(K*T**2)-q*q*divbeta/2
    full_force=-(tau*dq+q)/K-T*T*q*divbeta/2
    assert sp.simplify(entropy.subs(force,full_force)-q*q/(K*T*T))==0
    truncated=sp.simplify(entropy.subs(force,-(tau*dq+q)/K))
    assert sp.simplify(truncated-(q*q/(K*T*T)-q*q*divbeta/2))==0
    assert sp.simplify(truncated.subs(divbeta,3/(K*T*T))+q*q/(2*K*T*T))==0
    b,r,d,e,j,a,kr,kt,v=sp.symbols('b r d e j a kr kt v',real=True)
    # Heat equation divided by tau, a=h/tau; baryon equation eliminates
    # expansion from the derivative of beta=tau/(K*T^2).
    M=sp.Matrix([[1,0,0,0],[0,b,2*j,0],[0,0,1,1],[-j*kr/2,-j*kt/2,a,1]])
    N=sp.Matrix([[0,0,1,0],[0,0,e,1],[r,d,2*j,0],[0,a,0,0]])
    coeff=[b*(1-a)-j*j*kt,
        j*(2*(d-b)-2*a-(1-e)*kt/2+b*kr/2),
        -(b*r+d*e+a*(1-e-d))+j*j*kt,
        j*(2*a+r*kt/2-d*kr/2),r*a]
    polynomial=sum(value*v**(4-i) for i,value in enumerate(coeff))
    assert sp.simplify((N-v*M).det()-polynomial)==0
    assert sp.simplify(M.det()-coeff[0])==0
    w,h,Q,R,X,Y,V,TT=sp.symbols('w h Q R X Y V TT',real=True)
    # X=rho log rate, Y=conditional zero-vdot temperature log rate,
    # V=coefficient of vdot in the temperature log rate.
    den=w*(TT-h)+TT*Q*kt*V/2
    vdot=TT*(R-Q*(kr*X+kt*Y)/2)/den
    qdot=R-w*vdot;tlog=Y+V*vdot
    assert sp.simplify(TT*qdot+w*h*vdot-TT*Q*(kr*X+kt*tlog)/2)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        entropy_current='In c=1 units S^mu=s*n*u^mu+q^mu/T-[tau/(2*K*T^2)]*q_alpha*q^alpha*u^mu. Assume baryon/energy conservation, a differentiable Gibbs relation and this heat-only quadratic nonequilibrium entropy ansatz.',
        divergence='D_mu S^mu=-q.(grad T+T*a)/T^2-tau*q.Dq/(K*T^2)-(q^2/2)*D_mu[(tau/(K*T^2))*u^mu].',
        full_heat_law='tau*Delta*Dq+q=-K*(Delta*grad T+T*a)-(K*T^2/2)*q*div[(tau/(K*T^2))*u]. This gives div S=q^2/(K*T^2)>=0 for positive K,T within the declared model.',
        truncated_boundary='Omitting the last transport term leaves div S=q^2/(K*T^2)-(q^2/2)*div(beta*u). A declared local jet with div(beta*u)=3/(K*T^2) makes this negative. This is an algebraic counterexample to an unconditional entropy guarantee, not an observed violation in the stored star.',
        frozen_constant_tau='For constant tau, define kr=dlnK/dlnrho+1 and kt=dlnK/dlnT+2. At initial v=0 and fixed composition, dln(beta)/dt+N*theta=-kr*dlnrho/dt-kt*dlnT/dt.',
        initial_equations='For Q_F=Q: Q_dot+w*v_dot=R; tau*Q_dot+w*h*v_dot=(tau*Q/2)*(kr*rho_log_dot+kt*T_log_dot); T_log_dot=Y+V*v_dot. The extra term changes the initial tangent, so the old truncated tangent is not reused.',
        principal_matrix='Variables (delta ln rho,delta ln T,delta v,delta Q/w), M=[[1,0,0,0],[0,b,2j,0],[0,0,1,1],[-j*kr/2,-j*kt/2,a,1]], N=[[0,0,1,0],[0,0,e,1],[r,d,2j,0],[0,a,0,0]], a=h/tau; system M*y_t+c*N*y_x=lower-order terms.',
        scope='Formal entropy identity and local principal matrix of this chosen quadratic heat-only model. It neither calibrates the ansatz/coefficients nor supplies physical EOS derivative error, nonequilibrium validity, nonlinear well-posedness, viscosity/convection/reactions, atmosphere or a finite GR trajectory.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    fields=dict(np.load(old.char.OUT/'coefficients.npz'))['coefficients']
    initial=dict(np.load(g.OUT/'gr-heat-initial-constraints/initial-data.npz'))
    frame=dict(np.load(g.OUT/'gr-heat-fluid-frame/rate-coefficients.npz'))
    threshold=dict(np.load(g.OUT/'gr-heat-relaxation-boundary/thresholds.npz'))
    op=dict(np.load(g.OUT/'gr-opacity/evaluation.npz'))['values']
    # K=16*sigma*T^3/(3*rho*kappa_eff), on the bound opacity branch.
    kr=-op[:,1];kt=5-op[:,2]
    c=g.c.gr.C*100;w=threshold['enthalpy_density'];Q=initial['proper_heat_flux_cgs']/c;j=Q/w
    R=2*c*initial['lapse']*initial['radial_extrinsic_curvature_per_cm']*Q
    frozen=[[F(float(x)) for x in row] for row in fields];tau_star=F(0)
    for B,P,D,E,H in frozen:
        sound=P+D*E/B
        lamstar=max(F(1),(B+1-E-D-P)/(B*(1-sound)),(2*B+1-E-D)/(B*(2-sound)))
        tau_star=max(tau_star,H*lamstar)
    times=[2*tau_star,3*tau_star];save('constant-times.json',dict(classification='Counterexample candidate',
        exact_proper_seconds=list(map(str,times)),floating_display=list(map(float,times)),physical_calibration=False))
    brackets=[];rates=[];failures=[];minimum_den=None;maxres=F(0)
    for index,(B,P,D,E,H) in enumerate(frozen):
        J=F(float(j[index]));KR=F(float(kr[index]));KT=F(float(kt[index]))
        for model,TAU in enumerate(times):
            A=H/TAU
            coeff=[B*(1-A)-J*J*KT,J*(2*(D-B)-2*A-(1-E)*KT/2+B*KR/2),
                -(B*P+D*E+A*(1-E-D))+J*J*KT,J*(2*A+P*KT/2-D*KR/2),P*A]
            def poly(x):
                value=F(0)
                for a in coeff:value=value*x+a
                return value
            Mf=np.array([[1,0,0,0],[0,float(B),2*float(J),0],[0,0,1,1],[-float(J*KR/2),-float(J*KT/2),float(A),1]])
            Nf=np.array([[0,0,1,0],[0,0,float(E),1],[float(P),float(D),2*float(J),0],[0,float(A),0,0]])
            roots=np.linalg.eigvals(np.linalg.solve(Mf,Nf));assert np.all(roots.imag==0),(index,model,roots)
            centers=sorted(F(float(x)) for x in roots.real);found=None
            for power in plan['root_bracket_halfwidth_powers']:
                half=F(1,2**power);trial=[(x-half,x+half) for x in centers]
                if (coeff[0]>0 and trial[0][0]>-1 and trial[-1][1]<1 and
                    all(trial[i][1]<trial[i+1][0] for i in range(3)) and
                    all(poly(lo)*poly(hi)<0 for lo,hi in trial)):
                    found=trial;break
            if found is None:failures.append(dict(cell=index,model=model));continue
            brackets.append(dict(cell=index,model=model,bounds=[[str(lo),str(hi)] for lo,hi in found]))
            WW=F(float(w[index]));QQ=F(float(Q[index]));RR=F(float(R[index]))
            X=F(float(frame['density_rate'][index]));Y=F(float(frame['conditional_zero_vdot_logT_dot'][index]))
            V=F(float(frame['coefficient_vdot_in_logT_dot'][index]))
            den=WW*(TAU-H)+TAU*QQ*KT*V/2;assert den>0
            minimum_den=den/(WW*TAU) if minimum_den is None else min(minimum_den,den/(WW*TAU))
            vdot=TAU*(RR-QQ*(KR*X+KT*Y)/2)/den;qdot=RR-WW*vdot;Tdot=Y+V*vdot
            residual=max(abs(qdot+WW*vdot-RR),abs(TAU*qdot+WW*H*vdot-TAU*QQ*(KR*X+KT*Tdot)/2))
            maxres=max(maxres,residual)
            # The dropped term relative to Q, including lapse, is a useful
            # finite diagnostic even when it is smaller than unity.
            correction=TAU*(KR*X+KT*Tdot)/(2*F(float(initial['lapse'][index])))
            truncated_vdot=TAU*RR/(WW*(TAU-H))
            rates.append([index,model,float(vdot),float(qdot),float(Tdot),float(V*vdot),float(correction),float(truncated_vdot)])
    save('root-brackets.json',dict(classification='Proven',rows=brackets))
    values=np.array(rates);np.savez_compressed(OUT/'initial-rates.npz',values=values,kr=kr,kt=kt,Q=Q,j=j,momentum_forcing=R)
    result=dict(classification='Proven',completed=True,cells=len(fields),certified_cones=len(brackets),
        exact_frozen_cones_passed=not failures,failures=failures,exact_frozen_rate_residual=str(maxres),
        minimum_scaled_tangent_denominator_display=float(minimum_den),
        maximum_abs_vdot_per_coordinate_second=float(abs(values[:,2]).max()),
        maximum_abs_added_logT_rate_per_coordinate_second=float(abs(values[:,5]).max()),
        maximum_abs_dropped_transport_term_over_Q=float(abs(values[:,6]).max()),
        nonzero_initial_acceleration_cells=int(np.count_nonzero(values[values[:,1]==0,2])),
        caveat='Exact arithmetic applies only to stored coefficients, constructed times and forcing numbers. Formal entropy positivity assumes the differentiable Gibbs relation and the chosen quadratic heat-only entropy current. No physical closure calibration or native derivative error enclosure.',
        physical_EOS_certified=False,physical_transport_calibrated=False,
        nonequilibrium_finite_wavenumber_stability_certified=False,full_GR_evolution=False)
    save('result.json',result)
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('FULL QUADRATIC ENTROPY HEAT CLOSURE',result,flush=True)
    assert not failures and maxres==0,result;verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    result=json.loads((OUT/'result.json').read_text())
    assert result['exact_frozen_cones_passed'] and result['exact_frozen_rate_residual']=='0'
    print('PASS QUADRATIC ENTROPY HEAT CLOSURE SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
