"""Local coupled-fluid characteristics and finite-wave-number stability.

The stored EOS derivatives define a declared linear model. Exact arithmetic
certifies its coefficients, not native EOS errors or a physical heat time.
"""
import json,sys
from fractions import Fraction as F
import numpy as np
import sympy as sp
import gr_heat_relaxation_boundary as prior

g=prior.g;OUT=g.OUT/'gr-heat-characteristics'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def symbolic():
    b,r,d,e,h,tau,s,k,z=sp.symbols('b r d e h tau s k z',real=True)
    # x=delta ln rho_B, y=delta ln T, v=velocity/c, q=heat_flux/(c*w).
    # h=K*T/(w*c^2) in seconds. Spatial derivative below is c*d/dx.
    M=sp.Matrix([[1,0,0,0],[0,b,0,0],[0,0,1,1],[0,0,h,tau]])
    N=sp.Matrix([[0,0,1,0],[0,0,e,1],[r,d,0,0],[0,h,0,0]])
    J=sp.diag(0,0,0,1)
    dispersion=sp.expand((s*M+sp.I*k*N+J).det()/b)
    sound=r+d*e/b;D2=1-e-d+b*sound
    expected=(tau-h)*s**4+s**3+k*k*(tau*sound+h*(1-e-d)/b)*s*s+sound*k*k*s+h*r*k**4/b
    assert sp.simplify(dispersion-expected)==0
    a4,a3,a2,a1,a0=[sp.expand(dispersion).coeff(s,i) for i in [4,3,2,1,0]]
    assert sp.simplify(a3*a2-a4*a1-h*k*k*D2/b)==0
    assert sp.simplify(a3*a2*a1-a4*a1*a1-a3*a3*a0-h*k**4*(sound*D2-r)/b)==0
    square=(b*r-d*(1-d))**2/b
    assert sp.simplify((sound*D2-r).subs(e,d)-square)==0
    lam,beta=sp.symbols('lam beta',real=True)
    char=sp.simplify((N-beta*M).det().subs(tau,lam*h)/h)
    expected_char=b*(lam-1)*beta**4-(lam*(b*r+d*e)+1-e-d)*beta**2+r
    assert sp.simplify(char-expected_char)==0
    C=lam*(b*r+d*d)+1-2*d
    disc=C*C-4*b*(lam-1)*r
    positive=(lam*b*r-lam*d*d-1+2*d)**2+4*b*r*(lam*d-1)**2
    assert sp.simplify(disc-positive)==0
    assert sp.simplify(C-(lam*b*r+lam*(d-1/lam)**2+1-1/lam))==0
    A=b*(lam-1);Cgeneral=lam*b*sound+1-e-d
    light=(b+1-e-d-r)/(b*(1-sound))
    vertex=(2*b+1-e-d)/(b*(2-sound))
    assert sp.simplify(A-Cgeneral+r-b*(1-sound)*(lam-light))==0
    assert sp.simplify(2*A-Cgeneral-b*(2-sound)*(lam-vertex))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        variables='w=epsilon+P; b=rho_B*cv*T/w; r=(dP/dlnrho_B)_T/w; d=(dP/dlnT)_rho/w; e=(P-rho_B*(du/dlnrho_B)_T)/w; h=K*T/(w*c^2); lambda=tau/h; q=delta physical_heat_flux/(c*w). Coefficients are constant local equilibrium values. k=c*physical_wave_number.',
        equations='x_dot+c*v_x=0; b*y_dot+c*e*v_x+c*q_x=0; v_dot+q_dot+c*r*x_x+c*d*y_x=0; h*v_dot+tau*q_dot+q+c*h*y_x=0.',
        dispersion='(tau-h)*s^4+s^3+k^2*(tau*g+h*(1-e-d)/b)*s^2+g*k^2*s+h*r*k^4/b=0, g=r+d*e/b.',
        Routh_Hurwitz='For b,h,r>0, tau>h, g>0, 1-e-d+b*g>0, and g*(1-e-d+b*g)-r>0, all roots have Re(s)<0 for every nonzero real k. The two nontrivial Hurwitz determinants are h*k^2*(1-e-d+b*g)/b and h*k^4*(g*(1-e-d+b*g)-r)/b. At k=0 conserved roots remain at zero.',
        Maxwell_specialization='For an exact differentiable thermodynamic EOS, e=d. Then g*(1-2*d+b*g)-r=(b*r-d*(1-d))^2/b. Strict damping follows unless this square is zero. Equality is a degenerate case not covered by the strict Hurwitz claim.',
        characteristics='Squared speed/c^2 is a root of A*z^2-C*z+r=0, A=b*(lambda-1), C=lambda*(b*r+d*e)+1-e-d. If A,C,r>0, discriminant>=0, C<=2*A, and A-C+r>=0, both roots lie in (0,1].',
        Maxwell_hyperbolicity='For e=d,b*r>0,lambda>1, C=lambda*b*r+lambda*(d-1/lambda)^2+1-1/lambda>0, and discriminant=(lambda*b*r-lambda*d^2-1+2*d)^2+4*b*r*(lambda*d-1)^2>=0.',
        causal_trial_rule='When 0<g<1 choose lambda_*=max(1,(b+1-e-d-r)/(b*(1-g)),(2*b+1-e-d)/(b*(2-g))). Test lambda=2*lambda_* and 3*lambda_* with exact discriminant and all Hurwitz conditions; these are declared candidate coefficients, not physical calibrations.',
        scope='Linear constant-coefficient zero-flux local baryon-frame Cattaneo subsystem with isotropic stress and no viscosity, reactions or external fields. This is not a nonlinear GR, nonzero background heat-current, atmosphere, native-derivative-error or physical closure certificate. No identification with chi or a stellar relaxation time.'))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    prior.verify()
    paths=[g.ROOT/'verification/gr_heat_characteristics.py',prior.OUT/'manifest.json',prior.OUT/'thresholds.npz',
        g.OUT/'initial-state-17-4.npz',g.OUT/'gr-microphysics/auxiliaries.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e6aa3f1',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        Maxwell_finite_relative_tolerance=1e-8,trial_multipliers=[2,3],
        coefficient_contract='First compute b,r,d,e,h as binary64 from the bound arrays, retaining e from actual u_lnrho instead of silently imposing e=d. Thereafter all inequalities and trial lambda/tau values use exact rational arithmetic on those frozen coefficients.',
        physical_heat_time_calibrated=False,physical_EOS_certified=False,full_GR_evolution=False))
    symbolic()


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=prior.micro.inputs();a=aux['eos'];rho=np.exp(state['lnd'])
    threshold=dict(np.load(prior.OUT/'thresholds.npz'));w=threshold['enthalpy_density'];h=threshold['homogeneous_tau_min_seconds']
    b=rho*a[:,10]/w;r=a[:,1]*a[:,5]/w;d=a[:,1]*a[:,6]/w;e=(a[:,1]-rho*a[:,9])/w
    score=abs(d-e)/np.maximum(np.maximum(abs(d),abs(b)),np.finfo(float).tiny)
    assert np.all(np.isfinite(score)) and score.max()<plan['Maxwell_finite_relative_tolerance']
    fields=np.column_stack([b,r,d,e,h]);np.savez_compressed(OUT/'coefficients.npz',coefficients=fields,Maxwell_scores=score)
    failures=[];rows=[];minimum_H2=None;minimum_H3=None;maximum_g=F(0)
    for index,values in enumerate(fields):
        B,R,D,E,H=[F(float(v)) for v in values]
        assert B>0 and R>0 and H>0
        sound=R+D*E/B;H2=1-E-D+B*sound;H3=sound*H2-R
        maximum_g=max(maximum_g,sound)
        minimum_H2=H2 if minimum_H2 is None else min(minimum_H2,H2)
        minimum_H3=H3 if minimum_H3 is None else min(minimum_H3,H3)
        if not (0<sound<1 and H2>0 and H3>0):
            failures.append(dict(cell=index,reason='EOS or strict Hurwitz inequalities'));continue
        lamstar=max(F(1),(B+1-E-D-R)/(B*(1-sound)),(2*B+1-E-D)/(B*(2-sound)))
        row=[]
        for multiplier in plan['trial_multipliers']:
            lam=multiplier*lamstar;A=B*(lam-1);C=lam*B*sound+1-E-D
            discriminant=C*C-4*A*R
            valid=(lam>1 and A>0 and C>0 and discriminant>=0 and C<=2*A and A-C+R>=0)
            if not valid:failures.append(dict(cell=index,multiplier=multiplier,reason='characteristic inequalities'))
            # Only the displayed speed is floating point; its light-cone
            # membership is checked above using exact frozen fractions.
            Af,Cf,Rf=float(A),float(C),float(R)
            fast=(Cf+np.sqrt(float(discriminant)))/(2*Af);slow=Rf/(Af*fast)
            row.extend([float(H*lam),np.sqrt(slow),np.sqrt(fast)])
        rows.append(row)
    np.savez_compressed(OUT/'trial-times-speeds.npz',values=np.array(rows))
    trials=np.array(rows);result=dict(classification='Proven',completed=True,cells=len(fields),
        frozen_coefficient_exact_inequalities_passed=not failures,failures=failures,
        exact_scope='All inequalities use exact frozen binary64 coefficients. Both declared trial times pass the linear model tests at all saved cells if the pass flag is true. The maximum/minimum/range summaries below are floating displays, not outward interval bounds.',
        maximum_Maxwell_finite_score=float(score.max()),minimum_Hurwitz_second_factor=float(minimum_H2),
        minimum_Hurwitz_third_factor=float(minimum_H3),maximum_adiabatic_speed_squared_over_c2=float(maximum_g),
        tau_seconds_ranges=[list(map(float,[trials[:,j].min(),trials[:,j].max()])) for j in [0,3]] if len(trials) else [],
        maximum_characteristic_speed_over_c=[float(trials[:,j].max()) for j in [2,5]] if len(trials) else [],
        relaxation_nonuniqueness_demonstrated=not failures,physical_heat_time_calibrated=False,
        physical_EOS_certified=False,nonzero_heat_background_or_full_GR_evolution=False)
    save('result.json',result)
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('COUPLED HEAT CHARACTERISTICS',result,flush=True)
    assert not failures,result
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['frozen_coefficient_exact_inequalities_passed']
    print('PASS COUPLED HEAT CHARACTERISTICS SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
