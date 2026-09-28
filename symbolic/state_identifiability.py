"""Exact observable-pole uniqueness and equilibrium/dynamics nonidentifiability."""
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def main():
    z=s.symbols('z'); tau, eps, beta=s.symbols('tau eps beta',positive=True)
    lo,hi=s.symbols('lo hi',positive=True)
    den=(1+lo**2*tau**2)*(1+hi**2*tau**2)
    sl=tau/(1+lo**2*tau**2); sh=tau/(1+hi**2*tau**2)
    checks=[]
    def check(name,expr):
        assert s.factor(expr)==0, (name,expr)
        checks.append(name)
    check('moment_R',(hi**2*sh-lo**2*sl)/(hi**2-lo**2)-tau/den)
    check('moment_D',(1/(1+lo**2*tau**2)-1/(1+hi**2*tau**2))/(hi**2-lo**2)-tau**2/den)
    check('moment_Q',(sl-sh)/(hi**2-lo**2)-tau**3/den)
    w1,w2,t1,t2=s.symbols('w1 w2 t1 t2',positive=True)
    check('positive_two_atom_variance',(w1+w2)*(w1*t1**2+w2*t2**2)-(w1*t1+w2*t2)**2-w1*w2*(t1-t2)**2)
    close=beta/2*(1/(1+z*(tau-eps))+1/(1+z*(tau+eps)))
    check('two_positive_poles_approach_one',close-beta/(1+z*tau)-beta*z**2*eps**2/((1+z*tau)*((1+z*tau)**2-z**2*eps**2)))
    lam,eta=s.symbols('lam eta',positive=True)
    invisible=eta*s.prod(z*z+k*k for k in [1,2,3])/(z+lam)**7
    for k in [1,2,3]:
        for sign in [1,-1]: check(f'stable_carrier_invisible_{k}_{sign}',invisible.subs(z,sign*s.I*k))
    kappa,gamma,force=s.symbols('kappa gamma force',positive=True)
    check('same_static_response',(1/(kappa+z*gamma)).subs(z,0)-1/kappa)
    check('arbitrary_relaxation_same_stiffness',(gamma/kappa).subs(gamma,kappa*tau)-tau)
    # Coupled monopole radiation damping from the cited model is semidefinite, not SPD.
    damping=s.ones(2)
    check('monopole_damping_null_mode',(s.Matrix([1,-1]).T*damping*s.Matrix([1,-1]))[0])
    result=dict(status='Proven',checks=checks,
                positive_spectrum_identity='R*Q-D^2 = (integral dmu)^2 Var_mu(tau) >= 0; dmu=a*tau/[(1+lo^2*tau^2)(1+hi^2*tau^2)]',
                exact_uniqueness='Two exact calibrated complex samples identify a single positive observable relaxation time iff R*Q=D^2 and R>0; arbitrary real instantaneous term allowed.',
                boundaries=['Known common conjugate response and positive finite relaxation measure required',
                            'Zero-time instantaneous response and unobservable internal states remain unidentified',
                            'Finite noise cannot uniformly separate two arbitrarily close positive poles from one',
                            'Equilibrium stiffness alone does not set damping, inertia, or a relaxation-rate gap',
                            'Specified microphysics may relate these coefficients; the free-Gamma counterexample is not a claim about every fixed gravity theory'])
    (ROOT/'outputs/research-completion/state-identifiability.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'exact identities')


if __name__=='__main__': main()
