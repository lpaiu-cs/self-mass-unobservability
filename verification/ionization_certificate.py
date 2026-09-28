"""Potential-consistent electron quadrature and ideal ionization boundaries.

Proven conditional constructions; no certification of the nonideal FreeEOS
equilibrium or replacement of the running GR model.
"""
from fractions import Fraction as F
import json, sys
import sympy as sp
from mpmath import iv
import electron_thermo_certificate as e
import verify_electron_thermo as audit
from interval_records import interval_text

g=e.g;OUT=g.OUT/'ionization-certificate'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/ionization_certificate.py',e.OUT/'result.json',e.OUT/'audit-manifest.json',
        g.ROOT/'verification/interval_records.py']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='46f1a9c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        eta_interval=[-17,24],beta_interval=['0','0.006'],precision_decimal_digits=60,
        criterion='Prove positive electron potential curvature and quantify dual value, gradient and relative curvature errors; prove ideal-ion constrained curvature and monotone charge balance. Do not infer nonideal curvature from convergence.',
        full_EOS_certified=False,full_GR_evolution=False,
        sources=['https://freeeos.sourceforge.net/convergence.pdf',
            'https://freeeos.sourceforge.net/documentation.html']))


def symbolic():
    x,b=sp.symbols('x b',positive=True);q=sp.symbols('q',positive=True)
    dweight=sp.sqrt(x)*sp.sqrt(1+b*x/2)*(1+b*x)
    first=dweight*q*(1-q);second=q*(1-q)*sp.diff(first,q)
    assert sp.simplify(second-first*(1-2*q))==0
    # Exact constrained three-stage Hessian, checked by independent differentiation.
    p1,p2,N,eta=sp.symbols('p1 p2 N eta',positive=True)
    p0=1-p1-p2;electron=sp.Function('electron');z=sp.symbols('z')
    ne=N*(p1+2*p2)
    energy=N*(p0*sp.log(p0)+p1*sp.log(p1)+p2*sp.log(p2))+electron(ne)
    curvature=sp.diff(electron(z),z,2).subs(z,ne)
    expected=N*sp.Matrix([[1/p1+1/p0,1/p0],[1/p0,1/p2+1/p0]])+curvature*N*N*sp.Matrix([[1,2],[2,4]])
    assert sp.simplify(sp.hessian(energy,[p1,p2])-expected)==sp.zeros(2)
    # Charge statistics: derivative of the mean is minus its variance.
    a0,a1,a2,t=sp.symbols('a0 a1 a2 t',positive=True)
    partition=a0+a1*t+a2*t*t;mean=(a1*t+2*a2*t*t)/partition
    variance=(a1*t+4*a2*t*t)/partition-mean*mean
    assert sp.simplify(-t*sp.diff(mean,t)+variance)==0
    assert sp.simplify(variance-(a0*a1*t+4*a0*a2*t*t+a1*a2*t**3)/partition**2)==0
    # Value closeness of a function need not preserve convexity.
    eps=sp.Rational(1,1000);omega=100
    perturbed=eta*eta/2+eps*sp.cos(omega*eta)
    assert sp.diff(perturbed,eta,2).subs(eta,0)==-9
    # Finite quadrature does not automatically preserve integration by parts.
    nodes,weights=e.u.rule();defect=2-9*sum(w*z**8 for z,w in zip(nodes,weights))
    assert defect.a>0
    return dict(classification='Proven',passed=True,
        electron_log_curvature_lipschitz='|K_eta_eta_eta|=|D_eta_eta|<=D_eta=K_eta_eta; hence |d ln K_eta_eta/d eta|<=1.',
        ideal_ion_hessian='For fixed positive N_a and T, Phi=F/(kT)=sum_(a,j) N_a p_aj[ln p_aj+c_aj]+Phi_e(sum_(a,j) j N_a p_aj). On sum_j dp_aj=0, d^2 Phi=sum N_a dp_aj^2/p_aj+Phi_e_second*(sum j N_a dp_aj)^2 >= sum N_a dp_aj^2.',
        ideal_ion_conditions='0<p_aj<=1, finite fixed partition/ionization constants, ideal nondegenerate ions, positive ideal-electron compressibility. Zero-abundance elements are removed from the coordinates. No molecular or nonideal coupling is included in the scalar reduction.',
        ideal_equilibrium='p_aj=a_aj exp(-j eta)/sum_k a_ak exp(-k eta), a_aj>0. The charge balance n_e(eta)-sum_a N_a <j>_a has derivative dn_e/deta+sum_a N_a Var_a(j)>0. With at least one positive-abundance element permitting a charged stage, limits at eta=-infinity and +infinity have opposite signs; the unique ideal root exists. All-neutral-only inventories require a separate zero-electron boundary.',
        coordinate_boundary='The weighted norm sum N_a dp_aj^2 gives a dimensionless strong-convexity constant at least one at fixed external density, temperature and abundances. Parameter derivatives require the explicit coordinate dependence or fixed reference weights; no parameter Hessian bound is inferred from this fixed-state assertion.',
        nonideal_margin='If, in the same fixed weighted coordinates and entire feasible domain, the nonideal projected Hessian has lower bound -kappa*I with kappa<1, total strong convexity is at least 1-kappa. Local solver convergence supplies neither this bound nor a global minimum certificate.',
        value_only_positive_control=dict(uniform_value_perturbation=str(eps),curvature_at_zero=-9),
        finite_quadrature_boundary_defect=interval_text(defect),
        physical_EOS_certified=False)


def run():
    plan=json.loads((OUT/'plan.json').read_text());assert not (OUT/'result.json').exists()
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    audit.verify();iv.dps=plan['precision_decimal_digits'];save('symbolic.json',symbolic())
    prior=json.loads((e.OUT/'result.json').read_text());m=F(prior['inverse_transfer']['density_slope_lower_rational'])
    errors={(r['quantity'],r['eta_order'],r['beta_order']):e.audit.interval(r['error_upper'])[1]
        for r in prior['uniform_error_components']}
    rows=[]
    for name,first,second in [('density_quadrature',errors[('D',0,0)],errors[('D',1,0)]),
            ('potential_quadrature',errors[('K',1,0)],errors[('K',2,0)])]:
        delta=first/m;relative=second/m;assert 0<relative<1
        lo=iv.exp(-e.interval(delta))/(1+e.interval(relative))
        hi=iv.exp(e.interval(delta))/(1-e.interval(relative))
        relative_error=max(1-lo.a,hi.b-1)
        rows.append(dict(map=name,positive_quadrature_curvature_lower_rational=str(m-second),
            eta_error_upper_rational=str(delta),relative_reciprocal_slope_error_upper=interval_text(relative_error.b),
            reciprocal_slope_ratio_interval=[interval_text(lo.a),interval_text(hi.b)],
            display_only=dict(eta_error=float(delta),relative_inverse_derivative_error=float(relative_error.b))))
    save('result.json',dict(classification='Proven',passed=True,inverse_derivative_transfers=rows,
        restricted_Legendre_value_error_upper_rational=str(errors[('K',0,0)]),
        finite_pressure_density_identity_defect_upper_rational=str(errors[('K',1,0)]+errors[('D',0,0)]),
        potential_construction='At fixed beta use the SAME quadrature scalar Q_K(eta,beta). Define its density as dQ_K/deta, not an independently discretized D. The certified second eta derivative is positive. Define the reduced Helmholtz potential as sup_(eta in [-17,24]) [eta*y-Q_K]. The true restricted dual uses K in the same supremum.',
        value_proof='Suprema over the same compact interval differ by at most sup|K-Q_K| for every density target. This value bound includes endpoint maximizers. The gradient and inverse-curvature comparisons require both true and quadrature stationary roots to lie in the interval interior.',
        inverse_derivative_proof='The two roots differ by delta<=E_first/m. Positive true slopes obey exp(-delta)<=true_slope(eta_true)/true_slope(eta_quad)<=exp(delta). The quadrature slope differs at eta_quad by relative at most epsilon=E_second/m. Therefore exp(-delta)/(1+epsilon)<=quad_inverse_slope/true_inverse_slope<=exp(delta)/(1-epsilon). This compares derivatives at DIFFERENT roots for the SAME density target.',
        assumptions='Exact declared quadrature and interior matched density targets on the certified eta/beta domain; numerical root residual and evaluation roundoff are additional unless explicitly enclosed. beta=0 refers to reduced functions only.',
        ideal_ionization_certified_conditionally=True,nonideal_curvature_measured=False,
        native_EOS_replaced=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('PASS POTENTIAL AND IDEAL IONIZATION',[(r['map'],r['display_only']) for r in rows],flush=True)


def freeze():
    assert json.loads((OUT/'result.json').read_text())['passed'] and not (OUT/'manifest.json').exists()
    paths=[p for p in OUT.iterdir() if p.is_file()]+[g.ROOT/'verification/ionization_certificate.py']
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS POTENTIAL AND IDEAL IONIZATION BINDINGS',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
