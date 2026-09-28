"""Conditional ideal-electron thermodynamics from certified Fermi primitives.

This does not certify ionization, nonideal plasma physics or the native EOS.
"""
from fractions import Fraction as F
import json, math, sys
import sympy as sp
from mpmath import iv
import fermi_uniform as u
import verify_fermi_uniform as audit
from interval_records import interval_text, exact_endpoint

g=u.g;OUT=g.OUT/'electron-thermo-certificate'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def upper(value): return exact_endpoint(value._mpi_[1])
def lower(value): return exact_endpoint(value._mpi_[0])
def interval(value):
    if isinstance(value,F): return iv.mpf(value.numerator)/value.denominator
    return iv.mpf(value)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/electron_thermo_certificate.py',
        audit.OUT/'manifest.json',audit.OUT/'uniform-certificate.json',
        g.ROOT/'verification/interval_records.py',
        audit.OLD/'fermi_dirac.f90']+[audit.OUT/f'point-{i}.json' for i in range(3)]
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='b729d93',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        precision_decimal_digits=60,common_panels=8192,
        eta_interval=[-17,24],beta_interval=['0','0.006'],
        inverse_control_half_width=1,
        inverse_targets='Exact rational linear combinations of the previously stored native Fermi binary64 outputs, using the declared decimal beta. These are diagnostic density targets, not new observed densities.',
        criterion='Prove positive density slope, quadrature density slope, exact symbolic thermodynamic identities, and inverse existence/uniqueness inside each preregistered local bracket. Report the certified widths without a post-hoc accuracy threshold.',
        electron_model='Ideal relativistic electrons with two spin states and kinetic chemical potential eta; no positron or interaction contribution. Physical scaling requires beta>0. At beta=0 only reduced functions and their right limits are used.',
        full_EOS_certified=False,full_GR_evolution=False))


def slope_bound(a,b):
    # Restrict the positive derivative integral to a unit x interval.
    x=max(F(1),F(a));z=max(abs(x-F(b)),abs(x+1-F(a)))
    e=iv.exp(-interval(z))
    return iv.sqrt(interval(x))*e/(1+e)**2


def combine(fields,beta):
    def field(nu,k,l): return fields[(nu,k,l)]
    result={}
    for name,nu1,nu2,a,b in [('D','0.5','1.5',F(1),F(1)),
            ('K','1.5','2.5',F(2,3),F(1,3)),('H','1.5','2.5',F(1),F(1))]:
        for k,l in u.ORDERS:
            value=interval(a)*field(nu1,k,l)+interval(b)*beta*field(nu2,k,l)
            if l: value+=interval(b)*l*field(nu2,k,l-1)
            result[(name,k,l)]=value
    return result


def thermodynamics(values,eta,beta):
    D,De,Db=[values[('D',*order)] for order in [(0,0),(1,0),(0,1)]]
    K,Kb=[values[('K',*order)] for order in [(0,0),(0,1)]]
    H,He,Hb=[values[('H',*order)] for order in [(0,0),(1,0),(0,1)]]
    assert D.a>0 and De.a>0 and K.a>0
    A=(iv.mpf('1.5')*D+beta*Db)/De
    return dict(reduced_density=D,reduced_pressure=K,reduced_energy=H,
        energy_per_particle_over_kT=H/D,entropy_per_particle_over_k=(H+K)/D-eta,
        deta_dlnnumber_at_T=D/De,deta_dlnT_at_number=-A,
        electron_pressure_chi_number=D*D/(K*De),
        electron_pressure_chi_temperature=iv.mpf('2.5')+beta*Kb/K-D*A/K,
        electron_cv_per_particle_over_k=(iv.mpf('2.5')*H+beta*Hb-He*A)/D)


def symbolic():
    eta,beta,x=sp.symbols('eta beta x',real=True,positive=True)
    f=[sp.Function('f'+str(i))(eta,beta) for i in range(3)]
    for a,b,left,right in [(1,1,f[0],f[1]),(sp.Rational(2,3),sp.Rational(1,3),f[1],f[2]),(1,1,f[1],f[2])]:
        for k,l in u.ORDERS:
            formula=a*sp.diff(left,eta,k,beta,l)+b*beta*sp.diff(right,eta,k,beta,l)
            if l: formula+=b*l*sp.diff(right,eta,k,beta,l-1)
            assert sp.simplify(sp.diff(a*left+b*beta*right,eta,k,beta,l)-formula)==0
    w=1+beta*x/2
    d=x**sp.Rational(1,2)*sp.sqrt(w)*(1+beta*x)
    k=sp.Rational(2,3)*x**sp.Rational(3,2)*w**sp.Rational(3,2)
    h=x**sp.Rational(3,2)*sp.sqrt(w)*(1+beta*x)
    assert sp.simplify(sp.diff(k,x)-d)==0
    assert sp.simplify(sp.Rational(3,2)*k+beta*sp.diff(k,beta)-h)==0
    # Independent total differentiation at fixed physical number density.
    D,K,H=[sp.Function(n)(eta,beta) for n in ['D','K','H']]
    A=(sp.Rational(3,2)*D+beta*sp.diff(D,beta))/sp.diff(D,eta)
    def fixed_n(expr): return sp.diff(expr,beta)-A/beta*sp.diff(expr,eta)
    cv=fixed_n(beta*H/D)
    expected=(sp.Rational(5,2)*H+beta*sp.diff(H,beta)-sp.diff(H,eta)*A)/D
    assert sp.simplify(cv-expected)==0
    pressure=beta**sp.Rational(5,2)*K
    chi=sp.simplify(beta*fixed_n(pressure)/pressure)
    assert sp.simplify(chi-(sp.Rational(5,2)+beta*sp.diff(K,beta)/K-sp.diff(K,eta)*A/K))==0
    return dict(classification='Proven',passed=True,product_rule_checks=27,
        identities=['K_eta=D by integration by parts; boundary terms vanish',
            'H=(3/2)K+beta*K_beta pointwise in the integrand',
            'Pressure and energy temperature derivatives along beta^(3/2)*D=constant'],
        first_law='n*=sqrt(2)*beta^(3/2)*D, p*=sqrt(2)*beta^(5/2)*K, u*=sqrt(2)*beta^(5/2)*H; entropy per particle/k=(H+K)/D-eta. Physical constants multiply these normalized functions and have not been assigned experimental error bounds.')


def run():
    plan=json.loads((OUT/'plan.json').read_text());assert not (OUT/'result.json').exists()
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    audit.verify();iv.dps=plan['precision_decimal_digits'];save('symbolic.json',symbolic())
    cert=json.loads((audit.OUT/'uniform-certificate.json').read_text())
    errors={}
    for row in cert['rows']:
        ratio=F(row['panels'],plan['common_panels'])**8;assert ratio<=1
        # E_N <= ratio*E_old + tail is conservative even with outward records.
        bound=ratio*audit.interval(row['uniform_error_upper'])[1]+audit.interval(row['tail_upper'])[1]
        errors[(row['nu'],row['eta_order'],row['beta_order'])]=interval(bound)
    combined=combine(errors,iv.mpf(plan['beta_interval'][1]))
    error_rows=[dict(quantity=n,eta_order=k,beta_order=l,error_upper=interval_text(v.b)) for (n,k,l),v in combined.items()]
    slabs=[dict(eta=[a,a+1],density_slope_lower=interval_text(slope_bound(a,a+1).a)) for a in range(-17,24)]
    m=min(audit.interval(r['density_slope_lower'])[0] for r in slabs)
    e0=upper(combined[('D',0,0)]);e1=upper(combined[('D',1,0)])
    assert 0<e1<m
    # This uniform exact-quadrature inverse bound requires a bracketed target
    # and excludes floating evaluation error until separately enclosed.
    inverse=dict(density_slope_lower_rational=str(m),quadrature_density_slope_lower_rational=str(m-e1),
        eta_error_from_zero_quadrature_residual_upper_rational=str(e0/m),
        reciprocal_slope_error_at_same_eta_upper_rational=str(e1/(m*(m-e1))),
        general_residual_bound='|eta_true-eta_hat| <= (certified |Q_D(eta_hat,beta)-y| + E_D)/m, conditional on root bracketing. Evaluation roundoff must be included in the residual.',
        inverse_slope_scope='The reciprocal-slope error compares 1/D_eta and 1/Q_D_eta at the SAME eta, not derivatives evaluated at two different roots.')
    controls=[]
    for i in range(3):
        point=json.loads((audit.OUT/f'point-{i}.json').read_text());eta=iv.mpf(point['point'][0]);beta=iv.mpf(point['point'][1])
        fields={(r['nu'],r['eta_order'],r['beta_order']):iv.mpf(r['enclosure']) for r in point['rows']}
        values=combine(fields,beta);observed=thermodynamics(values,eta,beta)
        assert (values[('K',1,0)]-values[('D',0,0)]).a<=0<=(values[('K',1,0)]-values[('D',0,0)]).b
        assert (values[('H',0,0)]-iv.mpf('1.5')*values[('K',0,0)]-beta*values[('K',0,1)]).a<=0<=(values[('H',0,0)]-iv.mpf('1.5')*values[('K',0,0)]-beta*values[('K',0,1)]).b
        assert observed['electron_cv_per_particle_over_k'].a>0
        native={(r['nu'],r['eta_order'],r['beta_order']):F.from_float(r['native_value']) for r in point['rows']}
        target=native[('0.5',0,0)]+F(str(point['point'][1]))*native[('1.5',0,0)]
        D=values[('D',0,0)];residual=max(abs(lower(D)-target),abs(upper(D)-target))
        center=F(point['point'][0]);width=F(plan['inverse_control_half_width'])
        local_m=lower(slope_bound(center-width,center+width));delta=residual/local_m
        assert 0<=delta<width and target>0
        # At center +/- width, monotonicity moves D by at least m*width.
        # Thus a unique exact integral root exists and lies within +/- delta.
        controls.append(dict(point=point['point'],thermodynamics={k:interval_text(v) for k,v in observed.items()},
            target_reduced_density_rational=str(target),residual_upper_rational=str(residual),
            local_density_slope_lower_rational=str(local_m),eta_radius_upper_rational=str(delta),
            eta_enclosure_rational=[str(center-delta),str(center+delta)],
            eta_radius_display=float(delta),inverse_exists_and_is_unique=True,
            beta_zero_is_reduced_limit=point['point'][1]=='0'))
    save('result.json',dict(classification='Proven',passed=True,uniform_error_components=error_rows,
        monotonicity_slabs=slabs,inverse_transfer=inverse,point_controls=controls,
        display_only=dict(uniform_density_error=float(e0),uniform_density_derivative_error=float(e1),
            minimum_density_slope=float(m),uniform_quadrature_inverse_error=float(e0/m)),
        source_formula_file='outputs/direct-eos-gr33/fermi-uniform/fermi_dirac.f90',
        source_formula_sections='Scaling comments lines 28-55 and direct morder=21 branch; density, pressure, energy and beta product rules are reconstructed independently.',
        scope='Ideal-electron mathematical model, continuous exact-quadrature errors and finite interval thermodynamic/inverse controls. The local inverse derivative lower bounds apply on eta0 +/- 1, including outside the quadrature box, by a separate positive-integrand argument; no quadrature there is asserted.',
        native_full_EOS_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('PASS IDEAL ELECTRON TRANSFER',len(error_rows),'uniform components; inverse radii',
        [r['eta_radius_display'] for r in controls],'; quadrature inverse bound',float(e0/m),flush=True)


def freeze():
    assert json.loads((OUT/'result.json').read_text())['passed'] and not (OUT/'manifest.json').exists()
    paths=[p for p in OUT.iterdir() if p.is_file()]+[g.ROOT/'verification/electron_thermo_certificate.py']
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS IDEAL ELECTRON CERTIFICATE BINDINGS',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
