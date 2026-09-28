"""Conditional interval quadrature for the repaired nonrelativistic primitive."""
import ctypes, json, math, sys
import numpy as np
import sympy as sp
from mpmath import iv
import direct_ion_eos as d


OUT=d.OUT/'fermi-interval'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='7e0c32c',
        integral='F_nu^(k)(eta)=integral_0^infinity x^nu P_k(1/(1+exp(x-eta))) dx',
        cases=[dict(eta=eta,nu=nu,k=k) for eta in [1,4]
            for nu,orders in [('1.5',[0]),('0.5',range(4)),('-0.5',[3])]
            for k in orders],
        transform='x=t^2 removes the origin singularity: 2 t^(2 nu+1) P_k(q).',
        finite_interval_t=[0,8],panels=[64,128,256,512,1024],precision_decimal_digits=60,
        rule='Four-point Gauss-Legendre; exact algebraic nodes and weights enclosed with interval arithmetic. Bound the eighth t derivative on each panel and add the proved upper-tail bound at x=64.',
        maximum_enclosure_radius='1e-11',maximum_native_absolute_error='1e-9',
        scope='Registered nonrelativistic Fermi primitive values and eta derivatives at specified points; not all EOS states, beta derivatives, plasma model error or stellar propagation.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text());iv.dps=plan['precision_decimal_digits']
    t,q=sp.symbols('t q');polys=[q]
    for _ in range(5): polys.append(sp.expand(q*(1-q)*sp.diff(polys[-1],q)))
    # Validate the algebraic rule itself through its degree-seven moments.
    outer=iv.sqrt((3+2*iv.sqrt(iv.mpf(6)/5))/7)
    inner=iv.sqrt((3-2*iv.sqrt(iv.mpf(6)/5))/7)
    wo=(18-iv.sqrt(30))/36;wi=(18+iv.sqrt(30))/36
    nodes=[-outer,-inner,inner,outer];weights=[wo,wi,wi,wo]
    for k in range(8):
        result=sum(w*x**k for x,w in zip(nodes,weights));target=iv.mpf(2)/(k+1) if k%2==0 else iv.mpf(0)
        assert result.a<=target.a and result.b>=target.b,(k,result,target)
    constant=iv.mpf(math.factorial(4)**4)/(9*math.factorial(8)**3)
    monic=sp.Poly(sp.legendre(4,t),t).monic().as_expr()
    gamma=sp.integrate(monic**2,(t,-1,1))
    assert gamma/(math.factorial(8)*2**9)==sp.Rational(math.factorial(4)**4,9*math.factorial(8)**3)
    gamma_iv=iv.mpf(int(sp.numer(gamma)))/int(sp.denom(gamma))
    degree_eight_error=iv.mpf(2)/9-sum(w*x**8 for x,w in zip(nodes,weights))
    assert degree_eight_error.a>0 and degree_eight_error.a<=gamma_iv.a and degree_eight_error.b>=gamma_iv.b
    save('rule-control.json',dict(classification='Proven',passed=True,exactness_degree=7,
        degree_eight_positive_control=str(degree_eight_error),monic_norm=str(gamma),
        error_constant_identity=True,primary_reference='https://dlmf.nist.gov/3.5#v'))
    eos=d.EOS('full');call=eos.lib.__mod_fermi_dirac_MOD_fermi_dirac_ct
    call.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
    call.restype=ctypes.c_double;zero=ctypes.c_int(0);rows=[]
    for case in plan['cases']:
        eta=iv.mpf(case['eta']);nu=iv.mpf(case['nu']);k=case['k'];power=int(2*float(case['nu'])+1)
        coeff=[int(sp.Poly(polys[k],q).nth(j)) for j in range(k+2)]
        g=2*t**power*polys[k]
        for _ in range(8): g=sp.expand(sp.diff(g,t)-2*t*q*(1-q)*sp.diff(g,q))
        derivative=sp.lambdify((t,q),sp.horner(g,q),modules='math')
        def integrand(x):
            z=1/(1+iv.exp(x*x-eta));p=iv.mpf(0)
            for c in reversed(coeff): p=p*z+c
            return 2*x**power*p
        delta=iv.exp(eta-64)
        tail=sum(abs(c)*delta**(j-1) for j,c in enumerate(coeff) if j>0)*delta*iv.mpf(64)**nu/(1-max(float(case['nu']),0)/64)
        attempts=[];enclosure=None
        for n in plan['panels']:
            width=iv.mpf(8)/n;half=width/2;total=iv.mpf(0);error=iv.mpf(0)
            for j in range(n):
                a=width*j;b=width*(j+1);center=width*(iv.mpf(j)+iv.mpf('.5'))
                total+=half*sum(w*integrand(center+half*x) for x,w in zip(nodes,weights))
                interval=iv.mpf([a.a,b.b]);occupation=1/(1+iv.exp(interval**2-eta))
                bound=abs(derivative(interval,occupation)).b
                error+=constant*width**9*bound
            radius=error+tail
            enclosure=total+iv.mpf([-radius.b,radius.b])
            attempts.append(dict(panels=n,quadrature_error_bound=str(error),tail_bound=str(tail),enclosure=str(enclosure)))
            print('FERMI INTERVAL',case,n,'radius',str(radius),flush=True)
            if (enclosure.b-enclosure.a).b<=2*iv.mpf(plan['maximum_enclosure_radius']).a: break
        ct_order=0 if case['nu']=='1.5' else k+1 if case['nu']=='0.5' else 5
        divisor=1. if case['nu']=='1.5' else 1.5 if case['nu']=='0.5' else .75
        argument=ctypes.c_double(case['eta'])
        native=call(ctypes.byref(argument),ctypes.byref(ctypes.c_int(ct_order)),ctypes.byref(zero))/divisor
        native_bound=abs(iv.mpf(native)-enclosure).b
        passed=bool((enclosure.b-enclosure.a).b<=2*iv.mpf(plan['maximum_enclosure_radius']).a and
            native_bound<=iv.mpf(plan['maximum_native_absolute_error']).a)
        rows.append(dict(**case,attempts=attempts,native_value=native,native_absolute_error_bound=str(native_bound),passed=passed))
        save('result.json',dict(classification='Proven',rows=rows,complete=len(rows)==len(plan['cases']),
            passed=all(r['passed'] for r in rows),plan_sha256=d.e.c.sha(OUT/'plan.json'),
            library_sha256=json.loads((d.OUT/'full-integral-build.json').read_text())['sha256'],
            scope=plan['scope']))
        assert passed,case
    print('PASS FERMI INTERVAL',len(rows),'registered primitive cases',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
