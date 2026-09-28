"""Uniform conditional quadrature bounds for relativistic Fermi primitives.

Proven for a declared composite Gauss rule on a closed eta/beta box. The
existing adaptive native routine is compared only at the registered points.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import ctypes, json, math, shutil, subprocess, sys
import numpy as np
import sympy as sp
from mpmath import iv
import direct_eos_gr as g

OUT=g.OUT/'fermi-uniform';CACHE=g.CACHE/'fermi-uniform'
ORDERS=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2)]
CASES=[dict(nu=nu,eta_order=k,beta_order=l) for nu in ['0.5','1.5','2.5'] for k,l in ORDERS]
CASES += [dict(nu='-0.5',eta_order=3,beta_order=0)]
POLYS=[[0,1],[0,1,-1],[0,1,-3,2],[0,1,-7,12,-6]]


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir(exist_ok=True)
    paths=[g.ROOT/'verification/fermi_uniform.py',g.d.OUT/'full-integral-build.json',
        g.d.OUT/'fermi-interval/result.json',g.OUT/'opacity/electrons/result.json']
    for name in ['fermi_dirac_direct.f90','fermi_dirac.f90','fermi_dirac_ct.f90']:
        target=OUT/name;shutil.copy2(g.d.CACHE/'full-integral-source/src'/name,target);paths.append(target)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='cc6e137',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        eta_interval=[-17,24],beta_interval=['0','0.006'],cases=CASES,
        transformed_interval=[0,10],coarse_intervals=256,
        fine_panels=[256,512,1024,2048,4096,8192,16384,32768,65536],
        precision_decimal_digits=60,maximum_uniform_quadrature_error='1e-11',
        point_controls=[[-17,'0'],[0,'0'],[24,'0.006']],processes=3,
        native_point_error_relative_with_floor_one='1e-9',
        scope='Uniform value and specified mixed eta/beta derivative errors of a newly declared exact composite four-point Gauss rule, plus interval-arithmetic evaluation at three controls. No uniform error claim for the old adaptive binary, implicit EOS root, physical plasma or GR evolution.',
        physical_EOS_certified=False))


def rule():
    outer=iv.sqrt((3+2*iv.sqrt(iv.mpf(6)/5))/7)
    inner=iv.sqrt((3-2*iv.sqrt(iv.mpf(6)/5))/7)
    return [-outer,-inner,inner,outer],[(18-iv.sqrt(30))/36,(18+iv.sqrt(30))/36,
        (18+iv.sqrt(30))/36,(18-iv.sqrt(30))/36]


def tail(case,eta_max,beta_max,S):
    nu=iv.mpf(case['nu']);k=case['eta_order'];l=case['beta_order']
    delta=iv.exp(eta_max-S)
    coefficient=sum(abs(c)*delta**(j-1) for j,c in enumerate(POLYS[k]) if j>0)
    factor=[iv.mpf(1),iv.mpf(1)/4,iv.mpf(1)/16][l]
    power=nu+l
    if l==0:
        factor*=iv.sqrt(1+beta_max*S/2)/iv.sqrt(S);power+=iv.mpf('.5')
    return factor*coefficient*delta*S**power/(1-max(power,iv.mpf(0))/S)


def certify():
    plan=json.loads((OUT/'plan.json').read_text());iv.dps=plan['precision_decimal_digits']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    nodes,weights=rule();x=sp.symbols('x')
    monic=sp.Poly(sp.legendre(4,x),x).monic().as_expr()
    norm=sp.integrate(monic**2,(x,-1,1))
    exact=sp.Rational(math.factorial(4)**4,9*math.factorial(8)**3)
    assert norm/(math.factorial(8)*2**9)==exact
    for k in range(8):
        moment=sum(w*z**k for z,w in zip(nodes,weights));value=iv.mpf(2)/(k+1) if k%2==0 else iv.mpf(0)
        assert moment.a<=value.a and moment.b>=value.b
    positive=iv.mpf(2)/9-sum(w*z**8 for z,w in zip(nodes,weights))
    assert positive.a>0
    # Logistic eta-derivative and beta-derivative prefactors are exact.
    q,b=sp.symbols('q b');poly=q
    for k in range(4):
        assert sp.expand(poly-sum(c*q**j for j,c in enumerate(POLYS[k])))==0
        poly=sp.expand(q*(1-q)*sp.diff(poly,q))
    for l,factor in enumerate([sp.Rational(1),sp.Rational(1,4),-sp.Rational(1,16)]):
        assert sp.simplify(sp.diff(sp.sqrt(1+b*x/2),b,l)-factor*x**l*(1+b*x/2)**(sp.Rational(1,2)-l))==0
    constant=iv.mpf(int(exact.p))/int(exact.q)*math.factorial(8)
    U=iv.mpf(plan['transformed_interval'][1]);beta=iv.mpf(plan['beta_interval'][1]);eta=iv.mpf(plan['eta_interval'][1])
    count=plan['coarse_intervals'];width=U/count;totals=[iv.mpf(0) for _ in CASES]
    for j in range(count):
        a=width*j;b=width*(j+1);R=1/(2*(b+1));modulus=b+R
        imaginary=2*(b+R)*R;assert imaginary.b<iv.mpf('1.01').a
        lower_w=1-beta*R**2/2;upper_w=1+beta*modulus**2/2;assert lower_w.a>0
        real_square=max(iv.mpf(0),a-R)**2-R**2
        qmax=min(iv.mpf(1),iv.exp(eta-real_square))
        for i,case in enumerate(CASES):
            nu=Fraction(case['nu']);k=case['eta_order'];l=case['beta_order']
            power=int(2*nu+1+2*l);exponent=iv.mpf('.5')-l
            wbound=upper_w**exponent if l==0 else lower_w**exponent
            polynomial=sum(abs(c)*qmax**n for n,c in enumerate(POLYS[k]))
            M=[iv.mpf(2),iv.mpf('.5'),iv.mpf(1)/8][l]*modulus**power*wbound*polynomial
            totals[i]+=width*M/R**8
    rows=[]
    for case,total in zip(CASES,totals):
        upper_tail=tail(case,eta,beta,U**2);attempts=[]
        for n in plan['fine_panels']:
            error=constant*(U/n)**8*total+upper_tail
            attempts.append(dict(panels=n,error_upper=str(error.b)))
            if error.b<=iv.mpf(plan['maximum_uniform_quadrature_error']).a: break
        passed=bool(error.b<=iv.mpf(plan['maximum_uniform_quadrature_error']).a)
        rows.append(dict(**case,passed=passed,panels=n,uniform_error_upper=str(error.b),
            tail_upper=str(upper_tail.b),attempts=attempts))
    census=json.loads((g.OUT/'opacity/electrons/result.json').read_text())
    covered=bool(census['sampled_eta_range'][0]>=-17 and census['sampled_eta_range'][1]<=24
        and census['sampled_beta_range'][0]>=0 and census['sampled_beta_range'][1]<=.006)
    save('uniform-certificate.json',dict(classification='Proven',passed=all(r['passed'] for r in rows),rows=rows,
        complex_disk_proof='For t in [a,b], take R=1/(2(b+1)). Every disk has |Im z^2|<=2(b+R)R<=1<pi/2, so Re exp(z^2-eta)>0, and |q|<=min(1,exp(eta_max-max(0,a-R)^2+R^2)). Re(1+beta*z^2/2)>=1-beta_max*R^2/2>0; the principal fractional power is analytic. The displayed M bounds each transformed mixed derivative integrand on that disk. Cauchy gives |g^(8)|<=8!*M/R^8. Sum the four-point Gauss remainder over fine panels contained in the coarse cells.',
        tail_proof='For x>=S, q<=exp(eta_max-x). Bound its polynomial derivatives termwise. For beta order zero, sqrt(1+beta*x/2)<=sqrt(1+beta_max*S/2)*sqrt(x/S); for beta orders one and two the negative powers are <=1. Use integral_S^infinity x^p exp(-x) dx <= S^p exp(-S)/(1-max(p,0)/S).',
        original_saved_sample_domain_covered=covered,finite_sample_census_is_not_a_trajectory_enclosure=True,
        primary_references=['https://dlmf.nist.gov/1.9#iii','https://dlmf.nist.gov/3.5#v'],
        plan_sha256=g.c.sha(OUT/'plan.json'),old_native_uniform_error_certified=False,full_EOS_certified=False))
    print('UNIFORM FERMI BOUNDS',len(rows),all(r['passed'] for r in rows),'largest panels',max(r['panels'] for r in rows),'sample domain',covered,flush=True)
    assert all(r['passed'] for r in rows)


def build_bridge():
    source=OUT/'primitive_bridge.f90'
    source.write_text('''subroutine primitive(nu,eta,beta,answer) bind(C)
  use iso_c_binding
  use mod_fermi_dirac, only: fermi_dirac_direct
  implicit none
  real(c_double), value :: nu,eta,beta
  real(c_double), intent(out) :: answer(9)
  call fermi_dirac_direct(nu,eta,beta,answer)
end subroutine primitive
''')
    root=g.d.CACHE/'full-integral-build';module=next(root.rglob('mod_fermi_dirac.mod')).parent
    target=CACHE/'primitive.so';library=root/'src/libfree_eos_direct24_integral_full.so.1.0.0'
    expected=json.loads((g.d.OUT/'full-integral-build.json').read_text())['sha256'][str(library)]
    assert g.c.sha(library)==expected
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(module),str(source),
        '-L'+str(root/'src'),'-Wl,-rpath,'+str(root/'src'),'-lfree_eos_direct24_integral_full','-o',str(target)]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'bridge-build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr;assert g.c.sha(library)==expected
    save('bridge.json',dict(classification='Counterexample candidate',passed=True,command=command,
        source_sha256=g.c.sha(source),library_sha256=expected,bridge_sha256=g.c.sha(target),original_library_unchanged=True))


def point(index):
    plan=json.loads((OUT/'plan.json').read_text());cert=json.loads((OUT/'uniform-certificate.json').read_text())
    assert cert['passed'];iv.dps=plan['precision_decimal_digits'];eta_arg,beta_arg=plan['point_controls'][index]
    eta=iv.mpf(eta_arg);beta=iv.mpf(beta_arg);n=max(r['panels'] for r in cert['rows']);h=iv.mpf(10)/n
    nodes,weights=rule();total=[iv.mpf(0) for _ in CASES]
    for j in range(n):
        center=h*(iv.mpf(j)+iv.mpf('.5'))
        for z,weight in zip(nodes,weights):
            t=center+h*z/2;q=1/(1+iv.exp(t*t-eta));w=1+beta*t*t/2
            polynomials=[]
            for coeff in POLYS:
                value=iv.mpf(0)
                for c in reversed(coeff): value=value*q+c
                polynomials.append(value)
            for i,case in enumerate(CASES):
                l=case['beta_order'];power=int(2*Fraction(case['nu'])+1+2*l)
                value=[iv.mpf(2),iv.mpf('.5'),-iv.mpf(1)/8][l]*t**power*w**(iv.mpf('.5')-l)*polynomials[case['eta_order']]
                total[i]+=h/2*weight*value
        if j%2048==0: print('CERTIFIED FERMI POINT',index,j,'/',n,flush=True)
    lib=ctypes.CDLL(str(CACHE/'primitive.so'));call=lib.primitive
    call.argtypes=[ctypes.c_double]*3+[np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')];call.restype=None
    native={}
    for nu in ['-0.5','0.5','1.5','2.5']:
        a=np.full(9,np.nan);call(float(nu),float(eta_arg),float(beta_arg),a);assert np.all(np.isfinite(a));native[nu]=a
    rows=[]
    for value,case,bound in zip(total,CASES,cert['rows']):
        # The finer common panel count can retain each earlier, larger bound.
        radius=iv.mpf(bound['uniform_error_upper']).b
        interval=value+iv.mpf([-radius,radius]);actual=native[case['nu']][ORDERS.index((case['eta_order'],case['beta_order']))]
        error=abs(iv.mpf(float(actual))-interval).b;scale=max(iv.mpf(1),abs(iv.mpf(float(actual))))
        passed=bool(error<=iv.mpf(plan['native_point_error_relative_with_floor_one'])*scale)
        rows.append(dict(**case,enclosure=str(interval),native_value=float(actual),native_absolute_error_upper=str(error),passed=passed))
    result=dict(classification='Proven',point=[eta_arg,beta_arg],panels=n,rows=rows,
        passed=all(r['passed'] for r in rows),old_native_uniform_error_certified=False)
    save(f'point-{index}.json',result);assert result['passed'],index;return index


def points():
    with ProcessPoolExecutor(max_workers=3) as pool:
        for done in as_completed([pool.submit(point,i) for i in range(3)]): print('FERMI POINT COMPLETE',done.result(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
