"""Uniform value/mixed derivative bounds for a declared plasma-moment Gauss rule."""
from fractions import Fraction as F
import json, math, sys
import numpy as np
import sympy as sp
from mpmath import iv
import fermi_uniform as fermi
import gr_plasma_dispersion as dispersion
from interval_records import interval_text, exact_endpoint

g=dispersion.g;OUT=g.OUT/'gr-plasma-moment-certificate'
ORDERS=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def I(x): return iv.mpf(str(x))
def upper(x): return exact_endpoint(x._mpi_[1])
def lower(x): return exact_endpoint(x._mpi_[0])
def rational(x): return iv.mpf(x.numerator)/x.denominator
def box(s): return iv.mpf(s)


def prepare():
    assert not OUT.exists();OUT.mkdir();dispersion.verify()
    paths=[g.ROOT/'verification'/name for name in ['gr_plasma_moment_certificate.py','fermi_uniform.py','interval_records.py']]
    paths += [dispersion.OUT/'manifest.json',dispersion.plasma.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d5246ee',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        eta_interval=[-17,24],beta_interval=['0','0.006'],velocity_powers=list(range(26)),orders=ORDERS,
        normalization='Bounds apply to exp(-min(eta,0)) times the value or partial derivative of the UNNORMALIZED moment. The normalizing factor is held constant while differentiating; no derivative through min(eta,0). Common physical sqrt(2)*beta^(3/2) is factored out.',
        definition='J_j=integral_0^infinity sqrt(t)*sqrt(1+beta*t/2)*[1-(1+beta*t)^(-2)]^j*(f_e+f_p)dt, with f_e=1/(1+exp(t-eta)), f_p=1/(1+exp(t+eta+2/beta)). At beta=0 use the electron right limit and zero positron limit.',
        exact_rule='Composite four-point Gauss on s=sqrt(t) in [0,10], evaluating only the electron term. Include electron tail and all omitted positron derivatives in the error. This is a newly declared exact-real rule, not a certificate for the prior adaptive binary.',
        coarse_intervals=256,fine_panels=[256,512,1024,2048,4096,8192,16384,32768,65536,131072],
        eta_complex_radius='0.25',beta_complex_radius='0.001',precision_decimal_digits=60,
        uniform_moment_absolute_error='1e-14',uniform_kernel_absolute_error='1e-10',
        corner=[24,'0.006'],corner_components=[25,'longitudinal_tail_24'],
        controls='Verify symbolic parameter and disk derivative bounds. Independently validate exported decimal arithmetic using exact Fractions. Interval-evaluate the dominating high-eta/high-beta corner for both series-tail moments.',
        scope='Conditional continuous-domain moment value and derivative quadrature errors. Interval corner includes arithmetic/serialization. Does not certify all native binary values, all state derivatives of the dispersion root, physical higher-order plasma effects or actual GR evolution.'))


def gamma_half(n):
    ans=iv.sqrt(iv.pi)
    for j in range(n):ans*=I(2*j+1)/2
    return ans


def gamma_tail_integer(n,S):
    # e^S * integral_S^infinity t^n e^-t dt, exactly for integer n.
    return sum(I(math.factorial(n))/math.factorial(j)*S**j for j in range(n+1))


def omitted(j,l,beta,S):
    c=iv.sqrt(1+beta*S/2)/iv.sqrt(S)
    a,b=[(I(0),I(1)),(I('.25'),I(2*j)),(I(1)/16+j,I(4*j*j+2*j))][l]
    electron=iv.exp(24-S)*(a/iv.sqrt(S)+b*c)*gamma_tail_integer(l+1,S)
    h=iv.sqrt(beta/2)
    A0=gamma_half(1)+h
    A1=gamma_half(2)/4+2*j*(gamma_half(2)+2*h)
    A2=(I(1)/16+j)*gamma_half(3)+(4*j*j+2*j)*(gamma_half(3)+6*h)
    pair=iv.exp(-2/beta+34)*[A0,A1+2/beta**2*A0,
        A2+4/beta**2*A1+(4/beta**3+4/beta**4)*A0][l]
    return electron,pair


def tail_extra(beta,S):
    c=iv.sqrt(1+beta*S/2)/iv.sqrt(S)
    electron=iv.exp(24-S)*c*(gamma_tail_integer(1,S)+2/I(53)*(
        2*beta*gamma_tail_integer(2,S)+beta**2*gamma_tail_integer(3,S)))
    h=iv.sqrt(beta/2)
    pair=iv.exp(-2/beta+34)*(gamma_half(1)+h+2/I(53)*(
        2*beta*(gamma_half(2)+2*h)+beta**2*(gamma_half(3)+6*h)))
    return electron,pair


def certify():
    plan=json.loads((OUT/'plan.json').read_text());iv.dps=plan['precision_decimal_digits']
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    b,t,u,L,M,D=sp.symbols('b t u L M D',positive=True)
    w=1-(1+b*t)**-2;s=sp.sqrt(1+b*t/2)
    assert sp.simplify(sp.diff(w,b)-2*t/(1+b*t)**3)==0
    assert sp.simplify(sp.diff(w,b,2)+6*t*t/(1+b*t)**4)==0
    assert sp.simplify(sp.diff(s,b)-t/(4*s))==0
    assert sp.simplify(sp.diff(s,b,2)+t*t/(16*s**3))==0
    ratio=u*M*(2+u*M)/(D+u*L)**2
    assert sp.simplify(sp.diff(ratio,u)-2*M*(D+u*(D*M-L))/(D+u*L)**3)==0
    nodes,weights=fermi.rule()
    for k in range(8):
        val=sum(weight*node**k for node,weight in zip(nodes,weights));ref=I(2)/(k+1) if k%2==0 else I(0)
        assert lower(val)<=lower(ref) and upper(val)>=upper(ref)
    beta=I('.006');rb=I(plan['beta_complex_radius']);re=I(plan['eta_complex_radius'])
    count=plan['coarse_intervals'];width=I(10)/count;total=[I(0) for _ in range(27)]
    for j in range(count):
        a=width*j;b=width*(j+1);R=1/(2*(b+1));m=b+R;ms=m*m
        real_square=max(I(0),a-R)**2-R*R
        assert upper(2*m*R+re)<lower(iv.pi/2)
        gamma_lower=1+beta*real_square-rb*ms
        assert lower(1-rb*ms+min(I(0),beta*real_square))>0
        assert lower(1-rb*ms/2+min(I(0),beta*real_square/2))>0
        d=1-rb*(ms+real_square);slope=d+(beta+rb)*min(I(0),d*ms-real_square)
        assert lower(slope)>0
        W=(beta+rb)*ms*(2+(beta+rb)*ms)/gamma_lower**2
        qmax=min(max(I(1),iv.exp(re-real_square)),iv.exp(24+re-real_square))
        base=2*ms*iv.sqrt(1+(beta+rb)*ms/2)*qmax
        for n in range(26):total[n]+=width*base*W**n/R**8
        total[26]+=width*base*W**25*(1+2*(beta+rb)*ms*(2+(beta+rb)*ms)/53)/R**8
    constant=I(math.factorial(4)**4)/(9*math.factorial(8)**2)
    rows=[]
    for component in list(range(26))+['longitudinal_tail_24']:
        for k,l in (ORDERS if isinstance(component,int) else [(0,0)]):
            value=total[component if isinstance(component,int) else 26]*math.factorial(k)*math.factorial(l)/re**k/rb**l
            tail,pair=omitted(component,l,beta,I(100)) if isinstance(component,int) else tail_extra(beta,I(100))
            attempts=[]
            for panels in plan['fine_panels']:
                error=constant*(I(10)/panels)**8*value+tail+pair
                attempts.append(dict(panels=panels,error=interval_text(error)))
                if upper(error)<=F(plan['uniform_moment_absolute_error']):break
            rows.append(dict(component=component,eta_order=k,beta_order=l,panels=panels,error=interval_text(error),
                electron_tail=interval_text(tail),omitted_positron=interval_text(pair),attempts=attempts,
                passed=upper(error)<=F(plan['uniform_moment_absolute_error'])))
    common=max(row['panels'] for row in rows)
    for row in rows:
        previous=box(row['error']);omission=box(row['electron_tail'])+box(row['omitted_positron'])
        row['common_rule_error']=interval_text((I(row['panels'])/common)**8*(previous-omission)+omission)
    save('uniform-certificate.json',dict(classification='Proven',passed=all(r['passed'] for r in rows),rows=rows,common_panels=common,
        lower_normalized_Ip=interval_text(I(4)/(9*(1+iv.e))),
        proof='On each spatial disk |s-x|<=R=1/[2(b+1)], and independent parameter disks |delta eta|<=1/4, |delta beta|<=1/1000, the imaginary exponent stays below pi/2 and both gamma and the square-root argument have positive real parts. Bounds on the velocity rational function are monotone in the real beta center by the verified derivative. Multivariable Cauchy bounds the 8th spatial derivative and up to second parameter derivatives. The composite four-point Gauss remainder is summed over aligned coarse cells. Real tails use w<=1, |w_beta|<=2t, |w_betabeta|<=6t^2, logistic eta derivatives <=f, and exact integer exponential tails. Positron derivatives are bounded using exp(-2/beta) beta^-a monotonicity for 0<beta<=.006.',
        denominator_proof='For t in [0,1], exp(-min(eta,0))*f_e>=1/(1+e), sqrt(1+beta*t/2)>=1 and 1-w/3>=2/3; integrate sqrt(t) to obtain Ip>=4/[9(1+e)].',
        normalization_derivative_cusp_excluded=True,native_adaptive_binary_uniform_error_certified=False,physical_EOS_certified=False))
    print('PLASMA UNIFORM',len(rows),'components',all(r['passed'] for r in rows),'common panels',common,flush=True)
    assert all(r['passed'] for r in rows)


def corner():
    plan=json.loads((OUT/'plan.json').read_text());cert=json.loads((OUT/'uniform-certificate.json').read_text())
    assert cert['passed'];iv.dps=plan['precision_decimal_digits'];eta=I(24);beta=I('.006')
    selected=[next(r for r in cert['rows'] if r['component']==j and r['eta_order']==r['beta_order']==0) for j in plan['corner_components']]
    n=max(r['panels'] for r in selected);h=I(10)/n;nodes,weights=fermi.rule();totals=[I(0),I(0)]
    for j in range(n):
        center=h*(I(j)+I('.5'))
        for node,weight in zip(nodes,weights):
            s=center+h*node/2;t=s*s;x2=beta*t*(2+beta*t);w=x2/(1+beta*t)**2
            value=2*t*iv.sqrt(1+beta*t/2)*w**25/(1+iv.exp(t-eta))
            totals[0]+=h/2*weight*value;totals[1]+=h/2*weight*value*(1+2*x2/53)
        if j%2048==0:print('PLASMA INTERVAL CORNER',j,'/',n,flush=True)
    entries=[]
    for row,val in zip(selected,totals):
        omission=box(row['electron_tail'])+box(row['omitted_positron'])
        error=(I(row['panels'])/n)**8*(box(row['error'])-omission)+omission
        val_text=interval_text(val);error_text=interval_text(error)
        radius=box(error_text).b;enclosure=box(val_text)+iv.mpf([-radius,radius])
        entries.append(dict(component=row['component'],panels=n,quadrature=val_text,enclosure=interval_text(enclosure),error=error_text))
    D=box(cert['lower_normalized_Ip']).a;factor=1+iv.exp(-24)
    tailT=factor*box(entries[0]['enclosure']).b/(51*D)
    tailL=factor*box(entries[1]['enclosure']).b/D
    errs={r['component']:box(r['common_rule_error']).b for r in cert['rows'] if r['eta_order']==r['beta_order']==0}
    ED=errs[0]+errs[1]/3;assert lower(D-ED)>0
    coeff_errors=[I(0)]+[(errs[j]/(2*j+1)+errs[j+1]/(2*j+3)+I('1.5')*ED/(2*j+1))/(D-ED) for j in range(1,25)]
    KT=sum(coeff_errors)+tailT;KL=sum((2*j+1)*error for j,error in enumerate(coeff_errors))+tailL
    actual=np.load(dispersion.plasma.OUT/'stellar-fermi-plasma.npz')
    covered=bool(np.all((actual['eta']>=-17)&(actual['eta']<=24)&(actual['beta']>=0)&(actual['beta']<=.006)))
    save('result.json',dict(classification='Proven',completed=True,corner=entries,
        saved_5735_state_parameter_box_covered=covered,future_trajectory_enclosed=False,
        corner_domination='For all eta in [-17,24], exp(-min(eta,0))*f_e(t,eta)<=(1+exp(-24))*f_e(t,24). Positrons use the separate uniform omitted-positron bound, not this electron domination inequality.',
        corner_electron_tail_domination='Both positive tail integrands increase with beta. The upper corner enclosure includes a uniform omitted-positron bound; multiplying it by 1+exp(-24)>1 safely covers that bound.',
        transverse_uniform_tail=interval_text(tailT),longitudinal_uniform_tail=interval_text(tailL),
        coefficient_errors=[interval_text(x) for x in coeff_errors],transverse_uniform_kernel_error=interval_text(KT),
        longitudinal_uniform_kernel_error=interval_text(KL),
        passed=bool(upper(KT)<F(plan['uniform_kernel_absolute_error']) and upper(KL)<F(plan['uniform_kernel_absolute_error'])),
        exact_Gauss_rule=True,all_native_float_states_interval_evaluated=False,all_root_state_derivative_errors_certified=False,
        physical_EOS_certified=False,full_GR_evolution=False))
    audit();save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def audit():
    cert=json.loads((OUT/'uniform-certificate.json').read_text());result=json.loads((OUT/'result.json').read_text())
    def ends(s):
        a,b=map(F,s[1:-1].split(','));assert a<=b;return a,b
    for row in cert['rows']:
        assert ends(row['error'])[1]<=F('1e-14') and ends(row['common_rule_error'])[1]<=F('1e-14')
        for key in ['electron_tail','omitted_positron']:assert ends(row[key])[0]>=0
    for row in result['corner']:
        a,b=ends(row['quadrature']);lo,hi=ends(row['enclosure']);radius=ends(row['error'])[1]
        # The enclosure was formed after directed parsing of the exported operands.
        assert lo<=a-radius and hi>=b+radius
        assert max(a-radius-lo,hi-b-radius)<F('1e-55')
    for key in ['transverse_uniform_kernel_error','longitudinal_uniform_kernel_error']:
        assert ends(result[key])[1]<F('1e-10')
    save('audit.json',dict(classification='Proven',passed=True,uniform_components=len(cert['rows']),
        exact_decimal_endpoint_checks=True,scope='Exact-rational exported threshold/containment checks; proof uses directed mpmath interval arithmetic and the stated analytic remainder.'))
    print('PASS PLASMA UNIFORM AUDIT',len(cert['rows']),flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'audit.json').read_text())['passed']
    print('PASS conditional plasma moment certificate; native/physical EOS remain separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
