"""남은 레버의 대수 증명, 구간 근/미분 검사, 실제 변경량 감사."""
import json
import math
from pathlib import Path
import sys
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as s
from sympy.printing.pycode import MpmathPrinter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/remaining-levers15'
iv=mp.iv
iv.dps=70
mp.mp.dps=80


def endpoints(x):
    x=iv.mpf(x)
    return [math.nextafter(float(x.a),-math.inf),math.nextafter(float(x.b),math.inf)]


def exact_endpoints(x):
    def rational(t):
        sign,man,exp,_=t
        return str(Fraction((-1)**sign*man)*Fraction(2)**exp)
    return list(map(rational,iv.mpf(x)._mpi_))


def prove():
    # The old scalar term is broadcast along (1,1,1); it is not a vector.
    a=s.Matrix([1,2,3]);b=s.Matrix([2,0,1]);one=s.ones(3,1)
    R=s.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    defect=2*(R*a).dot(R*b)*one-R*(2*a.dot(b)*one)
    assert defect!=s.zeros(3,1) and R*(a+b)==R*a+R*b
    # mass * v^2 * position / c^2 is mass dipole; velocity has wrong units.
    expected=s.Matrix([1,1,0]) # M,L,T exponents
    correct=s.Matrix([1,0,0])+2*s.Matrix([0,1,-1])+s.Matrix([0,1,0])-2*s.Matrix([0,1,-1])
    wrong=correct-s.Matrix([0,1,0])+s.Matrix([0,1,-1])
    assert correct==expected and wrong!=expected
    m,M,f=s.symbols('m M f',positive=True)
    F=m**3-f*(m+M)**2
    fm=s.diff(F,m).subs(f,m**3/(m+M)**2).factor()
    dm_df=(m+M)**3/(m**2*(m+3*M));dm_dM=2*m/(m+3*M)
    assert s.simplify(fm*dm_df-(m+M)**2)==0
    assert s.simplify(fm*dm_dM-2*m**3/(m+M))==0
    # Kepler E-e sin E = mean anomaly: derivative is bounded below by 1-e.
    E,e,l=s.symbols('E e l',real=True)
    K=E-e*s.sin(E)-l
    assert s.simplify(s.diff(K,E)/(1-e*s.cos(E)))==1
    # u+D(u,theta)=a(theta), proper time = u-Einstein(u,theta).
    Du,Dt,at,Eu,Et=s.symbols('D_u D_theta a_theta E_u E_theta')
    up=(at-Dt)/(1+Du);tp=(1-Eu)*up-Et
    assert s.simplify(up+Du*up+Dt-at)==0
    # Exact linear-delay positive control including explicit parameter dependence.
    th=s.symbols('theta');aa=10+3*th;D=lambda u: s.Rational(1,4)*u+2*th
    u=(aa-2*th)/s.Rational(5,4);proper=s.Rational(9,10)*u-5*th
    assert s.diff(proper,th)==tp.subs({Du:s.Rational(1,4),Dt:2,at:3,Eu:s.Rational(1,10),Et:5})
    result=dict(classification='Proven',rotation_counterexample=list(map(str,defect)),mass_dipole_units=list(correct),old_term_units=list(wrong),positive_mass_root_derivative=str(fm),mass_derivative_f=str(dm_df),mass_derivative_center_mass=str(dm_dM),emission_parameter_derivative=str(up),proper_time_parameter_derivative=str(tp),linear_delay_control=True)
    (OUT/'symbolic-audit.json').write_text(json.dumps(result,indent=2,default=str)+'\n')


def root_box(f,M):
    """Positive root existence, uniqueness and interval Newton contraction."""
    assert f.a>0 and M.a>0, 'Positive mass-function domain required'
    fm=mp.mpf(f.mid._mpi_[0]);mm=mp.mpf(M.mid._mpi_[0])
    guess=mp.findroot(lambda x:x**3-fm*(x+mm)**2,(fm**(mp.mpf(1)/3)*mm**(mp.mpf(2)/3),mm+1))
    # A bracket is accepted only after directed interval sign checks.
    radius=abs(guess)*mp.mpf('1e-12')+mp.mpf('1e-30')
    lo=iv.mpf(str(guess-radius));hi=iv.mpf(str(guess+radius))
    fun=lambda x:x**3-f*(x+M)**2
    assert fun(lo).b<0 and fun(hi).a>0
    X=iv.mpf([lo.a,hi.b])
    for _ in range(12):
        mid=X.mid;df=3*X**2-2*f*(X+M)
        assert df.a>0
        Y=mid-fun(mid)/df
        X=iv.mpf([max(X.a,Y.a),min(X.b,Y.b)])
    assert X.a>0
    return X,(X+M)**3/(X**2*(X+3*M)),2*X/(X+3*M)


def masses():
    """Active parameter-set 6 mass map, including all ten contributing parameters."""
    p={line.split()[0]:line.split()[1] for line in (OUT/'parfile.txt').read_text().splitlines() if len(line.split())>1 and not line.startswith('#')}
    names=['apsini_i','period_i','absini_o','abcosi_o','period_o','masspar_p','delta_i','asini_extra1','acosi_extra1','P_extra1']
    vals=[iv.mpf(p[n])*(1+iv.mpf([-1,1])*iv.mpf(2)**-60) for n in names]
    assert all(vals[k].a>0 for k in [0,1,2,4,5,7,9])
    a,Pi,b,c,Po,q,delta,ae,ce,Pe=s.symbols('a Pi b c Po q delta ae ce Pe')
    xs=[a,Pi,b,c,Po,q,delta,ae,ce,Pe]
    G=s.Rational('6.67408e-11')*s.Rational('1.9884754153381438e30')
    C=s.Integer(299792458);day=s.Integer(86400)
    A=s.sqrt(b*b+c*c)*C
    ap=a*C*s.sqrt(b*b+c*c)/(b*s.cos(delta*s.pi/180)+c*s.sin(delta*s.pi/180))
    mb=(ap*(1+1/q))**3*(2*s.pi/(Pi*day))**2/G
    expressions=[mb/(1+q),q*mb/(1+q),4*s.pi**2*A**3/(Po**2*day**2*G),4*s.pi**2*(s.sqrt(ae*ae+ce*ce)*C)**3/(Pe**2*day**2*G)]
    module={'sin':iv.sin,'cos':iv.cos,'sqrt':iv.sqrt,'pi':iv.pi,'mpf':iv.mpf}
    # MpmathPrinter preserves rational constants; Python's p/q would round first.
    printer=MpmathPrinter({'fully_qualified_modules':False})
    calc=lambda expr:s.lambdify(xs,expr,modules=[module],printer=printer)(*vals)
    v=[calc(ex) for ex in expressions]
    assert all(x.a>0 for x in v)
    J=[[calc(s.diff(ex,x)) for x in xs] for ex in expressions]
    mo,dof,doM=root_box(v[2],v[0]+v[1])
    Jo=[dof*J[2][k]+doM*(J[0][k]+J[1][k]) for k in range(10)]
    mx,dxf,dxM=root_box(v[3],v[0]+v[1]+mo)
    Jx=[dxf*J[3][k]+dxM*(J[0][k]+J[1][k]+Jo[k]) for k in range(10)]
    result=dict(classification='Proven',scope='input_exact_rational이 인증 영역이다. input_boxes는 표시용 바깥 반올림이다. 정확한 실수 질량 함수와 해석적 미분이며 C++ 전체 초기화 반올림 인증은 아니다.',names=names,input_boxes=list(map(endpoints,vals)),input_exact_rational=list(map(exact_endpoints,vals)),mass_boxes=list(map(endpoints,[v[0],v[1],mo,mx])),mass_jacobian_boxes=[[endpoints(x) for x in row] for row in [J[0],J[1],Jo,Jx]],root_sign_and_interval_newton=True,all_28_initial_state_derivatives_certified=False)
    (OUT/'mass-map-certificate.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')


def compare():
    old=np.load(OUT/'original-initialization.npz');new=np.load(OUT/'corrected-initialization.npz')
    d=new['state']-old['state'];r=new['res']-old['res']
    rel=d-d[0] # Common COM recentering cannot remove relative changes.
    result=dict(classification='Proven',position_change_m=np.linalg.norm(d[:,:3],axis=1).tolist(),velocity_change_m_s=np.linalg.norm(d[:,3:],axis=1).tolist(),relative_to_pulsar_position_change_m=np.linalg.norm(rel[:,:3],axis=1).tolist(),residual_difference_rms_us=float(np.sqrt(np.mean(r*r))),residual_difference_max_us=float(np.max(abs(r))),corrected_residual_rms_us=float(np.sqrt(np.mean(new['res']**2))),original_residual_rms_us=float(np.sqrt(np.mean(old['res']**2))),refitted=False,full_physical_inference_complete=False)
    assert np.max(abs(rel))>0
    (OUT/'initialization-comparison.json').write_text(json.dumps(result,indent=2)+'\n')


def hex_fraction(token):
    sign=-1 if token.startswith('-') else 1
    mant,exponent=token.lstrip('+-').lower().split('p')
    whole,_,frac=mant.removeprefix('0x').partition('.')
    return sign*Fraction(int(whole+frac,16),16**len(frac))*Fraction(2)**int(exponent)


def readout():
    """Native call-point algebraic readouts and their interval partial derivatives."""
    def read(name):return list(map(hex_fraction,(OUT/name).read_text().split()))
    def rational(x):return s.Rational(x.numerator,x.denominator)
    def evaluate(names,values,expr,domains,native):
        symbols=s.symbols(' '.join(names))
        # Nonzero boxes deliberately include rounding of the sampled native inputs.
        boxes=[iv.mpf(v.numerator)/v.denominator+iv.mpf([-1,1])*(abs(iv.mpf(v.numerator)/v.denominator)*iv.mpf('1e-12')+iv.mpf('1e-25')) for v in values]
        module={'sin':iv.sin,'cos':iv.cos,'sqrt':iv.sqrt,'log':iv.log,'pi':iv.pi,'mpf':iv.mpf}
        printer=MpmathPrinter({'fully_qualified_modules':False})
        calc=lambda ex:s.lambdify(symbols,ex,modules=[module],printer=printer,cse=True)(*boxes)
        checked={key:calc(ex) for key,ex in domains.items()}
        assert all(x.a>0 for x in checked.values()),checked
        y=calc(expr);exact_native=iv.mpf(native.numerator)/native.denominator
        assert y.a<=exact_native.a<=exact_native.b<=y.b, (endpoints(y),float(native))
        jac=[calc(s.diff(expr,x)) for x in symbols]
        return dict(names=names,input_exact_rational=list(map(exact_endpoints,boxes)),output_days=endpoints(y),native_days=float(native),native_sample_contained=True,domain_lower_bounds={k:endpoints(v)[0] for k,v in checked.items()},partial_derivative_boxes=list(map(endpoints,jac)))
    data=read('geom-native.hex')
    assert len(data)==22
    _,bat,epoch,dist,dr,cday,year,c,rad=data[:9]
    vec=np.array(data[9:21],dtype=object).reshape(3,4)
    names=['bat','epoch','dist','dr',*[f'n{k}' for k in range(3)],*[f'pm{k}' for k in range(3)],*[f'o{k}' for k in range(3)],*[f'x{k}' for k in range(3)]]
    B,T,L,V,*rest=s.symbols(' '.join(names))
    n=s.Matrix(rest[:3]);pm=s.Matrix(rest[3:6]);obs=s.Matrix(rest[6:9]);x=s.Matrix(rest[9:12])
    dm=L*rational(year*c);dt=B-T;kp=pm*L*dt;kpa=V*rational(rad)*L*dt
    xp=x-x.dot(n)*n;op=obs-obs.dot(n)*n
    shk=kp.dot(kp)*rational(cday)/(2*dm)
    geo=x.dot(n)/rational(cday)+xp.dot(xp)/(2*dm*rational(cday))-xp.dot(op)/(dm*rational(cday))+kp.dot(xp)/dm+shk-kpa*rational(cday)/dm*shk
    geometry=evaluate(names,[bat,epoch,dist,dr,*vec[:,0],*vec[:,1],*vec[:,2],*vec[:,3]],geo,{'distance_m':dm},data[-1])
    data=read('nongem-native.hex');assert len(data)==30
    _,mp0,mi,mo,freq,c,twopi,day=data[:8]
    states=np.array(data[11:29],dtype=object).reshape(6,3)
    names=['mp','mi','mo','freq',*[f'n{k}' for k in range(3)],*[f'r{j}{k}' for j in range(3) for k in range(3)],*[f'v{j}{k}' for j in range(3) for k in range(3)]]
    mass=list(s.symbols(' '.join(names[:3])));ff=s.Symbol('freq');n=s.Matrix(s.symbols(' '.join(names[4:7])))
    rr=list(s.symbols(' '.join(names[7:16])));vv=list(s.symbols(' '.join(names[16:25])))
    r=[s.Matrix(rr[j*3:j*3+3]) for j in range(3)];v=[s.Matrix(vv[j*3:j*3+3]) for j in range(3)]
    angular=sum((mass[j]*r[j].cross(v[j]) for j in range(3)),s.zeros(3,1));cross=angular.cross(n)
    ab=v[0].dot(cross)*s.sqrt(angular.dot(angular))/(cross.dot(cross)*rational(twopi*c)*ff)
    args=[(s.sqrt((r[0]-r[j]).dot(r[0]-r[j]))-(r[0]-r[j]).dot(n))/rational(c) for j in [1,2]]
    sh=-2*s.Rational('4.92521372097374e-6')/rational(day)*(mass[1]*s.log(args[0])+mass[2]*s.log(args[1]))
    nongem=evaluate(names,[mp0,mi,mo,freq,*data[8:11],*states[:3].T.ravel(),*states[3:].T.ravel()],sh+ab,{'shapiro_inner_log_argument':args[0],'shapiro_outer_log_argument':args[1],'spin_cross_squared':cross.dot(cross),'spin_frequency':ff},data[-1])
    # Rigorous inversion control: nonlinear monotone F(u)=u+u^2/8-a,
    # derivative 1+u/4 >=1 on [0,1], known positive root for a=9/8 is 1.
    approx=iv.mpf(1)-iv.mpf(2)**-20
    residual=abs(approx+approx**2/8-iv.mpf(9)/8)
    error=1-approx
    assert error.b<=residual.b
    report=dict(classification='Proven',scope='추출된 실제 호출 지점의 대수 관측식과 구간 편미분. 적분·보간·Tempo2·매개변수 초기화의 전체 연쇄 인증은 아니다.',geometry=geometry,shapiro_aberration=nongem,nonlinear_inversion_error_control=True,einstein_integral_certified=False,extra_body_delay_included=False,full_timing_certificate=False)
    (OUT/'readout-certificate.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')


def kepler():
    p={line.split()[0]:line.split()[1] for line in (OUT/'parfile.txt').read_text().splitlines() if len(line.split())>1 and not line.startswith('#')}
    def par(n):return iv.mpf(p[n])*(1+iv.mpf([-1,1])*iv.mpf(2)**-60)
    C=iv.mpf(299792458);G=iv.mpf('6.67408e-11')*iv.mpf('1.9884754153381438e30')
    masses=json.loads((OUT/'mass-map-certificate.json').read_text())['mass_boxes']
    mp0,mi,mo,mx=[iv.mpf(pair) for pair in masses]
    aa=iv.sqrt(par('absini_o')**2+par('abcosi_o')**2)*C
    delta=par('delta_i')*iv.pi/180
    ap=par('apsini_i')*aa/(par('absini_o')*iv.cos(delta)+par('abcosi_o')*iv.sin(delta))
    ae=iv.sqrt(par('asini_extra1')**2+par('acosi_extra1')**2)*C
    result={}
    for label,eta,kappa,period,tasc,ar,m1,mc in [
      ('inner','eta_p','kappa_p','period_i','tasc_p',ap,mp0,mi),
      ('outer','eta_b','kappa_b','period_o','tasc_b',aa,mp0+mi,mo),
      ('extra','eta_extra1','kappa_extra1','P_extra1','tasc_extra1',ae,mp0+mi+mo,mx)]:
        e=iv.sqrt(par(eta)**2+par(kappa)**2);assert 0<e.a and e.b<1
        omega=iv.atan2(par(eta),par(kappa))
        mean=2*iv.pi*(par('treference')-par(tasc))/par(period)-omega
        fn=lambda E:E-e*iv.sin(E)-mean
        X=mean+iv.mpf([-1,1])*(e.b+iv.mpf('1e-40'))
        assert fn(iv.mpf(X.a)).b<0 and fn(iv.mpf(X.b)).a>0
        for _ in range(12):
            mid=X.mid;derivative=1-e*iv.cos(X);assert derivative.a>0
            new=mid-fn(mid)/derivative
            X=iv.mpf([max(X.a,new.a),min(X.b,new.b)])
        if label!='extra':
            total=m1+mc;nu=m1*mc/total**2;aRR=total/mc*ar
            eRR=e*(1+G*total/(aRR*C**2)*(4-iv.mpf(3)/2*nu))
            etheta=eRR*(1+G*m1*mc/(2*total*aRR*C**2))
        else:etheta=e
        q=etheta/(1+iv.sqrt(1-etheta**2));assert q.b<1
        terms=40;tail=2*q**(terms+1)/((terms+1)*(1-q));dtail=2*q**(terms+1)/(1-q)
        d=1-e*iv.cos(X)
        result[label]=dict(input_e_exact=exact_endpoints(e),input_mean_exact=exact_endpoints(mean),unwrapped_E=endpoints(X),dE_dmean=endpoints(1/d),dE_de=endpoints(iv.sin(X)/d),derivative_denominator_lower=endpoints(d)[0],fourier_terms=terms,angle_tail_upper=endpoints(tail)[1],derivative_tail_upper=endpoints(dtail)[1])
    (OUT/'kepler-certificate.json').write_text(json.dumps(dict(classification='Proven',scope='활성 세 궤도의 연속 Kepler 근과 미분, 1PN 각도 급수의 절대 나머지. 전체 초기 상태 미분 인증은 아님.',orbits=result,all_initialization_operations_certified=False),ensure_ascii=False,indent=2)+'\n')


def check():
    prove();masses();kepler();readout();compare()
    try:root_box(iv.mpf(-1),iv.mpf(1))
    except AssertionError:pass
    else:raise AssertionError('Negative mass-function must be rejected')
    assert hex_fraction('-0x1.8p+2')==Fraction(-6)
    assert hex_fraction('0xa.000000000000001p-5')!=Fraction(float.fromhex('0xa.000000000000001p-5'))
    # Outward reporting is independently checked with exact rational endpoints.
    for value in [iv.mpf(1)/3,iv.mpf([-1,1])*iv.mpf('1e-70')]:
        low,high=map(Fraction,exact_endpoints(value));a,b=map(Fraction,endpoints(value))
        assert a<=low<=high<=b
    result=dict(symbolic_identities=True,positive_mass_domains=True,three_kepler_domains=True,native_readout_samples_contained=True,outward_reporting_exact_fraction_checked=True,hexadecimal_input_not_rounded_to_binary64=True,full_timing_certificate=False)
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 대수·질량·Kepler·관측식·정확 입력·바깥 반올림 검사')


if __name__=='__main__':{'prove':prove,'masses':masses,'compare':compare,'readout':readout,'kepler':kepler,'check':check}[sys.argv[1]]()
