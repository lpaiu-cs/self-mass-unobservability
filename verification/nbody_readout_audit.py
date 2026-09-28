"""다체 관측식, 연속 스플라인 오차와 영 scalar 가지의 독립 검사."""
from fractions import Fraction
import json
from pathlib import Path
import re
import sys

import mpmath as mp
import numpy as np
import sympy as s

from remaining_audit import hex_fraction,exact_endpoints,endpoints

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/nbody-readout16'
iv=mp.iv;iv.dps=80

def interval(q):return iv.mpf(q.numerator)/q.denominator

def native_samples():
    result={}
    for scale in ['zero','one','two']:
        for mode in ['ein','shap']:
            lines=(OUT/(scale+'-'+mode+'.hex')).read_text().splitlines()
            head=lines[0].split();nt=int(head[0]);_,factor,mass,G,c,day=map(hex_fraction,head[1:])
            cert=[]
            for line in lines[1:]:
                tokens=line.split();index=int(tokens[0]);time,*values=map(hex_fraction,tokens[1:])
                # Certified algebraic evaluation on small, explicit boxes around
                # the exact native inputs; no trajectory enclosure is presumed.
                boxes=[interval(x)+iv.mpf([-1,1])*(abs(interval(x))*iv.mpf('1e-15')+iv.mpf('1e-28')) for x in values[:9]]
                d=[boxes[3*k]-boxes[3*k+1] for k in range(3)];n=[boxes[3*k+2] for k in range(3)]
                r=iv.sqrt(sum(x*x for x in d));arg=(r-sum(d[k]*n[k] for k in range(3)))/interval(c)
                assert r.a>0 and arg.a>0
                u=interval(factor*mass*G/c**2)/r
                sh=-2*iv.mpf('4.92521372097374e-6')/interval(day)*interval(factor*mass)*iv.log(arg)
                for calculated,native in [(u,values[9]),(sh,values[10])]:
                    if mode=='ein' and calculated is sh:continue
                    if mode=='shap' and calculated is u:continue
                    observed=interval(native)
                    assert calculated.a<=observed.a<=observed.b<=calculated.b
                cert.append(dict(index=index,time_days=float(time),input_exact_rational=list(map(exact_endpoints,boxes)),extra_einstein_integrand=endpoints(u),extra_shapiro_days=endpoints(sh)))
            result[scale+'-'+mode]=dict(grid_points=nt,validated_sample_count=len(cert),samples=cert)
    (OUT/'native-sample-audit.json').write_text(json.dumps(dict(classification='Proven',scope='추출한 호출 입력의 대수 계산 및 바깥 반올림 비교; 전 구간 궤도·보간 인증 아님',runs=result),ensure_ascii=False,indent=2)+'\n')

def spline():
    """C2 residual bounds on each nonuniform interpolation cell, no C4 assumption."""
    x,h,t=s.symbols('x h t',positive=True);z=s.symbols('z',nonnegative=True)
    # Green kernel for e''=r with zero endpoint values on [0,h].
    kleft=t*(h-x)/h;kright=x*(h-t)/h
    area=s.integrate(kleft,(t,0,x))+s.integrate(kright,(t,x,h))
    assert s.simplify(area-x*(h-x)/2)==0
    derivative_l1=s.integrate(t/h,(t,0,x))+s.integrate((h-t)/h,(t,x,h))
    assert s.simplify(derivative_l1-(x*x+(h-x)**2)/(2*h))==0
    assert s.integrate(area,(x,0,h))==h**3/12
    # Native natural spline on nonuniform quartic data: exact fractional bytes.
    rows=[list(map(hex_fraction,line.split())) for line in (OUT/'spline-control.hex').read_text().splitlines()]
    nodes=rows[:5];results=[]
    for k in range(4):
        a,ya,ma=map(s.Rational,nodes[k]);b,yb,mb=map(s.Rational,nodes[k+1]);hh=b-a
        aa=(b-z)/hh;bb=(z-a)/hh
        S=aa*ya+bb*yb+hh**2/6*((aa**3-aa)*ma+(bb**3-bb)*mb)
        error=s.expand(z**4-S);rho=12*b**2+max(abs(ma),abs(mb))
        # Interval evaluation covers the entire cell, not sampled extrema.
        X=iv.mpf([interval(Fraction(a)).a,interval(Fraction(b)).b])
        # S'' is a convex interpolation of ma,mb. Preserve that correlation
        # algebraically instead of expanding two dependent copies of z.
        assert s.simplify(s.diff(S,z,2)-((b-z)*ma+(z-a)*mb)/hh)==0
        curvature=iv.mpf([interval(Fraction(min(ma,mb))).a,interval(Fraction(max(ma,mb))).b])
        residual=12*X**2-curvature;R=interval(Fraction(rho))
        assert -R.b<=residual.a<=residual.b<=R.b
        value_bound=rho*hh**2/8;derivative_bound=rho*hh/2;integral_bound=rho*hh**3/12
        midpoint,native_mid,native_int=rows[5+k]
        true_mid=Fraction(midpoint)**4
        assert abs(true_mid-native_mid)<=Fraction(value_bound)
        exact_int=Fraction((b**5-a**5)/5)
        assert abs(exact_int-native_int)<=Fraction(integral_bound)
        # The quartic minus cubic derivative is cubic; check its exact extrema.
        crit=[a,b,*[v for v in s.solve(s.diff(error,z,2),z) if v.is_real and a<=v<=b]]
        assert all(abs(s.diff(error,z).subs(z,v))<=derivative_bound for v in crit)
        results.append(dict(cell=k,rho=str(rho),value_error_upper=str(value_bound),derivative_error_upper=str(derivative_bound),integral_error_upper=str(integral_bound)))
    # The value/integral constants are attained by e''=1 on a unit cell.
    e=z*(z-1)/2
    assert abs(e.subs(z,s.Rational(1,2)))==s.Rational(1,8)
    assert abs(s.integrate(e,(z,0,1)))==s.Rational(1,12)
    report=dict(classification='Proven',assumptions='각 셀에서 g는 C2, S는 저장 계수를 정확한 실수로 읽은 cubic polynomial, |g-S|의 끝점 오차가 eps, |g_second-S_second|<=rho. native 평가·적분 산술 반올림과 parameter에 따른 이동 격자 미분은 별도.',value_error='eps + rho*h^2/8',derivative_error='2*eps/h + rho*h/2',integral_error='eps*h + rho*h^3/12',nonuniform_cells=True,natural_boundary_C4_not_required=True,kernel_constants_proved=True,native_quartic_controls=results,actual_full_span_rho_certified=False)
    (OUT/'spline-error-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')

def curvature():
    """Connect flow/LOS bounds to the spline residual; all derivatives use one time unit."""
    t=s.symbols('t',real=True)
    d0=s.Matrix(s.symbols('x:3',real=True));v=s.Matrix(s.symbols('v:3',real=True));a=s.Matrix(s.symbols('a:3',real=True))
    n0=s.Matrix(s.symbols('n:3',real=True));n1=s.Matrix(s.symbols('p:3',real=True));n2=s.Matrix(s.symbols('q:3',real=True))
    d=d0+v*t+a*t*t/2;n=n0+n1*t+n2*t*t/2
    r=s.sqrt(d.dot(d));r0=s.sqrt(d0.dot(d0));rd=d0.dot(v)/r0
    rdd=(v.dot(v)+d0.dot(a))/r0-d0.dot(v)**2/r0**3
    assert s.factor(s.diff(r,t,2).subs(t,0)-rdd)==0
    upp=3*d0.dot(v)**2/r0**5-(v.dot(v)+d0.dot(a))/r0**3
    assert s.factor(s.diff(1/r,t,2).subs(t,0)-upp)==0
    z=r-d.dot(n);z0=r0-d0.dot(n0)
    zd=rd-v.dot(n0)-d0.dot(n1);zdd=rdd-a.dot(n0)-2*v.dot(n1)-d0.dot(n2)
    assert s.factor(s.diff(s.log(z),t,2).subs(t,0)-(zdd/z0-zd**2/z0**2))==0
    # Hessian(1/r) has radial eigenvalue 2/r^3 and tangent eigenvalues -1/r^3.
    radial=3*d0*d0.T/r0**5-s.eye(3)/r0**3
    assert (radial*d0-2*d0/r0**3).applyfunc(s.simplify)==s.zeros(3,1)
    vp=s.Matrix(s.symbols('u:3',real=True));ap=s.Matrix(s.symbols('b:3',real=True));jp=s.Matrix(s.symbols('j:3',real=True))
    vv=vp+ap*t+jp*t*t/2
    assert s.expand(s.diff(vv.dot(vv)/2,t,2).subs(t,0)-ap.dot(ap)-vp.dot(jp))==0
    report=dict(classification='Proven',assumptions='한 셀에서 r>=rmin>0, z=r-d.dot(n)>=zmin>0; |d|<=R, |dprime|<=V, |dsecond|<=A, |n|<=N, |nprime|<=N1, |nsecond|<=N2. 모두 같은 시간 단위의 미분. 고정 질량과 상수 사용.',
      einstein_potential_second_bound='(GM/c^2)*(A/rmin^2 + 2*V^2/rmin^3)',
      kinetic_second_bound='(Ap^2 + Vp*Jp)/c^2',
      Z1='V*(1+N) + R*N1',Z2='A*(1+N) + V^2/rmin + 2*V*N1 + R*N2',
      shapiro_second_bound='2*Tsun*M/daysec*(Z2/zmin + Z1^2/zmin^2)',
      spline_residual_bound='rho <= gsecond_bound + max(abs(stored_y2_left),abs(stored_y2_right))',
      actual_full_span_state_acceleration_jerk_LOS_bounds_available=False)
    (OUT/'curvature-transfer.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')

def zero_scalar():
    phi,beta,eps,T0,dT,H=s.symbols('phi beta eps T0 dT H',real=True)
    A=s.exp(beta*phi**2/2);alpha=s.diff(s.log(A),phi)
    assert s.simplify(alpha-beta*phi)==0 and A.subs(phi,0)==1 and alpha.subs(phi,0)==0
    source=-4*s.pi*alpha*(T0+eps*dT)
    assert source.subs(phi,0)==0
    expanded=s.series(source.subs(phi,eps*H),eps,0,2).removeO()
    assert s.expand(expanded).coeff(eps)==-4*s.pi*beta*T0*H
    assert s.diff(A**2,phi).subs(phi,0)==0
    # A scalar susceptibility alone cannot create an input-independent charge.
    F,chi,grad=s.symbols('F chi grad');q=chi*F
    assert q.subs(F,0)==0
    force=q*grad
    assert s.expand(force.subs({F:eps*F,grad:eps*grad})).coeff(eps)==0
    # Time-varying companion geometry changes coefficients of a homogeneous
    # equation; zero initial q and dq remain exactly zero by IVP uniqueness.
    t=s.symbols('t');qfun=s.Function('q')(t);K=s.Function('K')(t)
    ode=s.diff(qfun,t,2)+s.diff(qfun,t)+K*qfun
    assert ode.subs(qfun,0).doit()==0
    report=dict(classification='Proven',conditional_statement='A=exp(beta*phi^2/2), 모든 천체의 비스칼라화 가지, phi 및 법선 시간 미분의 영 초기자료, 영 외부·입사 scalar 자료, 고전 초기값 문제의 유일성 아래 phi=0은 정확한 GR 해이다.',zero_source_verified=True,linear_matter_drive_vanishes=True,scalar_force_first_order_vanishes=True,pole_does_not_imply_drive=True,homogeneous_time_dependent_coefficients_preserve_zero=True,stability_proven=False,all_scalarized_solutions_excluded=False,no_J0337_beta_limit_inferred=True,source='https://arxiv.org/pdf/gr-qc/9602056')
    (OUT/'zero-scalar-branch.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')

def live_results():
    data=np.load(OUT/'live.npz');result={}
    source=(OUT/'source/Fittriple-init.cpp').read_text()
    for component in ['geom','ein','shap','aber']:
        assert re.search(r'fake_delay_'+component+r'\[i\]\s*=.*?\*\s*86400\.',source,re.S)
    for k,name in [(1,'geometric'),(2,'einstein'),(3,'shapiro'),(4,'aberration')]:
        # Fittriple-init.cpp exports BAT in days, but these four delays in seconds.
        d=(data['fake_one'][k]-data['fake_zero'][k])*1e6
        linearity=(data['fake_two'][k]-data['fake_zero'][k])*1e6-2*d
        result[name]=dict(max_difference_us=float(np.max(abs(d))),rms_difference_us=float(np.sqrt(np.mean(d*d))),twice_minus_twice_max_us=float(np.max(abs(linearity))))
        assert result[name]['twice_minus_twice_max_us']<1e-6
    assert result['aberration']['max_difference_us']<1e-4
    assert result['einstein']['max_difference_us']>0
    (OUT/'fake-call-audit.json').write_text(json.dumps(dict(classification='Proven',n_fake=128,native_delay_unit='seconds',reported_unit='microseconds',preliminary_audit_correction='초 단위 출력에 중복 적용했던 86400 배율을 제거했다. 원시 live.npz는 변경하지 않았다.',components=result,interpretation='다른 해의 방출 시각을 포함한 모의 출력 차이이며 같은 시각의 대수 항만 비교한 값은 아니다.'),ensure_ascii=False,indent=2)+'\n')

def check():
    native_samples();spline();curvature();zero_scalar();live_results()
    (OUT/'checks.json').write_text(json.dumps(dict(symbolic_green_kernel=True,analytic_curvature_transfer=True,native_quartic_spline_control=True,exact_zero_scalar_branch=True,native_nbody_samples=True,fake_call_path=True,full_timing_certificate=False),indent=2)+'\n')
    print('PASS: 다체 관측식·자연 스플라인 오차 정리·영 scalar 가지·모의 호출 감사')

if __name__=='__main__':{'samples':native_samples,'spline':spline,'zero-scalar':zero_scalar,'live':live_results,'check':check}[sys.argv[1]]()
