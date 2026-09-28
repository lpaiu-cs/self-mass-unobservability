"""Joint regular-response domains, including the near-zero screening scale."""
from fractions import Fraction as F
import json,math,sys
import numpy as np
import sympy as sp
from mpmath import iv
import gr_response_regular_runner as regular
import gr_highq_self_certificate as highq

cusp=highq.cusp;ROOT=highq.ROOT;OUT=highq.OUT.parent/'gr-self-regular-domain';I=highq.I;low=highq.low;high=highq.high;sha=highq.sha
ORDERS=[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    regular.verify();highq.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_self_regular_domain.py',regular.OUT/'manifest.json',regular.OUT/'result.json',highq.OUT/'manifest.json',
        cusp.ROOTS/'result.json',cusp.THERMO/'states.npz',highq.moments.OUT/'inputs.npz',highq.sector.OUT/'result.json',
        highq.neutral.density.ionic.OUT/'inputs.npz',highq.neutral.density.ionic.OUT/'constants.json',
        cusp.OUT.parent/'gr-smooth-response-defined/rule.json']
    save('plan.json',dict(classification='Proven',bits=128,strip_in_sqrt_beta='1/16',eta_radius='1/100',tau_radius='1/1000',thermal_split=176,
        checkpoint='3d54106',bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},cells=3206,control_positions=[0,1603,3205],
        outer_quadrature_field_budget='2e-8',
        scope='Prove the joint domain and create three exact dyadic Gauss panel schedules only. No actual outer-node response evaluation, point-error budget or full finite self integral is certified here.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def dyadic_below(x):
    assert x>0;e=x.numerator.bit_length()-x.denominator.bit_length();value=F(2)**e
    while value>x:value/=2
    while value*2<=x:value*=2
    assert value<=x<2*value;return value


def run():
    p=bindings();iv.prec=p['bits'];delta=I(p['strip_in_sqrt_beta']);he=I(p['eta_radius']);ht=I(p['tau_radius']);Tc=I(p['thermal_split'])
    sigma=delta*iv.sqrt(I('.006'));kappa=1/iv.sqrt(1-sigma*sigma);plus=iv.sqrt(1+sigma*sigma);D=delta*delta*kappa/2
    E=iv.exp(ht)-1;cl=1-E*(1+delta*kappa);cu=1+E*(1+delta*kappa);Dl=(1+E)*D+E*delta*kappa/2;Du=E*delta*kappa/2
    phase=(1+E)*delta*kappa*iv.sqrt(2*Tc)+E*(Tc+D)+he+iv.atan2(sigma*kappa*kappa,I(1))+ht
    assert low(cl)>0 and high(phase)<low(iv.pi/2)
    L=iv.cos(phase)/plus*iv.exp(-ht-(cu-1)*Tc-Du-he)
    old=json.loads((highq.sector.OUT/'result.json').read_text());W=I(old['uniform_W0_lower']);C=I(old['uniform_real_response_upper_coefficient'])/cl**4
    tail=iv.exp(24+he+Dl-cl*Tc/2)*C;positive=L*W-(2*L+4*kappa*iv.exp(ht))*tail;assert low(positive)>0
    majorants=[2*math.factorial(n)*math.factorial(m)/(he**n*ht**m) for n,m in ORDERS]
    state,ions,roots,data,B,_=highq.context();rows=[];minimum=None;maximum=None;schedules=[]
    norm=I(1)/((2*16+1)*math.comb(32,16)**2)
    for i,root in enumerate(roots):
        assert root['cell']==int(state['cells'][i]);beta=I(data['beta'][i]);eta=cusp.interval(root['root']);s=delta*iv.sqrt(beta)
        band=iv.sqrt(beta)/(2*iv.sqrt(1+4*beta)*(1+iv.exp(2-I(low(eta)))))
        negative_tail=4*kappa*iv.exp(ht+I(high(eta))+he+Dl-cl*Tc/2)*C;M=L*band-negative_tail;assert low(M)>0
        qD=dyadic_below(min(low(s/4),low(iv.sqrt(B*M)/4)));assert high(I(qD)**2/B)<=low(M)/16
        coeff=[cusp.interval(root['derivatives'][name]) for name in highq.neutral.density.FIELDS]
        counts=[I(x)/I(a) for x,a in zip(ions['X'][i],ions['A'])];Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts)
        scale=I(state['scale'][i]);pref=2*B*Z2*scale/beta
        field_bounds=[pref*sum(abs(c)*m for c,m in zip(line,majorants))/scale for line in highq.neutral.matrix(coeff)]
        row=dict(position=i,cell=root['cell'],strip_radius_lower=str(low(s)),small_Q_response_lower=str(low(M)),screening_disk_radius=str(qD),
            physical_integrand_majorants=[str(high(x)) for x in field_bounds]);rows.append(row)
        minimum=qD if minimum is None else min(minimum,qD);maximum=qD if maximum is None else max(maximum,qD)
        if i not in p['control_positions']:continue
        a=F.from_float(float(state['scale'][i]))/2**24;b=4*F.from_float(float(data['pcut'][i]));length=b-a
        stack=[(a,b)];panels=[];error=I(0);field=I(max(map(high,field_bounds)))
        while stack:
            lo,hi=stack.pop();center=(lo+hi)/2;h=(hi-lo)/2
            additive=qD if center<=2*qD else min(low(s),center/2);radius=max(additive,center/256)
            if radius>h:
                remainder=I(hi-lo)*norm*(I(hi-lo)/I(radius-h))**32*field
                if high(remainder)<=F(p['outer_quadrature_field_budget'])*(hi-lo)/length:
                    panels.append([str(lo),str(hi)]);error+=remainder;continue
            assert h>0 and len(panels)<100000;stack.extend([(center,hi),(lo,center)])
        panels.sort(key=lambda x:F(x[0]));assert F(panels[0][0])==a and F(panels[-1][1])==b
        assert all(F(x[1])==F(y[0]) for x,y in zip(panels,panels[1:]));assert high(error)<F(p['outer_quadrature_field_budget'])
        schedules.append(dict(position=i,cell=root['cell'],panels=panels,nodes=16*len(panels),field_remainder_upper=str(high(error)),passed=True))
    x=sp.symbols('x',nonnegative=True);assert sp.simplify((sp.sqrt(x)-1/sp.sqrt(2))**2-(x+sp.Rational(1,2)-sp.sqrt(2*x)))==0
    assert F(1)+F(9,15)==F(8,5)<2
    save('result.json',dict(classification='Proven',passed=True,cells=len(rows),rows=rows,schedules=schedules,
        phase_upper=str(high(phase)),positive_response_lower_coefficient=str(low(positive)),derivative_majorants=[str(high(m)) for m in majorants],
        display_only=dict(phase_upper=float(high(phase)),positive_coefficient=float(low(positive)),screening_radius_range=[float(minimum),float(maximum)],
            schedule_nodes=[r['nodes'] for r in schedules]),
        tensor_identity='In the positive double integral set u=w^2 and v=p*w. The Jacobian cancels 1/2, giving S(Q)=integral_0^1 dw integral_0^infinity dv H(sqrt(v^2+(1-w^2)*Q^2)). H is evaluated through gamma=sqrt(1+v^2+(1-w^2)*Q^2), f and log(1+exp(eta-t)); no logarithmic cusp or polylogarithm remains. The u=0 edge has measure zero, and positivity/dominance justifies the substitution and order exchange.',
        joint_energy='Let the former geometric t be t_g, t_a-D<=Re t_g<=t_a and |Im t_g|<=delta*kappa*sqrt(2t_a). For |tau|<=h_tau, E=exp(h_tau)-1 bounds both |Re(exp(-tau))-1| and |Im(exp(-tau))|. Since sqrt(2t_a)<=t_a+1/2, c_l*t_a-D_l<=Re(exp(-tau)*t_g)<=c_u*t_a+D_u, with the encoded c_l,c_u,D_l,D_u. Inner phase adds (1+E)delta*kappa*sqrt(2Tc)+E(Tc+D)+h_eta, inverse-gamma phase and beta phase h_tau.',
        joint_inner='The occupation magnitude relative to the real reference is >=exp(-(c_u-1)Tc-D_u-h_eta). The V1 factor beta contributes exp(-h_tau) in magnitude and h_tau in phase. The saved common L therefore gives Re H_joint>=L H_real throughout t_a<=Tc.',
        joint_tail='Beyond Tc, compare with real eta=0,beta=2beta0/c_l. Its response coefficient is at most C_old/c_l^4. Real and complex omitted parts are bounded by [2L+4kappa exp(h_tau)]exp(24+h_eta+D_l-c_l Tc/2) times that coefficient. Subtracting them leaves a uniform positive Re S lower throughout the joint strip.',
        small_Q='For |Re Q|<=s=sqrt(beta0)/16, restrict the positive real double integral to p in [sqrt(beta0),2sqrt(beta0)] and u in [0,1]. Here t_a<=2, yielding band=sqrt(beta0)/(2sqrt(1+4beta0)(1+exp(2-eta_min))). Subtract the full complex outer tail to get M>0 uniformly over the eta/tau polydisk. Choose an exact downward dyadic qD<=min(s/4,sqrt(B*M)/4).',
        dielectric_domains='For real q0<=2qD, the disk |Q-q0|<=qD has |Re Q|<=3qD<s and Re(Q^2/B+S)>=15M/16. Hence |G|<=1+9/15=8/5. For q0>=2qD use radius min(s,q0/2): |arg Q|<=pi/6, Re S>0 and the phase geometry yields |G|<=sec(pi/3)=2. Both hold jointly in eta/tau. The older relative q0/256 joint disk can be used wherever it is larger, with the same upper bound 2.',
        derivatives='Cauchy in eta/tau gives uniform complex-Q bounds 2*n!*m!/(h_eta^n*h_tau^m). Each true neutral Hessian and ion inventory converts these to physical integrand majorants. On a Q panel, the exact 16-point Gauss remainder is width*norm16*(width/(R-half_width))^32 times that majorant. Three saved dyadic schedules have exact coverage and summed remainders below 2e-8, conditional on evaluating the actual nodes with separately bounded errors.',
        actual_outer_nodes_evaluated=False,finite_outer_integral_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['cells']==p['cells'] and F(r['positive_response_lower_coefficient'])>0
    for i,row in enumerate(r['rows']):assert row['position']==i and F(row['screening_disk_radius'])>0 and F(row['small_Q_response_lower'])>0
    for row in r['schedules']:assert row['passed'] and F(row['field_remainder_upper'])<F(p['outer_quadrature_field_budget']) and row['nodes']==16*len(row['panels'])
    print('PASS screening-scale joint G domains and three exact outer schedules; actual nodes/full finite integral remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
