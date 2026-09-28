"""Bound predictor-to-neutral-root error separately from numerical quadrature."""
from fractions import Fraction as F
import json, math, sys
import numpy as np
import sympy as sp
from mpmath import iv
import gr_polarization_thermodynamics as thermo
import verify_electron_density_export as export
from interval_records import interval_text

g=thermo.g;density=thermo.density;OUT=g.OUT/'gr-polarization-neutral-error'
I=density.I;high=density.high;low=density.low


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def matrix(coeff):
    en,et,enn,ent,ett=coeff
    return [[1,0,0,0,0,0],[1,-et,0,-1,0,0],[0,en,0,0,0,0],
        [0,-et,0,-1,0,0],[0,et-ett,-et**2,1,-2*et,-1],
        [0,ent,en*et,0,en,0],[0,en+enn,en**2,0,0,0]]


def symbolic():
    q=sp.symbols('q',real=True);poly=1-6*q+6*q*q
    assert sp.expand(1-poly-6*q*(1-q))==0
    assert sp.expand(1+poly-(6*(q-sp.Rational(1,2))**2+sp.Rational(1,2)))==0
    x,s=sp.symbols('x s',positive=True)
    assert sp.factor(1/(4*s)-x/(x+s)**2)==(s-x)**2/(4*s*(s+x)**2)
    raw=sp.symbols('J Je Jee Jt Jet Jtt');coeff=sp.symbols('en et enn ent ett')
    expected=thermo.physical_fields(np.array(raw,dtype=object),np.array(coeff,dtype=object),1)
    for row,want in zip(matrix(coeff),expected):assert sp.expand(sum(a*b for a,b in zip(row,raw))-want)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        occupation='For 0<=q<=1, |q_eta|, |q_etaeta|, |q_etaetaeta|<=q. The cubic logistic factor uses -1<=1-6q+6q^2<=1.',
        response='With H_j=X*M_j/(X+S)^2, the six raw-integral absolute bounds are [I0,I0,3I0,H1,3H1,H1+3H2]. Their eta derivatives are bounded by [I0,3I0,13I0,3H1,13H1,3H1+13H2]. Here H_j denotes its z integral. These follow by differentiating S/(X+S), positivity and M1^2<=S*M2.',
        lower_response='For 0<=Q<=a choose pL>=a and pH>pL. K_p>=p^2/sqrt(1+p^2) on p>=Q. Thus R>=q(pH,eta_min)*(pH^3-pL^3)/(3*sqrt(1+pH^2)*Sref)=rmin>0.',
        integrated='Let c=a^2/(B*Sref), A_j=int M_j(Q)dQ/(a*Sref). Then I0<=0.5*(6*A0/c)^(1/3). On z<=1 use X/(X+S)^2<=1/(4S), and on z>=1 use <=1/X: int H_j dz<=max(1/(4*rmin),1/c)*A_j.',
        error='For the exact matrix M(coeff), reference center c_eta and true neutral root e, |field(e,coeff_true)-field(c_eta,coeff_saved)| <= pref*sum_i[|M_true_i|*L_i*|e-c_eta|+|M_true_i-M_saved_i|*B_i]. Coefficients at the true root use their certified intervals. No numerical-interior or floating-operation error is included.',
        units='All seven errors are per ion in the same reference kBT/kB units as the thermodynamic numerical table. Constants and stored coordinate scales are exact declared inputs; alpha/pi is evaluated outward.'))


def prepare():
    export.verify();assert not OUT.exists();OUT.mkdir()
    assert json.loads((thermo.OUT/'numerical.json').read_text())['passed']
    paths=[g.ROOT/'verification/gr_polarization_neutral_error.py',
        g.ROOT/'verification/gr_polarization_thermodynamics.py',thermo.OUT/'plan.json',
        thermo.OUT/'numerical.json',thermo.OUT/'states.npz',export.OUT/'manifest.json',
        export.OUT/'result.json',density.OUT/'candidates.npz',density.ionic.OUT/'inputs.npz',density.ionic.OUT/'constants.json']
    save('plan.json',dict(classification='Proven',checkpoint='5e4a0c7',bits=128,cells=3206,
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        reference_budget='2e-7',budget_origin='Same absolute per-field budget as the fixed self-thermodynamic finite refinement. A failed budget is retained as a failure; mathematical bounds remain valid.',
        method='Bound analytic center-to-root displacement and inverse-derivative coefficient errors using a positive response lower bound on 0<=Q<=a and all-Q thermal moments. No new quadrature or physical EOS approximation.',
        scope='Exact declared ideal-electron RPA self term only. Does not certify the numerical interior, floating evaluations, full correlated EOS, GR evolution or observations.'))
    symbolic()


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return p


def run():
    plan=bindings();iv.prec=plan['bits'];a=thermo.inputs()
    roots=json.loads((export.OUT/'result.json').read_text())['rows']
    states=dict(np.load(thermo.OUT/'states.npz'));ions=dict(np.load(density.ionic.OUT/'inputs.npz'))
    const=json.loads((density.ionic.OUT/'constants.json').read_text())['native_binary64'];B=I(const['alpha'])/iv.pi
    records=[];maxima=[F(0)]*7;center_max=[F(0)]*7;coeff_max=[F(0)]*7
    for i,r in enumerate(roots):
        assert r['position']==i and r['cell']==int(a['cells'][i])
        elo,ehi=export.ends(r['root']);center=F.from_float(float(a['eta_center'][i]))
        emin,emax=I(min(elo,center)),I(max(ehi,center));delta=I(F(r['root_error_from_center_upper']))
        beta=I(a['beta'][i]);scale=I(a['scale'][i]);Sref=I(a['Sref'][i]);cut=I(max(F(0),high(emax)))
        def moment(k):
            return cut**(k+1)/(k+1)+iv.exp(emax-cut)*sum(I(math.factorial(k))/math.factorial(j)*cut**j for j in range(k+1))
        A=[iv.pi**2*beta/4*(moment(j)+2*beta*moment(j+1)+beta**2*moment(j+2))/(scale*Sref) for j in range(3)]
        pl=I(max(high(scale),high(iv.sqrt(beta))));ph=pl+iv.sqrt(beta);gamma=iv.sqrt(1+ph**2)
        qmin=1/(1+iv.exp(ph**2/(beta*(gamma+1))-emin))
        rmin=qmin*(ph**3-pl**3)/(3*gamma*Sref);assert low(rmin)>0
        c=scale**2/(B*Sref);factor=I(max(high(1/(4*rmin)),high(1/c)))
        H1,H2=factor*A[1],factor*A[2];J=(6*A[0]/c)**(I(1)/3)/2
        bounds=[J,J,3*J,H1,3*H1,H1+3*H2]
        lipschitz=[J,3*J,13*J,3*H1,13*H1,3*H1+13*H2]
        coeff=[density.read_interval(r['derivatives'][name][1:-1].split(',')) for name in density.FIELDS]
        actual=matrix(coeff);reference=matrix([I(v) for v in states['neutral_coefficients'][:,i]])
        counts=[I(x)/I(m) for x,m in zip(ions['X'][i],ions['A'])]
        Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts);pref=2*B*Z2*scale/beta
        shift=[pref*delta*sum(abs(I(v))*bound for v,bound in zip(row,lipschitz)) for row in actual]
        inverse=[pref*sum(abs(I(v)-I(w))*bound for v,w,bound in zip(row,ref,bounds)) for row,ref in zip(actual,reference)]
        total=[x+y for x,y in zip(shift,inverse)]
        maxima=[max(m,high(v)) for m,v in zip(maxima,total)]
        center_max=[max(m,high(v)) for m,v in zip(center_max,shift)]
        coeff_max=[max(m,high(v)) for m,v in zip(coeff_max,inverse)]
        records.append(dict(cell=int(a['cells'][i]),positive_response_lower=interval_text(rmin),
            field_error_upper={name:interval_text(v) for name,v in zip(thermo.FIELDS,total)}))
    save('result.json',dict(classification='Proven',cells=len(records),records=records,mathematical_bounds_valid=True,
        budget_passed=max(maxima)<F(plan['reference_budget']),maximum_exact=dict(zip(thermo.FIELDS,map(str,maxima))),
        display_only=dict(total=dict(zip(thermo.FIELDS,map(float,maxima))),center=dict(zip(thermo.FIELDS,map(float,center_max))),
            inverse_coefficients=dict(zip(thermo.FIELDS,map(float,coeff_max)))),
        interior_quadrature_certified=False,full_physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['mathematical_bounds_valid'] and r['cells']==p['cells']
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert r['budget_passed']==(max(map(F,r['maximum_exact'].values()))<F(p['reference_budget']))
    print('PASS analytic predictor displacement/neutral coefficient bounds; budget passed:',r['budget_passed'],r['display_only'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
