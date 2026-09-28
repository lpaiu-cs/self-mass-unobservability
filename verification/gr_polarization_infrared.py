"""Rigorous infrared endpoint of the same neutral self free energy."""
from fractions import Fraction as F
import json,sys
import numpy as np
import sympy as sp
from mpmath import iv
import gr_response_complex_runner as complex_domain
import gr_polarization_neutral_error as neutral

g=complex_domain.g;cusp=g.cusp;ROOT=g.ROOT;OUT=g.OUT.parent/'gr-polarization-infrared'
I=g.I;low=g.low;high=g.high;density=neutral.density;FIELDS=neutral.thermo.FIELDS
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    complex_domain.verify();assert not OUT.exists();OUT.mkdir()
    paths=[ROOT/'verification/gr_polarization_infrared.py',complex_domain.OUT/'manifest.json',complex_domain.OUT/'result.json',
        cusp.ROOTS/'manifest.json',cusp.ROOTS/'result.json',cusp.THERMO/'states.npz',density.ionic.OUT/'inputs.npz',density.ionic.OUT/'constants.json',
        ROOT/'verification/gr_polarization_neutral_error.py',ROOT/'verification/gr_polarization_thermodynamics.py']
    save('plan.json',dict(classification='Proven',checkpoint='e37efbd',bits=128,z_cut='1/16777216',field_budget='2e-7',
        bindings={p.relative_to(ROOT).as_posix():g.sha(p) for p in paths},
        target='Approximate the six raw self-integral contributions on 0<=Q<=scale*z_cut by [z_cut,0,0,0,0,0] after dividing the Q integral by the fixed scale. Propagate certified cubic remainder bounds through the true neutral Hessian and exact ion charge inventory into seven physical fields.',
        lower_response='For p>=Q, log((p+Q)/(p-Q))>=2Q/p gives K_p>=1/gamma. On p in [sqrt(beta),2sqrt(beta)], Q<=sqrt(beta) implies R=S/Sref>=sqrt(beta)/(sqrt(1+4beta)*Sref*(1+exp(2-eta_min)))=L>0.',
        endpoint='For X=Q^2/(B*Sref), 0<=1-G<=X/L. The five derivative magnitudes are <=X times the quotient-rule majorants with denominator L. Their integral divided by scale is bounded by Omega=scale^2*z_cut^3/(3B*Sref) times those constants. The value uses 1/L.',
        scope='Same full-ionization ideal-electron self term, exact declared inventory and alpha/pi, certified neutral eta and implicit derivative intervals. Only the infrared endpoint contribution; finite outer integration, full correlated EOS and GR/observations remain open.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def run():
    plan=bindings();iv.prec=plan['bits'];z=I(plan['z_cut']);state=dict(np.load(cusp.THERMO/'states.npz'));ions=dict(np.load(density.ionic.OUT/'inputs.npz'))
    roots=json.loads((cusp.ROOTS/'result.json').read_text())['rows'];complex_rows=json.loads((complex_domain.OUT/'result.json').read_text())['rows']
    constants=json.loads((density.ionic.OUT/'constants.json').read_text())['native_binary64'];B=I(constants['alpha'])/iv.pi
    p,Q,Gamma=sp.symbols('p Q Gamma',positive=True);assert sp.expand(p*p+(1+p*p-Q*Q)-1-(2*p*p-Q*Q))==0
    control=z-iv.atan2(z,I(1));assert 0<low(control) and high(control)<=high(z**3/3)
    rows=[];maximum=[F(0)]*7;minimum_L=None
    for i,root in enumerate(roots):
        assert root['cell']==complex_rows[i]['cell']==int(state['cells'][i])
        beta=I(state['beta'][i]);scale=I(state['scale'][i]);Sref=I(state['Sref'][i]);eta=cusp.interval(root['root']);assert high(scale*z)<=low(iv.sqrt(beta))
        L=iv.sqrt(beta)/(iv.sqrt(1+4*beta)*Sref*(1+iv.exp(2-I(low(eta)))));minimum_L=low(L) if minimum_L is None else min(minimum_L,low(L))
        m=list(map(I,complex_rows[i]['response_majorants']));C=[1/L,m[1]/L**2,m[2]/L**2+2*m[1]**2/L**3,m[3]/L**2,
            m[4]/L**2+2*m[1]*m[3]/L**3,m[5]/L**2+2*m[3]**2/L**3]
        omega=scale*scale*z**3/(3*B*Sref);raw=[omega*v for v in C]
        coeff=[cusp.interval(root['derivatives'][name]) for name in density.FIELDS]
        counts=[I(x)/I(a) for x,a in zip(ions['X'][i],ions['A'])];Z2=sum(x*I(k)**2 for x,k in zip(counts,ions['Z']))/sum(counts);pref=2*B*Z2*scale/beta
        errors=[pref*sum(abs(v)*bound for v,bound in zip(row,raw)) for row in neutral.matrix(coeff)]
        approx=[-pref*z,-pref*z,*[I(0)]*5];enclosures=[];points=[];scores=[]
        for a,error in zip(approx,errors):
            text=cusp.interval_text(a+cusp.symmetric(high(error)));lo,hi=cusp.endpoints(text);point=float((lo+hi)/2);score=max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))
            enclosures.append(text);points.append(point);scores.append(score)
        maximum=[max(old,new) for old,new in zip(maximum,scores)]
        rows.append(dict(cell=root['cell'],response_lower=str(low(L)),raw_error_upper=list(map(lambda v:str(high(v)),raw)),
            field_enclosures=enclosures,field_approximation=points,field_error_upper=list(map(str,scores)),passed=max(scores)<F(plan['field_budget'])))
    save('result.json',dict(classification='Proven',passed=all(r['passed'] for r in rows),cells=len(rows),fields=FIELDS,rows=rows,
        maximum_errors=dict(zip(FIELDS,map(str,maximum))),minimum_response_lower=str(minimum_L),display_only=dict(zip(FIELDS,map(float,maximum))),
        symbolic_passed=True,constant_response_control=dict(passed=True,actual_error=cusp.interval_text(control),cubic_bound=cusp.interval_text(z**3/3)),
        finite_outer_integral_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():g.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['symbolic_passed'] and r['constant_response_control']['passed'] and r['cells']==3206
    for row in r['rows']:
        assert row['passed'] and F(row['response_lower'])>0
        for text,p,bound in zip(row['field_enclosures'],row['field_approximation'],row['field_error_upper'],strict=True):
            lo,hi=cusp.endpoints(text);error=max(abs(lo-F.from_float(p)),abs(hi-F.from_float(p)));assert error==F(bound)<F(plan['field_budget'])
    print('PASS infrared self endpoint and seven neutral thermodynamic field enclosures; finite outer interval/full EOS remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
