"""Audit cusp quadrature and all product jets using independent cold convolution."""
from fractions import Fraction as F
from types import FunctionType
import json,sys
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad_vec
import gr_response_product_reference as previous

g=previous.g;ROOT=g.ROOT;OUT=g.OUT.parent/'gr-response-product-convolution';sha=g.sha
thermo=g.rule.highq.neutral.thermo


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    previous.bindings();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_product_convolution.py',ROOT/'verification/gr_finite_wavenumber_response.py',
        previous.OUT/'plan.json',previous.OUT/'manifest.json',previous.OUT/'result.json',
        previous.defined.OUT/'manifest.json',previous.defined.OUT/'result.json',g.cusp.THERMO/'states.npz',g.cusp.OUT/'inputs.npz',
        ROOT/'verification/gr_polarization_thermodynamics.py']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d13ca7c',digits=60,absolute_gate='2e-11',
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        target='Independently integrate the source cold response times minus the t derivative of each fixed-p occupation jet. Replay the four failed scalar calls with quad_vec full_output to expose convergence status. Preserve every earlier result and threshold.',
        scope='High-precision finite numerical controls of all 36 points and six components. These do not create a new rigorous quadrature certificate or repair the old thermodynamic table.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def weights(t,eta):
    q=1/(1+mp.exp(t-eta));v=1-q;w=q*v;a=1-2*q;b=1-6*q+6*q*q
    return [w,w*a,w*b,(t*a-1)*w,(t*b-a)*w,(t*t*b-3*t*a+1)*w]


def symbolic():
    t,e=sp.symbols('t e',real=True);q=1/(1+sp.exp(t-e));v=1-q;w=q*v;a=1-2*q;b=1-6*q+6*q*q
    jets=[q,w,w*a,t*w,t*w*a,(t*t*a-t)*w]
    expected=[w,w*a,w*b,(t*a-1)*w,(t*b-a)*w,(t*t*b-3*t*a+1)*w]
    for f,target in zip(jets,expected):assert sp.simplify(-sp.diff(f,t)-target)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        identity='For each fixed-p occupation jet f_j, integration by parts gives integral K_prime(p,Q) f_j dp = integral K(sqrt(beta*t*(2+beta*t)),Q) [-d f_j/dt] dt. K(0,Q)=0 and all thermal boundary terms vanish. This is a representation at the reference beta, not differentiation of the transformed integration coordinate.',
        weights='w=f(1-f), a=1-2f, b=1-6f+6f^2: [w,wa,wb,(ta-1)w,(tb-a)w,(t^2b-3ta+1)w]. Exact symbolic check of all six.'))


def run():
    p=bindings();symbolic();mp.mp.dps=p['digits'];state=dict(np.load(g.cusp.THERMO/'states.npz'));data=dict(np.load(g.cusp.OUT/'inputs.npz'))
    old=json.loads((previous.OUT/'result.json').read_text());products=json.loads((previous.defined.OUT/'result.json').read_text())['results']
    diagnostics=[]
    for failed in [r for r in old['checks'] if not r['passed']]:
        i,j=failed['position'],failed['z_index'];records=[]
        def audited(*args,**kwargs):
            value,error,info=quad_vec(*args,**kwargs,full_output=True)
            records.append(dict(estimated_error=float(error),success=bool(info.success),status=int(info.status),message=info.message,neval=int(info.neval),intervals=len(info.intervals)))
            return value,error
        fn=FunctionType(thermo.response.__code__,dict(thermo.__dict__,quad_vec=audited))
        fn(*[np.array([float(state[k][i])]) for k in ['eta','beta']],np.array([float(data['Q'][i,j])]),np.array([float(state['Sref'][i])]),2e-13)
        diagnostics.append(dict(position=i,z_index=j,parts=records))
    save('adaptive-diagnostics.json',dict(classification='Counterexample candidate',checks=diagnostics))
    cold=thermo.previous.BASE['cold'];checks=[]
    for result in products:
        i=result['position'];eta,beta,Sref=[mp.mpf(float(state[k][i])) for k in ['eta','beta','Sref']]
        for row in result['rows']:
            j=row['z_index'];Q=mp.mpf(float(data['Q'][i,j]));tq=Q*Q/(beta*(mp.sqrt(1+Q*Q)+1));top=max(mp.mpf(2),eta)
            points=sorted(set([mp.mpf(0),mp.mpf(1),top,top+8,top+40]+([tq] if 0<tq<top+80 else [])))+[mp.inf]
            values=[mp.quad(lambda t: cold(mp.sqrt(beta*t*(2+beta*t)),Q)*weights(t,eta)[k]/Sref,points) for k in range(6)]
            score=max(abs(value-mp.mpf(point)) for value,point in zip(values,row['approximation']))
            # The product midpoint has a certified 1e-15 bound; the independent finite reference gate is unchanged.
            checks.append(dict(position=i,z_index=j,reference=[str(x) for x in values],absolute_difference=str(score),passed=score<mp.mpf(p['absolute_gate'])))
            save('progress.json',dict(checks=checks));print('COLD CONVOLUTION',i,j,float(score),flush=True)
    save('result.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in checks),components=6*len(checks),checks=checks,
        maximum_absolute_difference=str(max(F(r['absolute_difference']) for r in checks))))
    save('manifest.json',dict(sha256={x.relative_to(ROOT).as_posix():sha(x) for x in OUT.iterdir() if x.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['components']==216
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    old=json.loads((previous.OUT/'result.json').read_text());assert not old['passed'] and sum(not x['passed'] for x in old['checks'])==4
    assert all(all(x['first_three_interval_overlaps']) for x in old['checks'])
    print('PASS 216 independent cold-convolution components; four original scalar failures retained:',float(F(r['maximum_absolute_difference'])),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
