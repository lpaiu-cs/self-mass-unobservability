"""Independent finite complex-ray checks of the analytic response majorants."""
from fractions import Fraction as F
import json,sys
import numpy as np
import mpmath as mp
import sympy as sp
import gr_response_complex_runner as certificate

g=certificate.g;cusp=g.cusp;ROOT=g.ROOT;OUT=g.OUT.parent/'gr-response-complex-controls'
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def M(x):
    q=F(x);return mp.mpf(q.numerator)/q.denominator


def prepare():
    certificate.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_complex_controls.py',certificate.OUT/'manifest.json',certificate.OUT/'result.json',cusp.OUT/'inputs.npz']
    save('plan.json',dict(classification='Counterexample candidate',digits=64,positions=[0,1603,3205],z_indices=[0,5,8],relative_imaginary_Q='1/128',
        bindings={p.relative_to(ROOT).as_posix():g.sha(p) for p in files},
        method='Integrate the six complex response fields along p=Q*u at the outer analytic circle point Q=q0*(1+i/128). Use real p_a=a*u, a=q0*sqrt(1-2r), to resolve the thermal support. Subtract B(1) from the logarithmic coefficient on u in [0,2], adding back its exact integral. Compare absolute values with the independent uniform analytic majorants.',
        boundary='Finite controls do not replace the analytic proof and do not certify an outer integral. These circle points are outside the smaller certified zero-free neighborhoods; no nonvanishing claim is made there.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def run():
    plan=bindings();mp.mp.dps=plan['digits'];data=dict(np.load(cusp.OUT/'inputs.npz'));rows=json.loads((certificate.OUT/'result.json').read_text())['rows'];checks=[]
    r=M(plan['relative_imaginary_Q']);assert sp.simplify(sp.I**2+1)==0
    for i in plan['positions']:
        for j in plan['z_indices']:
            eta=sum(map(M,cusp.endpoints(str(data['root_intervals'][i]))))/2;beta=mp.mpf(float(data['beta'][i]));Sref=mp.mpf(float(data['Sref'][i]))
            q0=mp.mpf(float(data['Q'][i,j]));Q=q0*(1+1j*r);a=q0*mp.sqrt(1-2*r);_,B0=cusp.kernel(mp.mpf(1),eta,beta,Q,Sref,mp);values=[]
            support=[mp.sqrt(beta)*v for v in map(mp.mpf,['.1','1','3','8','15'])]
            for k in range(6):
                def finite(p):
                    u=p/a;smooth,B=cusp.kernel(u,eta,beta,Q,Sref,mp)
                    return (smooth[k]-(B[k]-B0[k])*mp.log(abs(1-u)))/a if u!=1 else smooth[k]/a
                def infinite(p):
                    u=p/a;smooth,B=cusp.kernel(u,eta,beta,Q,Sref,mp);return (smooth[k]-B[k]*mp.log(u-1))/a
                left=sorted(set([mp.mpf(0),a]+[p for p in support if 0<p<a]));right=sorted(set([a,2*a]+[p for p in support if a<p<2*a]))
                end=sorted(set([2*a]+[p for p in support if p>2*a]))+[mp.inf]
                values.append(mp.quad(finite,left)+mp.quad(finite,right)+2*B0[k]+mp.quad(infinite,end))
            bounds=list(map(M,rows[i]['response_majorants']));ratios=[abs(v)/b for v,b in zip(values,bounds)];passed=max(ratios)<1
            checks.append(dict(cell=rows[i]['cell'],z_index=j,passed=passed,reference=[dict(real=str(v.real),imag=str(v.imag)) for v in values],maximum_fraction_of_majorant=str(max(ratios))))
            print('COMPLEX RAY CONTROL',rows[i]['cell'],j,passed,str(max(ratios)),flush=True)
    save('result.json',dict(classification='Counterexample candidate',passed=all(c['passed'] for c in checks),checks=checks,
        negative_control='For the admissible positive real-axis constant response S=1 and B*Sref=1, the denominator Q^2+1 vanishes at Q=i. Thus real-axis positivity alone is not a global complex zero-free claim. The constructive local disk condition is necessary to justify the later quotient bounds.'))
    assert len(checks)==9 and all(c['passed'] for c in checks)
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():g.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and len(r['checks'])==9 and all(c['passed'] for c in r['checks'])
    print('PASS 54 independent complex-ray response controls and explicit out-of-domain pole control',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
