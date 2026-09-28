"""Separate scalar-Q numerical controls from wide-range vector quadrature."""
from fractions import Fraction as F
import gzip,json,sys
import numpy as np
from mpmath import iv
import gr_response_product_defined as defined

g=defined.original;ROOT=g.ROOT;OUT=g.OUT.parent/'gr-response-product-reference';I=g.I;low=g.low;high=g.high;cusp=g.cusp;sha=g.sha

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    defined.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_product_reference.py',defined.OUT/'manifest.json',defined.OUT/'result.json',
        cusp.THERMO/'states.npz',cusp.OUT/'inputs.npz',cusp.ROOTS/'result.json',cusp.OUT.parent/'gr-smooth-response-defined/manifest.json',
        cusp.OUT.parent/'gr-smooth-response-defined/complete-response.jsonl.gz',ROOT/'verification/gr_polarization_thermodynamics.py']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='7bb5d3b',bits=256,scalar_reference_tolerance=2e-13,
        scalar_absolute_gate=2e-11,bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        target='The twelve-Q vector reference was close to its finite gate at one state. Re-evaluate each Q separately so occupied-region rescaling is selected for that Q. Also check necessary overlap of the first three certified response components after an analytic reference-center shift.',
        scope='New independent audit. Preserve the original vector control and gate, the product intervals, and all earlier numerical thermodynamic verdicts. A necessary interval-overlap check is not a substitute for either certificate.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def run():
    p=bindings();iv.prec=p['bits'];state=dict(np.load(cusp.THERMO/'states.npz'));data=dict(np.load(cusp.OUT/'inputs.npz'));roots=json.loads((cusp.ROOTS/'result.json').read_text())['rows']
    product=json.loads((defined.OUT/'result.json').read_text())['results'];positions={r['position'] for r in product};certified={}
    with gzip.open(cusp.OUT.parent/'gr-smooth-response-defined/complete-response.jsonl.gz','rt') as stream:
        for line in stream:
            row=json.loads(line)
            if row['position'] in positions:certified[(row['position'],row['z_index'])]=row
    checks=[];maximum=0.;batch_maximum=0.
    for result in product:
        i=result['position'];eta=float(state['eta'][i]);beta=float(state['beta'][i]);Sref=float(state['Sref'][i]);Qs=data['Q'][i]
        batch=g.rule.highq.neutral.thermo.response(np.full(12,eta),np.full(12,beta),Qs,np.full(12,Sref),p['scalar_reference_tolerance'])
        lo,hi=cusp.endpoints(roots[i]['root']);delta=max(abs(F.from_float(eta)-lo),abs(F.from_float(eta)-hi));factor=I(delta)*iv.exp(I(delta))
        for row in result['rows']:
            j=row['z_index'];reference=g.rule.highq.neutral.thermo.response(np.array([eta]),np.array([beta]),np.array([Qs[j]]),np.array([Sref]),p['scalar_reference_tolerance']).reshape(6)
            score=float(np.max(np.abs(reference-np.array(row['approximation']))));batch_score=float(np.max(np.abs(reference-batch[:,j])));maximum=max(maximum,score);batch_maximum=max(batch_maximum,batch_score)
            old=certified[(i,j)];value=cusp.interval(old['enclosures'][0]);shift=high(factor*I(high(value)))
            reference_ranges=[value*iv.exp(cusp.symmetric(delta)),cusp.interval(old['enclosures'][1])+cusp.symmetric(shift),cusp.interval(old['enclosures'][2])+cusp.symmetric(shift)]
            overlap=[]
            for text,expected in zip(row['enclosures'][:3],reference_ranges):
                a,b=cusp.endpoints(text);overlap.append(max(a,low(expected))<=min(b,high(expected)))
            checks.append(dict(position=i,z_index=j,scalar_reference=reference.tolist(),maximum_absolute_difference=score,batch_to_scalar_difference=batch_score,
                first_three_interval_overlaps=overlap,passed=score<p['scalar_absolute_gate'] and all(overlap)))
    save('result.json',dict(classification='Counterexample candidate',passed=all(x['passed'] for x in checks),checks=checks,
        maximum_scalar_difference=maximum,maximum_batch_to_scalar_difference=batch_maximum,
        shift_identity='For the positive original kernel, f_eta, f_etaeta and f_etaetaeta have absolute value <=f. If the stored eta center lies within distance delta of the certified root interval, f along the joining segment is <=exp(delta) times its root value. The value has the sharper multiplicative enclosure exp([-delta,delta])*R_root; the first and second eta derivatives move by at most delta*exp(delta)*R_root_upper. This supplies the three necessary interval-overlap checks.',
        reference_boundary='A vector adaptive quadrature uses one global occupied-region scale for all its entries. Scalar calls select that scale independently at each Q. The difference is a finite numerical reference diagnostic and does not alter the proven response intervals or retrospectively certify the old thermodynamic table.'))
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and len(r['checks'])==36
    print('PASS scalar-Q product reference and sharp eta-shift consistency:',r['maximum_scalar_difference'],r['maximum_batch_to_scalar_difference'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
