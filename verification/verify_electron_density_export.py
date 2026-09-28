"""Exact-rational containment audit of serialized neutral-root certificates."""
from fractions import Fraction as F
import json, sys
import numpy as np
import gr_electron_density_certificate as original

g=original.g;OUT=g.OUT/'gr-electron-density-export'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def ends(text):
    a,b=map(F,text[1:-1].split(','));assert a<=b;return a,b


def run():
    original.verify();assert not OUT.exists();OUT.mkdir()
    files=[g.ROOT/'verification/verify_electron_density_export.py',original.OUT/'manifest.json',
        original.OUT/'result.json',original.OUT/'candidates.npz',original.OUT/'plan.json']
    save('plan.json',dict(classification='Proven',bindings={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in files},
        issue='The mathematical raw interval bounds remain valid, but outward decimal rendering may make the exported interval wider than its raw width/error scalars. Audit the encoded interval itself and retain any original scalar inconsistency.',
        repair='Keep every exported root and derivative enclosure unchanged. Expand only the corresponding scalar width, center-distance and finite-comparison budgets to contain those exact rational decimal endpoints. No root, EOS, arithmetic precision or science threshold change.'))
    p=json.loads((original.OUT/'plan.json').read_text());r=json.loads((original.OUT/'result.json').read_text());a=dict(np.load(original.OUT/'candidates.npz'))
    rows=[];failures=0;maximum_padding=F(0);widths=[];scores=[]
    for old in r['rows']:
        i=old['position'];center=F.from_float(float(a['eta_center'][i]));lo,hi=ends(old['root'])
        assert center-F(1,2**40)<=lo<=hi<=center+F(1,2**40)
        m=a['momentum'][:,i];en=m[0]/m[1];et=-m[3]/m[1]
        reference=[en,et,(m[0]-m[2]*en**2)/m[1],-en*(m[2]*et+m[4])/m[1],-(m[5]+2*m[4]*et+m[2]*et**2)/m[1]]
        encoded_score=F(0)
        for name,v in zip(p['fields'],reference):
            x,y=ends(old['derivatives'][name]);v=F.from_float(float(v))
            encoded_score=max(encoded_score,max(abs(x-v),abs(y-v))/max(F(1),abs(v)))
        values=[hi-lo,max(abs(lo-center),abs(hi-center)),encoded_score]
        keys=['root_width','root_error_from_center_upper','finite_momentum_comparison_upper'];updated=dict(old);padding={}
        for key,value in zip(keys,values):
            prior=F(old[key]);bound=max(prior,value);delta=bound-prior;failures+=int(delta>0);maximum_padding=max(maximum_padding,delta)
            updated[key]=str(bound);padding[key]=str(delta)
        updated['export_padding']=padding;updated['passed']=F(updated['finite_momentum_comparison_upper'])<F(p['finite_comparison_gate'])
        rows.append(updated);widths.append(F(updated['root_width']));scores.append(F(updated['finite_momentum_comparison_upper']))
    assert failures>0
    save('result.json',dict(classification='Proven',passed=all(x['passed'] for x in rows),cells=len(rows),rows=rows,
        raw_export_budget_consistent=False,raw_scalar_inconsistencies=failures,maximum_padding_exact=str(maximum_padding),
        maximum_root_width_exact=str(max(widths)),maximum_finite_comparison_exact=str(max(scores)),
        display_only=dict(maximum_root_width=float(max(widths)),maximum_finite_comparison=float(max(scores)),maximum_padding=float(maximum_padding)),
        scope='Exact rational coherence of the encoded certificate. Original enclosures and mathematical gates are preserved; full physical EOS remains open.'))
    save('manifest.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in OUT.iterdir() if x.is_file()}))
    verify()


def verify():
    original.verify()
    for key in ['plan.json','manifest.json']:
        obj=json.loads((OUT/key).read_text())
        for rel,digest in obj['bindings' if key=='plan.json' else 'sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());old=json.loads((original.OUT/'result.json').read_text())
    assert r['passed'] and r['cells']==3206 and r['raw_scalar_inconsistencies']>0
    for before,after in zip(old['rows'],r['rows']):
        assert before['root']==after['root'] and before['derivatives']==after['derivatives']
        lo,hi=ends(after['root']);assert hi-lo<=F(after['root_width'])
    print('PASS exact exported neutral-root budgets; original literal-budget discrepancies preserved',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
