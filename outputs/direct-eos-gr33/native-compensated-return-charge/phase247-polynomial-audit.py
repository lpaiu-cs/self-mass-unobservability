"""Independently compare composed coefficients with the two-field readout."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.interpolate import PPoly
sys.path.insert(0,'verification')
import read_dense_returned_source as source
check=len(sys.argv)>1 and sys.argv[1]=='check'
folder=Path('native-dense-returned246-work')/('check' if check else 'full')
out=Path('native-compensated-charge247-work')/('check' if check else 'full')
rows=[];files=[Path(__file__),Path(source.__file__)]
for n in [64,128]:
    path=folder/f'gr/source-{n}.npz';d=dict(np.load(path));files.append(path)
    knots,co=source.coefficients(d);times=np.concatenate([a+(b-a)*np.array([1/6,1/2,5/6]) for a,b in zip(knots[:-1],knots[1:])])
    gt=d['geometry_times'];ids=np.clip(np.searchsorted(gt,times,side='right')-1,0,len(gt)-2);errors={}
    for k in source.KEYS:
        expected=source.base.prior.prior.poly(d['t'],d['state_coeff_'+k])(times).astype(source.LD)
        if k not in source.KEYS[-2:]:
            for name in ['u','lambda']:
                a,b=d['geometry_'+name+'_'+k];geometry=a[ids]+(times-gt[ids])[:,None]*b[ids]
                expected+=geometry*PPoly(d[name+'_coeff'],d['metric_times'])(times)/source.AMP
        actual=PPoly(np.asarray(co[k][::-1],float),knots)(times)
        measure=lambda v:np.sum(abs(v),axis=-1) if v.ndim>1 else abs(v)
        errors[k]=float(np.max(measure(actual-expected))/max(np.max(measure(expected)),source.LD('1e-290')))
    rows.append(dict(clock=n,independent_interior_evaluations=len(times),relative=errors))
passed=max(v for r in rows for v in r['relative'].values())<1e-12
result=dict(classification='Counterexample candidate',passed=passed,rows=rows,
    scope='Algebra of the declared returned-metric Hermite source on all composed intervals; not a uniform temporal interpolation or physical error certificate.',
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
(out/'polynomial-audit.json').write_text(json.dumps(result,indent=2)+'\n');assert passed,result
print(json.dumps(result))
