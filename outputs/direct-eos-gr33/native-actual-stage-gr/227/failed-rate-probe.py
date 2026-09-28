"""Check the declared piecewise-affine metric's exact integral, without a solve."""
from pathlib import Path
import hashlib,json
import numpy as np

out=Path('native-stage-metric227-work');out.mkdir(exist_ok=True)
metric=Path('native-dense-return225-work/metric/metric-128-g8.npz')
d=np.load(metric);times=d['t'];rows=[];paths=[metric]
for n in [64,128]:
    path=Path(f'native-interval-return226-work/sweep-1/photons/return-{n}.npz');paths.append(path)
    p=np.load(path);edges=p['actual_step_edges'];h=np.diff(edges);stage=p['joint_stage_times']
    assert np.array_equal(p['joint_stage_weights'].reshape(-1,2),h[:,None]*[.75,.25])
    j=np.clip(np.searchsorted(times,stage,side='left')-1,0,len(times)-2)
    result={}
    for key in ['delta_u','delta_lambda']:
        value=d[key].astype(np.longdouble);rates=np.diff(value,axis=0)/np.diff(times)[:,None]
        samples=rates[j].reshape(-1,2,value.shape[1]);numerical=h[:,None]*np.sum(samples*np.array([.75,.25])[None,:,None],axis=1)
        edge=np.array([np.interp(edges,times,v) for v in value.T]).T.astype(np.longdouble)
        exact=np.diff(edge,axis=0);defect=numerical-exact;scale=max(np.max(np.sum(abs(exact),axis=-1)),np.longdouble('1e-290'))
        result[key]=dict(maximum_step_relative=float(np.max(np.sum(abs(defect),axis=-1))/scale),
            whole_integral_relative=float(np.sum(abs(np.sum(defect,axis=0)))/max(np.sum(abs(edge[-1]-edge[0])),np.longdouble('1e-290'))),
            error_L1_by_step=np.sum(abs(defect),axis=-1).astype(float).tolist())
    rows.append(dict(clock=n,actual_steps=len(h),controls=result))
r=dict(classification='Counterexample candidate',rows=rows,
    scope='Exact primitive versus original Radau quadrature of the DECLARED piecewise-affine GR derivative on actual physical clocks. Demonstrates input quadrature aliasing only, not its full coupled error share.',
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},physical_steps=0,final_charge_conclusion='unadjudicated')
(out/'rate-alias.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
