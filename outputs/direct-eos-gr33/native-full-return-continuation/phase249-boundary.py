"""Preserve the endpoint derivative change that requires a short rollback."""
from pathlib import Path
import hashlib,json,numpy as np
root=Path('native-full-return249-work');root.mkdir(exist_ok=True)
paths=[Path(p)/'metric/metric-128-g8.npz' for p in ['native-stage-metric227-work','native-complete-return236-work']]
a,b=[np.load(p) for p in paths];ids=np.array([np.argmin(abs(b['t']-t)) for t in a['t']])
assert np.max(abs(b['t'][ids]-a['t']))<1e-18;rows=[]
for k in a.files:
    if (k.startswith('delta_') or k=='actual_delta_lambda_rate') and a[k].shape[0]==len(a['t']):
        v=b[k][ids];norm=max(np.max(abs(a[k])),1e-290)
        rows.append(dict(key=k,prefix=float(np.max(abs(v[:-1]-a[k][:-1]))/norm),old_final=float(np.max(abs(v[-1]-a[k][-1]))/norm)))
assert next(r for r in rows if r['key']=='actual_delta_lambda_rate')['old_final']>1e-3
r=dict(classification='Counterexample candidate',rows=rows,
    conclusion='The prior final source-polynomial derivative becomes an interior right-sided derivative after extension. Do not silently retain that formerly terminal stage. Restart one canonical interval earlier and re-evolve the overlap plus the last interval under the extended input.',
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*paths]},final_charge_conclusion='unadjudicated')
(root/'boundary-extension.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
