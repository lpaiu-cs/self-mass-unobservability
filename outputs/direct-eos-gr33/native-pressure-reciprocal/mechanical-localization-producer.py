from pathlib import Path
import json,time
import numpy as np
start=time.monotonic();root=Path('native-pressure-reciprocal175-work');LD=np.longdouble;AMP=LD(1e-26);rows=[]
for n in [64,128]:
    pairs=[];totals=[]
    for sweep in [0,1]:
        d=np.load(root/f'sweep-{sweep}/material/steps-{n}-reference-128.npz')
        p=np.load(root/f'sweep-{sweep}/photons/steps-{n}-reference-128.npz')
        ids=[np.argmin(abs(d['t']-t)) for t in p['t']]
        H=d['history_scaled'][ids,3].astype(LD)*AMP;collision=p['collision_transfer'][:,:,1].astype(LD)
        pairs.append(H-collision);totals.append(H)
    norm=lambda a:np.max(np.sum(abs(a),axis=-1))
    old,new=pairs;scale=norm(new);delta=norm(new-old)
    k,i=np.unravel_index(np.argmax(abs(new-old)),new.shape)
    rows.append(dict(clock=n,new_mechanical_H_L1=float(scale),old_mechanical_H_L1=float(norm(old)),
        change_L1=float(delta),change_over_new=float(delta/scale),total_H_L1=float(norm(totals[1])),
        subtraction_condition=float(norm(totals[1])/scale),float64_scale_over_change=float(np.finfo(float).eps*norm(totals[1])/delta),
        largest_cell=dict(k=int(k),cell=int(i),old=float(old[k,i]),new=float(new[k,i])),
        M_history_L1=np.sum(abs(new),axis=1).astype(float).tolist(),difference_history_L1=np.sum(abs(new-old),axis=1).astype(float).tolist()))
result=dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,no_evolution_steps=True)
(root/'mechanical-localization.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
