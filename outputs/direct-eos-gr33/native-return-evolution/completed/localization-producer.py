from pathlib import Path
import hashlib,json,time
import numpy as np
from evolve_same_solution_gr_return import AMP
start=time.monotonic();root=Path('native-return-evolution188-work');LD=np.longdouble
p=[dict(np.load(root/f'sweep-1/photons/return-{n}.npz')) for n in [64,128]]
g=[v['joint_native_rates_scaled'] for v in p];w=[v['joint_stage_weights'].astype(LD) for v in p]
integral=[np.sum(q[:,:,2]*h[:,None],axis=0,dtype=LD) for q,h in zip(g,w)]
floor=[v['material_floor_discard_scaled'][:,2] for v in p]
final=[v['conserved_material_history'][-1,0]/AMP for v in p]
delta=final[0]-final[1];native=integral[0]-integral[1];removed=floor[0]-floor[1]
norm=np.sum(abs(delta),dtype=LD);ids=np.argsort(abs(delta))[-12:][::-1]
rows=[]
for i in ids:
    rows.append(dict(cell=int(i),radius=float(p[0]['radius_E'][i]),difference=float(delta[i]),
        native_difference=float(native[i]),floor_difference=float(removed[i]),fraction=float(abs(delta[i])/norm)))
common=[]
for k,t in enumerate(p[0]['joint_stage_times']):
    j=int(np.argmin(abs(p[1]['joint_stage_times']-t)))
    if abs(p[1]['joint_stage_times'][j]-t)<1e-18:
        rate0=g[0][k,:,2];rate1=g[1][j,:,2]
        common.append(dict(time=float(t),rate_relative=float(np.sum(abs(rate0-rate1))/np.maximum(np.sum(abs(rate1)),LD('1e-290')))))
result=dict(classification='Counterexample candidate',scope='Saved failed return pair; zero new physical steps.',
    baryon_difference_L1_scaled=float(norm),native_difference_over_final_difference=float(np.sum(abs(native))/norm),
    floor_difference_over_final_difference=float(np.sum(abs(removed))/norm),
    ledger_relative=float(np.sum(abs(delta-native+removed))/norm),largest_cells=rows,common_stage_baryon_rate=common,
    clocks=[dict(times=v['joint_stage_times'].tolist(),B_rate_L1=np.sum(abs(q[:,:,2]),axis=1).astype(float).tolist()) for v,q in zip(p,g)],
    seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    original_time_gate=.02,original_baryon_relative=.05211566092658558,final_charge_conclusion='unadjudicated')
(root/'baryon-localization.json').write_text(json.dumps(result,indent=2)+'\n')
assert result['seconds']<15
assert result['ledger_relative']<1e-8
print(json.dumps(result))
