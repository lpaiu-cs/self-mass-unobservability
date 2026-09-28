from pathlib import Path
import json,resource
import numpy as np
import finish_returned_material_accuracy as producer

resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
AMP=producer.base.actual.AMP;LD=np.longdouble
allrows={};bycell={}
for n in [64,128]:
    p=np.load(f'native-material-accuracy257-work/sweep-1/photons/return-{n}.npz')
    state=p['material_history'][:,:,0]*p['material_energy_units']
    native=np.sum(p['joint_native_rates_scaled'][:,:,0],axis=1,dtype=LD)
    collision=np.sum(p['joint_collision_rates_scaled'][:,:,0],axis=1,dtype=LD)
    weights=p['joint_stage_weights'].astype(LD)
    floor=np.sum(p['material_floor_discard_history_scaled'][:,:,0],axis=1,dtype=LD)*AMP
    final=state[-1];bycell[n]=final
    nr=np.cumsum(native*weights,dtype=LD)*AMP;cr=np.cumsum(collision*weights,dtype=LD)*AMP
    ids=np.array([np.count_nonzero(p['joint_stage_times']<=t+1e-18)-1 for t in p['t']])
    nr=np.r_[0,nr[ids[1:]]];cr=np.r_[0,cr[ids[1:]]]
    observed=np.sum(state,axis=1,dtype=LD);relative=np.max(abs(observed-(nr+cr-floor)))/max(np.max(abs(observed)),LD('1e-290'))
    assert relative<1e-8,relative
    allrows[n]=np.array([observed,nr,cr,-floor]).T
    print(json.dumps(dict(clock=n,ledger_relative=float(relative),endpoint=[float(v) for v in allrows[n][-1]])))
diff=allrows[64]-allrows[128];delta=bycell[64]-bycell[128];ids=np.argsort(abs(delta))[::-1][:12]
print(json.dumps(dict(AMP=float(AMP),columns=['state','native','collision','minus_floor'],
    differences_by_canonical=diff.astype(float).tolist(),top_cells=[dict(cell=int(i),radius=float(p['radius_E'][i]),delta_erg=float(delta[i]),fraction=float(delta[i]/delta.sum())) for i in ids])))
out=Path('phase258-energy-owner.json');out.write_text(json.dumps(dict(classification='Counterexample candidate',columns=['state','native','collision','minus_floor'],
    difference_by_canonical=diff.astype(float).tolist(),endpoint_gas_difference_erg=float(delta.sum()),
    top_cells=[dict(cell=int(i),delta_erg=float(delta[i])) for i in ids],physical_steps=0),indent=2)+'\n')
