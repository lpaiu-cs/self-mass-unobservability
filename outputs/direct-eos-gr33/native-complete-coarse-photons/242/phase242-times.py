from pathlib import Path
import json,numpy as np
q=np.load('native-versioned-photon241-work/rejected-original-64.npz')
z=np.load('native-true-momentum238-work/sweep-1/photons/complete-64.npz')
paths=['native-flux-precision202-work/last-accepted-64.npz','native-momentum-continuation210-work/last-accepted-64.npz']
rows={}
for path in paths:
    p=np.load(path)
    rows[path]={k:dict(value=str(p[k]),dtype=str(p[k].dtype)) for k in ['next_time','next_step','macro_index','sub_index'] if k in p}
    rows[path]['edges_tail']=[str(v) for v in p['actual_edges'][-3:]]
rows['saved']={k:dict(value=str(q[k]),dtype=str(q[k].dtype)) for k in ['time','step_size','stage_times']}
rows['step_weights']=dict(value=str(z['joint_stage_weights'][230:232]),dtype=str(z['joint_stage_weights'].dtype))
rows['stage_times_exact']=[dict(value=str(t),ratio=[str(v) for v in t.as_integer_ratio()]) for t in q['stage_times']]
rows['edges_exact']=[dict(value=str(t),ratio=[str(v) for v in t.as_integer_ratio()]) for t in z['actual_step_edges'][115:117]]
defect=np.load('native-stored-time-photon242-work/actual-defect.npz')
original=np.load(paths[0])['g'];restored=defect['initial'][-original.size:].reshape(original.shape)
rows['initial_g261']={k:[str(v) for v in value[261]] for k,value in [('original',original),('restored',restored),('difference',restored-original)]}
rows['B_preimages_261']=[str(v) for v in q['gas'][:,261,2]]
print(json.dumps(rows,indent=2))
Path('native-stored-time-photon242-work/actual-time-types.json').write_text(json.dumps(rows,indent=2)+'\n')
