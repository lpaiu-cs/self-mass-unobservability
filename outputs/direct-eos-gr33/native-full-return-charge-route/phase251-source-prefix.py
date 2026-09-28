from pathlib import Path
import json,numpy as np
root=Path('native-full-charge251-work');rows=[]
for n in [64,128]:
 a=dict(np.load(f'native-full-captured244-work/gr/source-{n}.npz'));b=dict(np.load(f'native-complete-radau224-work/gr/source-{n}.npz'))
 count=len(b['t']);geometry=len(b['geometry_times']);diff=[]
 for k,v in b.items():
  w=a[k]
  if k=='t' or k in ['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','nonrest_stress_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg','metric_stress_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg']:w=w[:count]
  elif k.startswith('state_coeff_'):w=w[:,:count-1]
  elif k=='geometry_times':w=w[:geometry]
  elif k.startswith('geometry_coeff_'):w=w[:,:geometry-1]
  if not np.array_equal(w,v):
   ids=np.argwhere(w!=v);diff.append(dict(key=k,shape=v.shape,different=len(ids),first=ids[:4].tolist(),last=ids[-4:].tolist(),max_relative=float(np.max(abs(w-v))/max(np.max(abs(v)),1e-290))))
 rows.append(dict(clock=n,differences=diff,old_T=float(b['t'][-1]),new_same_T=float(a['t'][count-1]),old_geometry_times=b['geometry_times'].tolist(),new_geometry_times=a['geometry_times'].tolist()))
(root/'source-prefix-difference.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows))
