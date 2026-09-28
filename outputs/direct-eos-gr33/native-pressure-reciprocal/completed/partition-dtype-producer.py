from pathlib import Path
import numpy as np,json,time
w=Path('native-pressure-reciprocal175-work');start=time.monotonic()
with np.load(w/'sweep-0/material/steps-128-reference-128.npz') as d,np.load(w/'sweep-0/photons/steps-128-reference-128.npz') as p:
 ids=[int(np.argmin(abs(d['t']-t))) for t in p['t']]
 H=d['history_scaled'][ids,3].astype(np.longdouble);C=p['collision_transfer'][:,:,1]
 precise=H-C.astype(np.longdouble)/np.longdouble(1e-26)
 rounded=H-C/1e-26
 norm=max(np.max(abs(precise)),np.longdouble('1e-290'))
 r=dict(classification='Counterexample candidate',stored_collision_dtype=str(C.dtype),material_dtype=str(d['history_scaled'].dtype),check_conversion_relative=float(np.max(abs(precise-rounded))/norm),nominal_old_failure=7.131209897139719e-8,seconds=time.monotonic()-start)
 print(json.dumps(r));(w/'partition-dtype-evidence.json').write_text(json.dumps(r,indent=2)+'\n')