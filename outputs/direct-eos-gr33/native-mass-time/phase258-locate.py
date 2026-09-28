from pathlib import Path
import json,resource,time
import numpy as np

resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
root=Path('native-material-exterior257-work');a=dict(np.load(root/'mass-ledger-64.npz'));b=dict(np.load(root/'mass-ledger-128.npz'))
assert np.array_equal(a['t'],b['t'])
rows=[]
for i in range(17):
    t=b['t'][-1]*i/16;j=np.argmin(abs(b['t']-t));assert abs(b['t'][j]-t)<1e-18
    diff=a['terms'][j]-b['terms'][j]
    rows.append(dict(canonical=i,time=float(t),fine_mass_energy=float(b['total_erg'][j]),
        delta_mass_energy=float(sum(diff)),delta_gas_nonrest=float(diff[1]),
        relative=float(abs(sum(diff))/max(abs(b['total_erg'][j]),np.longdouble('1e-290')))))
print(json.dumps(dict(canonical_mass=rows)))
for n in [64,128]:
    p=np.load(f'native-material-accuracy257-work/sweep-1/photons/return-{n}.npz')
    keys=[k for k in p.files if any(v in k for v in ['material','ledger','units','scale','native','collision','mechanical'])]
    print(json.dumps(dict(clock=n,keys=keys,shapes={k:str(p[k].shape) for k in keys},times=p['t'].astype(float).tolist())))
