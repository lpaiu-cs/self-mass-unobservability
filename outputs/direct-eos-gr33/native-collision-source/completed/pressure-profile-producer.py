"""Counterexample candidate: read-only attribution of the saved pressure gap."""
from pathlib import Path
import resource,time,json
import numpy as np
import integrate_native_collision_source as s

resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
s.original.inf.incident.native.deadline(45)
start=time.monotonic();s.initialize(False);m=s.original.c.Response(128)
s.radau.OLD=s.OLD;s.radau.profile()
end=m.t[-1]/16;pmap=m.point(1)['pressure_map'];rows=[]
for phase,folder in [(168,s.OLD),(169,s.OUT)]:
    pair=[np.load(folder/f'sweep-1/photons/pilot-{n}.npz') for n in [64,128]]
    terms=[pmap*p['delta_material']*m.volume[:,None] for p in pair]
    pressure=[p['moments'][-1,6] for p in pair]
    norm=max(np.sum(abs(pressure[1]),dtype=s.LD),s.LD('1e-290'))
    discrepancy=abs(pressure[0]-pressure[1]);den=np.sum(discrepancy,dtype=s.LD)
    ids=np.argsort(discrepancy)[-8:][::-1]
    mapping=max(float(np.sum(abs(q.sum(1)-p),dtype=s.LD)/norm) for q,p in zip(terms,pressure))
    assert mapping<1e-12,mapping
    rows.append(dict(phase=phase,relative=float(den/norm),pressure_L1_erg=float(norm),mapping_relative=mapping,
        terms_L1_erg=[float(v) for v in np.sum(abs(terms[1]),axis=0,dtype=s.LD)],
        term_time_errors_over_pressure=[float(v/norm) for v in np.sum(abs(terms[0]-terms[1]),axis=0,dtype=s.LD)],
        cells=[dict(cell=int(i),error_fraction=float(discrepancy[i]/den),coarse=float(pressure[0][i]),fine=float(pressure[1][i]),
                    coarse_terms=terms[0][i].astype(float).tolist(),fine_terms=terms[1][i].astype(float).tolist()) for i in ids]))
result=dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,
    scope='Read-only decomposition of the same saved photon/gas prefixes; no new evolution or final charge.',
    bindings={str(p):s.sha(p) for p in [Path(__file__),Path(s.__file__)]+[f/f'sweep-1/photons/pilot-{n}.npz' for f in [s.OLD,s.OUT] for n in [64,128]]})
s.write(s.OUT/'pressure-profile.json',result);print(json.dumps(result),flush=True)
