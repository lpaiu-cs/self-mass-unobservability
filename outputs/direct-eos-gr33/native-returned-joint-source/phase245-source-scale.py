"""Compare components of the same accepted short solve, never add charges."""
from pathlib import Path
import hashlib,json
import numpy as np
high=Path('native-complete-radau224-check-work/gr/endpoint-128.npz')
low=Path('native-returned-source245-work/check/gr/endpoint-128.npz')
h,l=np.load(high),np.load(low);assert np.array_equal(h['t'],l['t'])
keys=['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','photon_energy_erg','metric_stress_erg']
value=dict(classification='Counterexample candidate',same_endpoint_times=True,
    horizon_seconds=float(l['t'][-1]),clock=128,
    maximum_spatial_L1_ratio={k:float(np.max(np.sum(abs(l[k]),axis=-1))/np.max(np.sum(abs(h[k]),axis=-1))) for k in keys},
    scope='Actual accepted227short high/low source component norms only; not a charge ratio, continuum error, contraction proof or whole-period result.',
    final_charge_conclusion='unadjudicated',full_goal_complete=False,
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [high,low,Path(__file__)]})
Path('native-returned-source245-work/check/component-scale.json').write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps(value))
