"""Phase273 part 1: history of the compact charge q(t) split into the instantaneous pulse x background part (direct
coupling, the static-coefficient-like response) and the dynamic state part (displacement), 4x solution.

Counterexample candidate (diagnostic). One unchanged setup per mask; the retarded field is evaluated at many times
(causal in the source, as the phase-267 readout does for t61..t63).
Usage: python3 phase273-history.py <readout folder> <mask: all|state|geometry> <n times>
"""
import json, shutil, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import read_full_captured_history as cap, extend_retarded_history as ext
out, spec, n = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
rch = cap.prior; KEYS = rch.KEYS
raw = dict(np.load(out/'gr/source-64.npz')); d = dict(np.load(out/'gr/field-source-64.npz')); M = float(d['M_cm']); T = float(d['t'][-1])
rr = dict(raw)
for key in KEYS:
    if spec == 'geometry': rr['state_coeff_' + key] = np.zeros_like(raw['state_coeff_' + key])
    if spec == 'state' and key not in KEYS[-2:]: rr['geometry_coeff_' + key] = np.zeros_like(raw['geometry_coeff_' + key])
cap.bind(rch.endpoint.initialize, OUT=out)(); setup = ext.previous.prior.Response.setup
folder = out/f'history-in-{spec}'; (folder/'gr').mkdir(parents=True); np.savez(folder/'gr/source-64.npz', **rr)
m = ext.Response(); ext.bind(setup, INPUT=folder)(m, d, 8)
times = np.linspace(T/n, T, n); free = ext.at(m, m.source, times); shutil.rmtree(folder)
q = [-float(v)/M for v in np.asarray(free[0])[:, -1]]
(out/f'history-{spec}.json').write_text(json.dumps(dict(mask=spec, times=times.tolist(), charge=q, T=T, D=float(raw['drive_duration'])), indent=1) + '\n')
print(spec, 'q(T/4) %.4e q(T/2) %.4e q(3T/4) %.4e q(T) %.4e' % (q[n//4 - 1], q[n//2 - 1], q[3*n//4 - 1], q[-1]))
