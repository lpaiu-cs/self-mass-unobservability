"""Phase266 zero-cost check: history of the compact charge q(t)=-U(t,obs)/M from the saved phase-248 primary fields.

Counterexample candidate (diagnostic; reads saved arrays only). Writes .phase266-charge-history.json in the runtime.
"""
import json
import numpy as np
C = 2.99792458e10
R = '/home/lpaiu/work/native-retained-tail-runtime/'
W = R + 'native-retarded-extension248-work/gr/'
src = np.load(W + 'source-128.npz'); M = float(src['M_cm']); edges = src['edges']
depth_edges = (edges[-1] - edges)/1e5
out = dict(classification='Counterexample candidate', scope='Saved phase-248 primary fields only; no new computation.',
           M_cm=M, interior_cell_widths_km=[float(v) for v in np.diff(edges)[:19]/1e5],
           interior_inner_edge_depths_km=[float(v) for v in depth_edges[:20]])
f = np.load(W + 'fields-128-g8.npz'); t = f['t']; q = -f['U'][:, -1]/M; T = t[-1]
rows = []
for frac in [1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 13/16, 7/8, 15/16, 31/32, 1.0]:
    i = int(np.argmin(abs(t - frac*T))); rows.append(dict(t_over_T=float(t[i]/T), q=float(q[i])))
lq = np.log(np.abs(q) + 1e-300); s = np.gradient(lq, t); idx = np.where(t > 0.5*T)[0]
out.update(fine_history=rows, sign_changes=int(np.sum(np.diff(np.sign(q[q != 0])) != 0)),
           last_half_log_slope_per_s=[float(s[idx].min()), float(s[idx].max())],
           visible_depth_at_T_km=float(C*T/2/1e5),
           e_folding_depth_km_range=[float(C/2/s[idx].max()/1e5), float(C/2/s[idx].min()/1e5)],
           endpoint={tag: -float(np.load(W + f'fields-{tag}.npz')['U'][-1, -1])/M for tag in ['128-g8', '64-g8', '128-g4']})
out['endpoint_64_vs_128_relative'] = (out['endpoint']['64-g8'] - out['endpoint']['128-g8'])/out['endpoint']['128-g8']
open(R + '.phase266-charge-history.json', 'w').write(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ['sign_changes', 'last_half_log_slope_per_s', 'visible_depth_at_T_km', 'e_folding_depth_km_range', 'endpoint_64_vs_128_relative']}))
