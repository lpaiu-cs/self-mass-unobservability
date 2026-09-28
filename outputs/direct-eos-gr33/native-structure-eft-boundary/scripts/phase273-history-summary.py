"""Phase273: summarize the charge history (all / state / geometry) of the 4x solution and the pulse at the surface."""
import json
import numpy as np
R = '/home/lpaiu/work/native-refined268-runtime/readout268-quad64-work/'
h = {m: json.load(open(R + f'history-{m}.json')) for m in ['all', 'state', 'geometry']}
t = np.array(h['all']['times']); D = h['all']['D']; T = h['all']['T']
qa, qs, qg = (np.array(h[m]['charge']) for m in ['all', 'state', 'geometry'])
p = t/D; f = np.where((p > 0) & (p < 1), (4*p*(1 - p))**4, 0.)  # incident pulse phase at the surface (x ~ 0)
closure = np.max(np.abs(qs + qg - qa))/np.max(np.abs(qa))
print('T %.6e D %.6e  closure max|state+geometry-all|/max|all| %.2e' % (T, D, closure))
print('max|geometry| %.3e at t=%.3e ; max|state| %.3e at t=%.3e' % (np.max(abs(qg)), t[np.argmax(abs(qg))], np.max(abs(qs)), t[np.argmax(abs(qs))]))
for i in list(range(7, len(t), 8)):
    share = qs[i]/qa[i] if qa[i] else float('nan')
    print('t %.4e (t/D %.3f, surface pulse %.3f)  all %+.4e  state %+.4e  geometry %+.4e  state share %.6f' % (t[i], p[i], f[i], qa[i], qs[i], qg[i], share))
json.dump(dict(T=T, D=D, closure=float(closure), max_geometry=float(np.max(abs(qg))), t_max_geometry=float(t[np.argmax(abs(qg))]),
               max_state=float(np.max(abs(qs))), t_max_state=float(t[np.argmax(abs(qs))]),
               state_share_at=[dict(t=float(t[i]), share=float(qs[i]/qa[i]) if qa[i] else None) for i in range(7, len(t), 8)]),
          open(R + 'history-summary.json', 'w'), indent=1)
