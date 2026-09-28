import json
from pathlib import Path
import numpy as np
new = Path('native-photon-boundary265-work'); old = Path('native-short-return259-work')
r = json.loads((new/'rejected-joint-stage.json').read_text())
print('rejected keys', list(r.keys())[:30])
print({k: (v if not isinstance(v, (list, dict)) else (str(v)[:300])) for k, v in r.items()})
z = np.load(new/'rejected-joint-stage.npz'); print('npz', {k: z[k].shape for k in z.files})
c = json.loads((old/'checks-64.json').read_text())
print('259 checks keys', list(c.keys()))
nw = c['newton']
print('259 newton entries', len(nw))
for i in range(min(4, len(nw))):
    print(i, [dict(relative=x['relative'], mat=max(x.get('material_relative', [0]))) for x in nw[i]][-3:])
# which material component peaks
for i in range(min(4, len(nw))):
    last = nw[i][-1]; m = last.get('material_relative')
    if m: print('259 step', i, 'last material_relative', ['%.3e' % v for v in m])
st = c.get('stages')
print('stages type', type(st), (len(st) if st else None))
if st: print(st[:3] if isinstance(st, list) else list(st)[:3])
