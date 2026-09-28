import json, time
from pathlib import Path
def load(p):
    try: return json.loads(Path(p).read_text())
    except Exception as e: return None
c = load('native-photon-boundary-hgate265-chain.json')
if c: print('chain', c['state'], c.get('step'), round(c['seconds']), c.get('error'))
w = Path('native-photon-boundary-hgate265-work')
z = load(w/'pipeline-status.json')
if z: print('pipeline', z['state'], z.get('action'), round(z.get('elapsed_seconds', 0)), z.get('error'))
for n in [64, 128]:
    f = load(w/f'capture-{n}.json')
    if f: print('capture', n, f)
for a in ['prepare', 'check', 'coarse', 'fine', 'audit']:
    x = load(w/f'{a}-receipt.json')
    if x: print(a, round(x['seconds'], 1), (x['error'] or 'ok')[:300])
r = load(w/'rejected-joint-stage.json')
if r: print('REJECTED', r['time'], [max(e['material_relative'][1], e['material_relative'][5]) for e in r['equations']][:4])
for name in ['charge', 'exterior']:
    d = Path(f'native-photon-boundary-hgate-{name}265-work')
    if d.exists():
        rec = sorted(str(q.relative_to(d)) for q in d.rglob('*receipt.json') if 'initialization' not in q.parts)
        bad = [q for q in rec if (load(d/q) or {}).get('error')]
        print(name, len(rec), 'receipts', rec[-4:], 'errors', bad)
print(time.strftime('%H:%M:%S'))
