import json, time
from pathlib import Path
def load(p):
    try: return json.loads(Path(p).read_text())
    except Exception: return None
c = load('native-photon-boundary-ld265-chain.json')
if c: print('chain', c['state'], c.get('steps'), round(c['seconds']), c.get('error'), [(d['step'], d['returncode']) for d in c.get('completed', [])])
w = Path('native-photon-boundary-ld265-work')
z = load(w/'pipeline-status.json')
if z: print('pipeline', z['state'], z.get('action'), round(z.get('elapsed_seconds', 0)), z.get('error'))
for n in [64, 128]:
    f = load(w/f'capture-{n}.json')
    if f: print('capture', n, f)
for a in ['prepare', 'boundary', 'check64', 'check128', 'coarse', 'fine', 'audit']:
    x = load(w/f'{a}-receipt.json')
    if x: print(a, round(x['seconds'], 1), (x['error'] or 'ok')[:300], 'ld_events', x.get('long_double_events'))
for n in [64, 128]:
    r = load(w/f'restart-regression-{n}.json')
    if r: print('restart', n, r['passed'], r['arrays'], r['replayed_stage_native_rates'], '%.2e' % r['max_saved_native_relative'])
for tag in ['coarse', 'fine']:
    e = load(w/f'long-double-{tag}.json')
    if e: print('LD', tag, [(v['actual_step'], v['accepted'], v['corrections'], round(v['seconds'])) for v in e['events']])
for stderr in ['.phase265-ld-chain-prepare.stderr.log', '.phase265-ld-chain-boundary.stderr.log', '.phase265-ld-chain-check64.stderr.log', '.phase265-ld-chain-check128.stderr.log']:
    p = Path(stderr)
    if p.exists() and p.stat().st_size: print(stderr, p.read_text()[-600:])
for name in ['charge', 'exterior']:
    d = Path(f'native-photon-boundary-ld-{name}265-work')
    if d.exists():
        rec = sorted(str(q.relative_to(d)) for q in d.rglob('*receipt.json') if 'initialization' not in q.parts)
        bad = [q for q in rec if (load(d/q) or {}).get('error')]
        print(name, len(rec), 'receipts', rec[-3:], 'errors', bad)
print(time.strftime('%H:%M:%S'))
