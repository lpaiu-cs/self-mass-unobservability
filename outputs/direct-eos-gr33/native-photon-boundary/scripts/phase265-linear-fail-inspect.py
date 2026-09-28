"""Read-only: the failed last coarse linear solve of the hgate run (refinement history)."""
import json
from pathlib import Path
w = Path('native-photon-boundary-hgate265-work')
r = json.loads((w/'right-calls-coarse.json').read_text())
calls, solves = r['calls'], r['solves']
print('solves', len(solves), 'calls', len(calls))
last = solves[-1]; print('last solve', {k: last[k] for k in last})
seg = calls[last['begin']:]
for i, c in enumerate(seg):
    hist = c.get('history') or []
    print(i, 'iters', c['iterations'], 'info', c['info'], 'true', c.get('true_relative'), 'sec', round(c.get('seconds', 0), 1), 'hist', [('%.2e' % v) if isinstance(v, float) else v for v in hist[:3]], '...', [('%.2e' % v) if isinstance(v, float) else v for v in hist[-3:]])
fl = w/'failed-linear.json'
if fl.exists(): print('failed-linear', fl.read_text()[:1200])
print('failed-linear.npz exists', (w/'failed-linear.npz').exists())
for name in ['last-short-exhaustion-coarse.json', 'capture-64.json', 'coarse-receipt.json']:
    p = w/name
    if p.exists(): print(name, p.read_text()[:600])
ph = sorted(p.name for p in (w/'sweep-1/photons').iterdir())
print('saved photons', ph)
