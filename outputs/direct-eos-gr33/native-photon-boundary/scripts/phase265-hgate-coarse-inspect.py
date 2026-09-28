"""Read-only: recent coarse linear solves and Newton behaviour of the hgate run."""
import json, time
from pathlib import Path
w = Path('native-photon-boundary-hgate265-work')
def load(p):
    try: return json.loads(Path(p).read_text())
    except Exception as e: return None
r = load(w/'right-calls-coarse.json')
if r:
    calls, solves = r['calls'], r['solves']
    print('coarse calls', len(calls), 'solves', len(solves), 'fallbacks', sum(s['original_fallback'] for s in solves))
    for s in solves[-12:]:
        seg = calls[s['begin']:s['end']]
        print(dict(begin=s['begin'], n=len(seg), seconds=round(s['seconds'], 1), fallback=s['original_fallback']),
              [(c['iterations'], c['info'], ('%.2e' % c['true_relative']) if isinstance(c.get('true_relative'), float) else c.get('true_relative')) for c in seg[-3:]])
    tot = sum(s['seconds'] for s in solves); print('total linear seconds', round(tot))
x = load(w/'last-short-exhaustion-coarse.json'); print('last-short-exhaustion-coarse', x)
f = load(w/'right-calls-fine.json')
if f: print('fine solves', len(f['solves']), 'fallbacks', sum(s['original_fallback'] for s in f['solves']), 'total linear s', round(sum(s['seconds'] for s in f['solves'])))
c = load(w/'capture-64.json'); print('capture-64', c)
st = Path(w/'coarse.stdout.log').read_text()[-1500:] if (w/'coarse.stdout.log').exists() else ''
print('coarse stdout tail:', st)
print(time.strftime('%H:%M:%S'))
