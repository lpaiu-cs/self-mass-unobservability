"""Read-only: late-step linear-solve cost in phase 259 (fine) versus the running hgate fine path."""
import json
from pathlib import Path
r = json.loads(Path('native-short-return259-work/right-calls.json').read_text())
s = r['solves']
print('259 fine solves', len(s), 'fallbacks', sum(x['original_fallback'] for x in s), 'total s', round(sum(x['seconds'] for x in s)))
tail = s[-30:]
print('259 last 30 solve seconds', [round(x['seconds']) for x in tail], 'fallback idx', [i for i, x in enumerate(s) if x['original_fallback']][:40])
h = json.loads(Path('native-photon-boundary-hgate265-work/right-calls-fine.json').read_text())['solves']
print('hgate fine solves', len(h), 'last 15 seconds', [round(x['seconds']) for x in h[-15:]])
