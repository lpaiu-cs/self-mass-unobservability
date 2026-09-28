import json
from pathlib import Path
w=Path('native-pressure-matter174-work')
for n in [64,128]:
 for folder in [Path('native-mixed-return162-work/sweep-1/material'),w/'sweep-1/material-flux-tangent']:
  p=folder/(f'steps-{n}-reference-128.json' if '162' in str(folder) else f'pilot-{n}.json')
  d=json.loads(p.read_text());print(str(p),{k:d[k] for k in ['steps','completed','substeps','raw_owner_calls','seconds','worker_wall_seconds'] if k in d})
print('spent',sum(json.loads(p.read_text())['seconds'] for p in w.glob('*-receipt.json'))+json.loads((w/'probe-localization.json').read_text())['seconds'])