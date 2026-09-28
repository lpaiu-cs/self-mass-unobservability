"""Phase 292: does the registered six-coefficient validation reproduce within one environment? Prints the a=0 inclusion count."""
import json, sys
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
sys.path.insert(0, str(root/'verification'))
import numpy as np
import simultaneous_inference as si
g = si.setup()
threshold = json.loads(si.CALPATH.read_text())['threshold']
y, rest = si.draws(g, 2026090913, 0.)
out = si.scan(g, y, rest)
print('numpy', np.__version__, 'hits', int(np.count_nonzero(out['statistic'] <= threshold)), 'eig0', float(g['eigenvalues'][0]), float(np.sum(g['q6'][:5, 0])))
