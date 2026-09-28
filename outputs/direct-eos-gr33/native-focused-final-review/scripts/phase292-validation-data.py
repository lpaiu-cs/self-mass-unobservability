"""Phase 292: keep the registered independent-validation simulations and replace only the data section.

`simultaneous_inference.py validate` recomputes both the registered simulations and the data section. The simulations do not
depend on the drive phases, but a rerun reproduces their counts only within Monte Carlo error (up to 22 of 8192 differ on the
current linear-algebra backend), so the registered rows of the frozen run are kept. The data section (omnibus statistic unchanged,
physical lag sections recomputed with the ELL1 periastron convention) comes from the rerun. The rerun rows are stored alongside.
"""
import json
from pathlib import Path
R = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/research-completion')
W = R/'withdrawn-periastron-convention'
registered = json.loads((W/'simultaneous-validation.json').read_text(encoding='utf-8'))
rerun = json.loads((R/'simultaneous-validation.json').read_text(encoding='utf-8'))
for k in ('seed', 'nsim', 'threshold', 'nominal_chi6', 'calibration_sha256', 'scope', 'status', 'signal_translation_control'):
    assert registered[k] == rerun[k], k
for k in ('a', 'scale', 'omnibus_null_statistic'):  # equal up to floating-point rounding (relative 1e-8)
    assert abs(registered['data'][k] - rerun['data'][k]) <= 1e-8*abs(registered['data'][k]), k
assert registered['data']['exceeds_threshold'] == rerun['data']['exceeds_threshold']
# The outputs keep the CRLF text format of the original files (Path.write_text default newline on Windows).
(W/'simultaneous-validation-rerun-phase292.json').write_text(json.dumps(rerun, indent=2, allow_nan=False) + '\n', encoding='utf-8')
out = dict(registered)
out['data'] = rerun['data']
out['data_revision'] = ('phase 292: physical lag sections recomputed with the ELL1 periastron convention (eta = e sin varpi); '
                        'simulation rows are the registered frozen run; a rerun with the same seeds on the current backend is kept in '
                        'withdrawn-periastron-convention/simultaneous-validation-rerun-phase292.json and differs only within Monte Carlo error')
(R/'simultaneous-validation.json').write_text(json.dumps(out, indent=2, allow_nan=False) + '\n', encoding='utf-8')
print('sections', [(s['tau'], round(s['minimum_statistic'], 4), s['empty']) for s in out['data']['physical_lag_sections']])
