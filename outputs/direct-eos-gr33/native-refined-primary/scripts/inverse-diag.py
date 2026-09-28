"""Which stored step-end stage states of seg-62 lack an exact coordinate preimage (restored_gas), and how far is it?"""
import os, sys, json
from pathlib import Path
from types import FunctionType
sys.path.insert(0, 'verification')
import numpy as np
work = Path('scratch-inverse267'); SRC = Path('native-true-momentum238-work')
if not work.exists():
    for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
        dst = work/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
    for part in ['sweep-1/photons', 'sweep-1/material']: (work/part).mkdir(parents=True)
import continue_true_momentum as ctm
import resume_coordinate_exact_photons as rc
import stabilize_full_interval_residual as stable; stable.OUT = work
FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=work))(False)
m = ctm.owner.Model(64); LD = np.longdouble
z = np.load('primary267-refined64-work/sweep-1/photons/seg-62.npz')
Q = z['joint_stage_conserved_scaled']; print('stages', len(Q), Q.dtype)
def search(q, limit):
    g = np.column_stack([(q[2]-m.kappa*q[0])/m.eu, q[3]/m.nu, q[0]/m.bu, q[1]/m.su]); worst = {}
    for component, coordinate in [(2, 0), (3, 1), (1, 3), (0, 2)]:
        for attempt in range(limit + 1):
            actual = m.conserved(g)[coordinate]; mask = actual != q[coordinate]
            if not np.any(mask): break
            if attempt == limit: worst[component] = np.flatnonzero(mask).tolist()[:8]; break
            target = np.where(actual[mask] < q[coordinate, mask], LD('inf'), LD('-inf'))
            g[mask, component] = np.nextafter(g[mask, component], target)
        else: pass
        if component in worst: return attempt, worst
    return 0, {}
rows = []
for i in range(1, len(Q), 2):
    a8, bad8 = search(Q[i].copy(), 8)
    if bad8:
        a64, bad64 = search(Q[i].copy(), 256)
        rows.append(dict(stage=i, fails_at_8=bad8, beyond_256=bool(bad64), cells_256=bad64))
print(json.dumps(rows, indent=1)[:3000])
import shutil; shutil.rmtree(work, ignore_errors=True)
