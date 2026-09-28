"""Phase267 step F: remaining inputs of the final driven model on the current grid.

Counterexample candidate. Run in the refined runtime. The immutable background alias and the Phase151 plan binding
point at the current grid's own undriven history. Inputs that the model reads but does not use (bitwise first-step
check, progress 3) become zero arrays of the current shape: 531/532-long axes become 539/540, other arrays and the
time knots are copied from the original. radius_E comes from the current Phase155 metric (the original files agree).
Usage: python3 .phase267-placeholders.py
"""
import hashlib, json, os, shutil
from pathlib import Path
import numpy as np
ORIG = Path('/home/lpaiu/work/native-retained-tail-runtime')
EV = Path('outputs/direct-eos-gr33/native-retained-completion/evolution')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
radius = np.load('native-incident-drive155-work/metric/corrected/metric-128-g8.npz')['radius_E']; n = len(radius)
old_radius = np.load(ORIG/'native-incident-drive155-work/metric/corrected/metric-128-g8.npz')['radius_E']
def zero(src, dst, keys=None):
    z = np.load(ORIG/src); out = {}
    for k in keys or z.files:
        v = z[k]
        if k == 'radius_E': assert np.array_equal(v, old_radius), src; out[k] = radius
        elif any(s in (531, 532) for s in v.shape): out[k] = np.zeros([{531: n, 532: n + 1}.get(s, s) for s in v.shape], v.dtype)
        else: out[k] = v
    Path(dst).parent.mkdir(parents=True, exist_ok=True); assert not Path(dst).exists(); np.savez_compressed(dst, **out)
    return {k: list(out[k].shape) for k in out}
made = {}
# Immutable undriven history: Phase150 alias is this grid's coupled-128; Phase151 plan binds the same SHA.
alias = Path('retained-native-return150-work/immutable-coupled-128.npz'); assert not alias.exists(); os.link(EV/'coupled-128.npz', alias)
Path('retained-motion-return151-work').mkdir()
Path('retained-motion-return151-work/plan.json').write_text(json.dumps(dict(classification='Counterexample candidate', phase267_refined_grid=True,
    bindings={(EV/'coupled-128.npz').as_posix(): sha(EV/'coupled-128.npz')}), indent=2) + '\n')
# Background template read by the undriven Capture constructor only.
Path('retained-metric-return152-work/fields').mkdir(parents=True); shutil.copyfile(EV/'source-128.npz', 'retained-metric-return152-work/fields/source-128.npz')
for n_ in [64, 128]:
    p = f'native-incident-drive155-work/photons-lift/steps-{n_}-reference-128.npz'; made[p] = zero(p, p, ['t', 'moments', 'collision_transfer'])
    for part, keys in [('photons', ['t', 'moments', 'collision_transfer']), ('material', ['t', 'history_scaled'])]:
        p = f'native-true-momentum238-work/sweep-0/{part}/steps-{n_}-reference-128.npz'; made[p] = zero(p, p, keys)
for p in ['native-incident-self-gr157-work/metric/metric-128-g8.npz', 'native-stage-collisions178-work/material-64.npz',
          'native-stage-collisions178-work/sweep-1/photons/pilot-64.npz']: made[p] = zero(p, p)
for name in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']:
    shutil.copyfile(ORIG/'native-true-momentum238-work'/name, Path('native-true-momentum238-work')/name)
Path('.phase267-placeholders.json').write_text(json.dumps(dict(cells=n, placeholders=made), indent=1) + '\n')
print(json.dumps(dict(cells=n, files=len(made))))
