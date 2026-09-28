"""Phase267 step D0: phase-150 work folder for the refined grid (plan binding of the undriven history, zero metric).

Counterexample candidate. The module's initialize() hashes EV/coupled-128.npz against plan.json and creates a
zero-additional-geometry metric by copying an earlier 531-cell metric with every delta_* set to zero. Here the same
file is made for the current grid: delta_* are zeros of the current shape, radius_E = re x R of the current source
template (identical construction verified on the original grid), t and the 17 outer boundary scalars are copied from
that earlier metric exactly as phase 150 did. Only phase-150 actions read this file; the driven model does not.
Usage: python3 .phase267-p150-setup.py <folder> <earlier metric npz (read-only source)> <R = radius_E/re of the original grid>
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
folder, earlier = Path(sys.argv[1]), Path(sys.argv[2])
EV = Path('outputs/direct-eos-gr33/native-retained-completion/evolution')
template = np.load('outputs/direct-eos-gr33/def-native-anisotropic-gr/source-128.npz')
assert not folder.exists()
for part in ['collision', 'hydro', 'photons', 'material', 'gr', 'metric/corrected']: (folder/part).mkdir(parents=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
(folder/'plan.json').write_text(json.dumps(dict(classification='Counterexample candidate', phase267_refined_grid=True,
    bindings={(EV/'coupled-128.npz').as_posix(): sha(EV/'coupled-128.npz')}), indent=2) + '\n')
z = dict(np.load(earlier)); cells = len(template['re']); R = float(sys.argv[3])
z2 = {}
for key, value in z.items():
    if key.startswith('delta_'):
        shape = list(value.shape); shape[-1] = cells + (1 if key == 'delta_nu_faces' else 0); z2[key] = np.zeros(shape, value.dtype)
    else: z2[key] = value
snap = np.load(EV/'coupled-128.npz')['snapshot_t']; assert np.max(abs(z['t'] - snap)) < 1e-18
z2['radius_E'] = template['re']*R; assert z2['radius_E'].shape == (cells,)
np.savez_compressed(folder/'metric/corrected/metric-128-g8.npz', **z2)
(folder/'metric/binding.json').write_text(json.dumps(dict(classification='Imported from prior work', path=str(earlier), sha256=sha(earlier),
    zero_additional_geometry=True, phase267_refined_grid_radius_E='re x R of the current source template', cells=cells), indent=2) + '\n')
print('cells', cells)
