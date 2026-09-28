"""Phase269: which phase-244 source components differ between the PL-off and PL/MHD 2x readouts at macro N?

Usage: python3 source-diff269.py <61|62|63|64>   (writes .phase269-source-diff-<N>.json in the PL-off runtime and prints)
"""
import json, sys
import numpy as np
n = sys.argv[1]
a = np.load(f'/home/lpaiu/work/native-refined267-runtime/readout267-refined{n}-work/gr/source-64.npz')
b = np.load(f'/home/lpaiu/work/native-eos269-runtime/readout269-eos{n}-work/gr/source-64.npz')
assert set(a.files) == set(b.files), sorted(set(a.files) ^ set(b.files))
rows = {}
for k in sorted(a.files):
    x, y = np.asarray(a[k]), np.asarray(b[k])
    if x.dtype.kind not in 'fc' or x.shape != y.shape: continue
    s = float(np.max(np.abs(x))) if x.size else 0.
    rows[k] = dict(max_abs_plmhd=s, max_diff_over_max=float(np.max(np.abs(y - x))/s) if s > 0 else 0., shape=list(x.shape))
    print('%-42s max|a| %.3e  max|b-a|/max|a| %.3e' % (k, s, rows[k]['max_diff_over_max']))
open(f'/home/lpaiu/work/native-eos269-runtime/.phase269-source-diff-{n}.json', 'w').write(json.dumps(dict(macro=int(n), components=rows), indent=1) + '\n')
