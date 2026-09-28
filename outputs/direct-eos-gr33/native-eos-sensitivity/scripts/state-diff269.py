"""Phase269 validity check: does the PL-off physics enter the evolved states? Compare PL-off (269) and PL/MHD (267) 2x runs.

Initial balanced state, undriven coupled history, and driven captures at matching stage indices (photon moments,
radial ports, collision rates): maximum relative differences per array.
"""
import glob, os, sys
import numpy as np
sys.stdout = open('/home/lpaiu/work/native-eos269-runtime/.phase269-state-diff.txt', 'w')  # recorded for publication
A = '/home/lpaiu/work/native-refined267-runtime/'; B = '/home/lpaiu/work/native-eos269-runtime/'
def rel(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != y.shape: return 'shape %s vs %s' % (x.shape, y.shape)
    scale = np.max(np.abs(x)) if x.size else 0.
    return '%.3e (max|a| %.2e)' % (np.max(np.abs(y - x))/scale if scale > 0 else 0., scale)
for rel_path in ['outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume/balanced-initial-state.npz',
                 'outputs/direct-eos-gr33/native-retained-completion/evolution/coupled-128.npz',
                 'outputs/direct-eos-gr33/native-retained-completion/evolution/source-128.npz',
                 'native-incident-drive155-work/fields/born-g8.npz']:
    a, b = np.load(A + rel_path), np.load(B + rel_path)
    print('==', rel_path)
    for k in a.files:
        if k in b.files and a[k].dtype.kind in 'fc': print('  ', k, rel(a[k], b[k]))
caps = sorted(glob.glob(A + 'primary267-refined64-work/captures/captured-64-*.npz'))
for p in [caps[0], caps[len(caps)//4], caps[len(caps)//2], caps[3*len(caps)//4], caps[-1]]:
    q = B + 'primary269-eos64-work/captures/' + os.path.basename(p)
    if not os.path.exists(q): print('missing', os.path.basename(p)); continue
    a, b = np.load(p), np.load(q)
    print('== capture', os.path.basename(p), 'time', float(a['time']), float(b['time']))
    for k in ['photon_moments', 'radial_ports', 'collision_rates']: print('  ', k, rel(a[k], b[k]))
