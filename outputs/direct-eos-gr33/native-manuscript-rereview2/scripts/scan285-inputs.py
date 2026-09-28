"""List, per result manifest, the bound scripts and result JSONs and the external input paths the scripts read."""
import json, re
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
names = ['native-quad-refined-primary', 'native-eos-sensitivity', 'native-charge-mechanism', 'native-structure-eft-boundary',
         'native-tidal-photosphere', 'native-final-closure', 'native-closure-transit', 'native-atmosphere-reconstruction']
pat = re.compile(r"""['"]([^'"\n]*(?:/home/lpaiu|/mnt/|\.npz|\.npy|work/|runtime)[^'"\n]*)['"]""")
for n in names:
    m = json.loads((root/f'outputs/direct-eos-gr33/{n}-manifest.json').read_text(encoding='utf-8'))
    files = list(m['sha256'])
    scripts = [f for f in files if f.endswith(('.py', '.sh'))]
    results = [f for f in files if f.endswith('.json') and '/scripts/' not in f and not f.endswith(('result.json', 'publication.json', '-manifest.json'))]
    print('==', n, f'({n}-manifest.json)', 'files', len(files))
    print('  scripts:', [Path(s).name for s in scripts])
    print('  results:', [f.split(n + '/')[-1] for f in results][:12], '...' if len(results) > 12 else '')
    ins = set()
    for s in scripts:
        for x in pat.findall((root/s).read_text(encoding='utf-8', errors='replace')): ins.add(x)
    print('  external inputs:', sorted(ins)[:14])
    other = {k: v for k, v in m.items() if k in ('inputs', 'input_sha256', 'external_inputs', 'runtime_inputs')}
    if other: print('  manifest input fields:', {k: (list(v)[:6] if isinstance(v, (dict, list)) else v) for k, v in other.items()})
