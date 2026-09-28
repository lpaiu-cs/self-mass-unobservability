"""Phase267 step A2: the phase-130 thermal/rate midpoint bank for any number of interior cells.

Counterexample candidate. def_native_refined_thermochemistry.bank() is reused with only its hard-coded 19-cell
loop generalized to the cells of the boundary-layer geometry it reads, its call cap scaled by 43/19 (resource
only), and its output folder given. In the original runtime (19 cells) it must reproduce the phase-130 bank.

Usage: python3 .phase268-thermal.py <output folder> [<audit json>]
"""
import inspect, json, os, sys, textwrap
from pathlib import Path
folder = Path(sys.argv[1]); audit = sys.argv[2] if len(sys.argv) > 2 else None
opened = {}
if audit:
    def hook(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
            e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
    sys.addaudithook(hook)
sys.path.insert(0, 'verification')
import numpy as np
import def_native_refined_thermochemistry as rt
source = textwrap.dedent(inspect.getsource(rt.bank))
for a, b in [("raw=np.zeros((19,2,4,21));rates=np.zeros((19,2,4,2,len(d['Einf'])));done=np.zeros((19,2,4),bool)",
              "N=len(d['r']);raw=np.zeros((N,2,4,21));rates=np.zeros((N,2,4,2,len(d['Einf'])));done=np.zeros((N,2,4),bool)"),
             ("for j in range(19):", "for j in range(N):"), ("n=chem.old.Native(cap=1000)", "n=chem.old.Native(cap=2500)"),
             ("reused_states=190", "reused_states=int(10*N)")]:
    assert source.count(a) == 1, a; source = source.replace(a, b)
folder.mkdir(parents=True, exist_ok=True); assert not (folder/'bank.npz').exists()  # phase 268: copied status files may precede the bank
ns = dict(vars(rt), OUT=folder); exec(compile(source, rt.__file__ + '#phase267', 'exec'), ns)
try:
    ns['bank']()
finally:
    if audit:
        root = os.path.abspath('.') + '/'
        rows = [dict(path=p.replace(root, ''), **e) for p, e in sorted(opened.items())
                if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(folder=str(folder), files=rows), indent=1) + '\n')
