"""Phase267 step B: the phase-120 finite-volume initial constraints on the current boundary-layer grid.

Counterexample candidate. def_native_initial_constraints.balanced() with FiniteVolumeData, exactly as its
finite() action, except that the original run's plan-binding hash check is removed (its bound files belong
to the original 19-cell run; provenance is recorded by phase 267) and the wall alarm is 120 s instead of 35 s.
In the original runtime it must reproduce the phase-120 finite-volume outputs.

Usage: python3 .phase267-initial.py <output folder> [<audit json>]
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
import def_native_initial_constraints as ic
source = textwrap.dedent(inspect.getsource(ic.balanced))
for a, b in [("    for p,h in json.loads((OUT/'balanced-plan.json').read_text())['bindings'].items():assert sha(p)==h,p\n", ""),
             ("signal.alarm(35)", "signal.alarm(120)")]:
    assert source.count(a) == 1, a; source = source.replace(a, b)
folder.mkdir(parents=True, exist_ok=False)
ns = dict(vars(ic), Data=ic.FiniteVolumeData, OUT=folder); exec(compile(source, ic.__file__ + '#phase267', 'exec'), ns)
try:
    ns['balanced']()
finally:
    if audit:
        root = os.path.abspath('.') + '/'
        rows = [dict(path=p.replace(root, ''), **e) for p, e in sorted(opened.items())
                if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(folder=str(folder), files=rows), indent=1) + '\n')
