"""Phase267 step D: phase-150 (def_retained_native_return) actions on the current grid, with an input audit.

Counterexample candidate. Runs initialize() and the named actions of the unchanged module in the folder given by
RETAINED_NATIVE_OUTPUT (the module reads it from the environment), without the original plan-binding checks.
Usage: RETAINED_NATIVE_OUTPUT=<folder> python3 .phase267-p150.py <audit json> <action> [<action> ...]
  actions: collision | bank | hydro_pilot | hydro_production (as in the module's execute())
"""
import json, os, sys, time
from pathlib import Path
audit = sys.argv[1]; actions = sys.argv[2:]
opened = {}
def hook(event, args):
    if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
        e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
        if any(w in p for w in TRACE) and 'stack' not in e:
            import traceback; e['stack'] = [f'{f.filename.split("/")[-1]}:{f.lineno} {f.name}' for f in traceback.extract_stack()[:-1]][-8:]
TRACE = [w for w in os.environ.get('PHASE267_TRACE', '').split(',') if w]
sys.addaudithook(hook)
sys.path.insert(0, 'verification')
import def_retained_native_return as r
rows = []; start = time.monotonic()
try:
    r.initialize()
    for a in actions:
        t0 = time.monotonic()
        if a == 'collision': r.collision_sources()
        elif a == 'bank': r.bank()
        elif a.startswith('hydro_'): r.hydro_parallel('pilot' in a)
        else: raise ValueError(a)
        rows.append(dict(action=a, seconds=time.monotonic() - t0))
        print(json.dumps(rows[-1]), flush=True)
finally:
    base = os.path.abspath('.') + '/'
    files = [dict(path=p.replace(base, ''), **e) for p, e in sorted(opened.items())
             if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
    Path(audit).write_text(json.dumps(dict(out=os.environ.get('RETAINED_NATIVE_OUTPUT'), actions=rows, seconds=time.monotonic() - start, files=files), indent=1) + '\n')
