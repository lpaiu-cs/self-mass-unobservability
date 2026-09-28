"""Phase267 step C: the undriven retained-tail path (phase 147/148 Capture) on the current grid, fine clock.

Counterexample candidate. evolve_native_retained_tail.Capture is used unchanged; only its evolution and GR
folders are given (phase 148 did the same for its continuation). resume=False starts from t=0.
Usage: python3 .phase267-undriven.py <root folder> <steps> [<count>] [<audit json>]
  writes <root>/evolution/{coupled-N,source-N,checkpoint-N,...} and <root>/gr/source-N.npz
"""
import json, os, sys, time
from pathlib import Path
from types import FunctionType
root = Path(sys.argv[1]); steps = int(sys.argv[2]); count = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != '-' else None
audit = sys.argv[4] if len(sys.argv) > 4 else None
opened = {}
if audit:
    def hook(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
            e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
    sys.addaudithook(hook)
sys.path.insert(0, 'verification')
import evolve_native_retained_tail as run
EV, GR = root/'evolution', root/'gr'; EV.mkdir(parents=True, exist_ok=False); GR.mkdir()
run.EV = EV; run.GR = GR
run.capture_run = FunctionType(run.capture_run.__code__, dict(run.capture_run.__globals__, OUT=EV))
start = time.monotonic(); row = None
try:
    model = run.Capture(); construct = time.monotonic() - start
    row = model.run(steps, count, False); row['construct_seconds'] = construct
    print(json.dumps({k: v for k, v in row.items() if not isinstance(v, (list, dict))}), flush=True)
finally:
    if audit:
        base = os.path.abspath('.') + '/'
        rows = [dict(path=p.replace(base, ''), **e) for p, e in sorted(opened.items())
                if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(root=str(root), steps=steps, count=count, seconds=time.monotonic() - start, files=rows), indent=1) + '\n')
