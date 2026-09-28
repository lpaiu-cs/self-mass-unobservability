"""Phase267 input audit: every file the final driven primary model reads while constructing and taking one t=0 step.

Diagnostic only. Builds the phase-238 final-arithmetic Model(64) in a scratch work directory, records every
opened file with Python audit hooks, runs one coarse step from t=0 (restart=None), and writes the list with
array shapes so grid-dependent inputs can be identified. No production directory is written.
"""
import os, sys, json, time
from pathlib import Path
opened = {}
def hook(event, args):
    if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.path.abspath(os.fsdecode(args[0])); mode = args[1] if len(args) > 1 else None
        write = isinstance(mode, str) and any(c in mode for c in 'wax+')
        e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if write else 'reads'] += 1
sys.addaudithook(hook)
sys.path.insert(0, 'verification')
import numpy as np
from types import FunctionType
SRC = Path('native-true-momentum238-work'); scratch = Path('scratch-audit267'); label = sys.argv[1] if len(sys.argv) > 1 else 'audit'
steps = int(sys.argv[2]) if len(sys.argv) > 2 else 1
report = Path(f'.phase267-audit-{label}.json'); result = dict(started=time.strftime('%H:%M:%S'))
assert not scratch.exists()
for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
t0 = time.monotonic(); import continue_true_momentum as ctm; result['import_seconds'] = time.monotonic() - t0
ctm.OUT = scratch
try:
    t0 = time.monotonic(); FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=scratch))(False)
    m = ctm.owner.Model(64); result['construct_seconds'] = time.monotonic() - t0
    result.update(n=int(m.n), q=int(getattr(m, 'q', -1)), nf=int(getattr(m, 'nf', -1)), bulk_n=int(m.model.bulk.n) if hasattr(m, 'model') else None,
                  run_globals_OUT=str(m.run.__globals__.get('OUT')))
    if steps:
        t0 = time.monotonic(); row = m.run(64, f'{label}-64', steps); result['run_seconds'] = time.monotonic() - t0
        result['row'] = {k: v for k, v in row.items() if not isinstance(v, (list, dict))}
except BaseException as exc:
    result['error'] = repr(exc)[:2000]
finally:
    files = []
    for p, e in sorted(opened.items()):
        if p.endswith('.py') or '/usr/lib/' in p or '/site-packages/' in p or '/dist-packages/' in p or p.startswith('/proc') or p.startswith('/dev'): continue
        info = dict(path=p, **e)
        if os.path.isfile(p):
            info['bytes'] = os.path.getsize(p)
            if p.endswith('.npz') and e['reads']:
                try:
                    with np.load(p, allow_pickle=False) as z: info['arrays'] = {k: list(z[k].shape) for k in z.files}
                except Exception as exc: info['arrays_error'] = repr(exc)[:200]
        files.append(info)
    result.update(files=files, finished=time.strftime('%H:%M:%S'))
    report.write_text(json.dumps(result, indent=1) + '\n')
    import shutil; shutil.rmtree(scratch, ignore_errors=True)
