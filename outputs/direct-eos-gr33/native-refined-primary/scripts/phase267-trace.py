"""Phase267 diagnostic: which arrays of each grid-dependent input the final driven model reads, and from where.

Diagnostic only. Same construction and first t=0 macro step as .phase267-audit.py (scratch work directory,
restart=None); records every NpzFile key access with the reading stack. No production directory is written.
Usage: python3 .phase267-trace.py <report json>
"""
import os, sys, json, time, traceback
from pathlib import Path
sys.path.insert(0, 'verification')
import numpy as np
try: from numpy.lib._npyio_impl import NpzFile
except ImportError: from numpy.lib.npyio import NpzFile
WATCH = ['fields/source-', 'native-incident-drive155-work', 'native-incident-self-gr157-work', 'native-stage-collisions178-work',
         'immutable-coupled', 'bank-128', 'evolution/source-128', 'scratch-audit267']
seen = {}
original = NpzFile.__getitem__; original_init = NpzFile.__init__
def init(self, *args, **kw):
    # Buffered owners load every .npz from BytesIO; recover the path from the calling frame's locals.
    name, f = None, sys._getframe(1)
    for _ in range(8):
        if f is None: break
        for var in ('path', 'p', 'file'):
            v = f.f_locals.get(var)
            if isinstance(v, (str, Path)) and str(v).endswith('.npz'): name = str(v); break
        if name: break
        f = f.f_back
    self._trace_name = name
    original_init(self, *args, **kw)
NpzFile.__init__ = init
def getitem(self, key):
    name = str(getattr(self, '_trace_name', None) or getattr(self.zip, 'filename', None) or '?')
    if any(w in name for w in WATCH):
        e = seen.setdefault(name, {})
        if key not in e:
            e[key] = [f'{f.filename.split("/")[-1]}:{f.lineno} {f.name}' for f in traceback.extract_stack()[:-1]
                      if 'site-packages' not in f.filename and '/usr/lib/' not in f.filename][-5:]
    return original(self, key)
NpzFile.__getitem__ = getitem
from types import FunctionType
SRC = Path('native-true-momentum238-work'); scratch = Path('scratch-audit267'); report = Path(sys.argv[1]); result = {}
assert not scratch.exists()
for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
    dst = scratch/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (scratch/part).mkdir(parents=True, exist_ok=True)
import continue_true_momentum as ctm
ctm.OUT = scratch
try:
    t0 = time.monotonic(); FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=scratch))(False)
    m = ctm.owner.Model(64); result['construct_seconds'] = time.monotonic() - t0
    result['after_construct'] = {k: sorted(v) for k, v in seen.items()}
    t0 = time.monotonic(); row = m.run(64, 'trace-64', 1); result['run_seconds'] = time.monotonic() - t0
except BaseException as exc:
    result['error'] = repr(exc)[:2000]
finally:
    result['reads'] = seen
    report.write_text(json.dumps(result, indent=1) + '\n')
    import shutil; shutil.rmtree(scratch, ignore_errors=True)
