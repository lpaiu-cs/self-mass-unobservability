"""Phase267 step G: the final driven primary from t=0 through the unchanged production evolve chain.

Counterexample candidate. Replays continue_true_momentum.coarse() (phase-238 arithmetic, stable operator, final
unseeded linear factory, stage guide updates and failure capture) with two changes: initialize() never loads the
phase-238 resume seed (seed=False), and the innermost Model.run(n,'complete-n',n,restart='interval-15-n') becomes
Model.run(n, <label>, <limit>, <restart>). Segments restart exactly from the saved restart_* state.
Every accepted stage is captured as in phases 235-238 (photon moments, radial ports, collision rates) for the readout.
Usage: python3 .phase267-driver.py <work> <label> <limit|-> <restart|-> [<audit json>] [<deadline seconds>]
"""
import json, os, resource, sys, time
from pathlib import Path
from types import FunctionType
work, label = Path(sys.argv[1]), sys.argv[2]
limit = None if sys.argv[3] == '-' else int(sys.argv[3]); restart = None if sys.argv[4] == '-' else sys.argv[4]
audit = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != '-' else None; deadline = int(sys.argv[6]) if len(sys.argv) > 6 else 14400
opened = {}
if audit:
    def hook(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
            e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
    sys.addaudithook(hook)
sys.path.insert(0, 'verification')
SRC = Path('native-true-momentum238-work')
if not work.exists():
    for p in list((SRC/'sweep-0').rglob('*.npz')) + [SRC/n for n in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
        dst = work/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
    for part in ['sweep-1/photons', 'sweep-1/material']: (work/part).mkdir(parents=True)
resource.setrlimit(resource.RLIMIT_AS, (12*1024**3, 12*1024**3))
import numpy as np
import continue_true_momentum as ctm
ctm.joint.previous.original.inf.incident.native.deadline(deadline)
class Done(Exception): pass
def initialize(seed=False):
    FunctionType(ctm.initialize.__code__, dict(ctm.initialize.__globals__, OUT=work))(False)
    constructor = ctm.owner.Model.__init__
    def construct(m, n):
        constructor(m, n); m.capture_clock = n
        def run(steps, *ignored, **ignored_kw):  # production asks for (n,'complete-n',n,'interval-15-n')
            raise Done(type(m).run(m, steps, label, limit, restart))
        m.run = run
    ctm.owner.Model.__init__ = construct
    # Stage capture of phases 235-238 (continue_integer_native.initialize with seed=True), for every accepted stage.
    boundary = ctm.owner.Model.boundary_ports; AMP = ctm.joint.AMP; (work/'captures').mkdir(exist_ok=True)
    if getattr(boundary, 'phase267_capture', False): return  # initialize() may run more than once
    def captured(m, t, x):
        value = boundary(m, t, x); indices = [i for i in range(len(m.stage_t)-2, len(m.stage_t)) if abs(m.stage_t[i]-t) < 1e-18]
        assert len(indices) == 1, ('Capture must be an actual accepted stage', t); i = indices[0]
        moments = np.array([np.sum(x*m.Eweight, axis=(1, 2)), np.sum(x*m.Eweight*m.model.bulk.mu2[None, :, None], axis=(1, 2)), np.sum(x*m.Nweight, axis=(1, 2))])*AMP
        np.savez_compressed(work/'captures'/f'captured-{m.capture_clock}-{i:03d}.npz', time=t, weight=m.stage_h[i], photon_moments=moments,
                            radial_ports=value*AMP, collision_rates=m.stage_collision[i])
        return value
    captured.phase267_capture = True; ctm.owner.Model.boundary_ports = captured
factory = lambda log, n: FunctionType(ctm.factory.__code__, dict(ctm.factory.__globals__, OUT=work))(log, n)
coarse = FunctionType(ctm.prior.coarse.__code__, dict(ctm.prior.coarse.__globals__, OUT=work, initialize=initialize, factory=factory))
start = time.monotonic(); result = dict(label=label, limit=limit, restart=restart)
try:
    coarse(); raise RuntimeError('evolve returned without running the hooked segment')
except Done as done:
    row = done.args[0]; result.update(passed=bool(row['passed']), row={k: v for k, v in row.items() if not isinstance(v, (list, dict))})
except BaseException as exc:
    result.update(passed=False, error=repr(exc)[:3000]); raise
finally:
    result['seconds'] = time.monotonic() - start; result['peak_RSS_bytes'] = 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (work/f'{label}-driver.json').write_text(json.dumps(result, indent=1) + '\n')
    if audit:
        base = os.path.abspath('.') + '/'
        files = [dict(path=p.replace(base, ''), **e) for p, e in sorted(opened.items())
                 if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(files=files), indent=1) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'row'}), flush=True)
sys.exit(0 if result.get('passed') else 3)
