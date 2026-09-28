"""Repair the grid-check summation, preserving its original tolerance."""
from pathlib import Path
import inspect,math,resource,time
import resolve_native_pressure_front as r
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));r.original.inf.incident.native.deadline(60)
start=time.monotonic();cpu=time.process_time();error=None
assert r.read(r.OUT/'check-receipt.json')['error'] and not (r.OUT/'grid-check.json').exists()
assert sum(r.read(p)['seconds'] for p in r.OUT.glob('*-receipt.json'))+60+220+30<r.TOTAL
r.write(r.OUT/'grid-sum-repair.json',dict(repair='Use compensated math.fsum for positive quadrature weights instead of naive repeated binary64 addition. Preserve the1e-18 absolute threshold, all physical rules and gates. No physical stage has run.',
    budget_seconds=60,total_action_seconds=r.TOTAL,
    bindings={str(p):r.sha(p) for p in [Path(__file__),Path(r.__file__),r.OUT/'check-receipt.json',r.OUT/'plan.json']}))
source=inspect.getsource(r.check);assert source.count('sum(weights)')==1
source=source.replace('sum(weights)','math.fsum(weights)');ns=dict(r.check.__globals__,math=math)
try:exec(compile(source,__file__,'exec'),ns);ns['check']()
except Exception as exc:error=repr(exc);raise
finally:r.write(r.OUT/'check-sum-receipt.json',dict(action='check_sum',seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
    peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=r.sha(__file__)))
