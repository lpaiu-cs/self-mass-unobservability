"""Run the boundary check in a fresh process; preserve non-reentrant init failure."""
import time
import propagate_returned_moments as adapter
b=adapter.base; start=time.monotonic(); error=None
assert 'AssertionError' in b.read(adapter.OUT/'check-receipt.json')['error']
assert b.read(adapter.OUT/'arithmetic-check.json')['passed']
try:
    b.Moments=adapter.Moments
    b.check()
except BaseException as exc:error=repr(exc); raise
finally:
    b.write(adapter.OUT/'check-repair-receipt.json',dict(seconds=time.monotonic()-start,error=error,
        reason='Legacy initializer mutates its imported function once. Run the independent boundary check in a fresh process, without a second initializer call.',
        source_sha256=b.sha(__file__),adapter_sha256=b.sha(adapter.__file__)))
