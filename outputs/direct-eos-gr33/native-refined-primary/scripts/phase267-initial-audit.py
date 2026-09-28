"""Phase267 step B2: the phase-120 finite-volume audit (native anchors) for the current grid.

Counterexample candidate. verify_native_initial_constraints.finite() reused verbatim except the native call
cap (100 -> 150, i.e. x27/19) and alarm (35 -> 60 s); resource only. Writes audit.json with native_anchors,
inventory-audit.npz and the corrected integrated-content controls into the given finite-volume folder.

Usage: python3 .phase267-initial-audit.py <finite-volume folder>
"""
import inspect, json, sys, textwrap
from pathlib import Path
from types import FunctionType, SimpleNamespace
folder = Path(sys.argv[1])
sys.path.insert(0, 'verification')
import numpy as np
import verify_native_initial_constraints as v
task = v.task
run_source = textwrap.dedent(inspect.getsource(v.run))
for a, b in [("native=task.previous.chem.old.Native(cap=100)", "native=task.previous.chem.old.Native(cap=150)"), ("signal.alarm(35)", "signal.alarm(60)")]:
    assert run_source.count(a) == 1, a; run_source = run_source.replace(a, b)
modified = SimpleNamespace(**dict(vars(task), OUT=folder, Data=task.FiniteVolumeData))
ns = dict(vars(v), task=modified, OUT=folder); exec(compile(run_source, v.__file__ + '#phase267', 'exec'), ns); ns['run']()
finite_source = textwrap.dedent(inspect.getsource(v.finite))
old = "    modified=SimpleNamespace(**dict(vars(task),OUT=OUT/'finite-volume',Data=task.FiniteVolumeData))\n    FunctionType(run.__code__,dict(globals(),task=modified,OUT=modified.OUT))()\n"
assert finite_source.count(old) == 1; finite_source = finite_source.replace(old, "")
ns2 = dict(vars(v), modified=modified, SimpleNamespace=SimpleNamespace, FunctionType=FunctionType)
exec(compile(finite_source, v.__file__ + '#phase267-finite', 'exec'), ns2); ns2['finite']()
