"""Counterexample candidate: localize the saved direct-Radau time difference."""
from pathlib import Path
import resource,time
import apply_native_direct_radau as r
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
r.original.inf.incident.native.deadline(20);start=time.monotonic()
r.radau.OUT=r.OUT;r.radau.OLD=Path('native-radau-transfer168-work');r.radau.profile()
r.write(r.OUT/'profile-receipt.json',dict(action='profile',seconds=time.monotonic()-start,
    source_sha256=r.sha(__file__),error=None,
    scope='Read only saved171and168prefix arrays; no new physical steps, model or EOS calls.'))
