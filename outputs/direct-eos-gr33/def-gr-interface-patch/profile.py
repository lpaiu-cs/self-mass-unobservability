"""Bounded timing of unchanged resolvents after the contrast budget stop."""
from pathlib import Path
import cProfile
import io
import pstats
import signal
import time
import def_gr_interface_patch as task

OUT=task.OUT;task.install();signal.alarm(60);start=time.monotonic()
task.write(OUT/'profile-plan.json',dict(classification='Counterexample candidate',
    claim='Identify repeated setup cost in the identical full resolvent after the four-contrast forecast exceeded the remaining budget.',
    budget=dict(hard_seconds=60,new_time_paths=0,new_EOS_calls=0),
    decision='Only algebraically identical caching/assembly changes may follow; verify all four full transfer readouts before new trajectories.',
    bindings={str(p):task.go.task.digest(p) for p in [Path(__file__),OUT/'contrasts-pilot.json',Path(task.go.__file__)]}))
p=task.go.Problem();prof=cProfile.Profile();prof.enable()
for k in [1,137,273,410,546,683,819,956,1092,1229,1365,1502,1638,1775,1911,2048]:p.transform(task.go.contour(k,12)[0])
prof.disable();buffer=io.StringIO();pstats.Stats(prof,stream=buffer).sort_stats('cumulative').print_stats(22)
(OUT/'profile.txt').write_text(buffer.getvalue());task.write(OUT/'profile-result.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-start,new_evolutions=0))
print(buffer.getvalue());signal.alarm(0)
