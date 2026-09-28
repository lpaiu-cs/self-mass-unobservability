"""Use a module-level multiprocessing entry point for the paired native probe."""
from pathlib import Path
import inspect
import def_reactive_paired_tangent as paired


def block(indices):
    return [paired.cell(int(i)) for i in indices]


ns=dict(paired.ns,block=block,OUT=paired.OUT/'worker-safe',__file__=__file__,paired=paired)
text=inspect.getsource(paired.old.run)
text=text.replace('files=[Path(__file__),','files=[Path(__file__),Path(paired.__file__),Path(paired.old.__file__),Path(paired.old.thermal.g.s.__file__),')
text=text.replace("claim='Compute", "preserved_failures=['Original EOS inverse pilot failed its unchanged energy gate.', 'Paired wrapper failed multiprocessing function identity before any probe result; use a named module-level worker.'],\n        claim='Compute")
exec(compile(text,__file__,'exec'),ns)

if __name__=='__main__':ns['run']()
