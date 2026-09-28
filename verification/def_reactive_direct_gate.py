"""Accept the original energy gate before chasing an unnecessary quarter gate."""
from pathlib import Path
import inspect
import def_reactive_paired_tangent as paired

cell_globals=dict(vars(paired))
source=inspect.getsource(paired.cell)
before="        best=old.thermal.g.s.enthalpy_inverse(query,lp,x,float(target),float(T+(-drest-loss)/cp),budget)"
after="""        guess=float(T+(-drest-loss)/cp)
        trial=query(1,lp,guess,x)
        defect=np.longdouble(trial[2])+np.longdouble(trial[1])/trial[0]-target
        if abs(defect)<=budget:
            best=(float(defect),guess,trial)
        else:
            best=old.thermal.g.s.enthalpy_inverse(query,lp,x,float(target),guess,budget)"""
assert source.count(before)==1
exec(compile(source.replace(before,after),__file__,'exec'),cell_globals)


def initialize():
    paired.initialize()
    cell_globals.update({k:getattr(paired,k) for k in ['data','sources','eos','rest','background']})


def block(indices):return [cell_globals['cell'](int(i)) for i in indices]


ns=dict(vars(paired.old),OUT=paired.OUT/'direct-gate',__file__=__file__,paired=paired,
    initialize=initialize,block=block)
text=inspect.getsource(paired.old.run)
text=text.replace('files=[Path(__file__),','files=[Path(__file__),Path(paired.__file__),Path(paired.old.__file__),Path(paired.old.thermal.g.s.__file__),')
text=text.replace("claim='Compute", "preserved_failures=['Original inverse and worker-identity failures retained.', 'Worker-safe pilot forecast 445.41 seconds exceeded the 200 second cap; full grid was not started.'],\n        efficiency_fix='Check the unchanged energy acceptance gate on the native initial guess before invoking the old inverse that targets a quarter of that gate. No acceptance tolerance change.',\n        claim='Compute")
exec(compile(text,__file__,'exec'),ns)

if __name__=='__main__':ns['run']()
