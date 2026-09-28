"""Two residual corrections test roundoff without changing meshes or gates."""
from pathlib import Path
from types import FunctionType
import inspect
import numpy as np
from scipy.sparse.linalg import splu
import def_free_surface_contrast as original

text=inspect.getsource(original.solve)
before="        y=spsolve(A.multiply((1/scale)[:,None]).tocsc(),rhs/scale)"
after="""        lu=splu(A.multiply((1/scale)[:,None]).tocsc())
        y=lu.solve(rhs/scale)
        extended=A.astype(np.clongdouble)
        for _ in range(2):
            residual_extended=rhs.astype(np.clongdouble)-extended@y.astype(np.clongdouble)
            y+=lu.solve(np.asarray(residual_extended/scale,complex))"""
assert text.count(before)==1
ns=dict(vars(original),OUT=original.repaired.OUT/'direct-contrast-refined',__file__=__file__,splu=splu)
exec(compile(text.replace(before,after),__file__,'exec'),ns)
main=inspect.getsource(original.main)
main=main.replace('files=[Path(__file__),','files=[Path(__file__),Path(original.__file__),')
main=main.replace("claim='Resolve", "precision_control='Two extended-precision residual corrections, same meshes, same equations and unchanged 2 percent gate.',\n        claim='Resolve")
ns['original']=original
exec(compile(main,__file__,'exec'),ns)

if __name__=='__main__':ns['main']()
