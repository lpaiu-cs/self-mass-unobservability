import json
from types import SimpleNamespace
import numpy as np
import gr_increment_structure as c

solver=object.__new__(c.Structure)
solver.mat=SimpleNamespace(B=7.)
solver.shells=np.zeros(1,dtype=np.longdouble)
solver.rhs=lambda x,y,i,B,outer:np.array([1e-18,2e-18,3e-18])*(-1 if outer else 1)
rows=[]
for sub in [4,8]:
    solver.sub=sub
    for outer in [False,True]:
        end=solver.step(0.,1.,np.ones(3),0,7.,outer)
        error=abs(float(solver.shells[0])/(14e-18)-1)
        assert error<1e-14 and end.dtype==np.longdouble
        assert np.all(end<1) if outer else np.all(end>1)
        rows.append(dict(subdivision=sub,outer=outer,shell_relative_error=error))
report=dict(classification='Proven',passed=True,rows=rows,
    scope='Manufactured constant RHS inner/outer sign and sub-binary64 shell increment control only.',
    source_sha256=c.g.c.sha(c.g.ROOT/'verification/gr_increment_structure.py'))
(c.OUT/'step-control.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS inner/outer RK increments retain sub-binary64 changes at both subdivisions')
