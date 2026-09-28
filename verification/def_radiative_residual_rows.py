"""Identify the exact equations limiting the full radiative joint solve."""
import signal
import time
import numpy as np
import def_orbital_radiative_transport as run

signal.alarm(90)
start = time.monotonic()
p = run.Problem()
source = run.direct_source[:run.direct_source.index('    perm=self.permutation;')]
source += '''    saved=np.load(OUT/'precision-state.npz'); answer=saved['answer']
    residual=bvec-block@answer
    scale=abs(bvec)+abs(block)@abs(answer)+1e-100
    error=abs(residual)/scale
    slots=np.argsort(error)[-12:][::-1]
    rows=[]
    for index in slots:
        position=np.argwhere(m.indices==index).tolist() if index<m.size else []
        row=block.getrow(index)
        rows.append(dict(index=int(index),position=position,error=float(error[index]),
            residual=[float(residual[index].real),float(residual[index].imag)],scale=float(scale[index]),
            terms=[[int(k),float(abs(v)),float(abs(answer[k]))] for k,v in zip(row.indices,row.data)]))
    write(OUT/'residual-rows.json',dict(classification='Counterexample candidate',rows=rows,
        seconds=time.monotonic()-START))
    print(rows,flush=True)
'''
ns = dict(run.namespace, START=start)
exec(compile(source, 'residual-rows-source.py', 'exec'), ns)
ns['solve'](p, 1, 'diagnosis')
