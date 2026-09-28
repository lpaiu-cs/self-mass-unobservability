"""Bounded precision diagnosis on the unchanged full-radius joint equations."""
import json
import signal
import time
import numpy as np
import def_orbital_radiative_transport as run


def main():
    signal.alarm(120)
    start = time.monotonic()
    p = run.Problem()
    source = (run.OUT/'second-direct-solver-source.py').read_text()
    anchor = 'for _ in range(4):answer+=invert(bvec-block@answer)'
    assert source.count(anchor) == 1
    source = source.replace(anchor, '''history=[]
    for iteration in range(24):
        answer+=invert(bvec-block@answer)
        defect=bvec-block@answer
        e=float(np.max(abs(defect)/(abs(bvec)+abs(block)@abs(answer)+1e-100)))
        history.append(e)
        print('REFINEMENT',iteration,e,flush=True)
        if e<1e-11:break
    np.savez_compressed(OUT/'precision-state.npz',answer=answer,bvec=bvec,energy_scale=self.energy_scale)
    write(OUT/'precision.json',dict(classification='Counterexample candidate',history=history,
        matrix_nnz=block.nnz,scaled_factor_nnz=factor.L.nnz+factor.U.nnz,
        elapsed_seconds=time.monotonic()-START,passed=e<1e-9))''')
    ns = dict(run.namespace, START=start)
    exec(compile(source, 'radiative-precision-source.py', 'exec'), ns)
    ns['solve'](p, 1, 'precision')


if __name__ == '__main__':
    main()
