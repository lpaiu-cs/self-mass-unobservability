"""Bounded comparison of the unchanged and updated heat matrix support."""
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
import def_background_transport_charge as run

signal.alarm(45)
resource.setrlimit(resource.RLIMIT_AS, (int(4e9), int(4e9)))
start = time.monotonic()
p = run.Problem()
old = [p.Gq.copy(), p.GEraw.copy()]
data = np.load(run.OUT/'p4-64-input.npz')
p.update(data)
rows = []
for name, a, b in zip(['Gq', 'GEraw'], old, [p.Gq, p.GEraw]):
    difference = b-a
    rows.append(dict(name=name, before_nonzeros=a.nnz, after_nonzeros=b.nnz,
        largest_before=float(max(abs(a.data))), largest_after=float(max(abs(b.data))),
        difference_nonzeros=difference.nnz))
# No-op must reproduce the original operator, not only its largest entries.
zeros = np.zeros(len(data['delta_lnT']))
p.update(dict(delta_lnT=zeros, delta_lnrho=zeros, opacity=p.model.heat.d['opacity'][:, 0]))
for row, a, b in zip(rows, old, [p.Gq, p.GEraw]):
    difference = b-a
    row['noop_max_absolute'] = float(max(abs(difference.data))) if difference.nnz else 0.
    row['noop_relative_max'] = row['noop_max_absolute']/row['largest_before']
    row['noop_nonzeros'] = b.nnz
    coo = b.tocoo()
    width = abs(coo.row-coo.col)
    far = width > 32
    row['far_nonzeros'] = int(sum(far))
    row['far_maximum'] = float(max(abs(coo.data[far]))) if any(far) else 0.
value = dict(classification='Counterexample candidate', rows=rows,
             seconds=time.monotonic()-start)
run.write(run.OUT/'sparsity.json', value)
print(value, flush=True)
