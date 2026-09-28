"""Repair serialized driver amplitude, preserving every physical source array."""
from pathlib import Path
import os
import numpy as np
import apply_driver_aware_radau_gr as r

out=r.OUT;assert not (out/'amplitude-fix.json').exists()
assert not r.read(out/'driver-polynomial-audit.json')['passed'];rows=[]
for n in [64,128]:
    p=out/'gr'/f'source-{n}.npz';old=p.with_name(f'source-{n}-wrong-amplitude.npz');assert not old.exists();os.link(p,old)
    d=dict(np.load(p));before=float(d['drive_amplitude']);assert before==r.AMP and before!=r.incident.ETA
    d['drive_amplitude']=np.array(r.incident.ETA);temp=p.with_name(f'source-{n}-fixed.npz');np.savez_compressed(temp,**d);os.replace(temp,p)
    a,b=dict(np.load(old)),dict(np.load(p));assert all(np.array_equal(a[k],b[k]) for k in a if k!='drive_amplitude')
    rows.append(dict(clock=n,before=before,actual_driver_ETA=r.incident.ETA,old_sha256=r.sha(old),corrected_sha256=r.sha(p),all_physical_arrays_bit_identical=True))
r.write(out/'amplitude-fix.json',dict(classification='Counterexample candidate',rows=rows,source_sha256=r.sha(__file__),
    change='Replace the response coordinate normalization stored by mistake as the driver amplitude with the unchanged live incident.ETA. No physical trajectory, source readout, pulse or acceptance gate changes. The wrong consumer metadata and failed direct-driver audit are preserved.'))
