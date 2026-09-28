"""Phase267 step E: Phase155 Born field and incident metric knots on the current grid.

Counterexample candidate. Runs the unchanged def_native_incident_drive.fields(reuse=False): the Driver builds the
Born field (born_return) for orders 4 and 8 and evaluates the incident metric at the 17 clock knots. Only the output
folder is redirected. In the original runtime every array must equal the Phase155 files bitwise.
Usage: python3 .phase267-born.py <new folder> [<reference folder to compare>]
"""
import json, sys
from pathlib import Path
sys.path.insert(0, 'verification')
import numpy as np
import def_native_incident_drive as drive
folder = Path(sys.argv[1]); ref = Path(sys.argv[2]) if len(sys.argv) > 2 else None
assert not folder.exists()
drive.OUT, drive.FIELDS, drive.METRIC = folder, folder/'fields', folder/'metric'
drive.FIELDS.mkdir(parents=True); (drive.METRIC/'corrected').mkdir(parents=True)
drive.fields(False)
if ref:
    report = {}
    for p in sorted(folder.rglob('*.npz')):
        a, b = np.load(p), np.load(ref/p.relative_to(folder))
        report[str(p.relative_to(folder))] = [k for k in a.files if not (a[k].shape == b[k].shape and np.array_equal(a[k], b[k]))] + \
            [f'missing:{k}' for k in b.files if k not in a.files]
    print(json.dumps(dict(differing=report), indent=1))
