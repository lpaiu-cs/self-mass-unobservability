"""Phase267 step C0: static source template (per-cell geometry and scalars) for the current grid.

Counterexample candidate. The undriven Capture copies every key of its template except the time histories it
overwrites or pops; only the static geometry and scalars pass through (and 'a' enters its energy balance).
They are built with the formulas of def_native_coupled_charge.sources() (radius/volume/edges from the bulk and
atmosphere, the model scalars) and the phase-121/122 Geometry (def_native_projected_evolution.Geometry) for
weight/delay/a/B/re. In the original runtime they must equal
the phase-122 template (def-native-anisotropic-gr/source-128.npz) bitwise.

Usage: python3 .phase267-template.py <output npz> [<reference npz to compare>]
"""
import sys
from pathlib import Path
sys.path.insert(0, 'verification')
import numpy as np
import def_native_stage_energy_history as seh
import def_native_projected_evolution as pe
out = Path(sys.argv[1]); ref = sys.argv[2] if len(sys.argv) > 2 else None
model = seh.flow.Coupled(); f, m, b = model.flow, model.m, model.bulk
radii = np.r_[b.d['r'], m.r]; volumes = np.r_[b.volume, 4*np.pi*m.RJ**2*m.vol]; edges = np.r_[b.d['edges'][:-1], m.rf]
w, delay, a, B, re = pe.Geometry(m)(radii - m.RJ)
template = dict(radius=radii, edges=edges, volume=volumes, weight=w, delay=delay, a=a, B=B, re=re, deep_cells=b.n, cx=f.eos.cx,
                M_cm=m.bg.M*m.R, K_cm=m.bg.K*m.R, inner_material_face_cm=m.rf[0], RJ=m.RJ)
if ref:
    z = np.load(ref); bad = [k for k in template if not (np.asarray(template[k]).shape == z[k].shape and np.array_equal(np.asarray(template[k]), z[k]))]
    print('static keys differing from reference:', bad)
    for k in bad: print('  ', k, np.asarray(template[k]).shape, z[k].shape, float(np.max(abs(np.asarray(template[k], float) - z[k].astype(float)))))
out.parent.mkdir(parents=True, exist_ok=True); assert not out.exists(); np.savez_compressed(out, **template)
print('cells', len(radii), 'deep', b.n)
