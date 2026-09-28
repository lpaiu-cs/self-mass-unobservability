"""Phase266 diagnostic: depth decomposition of the endpoint compact charge of the admitted primary GR field.

Counterexample candidate (diagnostic only; no production output, no new physical step). The accepted
charge 'high' is -U(T,observer)/M of the phase-248 primary field (fine clock, order 8). Its source is linear
in the per-cell state readout and in the per-cell incident-pulse x background-geometry terms, so masking
the source outside a depth band and propagating only to the endpoint T gives that band's share. The shares
must add up to the unmasked value. The potential return (-2 V U dx) is omitted: its endpoint change is
1.27e-83 in the same normalization (phase-248 row), 33 orders below the charge.

Usage: .phase251-reader-launch.py .phase266-depth-probe.py <label> <clock> <order> <band> [<band> ...]
  band = 'all' | 'boundary' | 'cells:a-b' (0-based cell indices, inclusive)
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import extend_retarded_history as ext
base = ext.base; KEYS = base.KEYS
SRC = Path('native-retarded-extension248-work')
label, n, q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); bands = sys.argv[4:]
work = Path(f'scratch-depth266-{label}')
out = Path(f'.phase266-depth-{label}.json')
result = dict(classification='Counterexample candidate', label=label, clock=n, order=q, started=time.strftime('%H:%M:%S'), rows=[])
def dump(): out.write_text(json.dumps(result, indent=1) + '\n')
assert not work.exists(), work
for p in list((SRC/'sweep-0').rglob('*.npz')) + list((SRC/'gr').glob('source-*.npz')) + [SRC/v for v in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
    dst = work/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (work/part).mkdir(parents=True)
try:
    s0 = time.monotonic(); ext.bind(base.endpoint.initialize, OUT=work)(); result['init_seconds'] = time.monotonic() - s0
    d = dict(np.load(SRC/f'gr/source-{n}.npz')); M = float(d['M_cm']); T = d['t'][-1]
    saved = np.load(SRC/f'gr/fields-{n}-g{q}.npz')
    result.update(saved_free_endpoint=-float(saved['direct_and_mass_stress_U'][-1, -1])/M, saved_endpoint=-float(saved['U'][-1, -1])/M,
                  cells=len(d['radius']), M_cm=M, T=float(T))
    depth = (d['edges'][-1] - d['edges'])/1e5
    def masked(keep, boundary):
        dd = dict(d)
        for key in KEYS:
            if key in KEYS[-2:]:
                if not boundary:
                    dd[key] = np.zeros_like(d[key]); dd['state_coeff_'+key] = np.zeros_like(d['state_coeff_'+key])
                continue
            dd[key] = np.where(keep[None, :], d[key], 0.)
            dd['state_coeff_'+key] = np.where(keep[None, None, :], d['state_coeff_'+key], 0)
            dd['geometry_coeff_'+key] = np.where(keep[None, None, :], d['geometry_coeff_'+key], 0)
        return dd
    for band in bands:
        keep = np.zeros(len(d['radius']), bool); boundary = False
        if band == 'all': keep[:] = True; boundary = True
        elif band == 'boundary': boundary = True
        else:
            a, b = map(int, band.split(':')[1].split('-')); keep[a:b+1] = True
        s0 = time.monotonic(); m = ext.Response(); m.setup(masked(keep, boundary), q); setup_s = time.monotonic() - s0
        s0 = time.monotonic(); free = ext.at(m, m.source, np.array([T])); prop_s = time.monotonic() - s0
        value = -float(free[0][-1, -1])/M
        cells = np.where(keep)[0]
        result['rows'].append(dict(band=band, boundary_terms=boundary, cells=[int(cells[0]), int(cells[-1])] if len(cells) else [],
            depth_km=[float(depth[cells[-1]+1]), float(depth[cells[0]])] if len(cells) else [], charge=value,
            setup_seconds=setup_s, propagate_seconds=prop_s))
        del m, free; dump()
except BaseException as exc:
    result['error'] = repr(exc)[:1500]; raise
finally:
    result['finished'] = time.strftime('%H:%M:%S'); dump()
    import shutil; shutil.rmtree(work, ignore_errors=True)
