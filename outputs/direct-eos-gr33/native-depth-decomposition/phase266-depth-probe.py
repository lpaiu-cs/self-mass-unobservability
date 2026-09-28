"""Phase266 diagnostic: depth decomposition of the endpoint compact charge of the admitted primary GR field.

Counterexample candidate (diagnostic only; no production output, no new physical step). The accepted
charge 'high' is -U(T,observer)/M of the phase-248 primary field. Its setup reads the phase-244 source
file (INPUT/gr/source-{n}.npz), not the array handed to it, so each band writes a masked copy of that file
into a scratch INPUT and binds the unchanged production setup to it. The source is linear in the per-cell
state readout and in the per-cell pulse x geometry terms; the band shares must add up to the unmasked
value. The potential return is omitted (2.5e-11 of the endpoint in the saved phase-248 field).

First version (preserved as .phase266-depth-b*-ineffective.json) masked the array handed to setup, which
setup does not read; every band returned the full value.

Usage: .phase251-reader-launch.py .phase266-depth-probe.py <label> <clock> <order> <band> [<band> ...]
  band = 'all' | 'boundary' | 'cells:a-b' (0-based cell indices, inclusive)
"""
import json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import extend_retarded_history as ext
base = ext.base; KEYS = base.KEYS
SRC = Path('native-retarded-extension248-work'); RAW = ext.INPUT  # native-full-captured244-work
setup = ext.previous.prior.Response.setup
label, n, q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); bands = sys.argv[4:]
work = Path(f'scratch-depth266-{label}')
out = Path(f'.phase266-depth-{label}.json')
result = dict(classification='Counterexample candidate', label=label, clock=n, order=q, raw_source=str(RAW/f'gr/source-{n}.npz'),
              started=time.strftime('%H:%M:%S'), rows=[])
def dump(): out.write_text(json.dumps(result, indent=1) + '\n')
assert not work.exists(), work
for p in list((SRC/'sweep-0').rglob('*.npz')) + list((SRC/'gr').glob('source-*.npz')) + [SRC/v for v in ['normalization.json', 'photon-conservation-plan.json', 'check-result.json']]:
    dst = work/p.relative_to(SRC); dst.parent.mkdir(parents=True, exist_ok=True); os.link(p, dst)
for part in ['sweep-1/photons', 'sweep-1/material']: (work/part).mkdir(parents=True)
try:
    s0 = time.monotonic(); ext.bind(base.endpoint.initialize, OUT=work)(); result['init_seconds'] = time.monotonic() - s0
    d = dict(np.load(SRC/f'gr/source-{n}.npz')); raw = dict(np.load(RAW/f'gr/source-{n}.npz'))
    M = float(d['M_cm']); T = d['t'][-1]
    saved = np.load(SRC/f'gr/fields-{n}-g{q}.npz')
    result.update(saved_free_endpoint=-float(saved['direct_and_mass_stress_U'][-1, -1])/M, saved_endpoint=-float(saved['U'][-1, -1])/M,
                  cells=len(d['radius']), M_cm=M, T=float(T))
    depth = (d['edges'][-1] - d['edges'])/1e5
    def masked(keep, boundary):
        rr = dict(raw)
        for key in KEYS:
            if key in KEYS[-2:]:
                if not boundary:
                    rr[key] = np.zeros_like(raw[key]); rr['state_coeff_'+key] = np.zeros_like(raw['state_coeff_'+key])
                continue
            rr[key] = np.where(keep[None, :], raw[key], 0)
            rr['state_coeff_'+key] = np.where(keep[None, None, :], raw['state_coeff_'+key], 0)
            rr['geometry_coeff_'+key] = np.where(keep[None, None, :], raw['geometry_coeff_'+key], 0)
        return rr
    for i, band in enumerate(bands):
        keep = np.zeros(len(d['radius']), bool); boundary = False
        if band in ('all', 'copy'): keep[:] = True; boundary = True
        elif band == 'boundary': boundary = True
        else:
            a, b = map(int, band.split(':')[1].split('-')); keep[a:b+1] = True
        folder = work/f'in-{i}'; (folder/'gr').mkdir(parents=True)
        if band == 'copy': os.link(RAW/f'gr/source-{n}.npz', folder/f'gr/source-{n}.npz')
        else: np.savez(folder/f'gr/source-{n}.npz', **masked(keep, boundary))
        s0 = time.monotonic(); m = ext.Response(); ext.bind(setup, INPUT=folder)(m, d, q); setup_s = time.monotonic() - s0
        s0 = time.monotonic(); free = ext.at(m, m.source, np.array([T])); prop_s = time.monotonic() - s0
        value = -float(free[0][-1, -1])/M
        cells = np.where(keep)[0]
        result['rows'].append(dict(band=band, boundary_terms=boundary, cells=[int(cells[0]), int(cells[-1])] if len(cells) else [],
            depth_km=[float(depth[cells[-1]+1]), float(depth[cells[0]])] if len(cells) else [], charge=value,
            setup_seconds=setup_s, propagate_seconds=prop_s))
        del m, free; shutil.rmtree(folder); dump()
except BaseException as exc:
    result['error'] = repr(exc)[:1500]; raise
finally:
    result['finished'] = time.strftime('%H:%M:%S'); dump()
    shutil.rmtree(work, ignore_errors=True)
