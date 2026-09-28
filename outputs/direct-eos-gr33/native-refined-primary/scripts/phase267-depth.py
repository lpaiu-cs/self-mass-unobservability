"""Phase267 diagnostic: depth decomposition of the endpoint compact charge of a phase267 readout (phase-266 method).

Counterexample candidate (diagnostic only). The source is linear in the per-cell state readout and the per-cell
pulse x geometry terms, so masking cells of the dense source (gr/source-64.npz) and re-running the unchanged
setup/propagator at t=T splits the endpoint value into depth bands; the bands must add up to the unmasked value.
Usage: python3 .phase267-depth.py <readout folder> <label> <band> [<band> ...]   band = all | boundary | cells:a-b
"""
import json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import read_full_captured_history as cap
import extend_retarded_history as ext
rch = cap.prior; KEYS = rch.KEYS
out, label, bands = Path(sys.argv[1]), sys.argv[2], sys.argv[3:]
report = out/f'depth-{label}.json'; result = dict(classification='Counterexample candidate', label=label, rows=[])
def dump(): report.write_text(json.dumps(result, indent=1) + '\n')
d = dict(np.load(out/'gr/field-source-64.npz')); raw = dict(np.load(out/'gr/source-64.npz')); M = float(d['M_cm']); T = d['t'][-1]
T = float(os.environ.get('PHASE267_TIME', T)); assert T <= d['t'][-1]  # evaluation time (retarded, causal)
depth = (d['edges'][-1] - d['edges'])/1e5; result.update(cells=len(d['radius']), M_cm=M, evaluation_time=float(T))
def masked(keep, boundary):
    rr = dict(raw)
    for key in KEYS:
        if key in KEYS[-2:]:
            if not boundary: rr[key] = np.zeros_like(raw[key]); rr['state_coeff_'+key] = np.zeros_like(raw['state_coeff_'+key])
            continue
        rr[key] = np.where(keep[None, :], raw[key], 0)
        rr['state_coeff_'+key] = np.where(keep[None, None, :], raw['state_coeff_'+key], 0)
        rr['geometry_coeff_'+key] = np.where(keep[None, None, :], raw['geometry_coeff_'+key], 0)
    return rr
cap.bind(rch.endpoint.initialize, OUT=out)()
setup = ext.previous.prior.Response.setup
for i, band in enumerate(bands):
    keep = np.zeros(len(d['radius']), bool); boundary = False
    if band == 'all': keep[:] = True; boundary = True
    elif band == 'boundary': boundary = True
    else: a, b = map(int, band.split(':')[1].split('-')); keep[a:b+1] = True
    folder = out/f'depth-in-{label}-{i}'; (folder/'gr').mkdir(parents=True)
    np.savez(folder/'gr/source-64.npz', **masked(keep, boundary))
    s0 = time.monotonic(); m = ext.Response(); ext.bind(setup, INPUT=folder)(m, d, 8)
    free = ext.at(m, m.source, np.array([T])); cells = np.where(keep)[0]
    result['rows'].append(dict(band=band, boundary_terms=boundary, cells=[int(cells[0]), int(cells[-1])] if len(cells) else [],
        depth_km=[float(depth[cells[-1]+1]), float(depth[cells[0]])] if len(cells) else [], charge=-float(free[0][-1, -1])/M, seconds=time.monotonic() - s0))
    del m, free; shutil.rmtree(folder); dump()
parts = [r['charge'] for r in result['rows'] if r['band'] != 'all']; total = [r['charge'] for r in result['rows'] if r['band'] == 'all']
if total: result.update(band_sum=sum(parts), band_sum_relative=(sum(parts) - total[0])/abs(total[0]))
dump(); print(json.dumps({k: v for k, v in result.items() if k != 'rows'}), flush=True)
