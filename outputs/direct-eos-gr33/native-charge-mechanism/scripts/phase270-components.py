"""Phase270 diagnostic: split the endpoint compact charge by source component and by part (state vs pulse x geometry).

Counterexample candidate (diagnostic only). The readout source of every key is (state readout polynomial from
state_coeff_<key>) + (incident pulse and Born field x geometry_coeff_<key>) (apply_driver_aware_radau_gr.coefficients),
and the retarded field is linear in it. Masking coefficient arrays of the stored source (gr/source-64.npz) and
re-running the unchanged setup/propagator at t=T (as the phase-266/267 depth tool does for cells) splits the endpoint
value exactly; the parts must add up to the unmasked value.
Specs: all | state | geometry | state:<key> | geometry:<key> | key:<key>
Usage: python3 phase270-components.py <readout folder> <label> <spec> [<spec> ...]
"""
import json, os, shutil, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, 'verification')
import read_full_captured_history as cap
import extend_retarded_history as ext
rch = cap.prior; KEYS = rch.KEYS; BOUNDARY = KEYS[-2:]
out, label, specs = Path(sys.argv[1]), sys.argv[2], sys.argv[3:]
report = out/f'components-{label}.json'; result = dict(classification='Counterexample candidate', label=label, keys=KEYS, rows=[])
def dump(): report.write_text(json.dumps(result, indent=1) + '\n')
d = dict(np.load(out/'gr/field-source-64.npz')); raw = dict(np.load(out/'gr/source-64.npz')); M = float(d['M_cm']); T = float(d['t'][-1])
def masked(spec):
    part, _, which = spec.partition(':')
    rr = dict(raw)
    for key in KEYS:
        keep_state = spec == 'all' or (part == 'state' and which in ('', key)) or (part == 'key' and which == key)
        keep_geometry = key not in BOUNDARY and (spec == 'all' or (part == 'geometry' and which in ('', key)) or (part == 'key' and which == key))
        if not keep_state: rr['state_coeff_' + key] = np.zeros_like(raw['state_coeff_' + key]); rr[key] = np.zeros_like(raw[key])
        if key not in BOUNDARY and not keep_geometry: rr['geometry_coeff_' + key] = np.zeros_like(raw['geometry_coeff_' + key])
    return rr
cap.bind(rch.endpoint.initialize, OUT=out)()
setup = ext.previous.prior.Response.setup
for i, spec in enumerate(specs):
    folder = out/f'components-in-{label}-{i}'; (folder/'gr').mkdir(parents=True)
    np.savez(folder/'gr/source-64.npz', **masked(spec))
    s0 = time.monotonic(); m = ext.Response(); ext.bind(setup, INPUT=folder)(m, d, 8)
    free = ext.at(m, m.source, np.array([T]))
    result['rows'].append(dict(spec=spec, charge=-float(free[0][-1, -1])/M, seconds=time.monotonic() - s0))
    print(spec, result['rows'][-1]['charge'], flush=True)
    del m, free; shutil.rmtree(folder); dump()
rows = {r['spec']: r['charge'] for r in result['rows']}
if 'all' in rows:
    for group in ['state', 'geometry']:
        if group in rows: result[f'{group}_share'] = rows[group]/rows['all']
    if 'state' in rows and 'geometry' in rows: result['state_plus_geometry_relative'] = (rows['state'] + rows['geometry'] - rows['all'])/abs(rows['all'])
    keys = [k for k in rows if k.startswith('key:')]
    if keys: result['key_sum_relative'] = (sum(rows[k] for k in keys) - rows['all'])/abs(rows['all'])
dump(); print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'keys')}), flush=True)
