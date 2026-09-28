"""Scratch end-to-end test of read_photon_boundary_exterior.run (no production output).

Stand-ins: the first attempt's shifted metric (native-photon-boundary265-work) for the
low launch term and the 259 frozen-exterior reading for the rows. Only code paths,
controls and magnitudes are checked; nothing here is a result.
"""
import inspect, json, shutil, sys, time
from pathlib import Path
sys.path.insert(0, 'verification')
import read_photon_boundary_exterior as r
scratch = Path('scratch-complete265-test')
if scratch.exists(): shutil.rmtree(scratch)
r.ACTUAL = Path('native-photon-boundary265-work')
r.FROZEN = Path('native-short-return-exterior259-work/full')
r.OUT = scratch/'complete'; r.p.OUT = r.photons.OUT = r.OUT
src = inspect.getsource(r.run)
for a, b in [("assert audit['passed'] and fro['passed'] and actual['passed'] and actual['photon_geometric_boundary_applied']", "assert audit['passed'] and fro['passed']"),
             ("actual=read(ACTUAL/'result.json')", "actual={}"),
             ("        ACTUAL/'result.json',GEOM/'result.json'", "        GEOM/'result.json'")]:
    assert src.count(a) == 1, a; src = src.replace(a, b)
ns = dict(r.run.__globals__, ACTUAL=r.ACTUAL, FROZEN=r.FROZEN, OUT=r.OUT)
exec(compile(src, 'scratch-run', 'exec'), ns)
start = time.monotonic()
try:
    ns['run']()
except AssertionError as exc:
    print('ASSERTION', str(exc)[:600])
res = json.loads((r.OUT/'result.json').read_text())
print('passed', res['passed'], 'seconds', round(time.monotonic() - start))
print('controls', json.dumps(res['controls']))
print('representation max', max(max(v) for v in res['representation_controls'].values()), {k: ['%.2e' % x for x in v] for k, v in res['representation_controls'].items() if '128-a8-r8' in k})
print('cross', res['cross_checks'])
for c in res['components']: print('component', {k: (v[:22] if isinstance(v, str) else v) for k, v in c.items()})
for row in res['rows']:
    if (row['clock'], row['angular'], row['radial']) == (128, 8, 8):
        print(row['component'], {k: (row[k][:14] if isinstance(row[k], str) else row[k]) for k in ['frozen_exterior', 'geometric_exterior', 'geometric_arrived_energy_erg', 'geometric_exited_energy_erg', 'kappa', 'epsilon', 'relative_change_from_frozen_only']})
shutil.rmtree(scratch)
print('scratch removed', not scratch.exists())
