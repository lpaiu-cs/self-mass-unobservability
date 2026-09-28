"""Scratch dry run: patched gate names resolve in the live stage and linear-solve namespaces (no physical step)."""
import inspect, json, shutil, sys
from pathlib import Path
sys.path.insert(0, 'verification')
import apply_photon_geometric_boundary_hgate as h
first = h.first
scratch = Path('scratch-hgate265-dryrun')
assert not scratch.exists()
first.OUT = h.OUT = scratch/'work'; first.CHARGE = scratch/'charge'; first.EXTERIOR = scratch/'exterior'
scratch.mkdir()
try:
    print('patched', h.material.initialize.hgate, 'gate strings', h.PATCHED[0].count('gas_gate('), 'newton', 'range(12)' in h.PATCHED[0])
    print('gas_gate first proposal of 265b ->', h.gas_gate([1.07898e-16, 1.00394e-13, 2.69529e-17, 4.01875e-20, 8.59874e-17, 1.93826e-14, 4.32136e-17, 4.95527e-20]))
    print('gas_gate H 2.1e-13 ->', h.gas_gate([0, 2.1e-13, 0, 0, 0, 0, 0, 0]), '; Etilde 1.1e-13 ->', h.gas_gate([1.1e-13, 0, 0, 0, 0, 0, 0, 0]), '; stage1 S 1.1e-13 ->', h.gas_gate([0, 0, 0, 0, 0, 0, 0, 1.1e-13]))
    h.prepare()
    plan = json.loads((h.OUT/'plan.json').read_text())
    print('plan gate', plan['internal_material_gate'], 'H attempts', {k: [round(x/1e-13, 5) for x in v] for k, v in plan['attempt_H_material_relative'].items()})
    Model = first.initialize()
    stages = Model.run.__globals__['stages']
    print('stage globals gas_gate', 'gas_gate' in stages.__globals__, 'names', 'gas_gate' in stages.__code__.co_names, 'consts 12', 12 in stages.__code__.co_consts, 'consts 11', 11 in stages.__code__.co_consts)
    wrapper = stages.__globals__['solve']; cv = inspect.getclosurevars(wrapper).nonlocals
    full = cv['full']
    print('full solve globals gas_gate', 'gas_gate' in full.__globals__, 'names', 'gas_gate' in full.__code__.co_names, 'refinement consts', [c for c in full.__code__.co_consts if isinstance(c, int) and c in (4, 11, 12)])
    print('short solve present', 'short' in cv, 'OUT', stages.__globals__.get('OUT'))
finally:
    shutil.rmtree(scratch)
    print('scratch removed', not scratch.exists())
