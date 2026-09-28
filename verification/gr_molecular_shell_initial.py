"""Conservative finite-volume initialization from retained GR shell increments.

Counterexample candidate: a separately defined restriction, not completion of
the interrupted 8/16-node quadrature or a continuous EOS certificate.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType
import argparse
import inspect
import json
import time

import numpy as np
import gr_molecular_conservative_initial as original

ROOT, ld, e = original.ROOT, original.ld, original.e
OUT = e.g.OUT/'gr-molecular-shell-initial'


def save(name, value):
    e.write(OUT/name, value)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel, value in plan['bindings'].items():
        assert e.digest(ROOT/rel) == value, rel
    original.bindings()
    return plan


def inputs():
    data = dict(np.load(original.molecular.OUT/'molecular-state-17-4.npz'))
    moments = dict(np.load(OUT/'moments.npz'))
    baryons, energy = moments['baryons'][::-1], moments['energy'][::-1]
    rf, r = data['radius_faces_m'][::-1].astype(ld)*100, data['r_mid_m'][::-1].astype(ld)*100
    volume = 4*ld(str(np.pi))/3*np.diff(rf**3)
    E = energy/volume
    mf = np.r_[ld(0), np.cumsum(e.GRAV*energy)]
    mass = mf[:-1]+e.GRAV*E*4*np.pi/3*(r**3-rf[:-1]**3)
    a = 1/np.sqrt(1-2*mass/r)
    lr = np.log(baryons/(a*volume))
    X = data['X'][::-1].astype(ld)
    rest = (X/e.g.c.A)@e.g.c.W*e.C**2
    u = E/np.exp(lr)-rest
    return data, list(zip(np.arange(len(r))[::-1], lr, data['lnT'][::-1], u, X))


def prepare():
    assert not OUT.exists()
    prior = original.bindings()
    branch = json.loads((original.OUT/'branches.json').read_text())
    assert branch['passed'] and branch['path_sha256'] == e.digest(original.OUT/'path-4.npz')
    refined = json.loads((original.molecular.OUT/'result.json').read_text())
    assert refined['completed'] and refined['all_interfaces_passed'] and refined['finite_face_refinement_passed']
    data = dict(np.load(original.molecular.OUT/'molecular-state-17-4.npz'))
    shell = np.load(original.OUT/'path-4.npz')['shell_mass_geom_m']
    assert np.all(shell > 0) and len(shell) == len(data['dm'])
    OUT.mkdir()
    np.savez_compressed(OUT/'moments.npz', baryons=data['dm'].astype(ld), energy=shell*100/e.GRAV)
    files = [Path(__file__), Path(original.__file__), original.OUT/'plan.json',
        original.OUT/'branches.json', original.OUT/'path-4.npz',
        original.molecular.OUT/'molecular-state-17-4.npz', original.molecular.OUT/'result.json',
        original.molecular.OUT/'manifest.json', OUT/'moments.npz']
    save('plan.json', dict(classification='Counterexample candidate',
        bindings={p.relative_to(ROOT).as_posix():e.digest(p) for p in files},
        cells=5735, pilot_cells=[0,1175,1176,1792,1920,2972,5734],
        restriction_gates=prior['restriction_gates'],
        method='Use each original baryon inventory and each independently retained RK4 shell mass increment as the finite-volume moments. Reuse the original conservative primitive inverse, mass/lapse reconstruction and fresh native molecular EOS/opacity diffusion initialization unchanged. No global mass fit, requested target substitution, shell rescaling or entropy fitting.',
        distinction='This is a separately frozen finite-volume initial-data definition. It does not finish or pass the original 8/16-node quadrature. Subcell entropy and intermediate spatial error are not certified. Existing quadrature blocks and failures remain intact.',
        pilot='Seven actual shell targets, including both sides of the molecular join and the interrupted region. Reuse their primitive solutions in the full-grid assembly. Stop on any original restriction failure.',
        budget=dict(workers=4,blas_threads=1,gpu=False,pilot_timeout_seconds=90,
            full_timeout_seconds=600,maximum_full_runs=1,
            authorization_gate='Extrapolated pilot CPU cost divided by four must be below 360 seconds; retain a 600-second hard cap. Failure or time cap stops this method without automatic refinement or relaxed gates.'),
        physical_EOS_certified=False, full_GR_evolution=False,
        native_subcell_quadrature_completed=False, continuous_errors_certified=False))
    data, rows = inputs()
    assert min(row[3] for row in rows) > 0
    print('PREPARED retained-shell molecular restriction',len(rows),flush=True)


def pilot():
    plan = bindings(); assert not (OUT/'pilot.json').exists()
    original.worker_init(); _, tasks = inputs()
    chosen = [r for r in tasks if r[0] in plan['pilot_cells']]
    began = time.monotonic(); rows = []
    for task in chosen:
        started = time.monotonic()
        value = original.restrict_cell(task)
        rows.append(dict(cell=int(task[0]),logT=str(value[0]),aux=[float(v) for v in value[1]],
            residual=value[2],calls=value[3],seconds=time.monotonic()-started))
    predicted = sum(row['seconds'] for row in rows)/len(rows)*plan['cells']/plan['budget']['workers']
    result = dict(classification='Counterexample candidate',passed=predicted<360,
        rows=rows,seconds=time.monotonic()-began,
        estimated_full_wall_seconds=predicted,
        estimate_limit='Seven states, not a measured distribution over the remaining grid; runtime cap is independent of this estimate.')
    save('pilot.json',result); assert result['passed'],result
    print('PASS native shell restriction pilot',json.dumps(result),flush=True)


class ReusePilot:
    def __init__(self, pool):
        self.pool = pool
        self.saved = {r['cell']:(ld(r['logT']),np.asarray(r['aux']),r['residual'],r['calls'])
            for r in json.loads((OUT/'pilot.json').read_text())['rows']}

    def map(self, function, tasks, chunksize=4):
        tasks = list(tasks)
        missing = [t for t in tasks if int(t[0]) not in self.saved]
        fresh = self.pool.map(function, missing, chunksize=chunksize)
        for task in tasks:
            index = int(task[0])
            yield self.saved[index] if index in self.saved else next(fresh)


def assemble(pool):
    # Reuse the entire native restriction, metric, opacity and heat initialization.
    # Only the frozen input moments change; no post-result gate is replaced.
    source = inspect.getsource(original.assemble)
    start = source.index('    baryons, energy = np.zeros')
    end = source.index('    baryons, energy = baryons[::-1]',start)
    source = source[:start]+"    moments = np.load(OUT/'moments.npz')\n    baryons, energy = moments['baryons'].copy(), moments['energy'].copy()\n"+source[end:]
    namespace = dict(original.assemble.__globals__, OUT=OUT, bindings=bindings, save=save)
    exec(compile(source,str(Path(__file__)),'exec'),namespace)
    namespace['assemble'](pool)


def run():
    plan = bindings(); pilot_result = json.loads((OUT/'pilot.json').read_text())
    assert pilot_result['passed'] and not (OUT/'initial.npz').exists()
    began = time.monotonic()
    try:
        with ProcessPoolExecutor(max_workers=plan['budget']['workers'],initializer=original.worker_init) as pool:
            assemble(ReusePilot(pool))
    except Exception as error:
        save('failure.json',dict(classification='Counterexample candidate',error=repr(error)))
        raise
    save('runtime.json',dict(seconds=time.monotonic()-began,workers=plan['budget']['workers'],
        reused_pilot_cells=plan['pilot_cells']))
    # finalize only after every result is closed; resource logs live outside OUT.
    save('initial-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():e.digest(p)
        for p in OUT.iterdir() if p.is_file() and p.name!='initial-manifest.json'}))


def initialize(pool):
    return FunctionType(original.initialize.__code__,dict(original.initialize.__globals__,
        OUT=OUT,bindings=bindings))(pool)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('action',choices=['prepare','pilot','run'])
    globals()[parser.parse_args().action]()
