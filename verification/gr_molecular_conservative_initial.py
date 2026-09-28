"""Counterexample candidate: conservative initial data for the molecular EOS.

Reconstruct local mass increments, integrate the new material reference, and
invert its baryon/coordinate-energy moments with the same molecular EOS.  The
existing live evolution and the molecular shooting run are separate inputs.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, MethodType, SimpleNamespace
import argparse
import inspect
import json
import shutil
import subprocess

import numpy as np
import gr_coupled_evolution as e
import gr_conservative_evolution as restriction
import gr_conservative_two_carrier as conservative
import gr_full_subcell_reference as quadrature
import gr_increment_structure as increments
import gr_molecular_structure_runner as molecular_runner
import gr_two_carrier_evolution as two

ROOT, ld = e.ROOT, e.ld
OUT = e.g.OUT/'gr-molecular-conservative-initial'
molecular = molecular_runner.module((molecular_runner.OUT/'candidate.py').read_text())
write, digest = e.write, e.digest


def save(name, value):
    write(OUT/name, value)


def worker_init():
    e.EOS = molecular.model.EOS()
    e.OPACITY = e.opacity_tables.Opacity()


def restrict_cell(row):
    fn = restriction.restrict_cell
    return FunctionType(fn.__code__, dict(fn.__globals__,
        evolution=SimpleNamespace(material=two.material)))(row)


class Structure(molecular.g.c.Structure):
    step = increments.Structure.step

    def __init__(self, data, subdivision, folder):
        fn = molecular.g.c.Structure.__init__
        FunctionType(fn.__code__, dict(fn.__globals__, OUT=molecular.OUT,
            EOS=molecular.model.EOS))(self, 'molecular', data, 17, subdivision)
        self.shells = np.zeros(len(self.mat.lp), dtype=ld)
        self.stats = dict(calls=0, evaluations=0, maximum_score=0., label=folder.name)
        self.inverse = molecular.inverse(self.stats, folder)
        fn = molecular.g.c.Structure.state
        self.state = MethodType(FunctionType(fn.__code__, dict(fn.__globals__,
            be=SimpleNamespace(invert=self.inverse))), self)


def selfcheck():
    """Exercise the actual native temperature inverse across the new EOS grid."""
    worker_init()
    data, _, _ = molecular.inputs()
    errors = []
    for i in [0, 1175, 1176, 2972, 5734]:
        lr, lt, x = data['lnd'][i], data['lnT'][i], data['X'][i]
        expected = two.material((lr, lt, x))
        actual, aux, residual, calls = restrict_cell((i, lr, ld(lt)+ld('.01'), ld(expected[2]), x))
        assert len(aux) == 30 and abs(actual-ld(lt)) < ld('2e-8'), (i, actual, lt)
        assert abs(residual) < 1e-8 and np.all(np.isfinite(aux))
        errors.append(dict(cell=i, lnT_error=float(abs(actual-lt)),
                           thermal_residual=residual, native_calls=calls))
    assert e.symbolic()['passed']
    print('MOLECULAR CONSERVATIVE INVERSE', json.dumps(errors), flush=True)
    return errors


def prepare(workers):
    assert not OUT.exists(), 'Preserve every started initial-data calculation.'
    assert molecular_runner.source()[0] == (molecular.OUT/'candidate.py').read_text()
    molecular.bindings()
    molecular.precision.verify()
    molecular.model.verify()
    report = json.loads((molecular.OUT/'GR-4.json').read_text())
    assert report['interface_max'] < 1e-8
    table = json.loads((molecular.OUT/'table-result.json').read_text())
    assert table['completed'] and table['output_sha256'] == digest(molecular.OUT/'molecular-adiabats-17.npz')
    checks = selfcheck()
    OUT.mkdir()
    source = molecular.reference.OUT/'reference-state.npz'
    shutil.copy2(source, OUT/'reference-state.npz')
    paths = [Path(__file__), Path(e.__file__), Path(restriction.__file__),
        Path(conservative.__file__), Path(quadrature.__file__), Path(increments.__file__),
        Path(two.__file__), Path(e.opacity_tables.__file__),
        ROOT/'verification/conservative_star.py', ROOT/'verification/gr_subcell_precision_fallback.py',
        molecular.OUT/'plan.json', molecular.OUT/'candidate.py', molecular.OUT/'table-result.json',
        molecular.OUT/'molecular-adiabats-17.npz', molecular.OUT/'molecular-state-17-4.npz',
        molecular.OUT/'GR-4.json', source, OUT/'reference-state.npz',
        molecular.model.OUT/'manifest.json', molecular.precision.OUT/'manifest.json']
    plan = dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        bindings={p.relative_to(ROOT).as_posix():digest(p) for p in paths},
        runtime_sha256={str(p):digest(p) for p in [molecular.model.LIB, molecular.precision.BRIDGE]},
        cells=5735, subdivision=4, nodes=[8, 16], workers=workers, block_size=128,
        rtol=2e-10, atol=1e-12, endpoint_tolerance=1e-8,
        volume_relative_tolerance=1e-5, shell_mass_relative_tolerance=1e-5,
        inventory_relative_tolerance=1e-10, finite_quadrature_relative_tolerance=1e-5,
        interface_tolerance=1e-8, saved_face_relative_tolerance=1e-8,
        restriction_gates=dict(maximum_cell_baryon_relative_defect=1e-14,
            maximum_cell_energy_relative_defect=1e-13, maximum_native_thermal_residual=1e-8,
            maximum_original_baryon_relative_defect=1e-11, maximum_original_mass_relative_defect=2e-12),
        native_inverse_controls=checks,
        geometry='Reuse the matched molecular GR-4 shooting parameters without fitting or tuning baryons. Reuse the retained-increment RK4 and the original one-metre central and surface analytic seeds. Save local shell mass changes independently of cumulative face subtraction.',
        moments='Reuse positive 8/16-node material-coordinate DOP853 quadratures with molecular tables, native molecular roots and matching extended-arithmetic fallback. Baryon weights are the given coordinate measure. Compare native energy and volume against independently retained shell increments and face geometry.',
        initialization='Use the new 16-node baryon and coordinate-energy moments. Reconstruct the constant-energy mass constraint, set rho=B/(a V), and solve native molecular u at that rho. Recompute the existing baryon-face diffusion law and its two-carrier split at the restricted state with fresh molecular EOS/opacity. No old-model cell moments, auxiliaries or luminosity are imported.',
        limits='Finite initial data for the same conditional LTE two-carrier model and reflecting wall, not an isentropic remap, physical EOS/atmosphere certificate, continuum bound or completed evolution. The pending molecular GR-8 comparison must pass before a new evolution may start.')
    save('plan.json', plan)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel, value in plan['bindings'].items():
        assert digest(ROOT/rel) == value, rel
    for name, value in plan['runtime_sha256'].items():
        assert digest(Path(name)) == value, name
    return plan


def branches():
    plan = bindings()
    assert not (OUT/'path-4.npz').exists() and not (OUT/'branch-roots').exists()
    data = dict(np.load(OUT/'reference-state.npz'))
    saved = dict(np.load(molecular.OUT/'molecular-state-17-4.npz'))
    parameters = json.loads((molecular.OUT/'GR-4.json').read_text())['parameters']
    solver = Structure(data, 4, OUT/'branch-roots')
    m, B = solver.mat, solver.mat.B
    error, inner, outer = solver.branches(parameters, record=True)
    pc, rs, ms = parameters
    rs, ms = m.R*np.exp(rs), molecular.g.c.gr.TARGET*np.exp(ms)
    p, energy, baryon = solver.state(m.lp[0], 0)
    f = 1-2*ms/rs
    massB = energy/baryon*np.sqrt(f)
    lpB = -(energy+p)*(ms+4*np.pi*rs**3*p)/(4*np.pi*rs**4*baryon*np.sqrt(f)*p)
    w0 = min(m.dm[0]*1e-6, 1e-8/abs(lpB*B))
    solver.shells[0] += ld(massB)*ld(B)*ld(w0)
    solver.shells[-1] += ld(inner[0, 2])*ld(B)
    assert np.all(solver.shells > 0) and np.all(np.isfinite(solver.shells))
    n = len(m.lp)
    faces = np.zeros((n+1, 3), dtype=ld)
    faces[0] = [rs/m.R, ms/B, m.lp[0]]
    faces[1:m.split+1] = outer[1:, 1:]
    faces[m.split:n] = inner[1:, 1:][::-1]
    faces[n] = [0, 0, pc]
    radius, mass = faces[:, 0]*ld(m.R), faces[:, 1]*ld(B)
    face_gap = max(float(abs(radius-saved['radius_faces_m']).max()/rs),
                   float(abs(mass-saved['mass_faces_geom']).max()/ms))
    np.savez_compressed(OUT/'path-4.npz', faces=faces, radius_m=radius,
        mass_geom_m=mass, shell_mass_geom_m=solver.shells, inner=inner, outer=outer)
    report = dict(classification='Counterexample candidate', interface_max=float(abs(error).max()),
        saved_face_relative_difference=face_gap, root_statistics=solver.stats,
        shell_integral_relative_boundary_difference=float(solver.shells.sum()/ld(ms)-1))
    report['passed'] = report['interface_max'] < plan['interface_tolerance'] and face_gap < plan['saved_face_relative_tolerance']
    report['path_sha256'] = digest(OUT/'path-4.npz')
    save('branches.json', report)
    assert report['passed'], report
    print('MOLECULAR CONSERVATIVE BRANCHES', json.dumps(report), flush=True)


def block(job):
    label, _ = job
    report = json.loads((OUT/'branches.json').read_text())
    assert report['passed'] and report['path_sha256'] == digest(OUT/'path-4.npz')
    # Reuse the entire existing node/seed/weight algorithm. Only its declared
    # reference, geometry provider, native inverse and output bindings differ.
    source = inspect.getsource(quadrature.block)
    before = """    inverse=FunctionType(g.audit.strict_invert.__code__,dict(g.audit.strict_invert.__globals__,
        ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))"""
    after = "    inverse=molecular.inverse(stats,OUT/('node-roots-'+label))"
    assert source.count(before) == 1
    source = source.replace(before, after)
    namespace = dict(quadrature.block.__globals__, OUT=OUT, save=save, bindings=bindings,
        molecular=molecular, g=SimpleNamespace(OUT=OUT, ROOT=ROOT, c=molecular.g.c),
        structure=SimpleNamespace(OUT=OUT, Structure=lambda data, sub:
            Structure(data, sub, OUT/('geometry-roots-'+label))))
    exec(compile(source, str(Path(__file__)), 'exec'), namespace)
    result = namespace['block'](job)
    assert result['all_passed'], ('Molecular quadrature failure retained', label)
    return result


def make_star(pool, saved):
    star = object.__new__(two.TwoCarrierStar)
    star.pool, star.material_cache = pool, {}
    for key in ['rf', 'r', 'volume', 'indices', 'base', 'qscale']:
        setattr(star, key, saved[key].copy())
    star.n = len(star.r)
    star.source = molecular.OUT/'molecular-state-17-4.npz'
    star.area_difference_over_volume = 4*np.pi*np.diff(star.rf**2)/star.volume
    star.original_baryons = ld(saved['original_baryons'])
    star.original_mass = ld(saved['original_mass'])
    if 'aux' in saved:
        star.initial_aux = saved['aux'].copy()
        star.material_cache = {e.material_key(row):value for row, value in zip(
            zip(star.base[:, 0], star.base[:, 1], star.base[:, 5:]), star.initial_aux)}
    return star


def assemble(pool):
    plan = bindings()
    original = dict(np.load(molecular.OUT/'molecular-state-17-4.npz'))
    count = plan['cells']
    baryons, energy = np.zeros((2, count), dtype=ld)
    seen = set()
    for path in sorted(OUT.glob('block-*-nodes-16.npz')):
        data = np.load(path)
        cells = data['cells']
        assert not seen.intersection(cells)
        seen.update(cells)
        weights, aux = data['coordinate_weights_cm3'].astype(ld), data['eos'].astype(ld)
        baryons[cells] = np.sum(aux[:, :, 0]*data['metric_a']*weights, axis=1)
        energy[cells] = np.sum(aux[:, :, 0]*(data['C_X'].astype(ld)[:, None]*e.C**2+aux[:, :, 2])*weights, axis=1)
    assert seen == set(range(count))
    baryons, energy = baryons[::-1], energy[::-1]
    rf, r = original['radius_faces_m'][::-1].astype(ld)*100, original['r_mid_m'][::-1].astype(ld)*100
    volume = 4*ld(str(np.pi))/3*np.diff(rf**3)
    assert rf[0] == 0 and np.all(np.diff(rf) > 0) and np.all((r > rf[:-1]) & (r < rf[1:]))
    mean_energy = energy/volume
    mf = np.r_[ld(0), np.cumsum(e.GRAV*energy)]
    mass = mf[:-1]+e.GRAV*mean_energy*4*np.pi/3*(r**3-rf[:-1]**3)
    a = 1/np.sqrt(1-2*mass/r)
    lr = np.log(baryons/(a*volume))
    x = original['X'][::-1].astype(ld)
    rest = (x/e.g.c.A)@e.g.c.W*e.C**2
    target_u = mean_energy/np.exp(lr)-rest
    rows = list(pool.map(restrict_cell, zip(np.arange(count)[::-1], lr,
        original['lnT'][::-1], target_u, x), chunksize=4))
    aux = np.asarray([row[1] for row in rows], dtype=ld)
    base = np.column_stack([lr, [row[0] for row in rows], np.zeros((count, 3)), x]).astype(ld)
    saved = dict(rf=rf, r=r, volume=volume, indices=np.arange(count)[::-1], base=base,
        qscale=np.exp(lr)*(rest+aux[:, 2])+aux[:, 1], aux=aux,
        original_baryons=np.sum(original['dm'].astype(ld)), original_mass=ld(original['mass_faces_geom'][0])*100)
    star = make_star(pool, saved)
    z = star.state(star.base)
    nu = np.log(z['N'])
    nu_faces = star.faces(nu)
    nu_faces[-1] = ld('.5')*np.log1p(-2*z['mf'][-1]/rf[-1])
    transport = dict(dm=baryons[::-1], lnT=star.base[::-1, 1], nu=nu[::-1],
        nu_faces=nu_faces[::-1], radius_faces_m=rf[::-1]/100)
    flux = molecular.g.s.baryon_face_diffusion(transport, aux[::-1, 21])
    lum = np.r_[ld(0), flux[::-1], ld(0)]
    Q = (lum[:-1]+lum[1:])/2/(4*np.pi*r*r*z['N']**2*e.C)
    star.base[:, 3] = Q/star.qscale
    star.base[:, 4] = star.base[:, 3]*aux[:, 21]/aux[:, 24]
    actual = star.state(star.base)
    actual_b, actual_e = actual['a']*actual['D']*volume, actual['E']*volume
    report = dict(classification='Counterexample candidate', cells=count,
        maximum_cell_baryon_relative_defect=float(abs(actual_b/baryons-1).max()),
        maximum_cell_energy_relative_defect=float(abs(actual_e/energy-1).max()),
        maximum_native_thermal_residual=max(abs(row[2]) for row in rows),
        maximum_original_baryon_relative_defect=float(abs(actual_b.sum()/saved['original_baryons']-1)),
        maximum_original_mass_relative_defect=float(abs(e.GRAV*actual_e.sum()/saved['original_mass']-1)),
        maximum_native_calls_per_cell=max(row[3] for row in rows),
        physical_EOS_certified=False, full_GR_evolution=False)
    report['passed'] = all(report[k] < value for k, value in plan['restriction_gates'].items())
    save('restriction.json', report)
    assert report['passed'], report
    saved.update(base=star.base, baryons=baryons, energy=energy, initial_face_Linf=lum,
        restriction_residual=np.asarray([row[2] for row in rows]))
    np.savez_compressed(OUT/'initial.npz', **saved)
    save('initial-manifest.json', dict(sha256={p.relative_to(ROOT).as_posix():digest(p)
        for p in OUT.rglob('*') if p.is_file() and p.name != 'initial-manifest.json'}))
    print('MOLECULAR CONSERVATIVE INITIAL DATA', json.dumps(report), flush=True)


def initialize(pool):
    bindings()
    for rel, value in json.loads((OUT/'initial-manifest.json').read_text())['sha256'].items():
        assert digest(ROOT/rel) == value, rel
    assert json.loads((OUT/'restriction.json').read_text())['passed']
    finer = json.loads((molecular.OUT/'result.json').read_text())
    assert finer['completed'] and finer['all_interfaces_passed'] and finer['finite_face_refinement_passed']
    molecular.verify()
    star = make_star(pool, dict(np.load(OUT/'initial.npz')))
    fn = conservative.initialize
    star = FunctionType(fn.__code__, dict(fn.__globals__,
        parent=SimpleNamespace(initialize=lambda _:star)))(pool)
    star.operator_step, star.operator_time = 0, ld(0)
    return star


def run():
    plan = bindings()
    try:
        branches()
        jobs = [(f'block-{start:04}', list(range(start, min(start+plan['block_size'], plan['cells']))))
                for start in range(0, plan['cells'], plan['block_size'])]
        with ProcessPoolExecutor(max_workers=plan['workers'], initializer=worker_init) as pool:
            completed = 0
            for result in pool.map(block, jobs):
                completed += len(result['rows'])
                save('progress.json', dict(classification='Counterexample candidate', completed_cells=completed))
            assemble(pool)
    except Exception as error:
        save('failure.json', dict(classification='Counterexample candidate', error=repr(error)))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['selfcheck', 'prepare', 'run'])
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    assert 1 <= args.workers <= 15
    if args.action == 'prepare':
        prepare(args.workers)
    else:
        globals()[args.action]()
