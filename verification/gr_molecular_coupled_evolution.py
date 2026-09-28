"""Counterexample candidate: evolve the conservatively initialized molecular EOS.

Reuse the frozen conservative equations, native composition tangent, time mesh
and gates.  Every path starts from the new initial state; old accepted time
states are never imported.  Physical EOS, atmosphere and observational closure
remain separate requirements.
"""
from types import SimpleNamespace
import gr_coupled_evolution as original_e
import gr_conservative_composition_tangent as method
import gr_molecular_conservative_initial as initial
from gr_conservative_composition_tangent import *

OUT = original_e.g.OUT/'gr-molecular-coupled-evolution'
REFINEMENT, PREFIX = None, 0
_original_material = original_e.material


def guarded_material(row):
    assert type(original_e.EOS) is initial.molecular.model.EOS, 'Wrong EOS in a molecular evolution worker'
    return _original_material(row)


def worker_init():
    initial.worker_init()
    # Both native state evaluations and composition probes call two.material,
    # which in turn uses this provider in each separate worker process.
    original_e.material = guarded_material


e = SimpleNamespace(**dict(vars(original_e), worker_init=worker_init))


def initialize(pool):
    if pool is None:
        pool = parent.parent.old.CachedOnly()
    return initial.initialize(pool)


_plan_bindings = FunctionType(prior.bindings.__code__, globals())


def bindings():
    plan = _plan_bindings()
    template = method.bindings()
    initial.bindings()
    for key in ['coordinate_edges_seconds', 'refinements', 'conduction_time',
                'finite_conservation_gates', 'time_refinement_gate']:
        assert plan[key] == template[key], key
    assert plan['imported_prefix_steps'] == 0
    assert plan['nonlinear_absolute_tolerances'] == ATOL.astype(float).tolist()
    return plan


def require_inputs():
    paths = [initial.OUT/'initial-manifest.json', initial.OUT/'restriction.json',
             initial.molecular.OUT/'manifest.json', initial.molecular.OUT/'result.json']
    missing = [str(p.relative_to(e.ROOT)) for p in paths if not p.is_file()]
    if missing:
        raise RuntimeError(('Molecular conservative inputs are still incomplete', missing))
    initial.bindings()
    for rel, digest in json.loads(paths[0].read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    assert json.loads(paths[1].read_text())['passed']
    result = json.loads(paths[3].read_text())
    assert result['completed'] and result['all_interfaces_passed'] and result['finite_face_refinement_passed']
    initial.molecular.verify()
    return paths


def prepare():
    assert not OUT.exists(), 'Preserve every previously started molecular path.'
    inputs = require_inputs()
    template = method.bindings()
    source_plan = initial.bindings()
    assert prior.symbolic()['passed']
    star = initialize(None)
    assert star.n == source_plan['cells']
    state = star.evaluate(np.zeros_like(star.base))
    assert np.all(state['dU'] == 0) and np.all(state['dBX'] == 0)
    for key in ['m', 'mf', 'a', 'N', 'Q', 'aux']:
        assert np.array_equal(state[key], star.reference[key]), key
    zero_flux = np.zeros(star.n+1, dtype=ld)
    initial_budget = budget(star, state, zero_flux, zero_flux)[0]
    initial_cone = parent.cones(state)
    assert budget_passed(initial_budget, template)
    assert initial_cone['sampled_cone_inside_light_cone']
    sources = [Path(__file__), Path(initial.__file__), method.OUT/'plan.json',
               initial.OUT/'plan.json', *inputs]
    plan = dict(classification='Counterexample candidate',
        checkpoint=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        cells=source_plan['cells'],
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in sources},
        runtime_sha256=dict(template['runtime_sha256'], **source_plan['runtime_sha256']),
        candidate='The molecular EOS, its native conservative initial data and recomputed two-carrier heat flux in the existing conservative GR evolution.',
        method=template['method'], correction=template['correction'],
        coordinate_edges_seconds=template['coordinate_edges_seconds'],
        refinements=template['refinements'], conduction_time=template['conduction_time'],
        nonlinear_absolute_tolerances=ATOL.astype(float).tolist(), nonlinear_relative_tolerance=0,
        maximum_stage_iterations=24, finite_conservation_gates=template['finite_conservation_gates'],
        initial_state_check=dict(classification='Counterexample candidate', cells=star.n,
            reference_arrays_equal=True, zero_conservative_increments=True,
            budget=initial_budget, sampled_cone=initial_cone, passed=True),
        time_refinement_gate=template['time_refinement_gate'],
        nonlinear_operator=template['nonlinear_operator'], composition_probe=template['composition_probe'],
        iteration_limitations=template['iteration_limitations'],
        imported_prefix_steps=0,
        initialization='Every refinement starts from gr-molecular-conservative-initial/initial.npz after its full-grid restriction and molecular GR-8 gates pass. No old-model native auxiliary arrays or accepted time states are imported.',
        worker_provider='Initialize the molecular native EOS and original opacity tables in every worker. Guard the common material entry point used by both state evaluation and composition probes against any other EOS class. Provider caches belong to the new star and new workers.',
        historical_bindings='Revalidate the frozen parent method and initial-data bindings as transitive provenance. Historical parent control artifacts are provenance, not starting states or new-model operator checks.',
        limitations=template['limitations'], physical_EOS_certified=False,
        nuclear_reactions_included=False, observational_closure=False)
    OUT.mkdir()
    e.write(OUT/'plan.json', plan)
    bindings()
    print('PREPARED MOLECULAR CONSERVATIVE GR EVOLUTION; NO IMPORTED TIME STATES', flush=True)


stage = FunctionType(method.stage.__code__, globals())
_native_run = FunctionType(prior.run.__code__, globals())
completed = FunctionType(prior.completed.__code__, globals())
compare = FunctionType(prior.compare.__code__, globals())


def run(refinement, workers):
    global REFINEMENT, PREFIX
    REFINEMENT, PREFIX = refinement, 0
    require_inputs()
    return _native_run(refinement, workers)


def wrong_provider_rejected(row):
    selected = original_e.EOS
    try:
        original_e.EOS = original_e.g.EOS()
        try:
            parent.parent.material(row)
        except AssertionError as error:
            assert str(error) == 'Wrong EOS in a molecular evolution worker'
            return True
        return False
    finally:
        original_e.EOS = selected


def selfcheck():
    """Use the actual process-pool/material route, including a wrong-EOS control."""
    data, _, _ = initial.molecular.inputs()
    cells = [0, 1175, 1176, 2972, 5734]
    rows = [(data['lnd'][i], data['lnT'][i], data['X'][i]) for i in cells]
    native = initial.molecular.model.EOS()
    expected = [native(2, float(lr), float(lt), x) for lr, lt, x in rows]
    with ProcessPoolExecutor(max_workers=1, initializer=worker_init) as pool:
        actual = list(pool.map(parent.parent.material, rows))
        for i, a, b in zip(cells, actual, expected):
            assert a.shape == (30,) and np.array_equal(a[:21], b), i
        assert pool.submit(wrong_provider_rejected, rows[0]).result()
        restored = pool.submit(parent.parent.material, rows[0]).result()
        assert np.array_equal(restored, actual[0])
    assert _native_run.__globals__['e'].worker_init is worker_init
    assert stage.__globals__['PREFIX'] == 0
    assert prior.symbolic()['passed']
    print('PASS actual molecular worker/material route, five native vectors, wrong-provider rejection and restoration, and conservative symbolic checks', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['selfcheck', 'prepare', 'run', 'completed', 'compare', 'chain'])
    parser.add_argument('--refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--workers', type=int, default=15)
    args = parser.parse_args()
    assert 1 <= args.workers <= 15
    if args.command == 'run':
        run(args.refinement, args.workers)
    elif args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1], indent=2))
    elif args.command == 'chain':
        for refinement in [1, 2, 4]:
            run(refinement, args.workers)
        compare()
    else:
        globals()[args.command]()
