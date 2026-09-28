"""Counterexample candidate: resume frozen calculations after the WSL restart.

Keep the interrupted directories unchanged. Replay accepted GR states through
the original native residual and budget code; recompute only unfinished work.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import os
import shutil

import numpy as np
import gr_conservative_composition_tangent as gr
import gr_molecular_conservative_initial as initial

ROOT = gr.e.ROOT
GR_SOURCE, INITIAL_SOURCE = gr.OUT, initial.OUT
GR_OUT = GR_SOURCE.with_name('gr-conservative-composition-recovery-20260914')
INITIAL_OUT = INITIAL_SOURCE.with_name('gr-molecular-initial-recovery-20260914')
NATIVE_STAGE, NATIVE_SAVE = gr.stage, gr.save_state
STATE_KEYS = ['m', 'mf', 'a', 'N', 'Q', 'aux', 'dU']


def copy_files(source, target, paths):
    """Only immutable arrays share disk storage; reports are independent copies."""
    for path in paths:
        destination = target/path.relative_to(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        assert not destination.exists(), destination
        if path.suffix == '.npz':
            os.link(path, destination)
        else:
            shutil.copy2(path, destination)


def receipt(target, source, inputs, **fields):
    gr.e.write(target/'recovery.json', dict(
        classification='Counterexample candidate', source=str(source.relative_to(ROOT)),
        source_sha256={str(p.relative_to(ROOT)):gr.e.digest(p) for p in inputs},
        recovery_source_sha256=gr.e.digest(Path(__file__)),
        prior_boot_id='9b968b4f-ecbd-483b-ada9-c1c871b6ab20',
        preparation_boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
        reason='The original process handles are absent, the WSL boot identifier changed, and no original scientific process remains.',
        physical_EOS_certified=False, continuous_errors_certified=False, **fields))


def verify_receipt(target):
    report = json.loads((target/'recovery.json').read_text())
    assert report['recovery_source_sha256'] == gr.e.digest(Path(__file__))
    for rel, digest in report['source_sha256'].items():
        assert gr.e.digest(ROOT/rel) == digest, rel
    return report


def iteration_rows():
    rows = {}
    for line in (GR_SOURCE/'path-4/iterations.jsonl').read_text().splitlines():
        row = json.loads(line)
        rows.setdefault(row['step'], []).append(row)
    return rows


def same_arrays(saved, arrays):
    assert set(saved.files) == set(arrays)
    for key, value in arrays.items():
        assert np.array_equal(saved[key], value), key


def restored_state(star, step):
    with np.load(GR_SOURCE/f'path-4/step-{step:04d}.npz') as saved:
        delta = saved['delta'].copy()
        y = star.base+delta
        star.material_cache = {gr.e.material_key(row):aux for row, aux in zip(
            zip(y[:, 0], y[:, 1], y[:, 5:]), saved['aux'])}
        state = star.evaluate(delta)
        for key in STATE_KEYS:
            assert np.array_equal(state[key], saved[key]), (step, key)
    return delta, state


def prepare_gr():
    gr.bindings()
    step = json.loads((GR_SOURCE/'path-4/progress.json').read_text())['step']
    assert step == 90 and not (GR_SOURCE/'path-4/failure.json').exists()
    rows = iteration_rows()
    for n in range(1, step+1):
        assert [r['iteration'] for r in rows[n]] == list(range(len(rows[n])))
        assert len(rows[n]) <= 24 and np.isfinite(rows[n][-1]['residual_norm'])
        assert rows[n][-1]['residual_norm'] <= 1
    assert not GR_OUT.exists()
    GR_OUT.mkdir()
    copied = [GR_SOURCE/'plan.json']
    for refinement in [1, 2]:
        folder = GR_SOURCE/f'path-{refinement}'
        manifest = json.loads((folder/'manifest.json').read_text())
        assert json.loads((folder/'result.json').read_text())['completed']
        for name, digest in manifest['sha256'].items():
            assert gr.e.digest(folder/name) == digest, name
        copied += [p for p in folder.iterdir() if p.is_file()]
    copy_files(GR_SOURCE, GR_OUT, copied)
    inputs = copied+[GR_SOURCE/f'path-4/step-{n:04d}.npz' for n in range(step+1)]
    inputs += [GR_SOURCE/'path-4/iterations.jsonl', GR_SOURCE/'path-4/progress.json']
    receipt(GR_OUT, GR_SOURCE, inputs, accepted_steps=step,
        replay='Replay states 0 through 90 with cached native auxiliaries, exact state arrays, all 31 native BDF residuals, original flux accumulation, conservation and cone gates. Import only accepted iteration records. Stage 91 starts from accepted step 90, not its interrupted iterate. The full original 1/2/4 comparison runs after step 156.')


def replay_stage(star, previous, older, h, coefficients, log):
    step = star.operator_step+1
    if step > RECOVERY['accepted_steps']:
        return NATIVE_STAGE(star, previous, older, h, coefficients, log)
    star.operator_step = step
    star.operator_time += h
    delta, state = restored_state(star, step)
    with np.load(GR_SOURCE/f'path-4/step-{step:04d}.npz') as saved:
        assert saved['time_seconds'] == star.operator_time
    native, state = gr.residual(star, delta, previous, older, h, coefficients)
    assert np.max(abs(native)/gr.ATOL) <= 1, step
    for row in ITERATIONS[step]:
        log({key:value for key, value in row.items() if key != 'step'})
    return delta, state


def replay_save(path, delta, state, t, energy, baryon, energy_defect, baryon_defect):
    step = int(path.stem.split('-')[-1]) if path.stem.startswith('step-') else None
    if step is not None and step <= RECOVERY['accepted_steps']:
        source = GR_SOURCE/path.parent.name/path.name
        with np.load(source) as saved:
            same_arrays(saved, dict(delta=delta, **{k:state[k] for k in STATE_KEYS},
                time_seconds=t, integrated_energy_flux=energy, integrated_baryon_flux=baryon,
                normalized_energy_defect=energy_defect, normalized_baryon_defect=baryon_defect))
        assert not path.exists()
        os.link(source, path)
    else:
        NATIVE_SAVE(path, delta, state, t, energy, baryon, energy_defect, baryon_defect)


def run_gr():
    global RECOVERY, ITERATIONS
    RECOVERY, ITERATIONS = verify_receipt(GR_OUT), iteration_rows()
    assert not (GR_OUT/'path-4').exists()
    gr.OUT, gr.stage, gr.save_state = GR_OUT, replay_stage, replay_save
    gr.run(4, 15)
    gr.compare()


def prepare_initial():
    plan = initial.bindings()
    assert not (INITIAL_SOURCE/'failure.json').exists()
    done = json.loads((INITIAL_SOURCE/'progress.json').read_text())['completed_cells']
    assert done == 1792 and not INITIAL_OUT.exists()
    files = [INITIAL_SOURCE/name for name in ['plan.json', 'reference-state.npz', 'path-4.npz', 'branches.json']]
    for start in range(0, done, plan['block_size']):
        label = f'block-{start:04d}'
        report = json.loads((INITIAL_SOURCE/(label+'.json')).read_text())
        assert report['all_passed'] and [r['cell'] for r in report['rows']] == list(range(start, start+plan['block_size']))
        files += [INITIAL_SOURCE/(label+'.json')]
        files += [INITIAL_SOURCE/f'{label}-nodes-{n}.npz' for n in plan['nodes']]
        for prefix in ['node-roots-', 'geometry-roots-']:
            files += [p for p in (INITIAL_SOURCE/(prefix+label)).rglob('*') if p.is_file()]
    INITIAL_OUT.mkdir()
    copy_files(INITIAL_SOURCE, INITIAL_OUT, files)
    receipt(INITIAL_OUT, INITIAL_SOURCE, files, completed_cells=done,
        resume='Reuse complete blocks and their native roots. Recompute the entire interrupted 1792-1919 block in this separate directory; preserve all original partial roots. Original molecular block, assembly, nodes, tolerances and scientific gates are unchanged.')
    gr.e.write(INITIAL_OUT/'progress.json', dict(classification='Counterexample candidate', completed_cells=done))


def run_initial():
    report = verify_receipt(INITIAL_OUT)
    initial.OUT = INITIAL_OUT
    plan = initial.bindings()
    with (INITIAL_OUT/'run-started.json').open('x') as stream:
        json.dump(dict(pid=os.getpid(), boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip()), stream)
    done = report['completed_cells']
    jobs = [(f'block-{n:04d}', list(range(n, min(n+plan['block_size'], plan['cells']))))
            for n in range(done, plan['cells'], plan['block_size'])]
    try:
        with ProcessPoolExecutor(max_workers=plan['workers'], initializer=initial.worker_init) as pool:
            for result in pool.map(initial.block, jobs):
                done += len(result['rows'])
                initial.save('progress.json', dict(classification='Counterexample candidate', completed_cells=done))
            initial.assemble(pool)
    except Exception as error:
        initial.save('failure.json', dict(classification='Counterexample candidate', error=repr(error)))
        raise


def selfcheck():
    gr.bindings()
    initial.bindings()
    assert gr.prior.symbolic()['passed']
    star = gr.prior.initialize(None)
    older, previous, (delta, state) = [restored_state(star, n) for n in [88, 89, 90]]
    times = gr.parent.wall.prior.time_nodes(gr.bindings(), 4)
    h = times[90]-times[89]
    native, _ = gr.residual(star, delta, previous, older, h, gr.weights(h, times[89]-times[88]))
    norm = float(np.max(abs(native)/gr.ATOL))
    assert norm <= 1
    with np.load(GR_SOURCE/'path-4/step-0090.npz') as saved:
        arrays = {k:saved[k].copy() for k in saved.files}
        same_arrays(saved, arrays)
        arrays['integrated_energy_flux'][1] = np.nextafter(
            arrays['integrated_energy_flux'][1], np.longdouble(np.inf))
        try:
            same_arrays(saved, arrays)
        except AssertionError:
            pass
        else:
            raise AssertionError('Corrupt flux history was not rejected')
    print('PASS native step 90 replay, altered-history rejection and symbolic checks', norm, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['selfcheck', 'prepare_gr', 'run_gr', 'prepare_initial', 'run_initial'])
    globals()[parser.parse_args().action]()
