"""Resume the interrupted fixed EOS population in a separate output directory.

No new states, tolerances or approximation are introduced. Completed native
results are verified in place; the four interrupted original folders stay intact.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import FunctionType
import json
import sys

import gr_outer_product_adaptive_table as previous

ROOT, sha = previous.ROOT, previous.sha
OUT = previous.OUT.parent/'gr-outer-product-recovery'
CACHE = previous.CACHE.parent/'outer-product-recovery'


def save(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2)+'\n')


def prepare():
    from gr_outer_product_transition import snapshot
    assert not OUT.exists() and not CACHE.exists()
    assert not any(len(r['command']) > 1 and Path(r['command'][0]).name == 'python3'
                   and r['command'][1] == 'verification/gr_outer_product_adaptive_table.py'
                   for r in snapshot().values())
    plan = previous.bindings()
    imported = previous.check_inventory()
    completed = {r['position']: dict(r, origin='original') for r in imported['records']}
    interrupted = []
    files = [Path(__file__), previous.OUT/'plan.json', previous.OUT/'inventory.json', previous.OUT/'progress.json']
    for folder in sorted(previous.OUT.glob('cell-*')):
        position = int(folder.name.split('-')[1])
        if (folder/'manifest.json').exists():
            result = previous.worker().verify_cell(position)
            assert result['passed'] and position not in completed
            completed[position] = dict(position=position, origin='adaptive',
                                        path=(folder/'manifest.json').relative_to(ROOT).as_posix())
            files.append(folder/'manifest.json')
        else:
            interrupted.append(position)
            files.extend(p for p in folder.iterdir() if p.is_file())
    assert set(completed) <= set(plan['positions'])
    assert not json.loads((previous.OUT/'progress.json').read_text())['failures']
    OUT.mkdir()
    CACHE.mkdir()
    (OUT/'candidate.py').write_bytes((previous.OUT/'candidate.py').read_bytes())
    files.append(OUT/'candidate.py')
    plan.update(bindings=dict(plan['bindings'], **{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        recovery='Old process identities no longer exist. Retain all prior files. Verify and reuse completed original/adaptive results in place, compute only the remaining fixed positions with the identical adaptive evaluator in a new output/cache directory.',
        completed_records=list(completed.values()), interrupted_positions_preserved=interrupted,
        remaining_positions=sorted(set(plan['positions'])-set(completed)), processes=4)
    save('plan.json', plan)
    print('EOS RECOVERY PREPARED', len(completed), 'reused;', len(plan['remaining_positions']), 'remaining', flush=True)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel, h in plan['bindings'].items(): assert sha(ROOT/rel) == h, rel
    return plan


def worker():
    return FunctionType(previous.worker.__code__, dict(vars(previous), OUT=OUT, CACHE=CACHE, bindings=bindings))()


def evaluate(position): return worker().evaluate(position)


def run():
    plan = bindings()
    done, failures = [], []
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        jobs = {pool.submit(evaluate, i):i for i in plan['remaining_positions']}
        for future in as_completed(jobs):
            i = jobs[future]
            try:
                assert future.result()['passed']
                done.append(i)
                print('RECOVERED EOS', i, len(done), flush=True)
            except Exception as error:
                failures.append(dict(position=i, error=repr(error)))
                print('RECOVERED EOS FAILURE', i, repr(error), flush=True)
            save('progress.json', dict(classification='Counterexample candidate',
                 completed_new_positions=sorted(done), failures=failures))
    assert not failures


def all_results():
    plan = bindings()
    old = {r['position']:r for r in plan['completed_records']}
    results = []
    for i in plan['positions']:
        record = old.get(i)
        if record is None: result = worker().verify_cell(i)
        elif record['origin'] == 'adaptive': result = previous.worker().verify_cell(i)
        elif record['kind'] == 'original_pilot': result = json.loads((ROOT/record['path']).read_text())
        else: result = previous.original.verify_cell(i)
        assert result['position'] == i and result['passed']
        results.append(result)
    assert len({r['cell'] for r in results}) == 3206
    return results


def finalize():
    plan = bindings()
    previous.original.H_table.verify()
    rows = all_results()
    save('result.json', dict(classification='Proven', passed=True, cells=3206, rows=rows,
        maximum_field_error_upper=str(max(previous.F(e) for r in rows for e in r['field_error_upper'])),
        physical_EOS_certified=False))
    files = [OUT/'plan.json', OUT/'result.json']+list(OUT.glob('cell-*/manifest.json'))
    save('manifest.json', dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}))
    print('RECOVERED FIXED EOS POPULATION COMPLETE', len(rows), flush=True)


if __name__ == '__main__': globals()[sys.argv[1]]()
