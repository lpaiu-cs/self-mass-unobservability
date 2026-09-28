"""Resume cancelled blocks separately; retain original failures and exact coverage."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import FunctionType
import json, shutil, sys, traceback
import numpy as np
import gr_subcell_precision_fallback as fallback

original=fallback.original;g=fallback.g;OUT=g.OUT/'gr-subcell-precision-resume'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();fallback.bindings()
    completed=sorted(p for p in original.OUT.glob('block-*.json') if '-failure' not in p.name)
    starts=[int(p.stem.split('-')[1]) for p in completed]
    assert starts==[i for i in range(0,3968,128) if i!=2688],starts
    paths=[original.OUT/'plan.json',original.OUT/'block-2688-failure.json',
        fallback.OUT/'plan.json',fallback.OUT/'candidate-block.py',
        g.ROOT/'verification/gr_subcell_precision_resume.py',g.ROOT/'verification/gr_subcell_precision_fallback.py']
    for path in completed:
        row=json.loads(path.read_text());assert row['all_passed'] and len(row['rows'])==128
        paths.extend([path]+[original.OUT/f'{path.stem}-nodes-{n}.npz' for n in [8,16]])
    shutil.copy2(g.ROOT/'outputs/gr-full-subcell33-run.log',OUT/'original-run-failure.log')
    paths.append(OUT/'original-run-failure.log')
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='cf062eb',retained_starts=starts,remaining_starts=list(range(3968,5735,128)),
        intervention='Continue only the 14 unstarted blocks cancelled by the original process-pool failure. Reuse the separately frozen precision fallback algorithm at unchanged physical formulas, geometry, nodes, root and quadrature gates. Per-block directories isolate failure records. Collect every worker result even if one fails.',
        aggregate='An explicit new mixed-precision candidate combines the 30 original successful blocks, separately recomputed block 2688, and these 14 new blocks. The original full-grid run remains failed. No continuous evaluator-switch, native, physical EOS or evolution error certificate.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    save('plan.json',plan)


def bindings():
    return FunctionType(original.bindings.__code__,dict(original.bindings.__globals__,OUT=OUT))()


def block(start):
    plan=bindings();folder=OUT/f'block-{start:04}';folder.mkdir()
    local_save=FunctionType(save.__code__,dict(save.__globals__,OUT=folder))
    inverse=FunctionType(fallback.make_inverse.__code__,dict(fallback.make_inverse.__globals__,OUT=folder,save=local_save))
    namespace=dict(original.block.__globals__,OUT=folder,save=local_save,bindings=bindings,make_inverse=inverse)
    source=fallback.OUT/'candidate-block.py'
    exec(compile(source.read_text(),str(source),'exec'),namespace)
    return namespace['block']((folder.name,list(range(start,min(start+plan['block_size'],plan['cells'])))))


def run():
    plan=bindings();assert not any(OUT.glob('block-*'))
    records=[];errors=[]
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        futures={pool.submit(block,i):i for i in plan['remaining_starts']}
        for future in as_completed(futures):
            start=futures[future]
            try:records.append(future.result())
            except Exception:
                error=dict(start=start,traceback=traceback.format_exc());errors.append(error)
                save(f'worker-{start:04}-failure.json',error)
            save('progress.json',dict(classification='Counterexample candidate',
                completed_cells=sum(len(r['rows']) for r in records),worker_errors=errors))
    save('run-result.json',dict(classification='Counterexample candidate',completed=True,
        completed_cells=sum(len(r['rows']) for r in records),worker_errors=errors,
        finite_passed=not errors and all(r['all_passed'] for r in records)))


def aggregate():
    plan=bindings();fallback.verify();assert not (OUT/'manifest.json').exists()
    assert json.loads((OUT/'run-result.json').read_text())['completed']
    paths=[original.OUT/f'block-{i:04}.json' for i in plan['retained_starts']]
    paths+=[fallback.OUT/'block-2688.json']
    paths+=[OUT/f'block-{i:04}'/f'block-{i:04}.json' for i in plan['remaining_starts']]
    paths.sort(key=lambda p:int(p.stem.split('-')[1]));rows=[];sources=[];bound=[]
    for path in paths:
        part=json.loads(path.read_text());cells=[r['cell'] for r in part['rows']];rows.extend(part['rows'])
        sources.append(dict(path=path.relative_to(g.ROOT).as_posix(),cells=cells,
            accepted_root_maximum_score=part['entropy_roots']['maximum_score'],
            extended_fallbacks=part['entropy_roots'].get('extended_fallbacks',0)))
        assert part['entropy_roots']['maximum_score']<=1
        bound.append(path)
        for number in plan['nodes']:
            nodepath=path.with_name(path.stem+f'-nodes-{number}.npz');a=dict(np.load(nodepath))
            assert list(a['cells'])==cells and np.all(a['coordinate_weights_cm3']>0)
            assert np.all(a['proper_weights_cm3']>0) and np.all(a['baryon_weights_g']>0)
            bound.append(nodepath)
    assert [r['cell'] for r in rows]==list(range(plan['cells']))
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(rows),sources=sources,
        all_passed=all(r['passed'] for r in rows),failed_cells=[r for r in rows if not r['passed']],
        maximum_finite_quadrature_difference=max(r['finite_quadrature_relative_difference'] for r in rows),
        maximum_shell_mass_difference=max(abs(q['native_shell_mass_relative_difference']) for r in rows for q in r['quadratures']),
        maximum_volume_difference=max(abs(q['native_coordinate_volume_relative_difference']) for r in rows for q in r['quadratures']),
        original_full_grid_run_passed=False,continuous_errors_certified=False,physical_EOS_certified=False,
        full_GR_evolution=False))
    bound.extend(p for p in OUT.rglob('*') if p.is_file())
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in bound}))
    verify()


def verify():
    plan=bindings();fallback.verify()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text())
    assert result['completed'] and result['cells']==plan['cells']
    assert [i for row in result['sources'] for i in row['cells']]==list(range(plan['cells']))
    print('PASS complete mixed-precision subcell coverage and bindings; original run failure retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
