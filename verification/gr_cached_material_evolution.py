"""Counterexample candidate: exact native cache and reuse of completed paths.

The failed 8/16/32 verdict remains frozen. A new 16/32/64 comparison copies
and revalidates the identical 16/32 histories and computes only the 64 path.
"""
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType
import argparse
import json
import os
import shutil
import time

import numpy as np
import gr_compatible_equilibrium_evolution as original
from gr_compatible_equilibrium_evolution import *

PARENT=BASE/'pilot'
CHECK=BASE/'cache-check'
OUT=BASE/'pilot-cached-16-32-64'


class CachedMaterialPool:
    def __init__(self,pool,limit=8*5735):
        self.pool,self.limit=pool,limit
        self.cache=OrderedDict()
        self.evaluations=self.hits=0

    def map(self,function,rows,**kwargs):
        assert function is original.parent.parent.material
        rows=list(rows)
        keys=[e.material_key(row) for row in rows]
        missing={key:row for key,row in zip(keys,rows) if key not in self.cache}
        self.evaluations+=len(missing)
        self.hits+=len(keys)-len(missing)
        for key,value in zip(missing,self.pool.map(function,missing.values(),**kwargs)):
            self.cache[key]=np.asarray(value).copy()
        result=[]
        for key in keys:
            result.append(self.cache[key].copy())
            self.cache.move_to_end(key)
        while len(self.cache)>self.limit:
            self.cache.popitem(last=False)
        return result


def initialize(pool):
    return original.initialize(None if pool is None else CachedMaterialPool(pool))


bindings=FunctionType(prior.bindings.__code__,globals())
run=FunctionType(prior.run.__code__,globals())
completed=FunctionType(prior.completed.__code__,globals())
compare=FunctionType(prior.compare.__code__,globals())


def check(workers=15):
    from gr_execution_acceleration import state
    assert not CHECK.exists()
    CHECK.mkdir()
    original.OUT=PARENT
    plan=original.bindings()
    times=parent.wall.prior.time_nodes(plan,2)
    with ProcessPoolExecutor(max_workers=workers,initializer=e.worker_init) as pool:
        cached=CachedMaterialPool(pool)
        star=original.initialize(cached)
        older=state(star,PARENT/'path-2',6)
        original.composition_response(star,*older)
        prior_cache=cached.cache.copy()
        previous=state(star,PARENT/'path-2',7)
        cached.evaluations=cached.hits=0
        began=time.perf_counter()
        reused=original.composition_response(star,*previous)
        cached_seconds=time.perf_counter()-began
        counts=dict(native_evaluations=cached.evaluations,cache_hits=cached.hits)
        star.pool=pool
        began=time.perf_counter()
        fresh=original.composition_response(star,*previous)
        fresh_seconds=time.perf_counter()-began
        assert np.array_equal(reused['coefficients'],fresh['coefficients'])
        assert np.array_equal(reused['inverse'],fresh['inverse'])
        # Restore only preceding-stage probes: do not time an artificially
        # precomputed current stage as a real steady-state speed improvement.
        star.pool=cached
        cached.cache=prior_cache
        cached.evaluations=cached.hits=0
        records=[]
        def log(row):
            records.append(row)
            print('CACHED NATIVE STAGE',row['iteration'],row['residual_norm'],flush=True)
        h,old_h=times[8]-times[7],times[7]-times[6]
        began=time.perf_counter()
        delta,z=original.stage(star,previous,older,h,weights(h,old_h),log)
        stage_seconds=time.perf_counter()-began
        with np.load(PARENT/'path-2/step-0008.npz') as cp:
            assert np.array_equal(delta,cp['delta'])
            for name in ['m','mf','a','N','Q','aux','dU']:
                assert np.array_equal(z[name],cp[name]),name
        result=dict(classification='Counterexample candidate',passed=True,
            source_sha256=e.digest(Path(__file__)),original_source_sha256=e.digest(Path(original.__file__)),
            counts=counts,composition_seconds=dict(cached=cached_seconds,fresh=fresh_seconds),
            complete_native_stage_seconds=stage_seconds,native_stage_evaluations=cached.evaluations,
            native_stage_hits=cached.hits,final_score=records[-1]['residual_norm'],
            native_coefficients_exact=True,complete_native_saved_state_exact=True,
            scope='One fixed actual middle stage and its preceding probe cache. Warm-cache timing; priming costs are not included. No approximate keys, EOS surrogate, changed arithmetic precision, equations, histories or acceptance tolerances. Not a full-trajectory speedup.')
        e.write(CHECK/'result.json',result)
        print('EXACT CACHE CHECK',json.dumps(result),flush=True)


def prepare():
    assert not OUT.exists()
    checked=json.loads((CHECK/'result.json').read_text())
    assert checked['passed'] and checked['source_sha256']==e.digest(Path(__file__))
    assert checked['original_source_sha256']==e.digest(Path(original.__file__))
    original.OUT=PARENT
    template=original.bindings()
    for relative,digest in json.loads((PARENT/'time-refinement-manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/relative)==digest,relative
    assert not json.loads((PARENT/'time-refinement.json').read_text())['passed']
    # Revalidate both histories against their original native equations first.
    imported=[original.completed(r)[1] for r in [2,4]]
    edges=parent.wall.prior.time_nodes(template,2)
    plan=dict(template,phase='additional native temporal refinement',
        coordinate_edges_seconds=[str(t) for t in edges],
        candidate='Same v2 spatial equations with exact native material memoization; additional 16/32/64 time comparison.',
        imported_paths=[dict(new_refinement=r,source_refinement=2*r) for r in [1,2]],
        imported_path_validation=imported,
        original_failed_pilot_verdict_preserved=True,
        reuse_scope='The completed 16/32-step paths solve these identical space/time equations. Copy them unchanged and compute only 64 steps. The earlier different-spatial-operator 39/78/156 paths are not reused as accepted states.',
        cache_scope='In-memory LRU, at most 45880 exact binary64 native argument tuples including all composition and opacity arguments. Actual EOS/opacity functions unchanged; copies returned to protect cached values.',
        new_time_refinement_steps=[16,32,64],additional_steps_to_compute=64,
        predeclared_new_gate=template['time_refinement_gate'])
    plan['bindings']=dict(template['bindings'])
    for p in [Path(__file__),CHECK/'result.json',PARENT/'time-refinement-manifest.json']:
        plan['bindings'][p.relative_to(e.ROOT).as_posix()]=e.digest(p)
    OUT.mkdir()
    e.write(OUT/'plan.json',plan)
    for r in [1,2]:
        source=PARENT/f'path-{2*r}'
        manifest=json.loads((source/'manifest.json').read_text())
        target=OUT/f'path-{r}'
        target.mkdir()
        for name,digest in manifest['sha256'].items():
            assert e.digest(source/name)==digest,name
            shutil.copy2(source/name,target/name)
        manifest.update(plan_sha256=e.digest(OUT/'plan.json'),
            imported_from=source.relative_to(e.ROOT).as_posix(),source_manifest_sha256=e.digest(source/'manifest.json'))
        e.write(target/'manifest.json',manifest)
    bindings()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=['check','prepare','chain','compare'])
    parser.add_argument('--workers',type=int,default=15)
    args=parser.parse_args()
    assert set(os.sched_getaffinity(0))<=set(range(16)) and 1<=args.workers<=16
    if args.command=='check':check(args.workers)
    elif args.command=='prepare':prepare()
    elif args.command=='chain':
        run(4,args.workers)
        compare()
    else:compare()
