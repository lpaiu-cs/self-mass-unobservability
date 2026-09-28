"""Counterexample candidate: time a saved native stage and reuse only guesses.

No saved trajectory becomes accepted history. Both arms solve the same new
31-equation stage, from the same two previously accepted states and tolerances.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType
import inspect
import json
import os
import time

import numpy as np
import gr_compatible_equilibrium_evolution as model

OUT=model.BASE/'acceleration'


class TimedPool:
    def __init__(self,pool):
        self.pool=pool
        self.seconds=0.
        self.calls=0

    def map(self,function,rows,**kwargs):
        rows=list(rows)
        began=time.perf_counter()
        result=list(self.pool.map(function,rows,**kwargs))
        self.seconds+=time.perf_counter()-began
        self.calls+=len(rows)
        return result


def state(star,folder,step):
    with np.load(folder/f'step-{step:04d}.npz') as cp:
        delta=cp['delta'].copy()
        y=star.base+delta
        star.material_cache.update({model.e.material_key(row):aux.copy() for row,aux in
            zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])})
        return delta,star.evaluate(delta)


def run(workers=15):
    assert not OUT.exists()
    OUT.mkdir()
    model.OUT=model.BASE/'pilot'
    plan=model.bindings()
    times=model.parent.wall.prior.time_nodes(plan,2)
    # Fixed middle stage before inspecting any seed scores or wall times.
    step=8
    h,old_h=times[step]-times[step-1],times[step-1]-times[step-2]
    coefficients=model.weights(h,old_h)
    metadata=dict(classification='Counterexample candidate',refinement=2,step=step,
        workers=workers,affinity=sorted(os.sched_getaffinity(0)),same_time_history=True,
        model_sha256=model.e.digest(Path(model.__file__)),source_sha256=model.e.digest(Path(__file__)),
        pilot_plan_sha256=model.e.digest(model.OUT/'plan.json'),production_changed=False)
    model.e.write(OUT/'plan.json',metadata)
    results=[]
    with ProcessPoolExecutor(max_workers=workers,initializer=model.e.worker_init) as raw_pool:
        for arm in ['previous_state','coarse_saved_guess','history_extrapolation']:
            pool=TimedPool(raw_pool)
            star=model.initialize(pool)
            older=state(star,model.OUT/'path-2',step-2)
            previous=state(star,model.OUT/'path-2',step-1)
            if arm=='coarse_saved_guess':
                guess=state(star,model.OUT/'path-1',step//2)[0]
            elif arm=='history_extrapolation':
                guess=previous[0]+h/old_h*(previous[0]-older[0])
            else:
                guess=previous[0].copy()
            # Extrapolating composition can leave its simplex. Species are
            # projected by the existing conservative solve before any use.
            guess[:,5:]=model.species(star,guess,previous,older,h,coefficients,previous[1])
            timings={}
            def timed(name,fn):
                def wrapped(*args,**kwargs):
                    began=time.perf_counter()
                    try:return fn(*args,**kwargs)
                    finally:timings[name]=timings.get(name,0.)+time.perf_counter()-began
                return wrapped
            context=dict(vars(model),REFINEMENT=2,PREFIX=0)
            for name in ['composition_response','jacobian','splu','residual']:
                context[name]=timed(name,getattr(model,name))
            source=inspect.getsource(model.method.stage)
            line='delta, history, factor = previous_state[0].copy(), [], None'
            assert source.count(line)==1
            if arm!='previous_state':
                source=source.replace(line,'delta, history, factor = BENCHMARK_GUESS.copy(), [], None')
            context['BENCHMARK_GUESS']=guess
            exec(compile(source,str(Path(__file__))+'::same_native_stage','exec'),context)
            log=[]
            def record(row):
                log.append(row)
                print('ACCELERATION STAGE',arm,row['iteration'],row['residual_norm'],flush=True)
            pool.seconds=pool.calls=0
            began=time.perf_counter()
            delta,z=context['stage'](star,previous,older,h,coefficients,record)
            elapsed=time.perf_counter()-began
            native,_=model.residual(star,delta,previous,older,h,coefficients)
            assert np.max(abs(native)/model.ATOL)<=1
            saved=np.load(model.OUT/'path-2'/f'step-{step:04d}.npz')['delta']
            result=dict(arm=arm,elapsed_seconds=elapsed,timings=timings,
                native_material_seconds=pool.seconds,native_material_evaluations=pool.calls,
                iteration_count=len(log)-1,initial_score=log[0]['residual_norm'],
                final_score=float(np.max(abs(native)/model.ATOL)),
                maximum_difference_from_saved=np.max(abs(delta[:,:5]-saved[:,:5]),axis=0).astype(float).tolist(),
                notes='One fixed saved stage; timings are not a full-run speedup or an independent time-convergence result.')
            results.append(result)
            np.savez_compressed(OUT/f'{arm}.npz',delta=delta)
            model.e.write(OUT/f'{arm}.json',dict(classification='Counterexample candidate',**result,iterations=log))
            print('ACCELERATION RESULT',json.dumps(result),flush=True)
    model.e.write(OUT/'result.json',dict(**metadata,results=results))
    model.e.write(OUT/'manifest.json',dict(sha256={p.relative_to(model.e.ROOT).as_posix():model.e.digest(p)
        for p in [Path(__file__),*OUT.iterdir()] if p.is_file()}))


if __name__=='__main__':
    assert set(os.sched_getaffinity(0))<=set(range(16))
    run()
