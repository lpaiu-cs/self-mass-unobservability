"""Counterexample candidate: continue native GR with a smaller future time grid."""
import argparse
import json
import os
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_compatible_full_duration as full
import gr_step42_subdivision as split

m,e,ld=full.cached,full.e,full.ld
OUT=m.BASE/'production-subdivided'
context=dict(full.context,OUT=OUT)
bindings=FunctionType(m.prior.bindings.__code__,context)
context['bindings']=bindings
completed=FunctionType(m.prior.completed.__code__,context)
context['completed']=completed
local=dict(vars(full),OUT=OUT,bindings=bindings,context=context)
check=FunctionType(full.check.__code__,local)
run=FunctionType(full.run.__code__,local)
common_times=FunctionType(full.common_times.__code__,local)
compare=FunctionType(full.compare.__code__,dict(local,common_times=common_times))


def prepare():
    assert not OUT.exists()
    old=full.bindings()
    result=json.loads((split.OUT/'result.json').read_text())
    assert result['passed'] and [r['step'] for r in result['rows']]==[42,43]
    for rel,digest in json.loads((split.OUT/'manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel)==digest,rel
    original=full.time_nodes(old,1)
    tail=np.concatenate([np.linspace(a,b,3,dtype=ld)[1:] for a,b in zip(original[41:-1],original[42:])])
    coarse=np.r_[original[:42],tail]
    assert len(coarse)==100 and coarse[-1]==original[-1]
    grids={'1':list(map(str,coarse))}
    for r in (2,4):
        prefix=full.time_nodes(old,r)[:32*r+1]
        tail=np.concatenate([np.linspace(a,b,r+1,dtype=ld)[1:] for a,b in zip(coarse[32:-1],coarse[33:])])
        times=np.r_[prefix,tail]
        assert len(times)==99*r+1 and times[-1]==coarse[-1]
        assert np.max(np.diff(times)[1:]/np.diff(times)[:-1])<=ld('2.000000000001')
        grids[str(r)]=list(map(str,times))
    OUT.mkdir()
    prefix=OUT/'accepted-prefix-1'
    prefix.mkdir()
    for step in range(44):
        source=(full.OUT/'path-1' if step<=41 else split.OUT)/f'step-{step:04d}.npz'
        with np.load(source) as cp:
            assert cp['time_seconds']==coarse[step],step
        os.link(source,prefix/source.name)
    rows=[json.loads(s) for s in (full.OUT/'path-1/iterations.jsonl').read_text().splitlines() if json.loads(s)['step']<=41]
    rows += [json.loads(s) for s in (split.OUT/'iterations.jsonl').read_text().splitlines()]
    for step in range(1,44):
        records=[r for r in rows if r['step']==step]
        assert [r['iteration'] for r in records]==list(range(len(records)))
        assert len(records)<=24 and records[-1]['residual_norm']<=1
    (prefix/'iterations.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    e.write(prefix/'manifest.json',dict(classification='Counterexample candidate',accepted_prefix_steps=43,
        source_plan_sha256=e.digest(full.OUT/'plan.json'),subdivision_manifest_sha256=e.digest(split.OUT/'manifest.json'),
        sha256={p.name:e.digest(p) for p in prefix.iterdir() if p.is_file()}))
    plan=dict(old,phase='Full original duration with halved future intervals after step-42 failure',
        candidate='Same native GR solver and equations; preserve accepted coarse history then halve each remaining original interval.',
        coordinate_edges_seconds=grids['1'],exact_time_nodes_seconds=grids,new_time_refinement_steps=[99,198,396],
        imported_prefix='Native replay of accepted 43/64/128 prefixes, preserving both BDF states and both accumulated flux histories.',
        continuation_grid='Retain every executed pilot clock and original coarse steps 0..41. Halve each remaining original coarse interval, then split new intervals by 1/2/4.',
        additional_steps_to_compute=458,imported_prefix_steps_by_refinement={'1':43,'2':64,'4':128},
        failed_full_plan= (full.OUT/'plan.json').relative_to(e.ROOT).as_posix(),
        maximum_future_coarse_step_seconds=str(np.diff(coarse)[41:].max()),
        original_native_solver_unchanged=True,failed_iteration_repairs_used=False)
    plan['imported_paths']=dict(old['imported_paths'])
    plan['imported_paths']['1']=dict(steps=43,directory=prefix.relative_to(e.ROOT).as_posix(),refinement=1,
        source_plans=[(full.OUT/'plan.json').relative_to(e.ROOT).as_posix(),(split.OUT/'plan.json').relative_to(e.ROOT).as_posix()])
    plan['bindings']=dict(old['bindings'])
    for p in [Path(__file__),split.OUT/'manifest.json',prefix/'manifest.json',full.OUT/'path-1/failure.json']:
        plan['bindings'][p.relative_to(e.ROOT).as_posix()]=e.digest(p)
    e.write(OUT/'plan.json',plan)
    check()
    print('SUBDIVIDED FULL DURATION PREPARED',plan['new_time_refinement_steps'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=['prepare','chain','compare'])
    parser.add_argument('--workers',type=int,default=15)
    args=parser.parse_args()
    assert set(os.sched_getaffinity(0))<=set(range(16)) and 1<=args.workers<=16
    if args.command=='prepare':
        prepare()
    elif args.command=='chain':
        assert json.loads((OUT/'restart-check.json').read_text())['passed']
        for r in (1,2,4):
            run(r,args.workers)
        compare()
    else:
        compare()
