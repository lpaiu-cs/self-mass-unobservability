"""Counterexample candidate: apply the repaired arithmetic to real evolution."""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import stabilize_full_interval_residual as stable

OUT=Path('native-stable-full194-work');ORIGIN=stable.OUT
old=stable.old;prior=old.prior.prior;owner,joint,LD=stable.owner,stable.joint,stable.LD
read,write,sha=stable.read,stable.write,stable.sha
CAPS=dict(prepare=25,coarse=500,fine=800,audit=30)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert read(ORIGIN/'consistent/stable-check-result.json')['passed']
    assert read(ORIGIN/'consistent_check-receipt.json')['error'] is None
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    inputs=list((ORIGIN/'sweep-0').rglob('*.npz'))+[ORIGIN/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [old.OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for src in inputs:
        rel=src.relative_to(ORIGIN if src.is_relative_to(ORIGIN) else old.OLD);dst=OUT/rel
        os.link(src,dst);reused[str(dst)]=sha(src)
    files=inputs+[ORIGIN/'consistent/stable-check-result.json',ORIGIN/'consistent_check-receipt.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='3035a0bdf',
        claim='Resolve the actual terminal coupled interval with consistent baryon arithmetic, retaining the accepted15/16history and all original physical gates.',
        evidence='Identical failing RHS and guess now pass at2.60469e-15 true linear vector and physical moments below6.67e-18. Exact polynomial arithmetic is used on all Krylov/residual inputs; original4refinements remain.',
        method='Change only the linear solver to the checked right-preconditioned and consistent baryon-row evaluation. Original actual nonlinear native stage equation, three Newton proposals, sources, amplitude, fronts, floor map, Radau times, ledger and physical tolerances remain.',
        decision='Continue coarse60/64to64 first. Only a passing actual coarse path admits fine120/128to128. Require exact saved prefixes, actual stage/energy/packet/port gates and all10time channels below2percent before accepting the full period.',
        reuse='No accepted15/16prefix replay. Only the remaining8coarse and16fine substeps. On any stage failure serialize its last accepted inputs and all pre-stage histories, instead of discarding new accepted work.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02,max_Newton=3,max_inner_solves=4),
        forecast=dict(coarse_seconds=[220,500],fine_seconds=[300,800],
            basis='Previous interval8/16steps cost151/274s including setup. The difficult saved system now takes25.33s Krylov plus1.23s exact-row work instead of repeated failed corrections. Later behavior is unmeasured; about2.5x prior path cost plus overhead is allowed.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Original gate,4refinement/3Newton limit,500/800s cap or prefix mismatch. No extra grid,period,solver ladder or automatically repeated failed path.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))


def evolve(n):
    if n==128:assert read(OUT/'path-64.json')['passed']
    stable.OUT=OUT;stable.initialize();logs=[];arithmetic=[];solves=[]
    solver=old.right_solver(logs)
    def solve(m,op,P,rhs,guess):
        first=len(logs);correct=stable.stable_operator(m,op,arithmetic,True)
        try:return solver(m,correct,P,rhs,guess)
        except BaseException as exc:
            tb=exc.__traceback__
            while tb and tb.tb_frame.f_code!=joint.solve.__code__:tb=tb.tb_next
            if tb:
                sol=tb.tb_frame.f_locals['sol'];np.savez_compressed(OUT/f'failed-linear-{n}.npz',rhs=rhs,guess=guess,solution=sol,residual=rhs-correct.matvec(sol))
            raise
        finally:
            solves.append(dict(begin=first,end=len(logs)))
            write(OUT/f'linear-{n}.json',dict(classification='Counterexample candidate',calls=logs,solves=solves,decimal_seconds=sum(arithmetic)))
    stage=owner.Model.run.__globals__['stages'];actual=FunctionType(stage.__code__,dict(stage.__globals__,solve=solve))
    def saved_stage(m,t,h,x,g,lus):
        state=prior.prior.state(m);state['material_floor_discard_scaled']=m.floor_discard.copy();state['restart_guide']=m.guide_g.copy()
        angular_t=list(m.angular_times);angular=list(m.angular)
        try:
            result=actual(m,t,h,x,g,lus)
            write(OUT/f'stage-progress-{n}.json',dict(last_stage_time=float(t+h),solved_actual_steps=len(m.newton_iterations),
                original_stage_relative=m.newton_iterations[-1][-1]['relative'],physical_stage_relative=m.newton_iterations[-1][-1]['moments']))
            return result
        except BaseException as exc:
            frame=sys._getframe(1).f_locals
            keys=['x','g','ledger','escape','impulse','ports','transfer','times','records','port_history','photon_history','gas_history','transfer_history','stage_weights','actual_edges']
            payload={k:frame[k] for k in keys};payload.update(state)
            np.savez_compressed(OUT/f'last-accepted-{n}.npz',**payload,angular_times=angular_t,angular=angular,next_time=t,next_step=h,macro_index=frame['k'],sub_index=frame['sub'])
            write(OUT/f'failure-{n}.json',dict(classification='Counterexample candidate',error=repr(exc),time=float(t),step=float(h),
                actual_accepted_steps=len(frame['actual_edges'])-1,macro_index=frame['k'],sub_index=frame['sub'],
                last_accepted_inputs_saved=True,serialized_substep_restart_validated=False,final_charge_conclusion='unadjudicated'))
            raise
    run=owner.Model.run;owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=saved_stage),argdefs=run.__defaults__)
    m=owner.Model(n);row=m.run(n,f'complete-{n}',n,restart=f'interval-15-{n}')
    report,_,_,_=prior.verify(OUT/f'sweep-1/photons/complete-{n}.npz',OUT/f'sweep-1/photons/interval-15-{n}.npz')
    row.update(audit=report,physical_equation_and_gates_unchanged=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated')
    write(OUT/f'path-{n}.json',row);assert row['passed']


def audit():
    arrays=[];clocks=[];rows=[read(OUT/f'path-{n}.json') for n in [64,128]]
    for n in [64,128]:
        report,ph,gas,clock=prior.verify(OUT/f'sweep-1/photons/complete-{n}.npz',OUT/f'sweep-1/photons/interval-15-{n}.npz')
        arrays.append((ph,gas));clocks.append(clock)
    assert np.array_equal(*clocks)
    relative=joint.previous.run.c.relative;errors=relative(arrays[0][0],arrays[1][0])+relative(arrays[0][1],arrays[1][1])
    assert read(prior.OLD/'time-result.json')['passed']
    for k in range(2,16):assert read(prior.OUT/f'comparison-{k:02d}.json')['passed']
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,rows=rows,time_relative=errors,
        full_horizon_seconds=float(clocks[0][-1]),full_horizon_photon_material_completed=True,accepted15of16preserved=True,
        original_failures_preserved=True,physical_equation_and_gates_unchanged=True,maximum_inner_solves=4,
        self_GR_return_closed=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);assert result['passed']


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                actual=OUT/'prepared-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(actual)==h,p
            write(OUT/'predispatch-amendment.json',dict(classification='Proven',
                scope='Before first evolution, bind the accepted first interval to184time-result.json;185created comparison files02..15 only. Add stage progress serialization without changing equations.',
                prepared_source_sha256=sha(OUT/'prepared-producer.py'),executed_source_sha256=sha(__file__)))
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
