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
        m.guide_g=g.copy()
        state=prior.prior.state(m);state['material_floor_discard_scaled']=m.floor_discard.copy();state['restart_guide']=m.guide_g.copy()
        angular_t=list(m.angular_times);angular=list(m.angular)
        frame=sys._getframe(1).f_locals
        keys=['x','g','ledger','escape','impulse','ports','transfer','times','records','port_history','photon_history','gas_history','transfer_history','stage_weights','actual_edges']
        payload={k:frame[k] for k in keys};payload.update(state)
        np.savez_compressed(OUT/f'last-accepted-{n}.npz',**payload,angular_times=angular_t,angular=angular,next_time=t,next_step=h,macro_index=frame['k'],sub_index=frame['sub'])
        try:
            result=actual(m,t,h,x,g,lus)
            write(OUT/f'stage-progress-{n}.json',dict(last_stage_time=float(t+h),solved_actual_steps=len(m.newton_iterations),
                original_stage_relative=m.newton_iterations[-1][-1]['relative'],physical_stage_relative=m.newton_iterations[-1][-1]['moments']))
            return result
        except BaseException as exc:
            write(OUT/f'failure-{n}.json',dict(classification='Counterexample candidate',error=repr(exc),time=float(t),step=float(h),
                actual_accepted_steps=len(frame['actual_edges'])-1,macro_index=frame['k'],sub_index=frame['sub'],
                last_accepted_inputs_saved=True,serialized_substep_restart_validated=False,final_charge_conclusion='unadjudicated'))
            raise
    run=owner.Model.run;owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=saved_stage),argdefs=run.__defaults__)
    m=owner.Model(n);row=m.run(n,f'complete-{n}',n,restart=f'interval-15-{n}')
    report,_,_,_=prior.verify(OUT/f'sweep-1/photons/complete-{n}.npz',OUT/f'sweep-1/photons/interval-15-{n}.npz')
    row.update(audit=report,physical_equation_and_gates_unchanged=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated')
    write(OUT/f'path-{n}.json',row);assert row['passed']
