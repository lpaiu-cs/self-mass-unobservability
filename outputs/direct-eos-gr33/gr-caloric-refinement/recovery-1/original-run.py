def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in {**plan['bindings'],**plan['diagnostic_binding']}.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=old.micro.inputs();base=aux['eos'];dm=state['dm'];mass=dm*np.exp(state['nu']);n=len(dm)
    p=dict(np.load(g.OUT/'gr-opacity/new-GR-captured.npz'))['parameters'].copy();model=old.opacity.Opacity()
    duration=json.loads((g.OUT/'gr-nonlinear-thermal/duration.json').read_text())['coordinate_seconds']
    records=[];previous=None;ld=np.longdouble
    with ProcessPoolExecutor(max_workers=plan['processes'],initializer=initialize) as pool:
        def evaluate(shift,delta,number=4):
            pieces=list(pool.map(samples,[(i,shift[i:i+128],delta[i:i+128],number) for i in range(0,n,128)]))
            integral=np.concatenate([r[0] for r in pieces])
            if number==2:return integral
            a=np.concatenate([r[1] for r in pieces]);par=p.copy()
            par[:,4]=np.asarray(state['lnT'].astype(ld)+shift+delta,float)/np.log(10)
            op=np.array([model(row) for row in par]);L,fo,fi=fluxes(state['lnT'],shift+delta,state['nu'],dm,state['radius_faces_m'],state['nu_faces'],op)
            return integral,a,op,L,fo,fi,np.concatenate([r[2] for r in pieces])
        for steps in plan['step_counts']:
            shift=np.zeros(n,dtype=ld);old_a=base.copy();total_energy=np.zeros(n,dtype=ld);exchanged=ld(0);history=[];dt=duration/steps
            for step in range(steps):
                delta=np.zeros(n);evaluated=evaluate(shift,delta)
                for iteration in range(plan['max_Newton_iterations']):
                    U,a,op,L,fo,fi,S=evaluated;div=np.r_[L,ld(0)]-np.r_[ld(0),L]
                    energy=mass.astype(ld)*U;target=ld(dt)*div;residual=energy-target
                    scale=mass.astype(ld)*a[:,10];norm=float(abs(residual/scale).max())
                    exchange=abs(target).sum();global_error=float(abs(residual.sum())/max(ld(1),exchange))
                    print('CALORIC INCREMENT',steps,step,iteration,norm,global_error,flush=True)
                    if norm<=plan['local_energy_residual_scaled_tolerance'] and global_error<=plan['global_energy_relative_to_exchange_tolerance']:break
                    band=np.zeros((3,n));band[1]=np.asarray((scale-dt*(np.r_[fo,ld(0)]-np.r_[ld(0),fi]))/scale,float)
                    band[0,1:]=np.asarray(-dt*fi/scale[:-1],float);band[2,:-1]=np.asarray(dt*fo/scale[1:],float)
                    update=solve_banded((1,1),band,np.asarray(-residual/scale,float))
                    update*=min(1.,plan['temperature_step_cap']/max(float(abs(update).max()),1e-300))
                    # The merit includes the global gate, so convergence of
                    # large local energy scales cannot mask a lost increment.
                    merit=max(norm/plan['local_energy_residual_scaled_tolerance'],global_error/plan['global_energy_relative_to_exchange_tolerance'])
                    for backtrack in range(plan['max_backtracks']):
                        proposed=delta+update*(.5**backtrack);trial=evaluate(shift,proposed)
                        tr=mass.astype(ld)*trial[0]-ld(dt)*(np.r_[trial[3],ld(0)]-np.r_[ld(0),trial[3]])
                        tn=float(abs(tr/(mass.astype(ld)*trial[1][:,10])).max())
                        tg=float(abs(tr.sum())/max(ld(1),abs(ld(dt)*(np.r_[trial[3],ld(0)]-np.r_[ld(0),trial[3]])).sum()))
                        if max(tn/plan['local_energy_residual_scaled_tolerance'],tg/plan['global_energy_relative_to_exchange_tolerance'])<merit:
                            delta=proposed;evaluated=trial;break
                    else:raise AssertionError(('Caloric line search',steps,step,norm,global_error,tn,tg))
                else:raise AssertionError(('Caloric root',steps,step,norm,global_error))
                lower=evaluate(shift,delta,2);qscore=float(np.sum(mass.astype(ld)*abs(U-lower))/max(exchange,ld(1)))
                native=a[:,2].astype(ld)-old_a[:,2].astype(ld)
                budget=32*(np.spacing(abs(a[:,2]))+np.spacing(abs(old_a[:,2]))).astype(ld)
                unresolved=abs(native-U)<=budget
                row=dict(step=step,iterations=iteration,local_scaled_energy_residual=norm,
                    global_energy_relative_to_exchange=global_error,finite_two_four_point_difference_relative_to_exchange=qscore,
                    finite_quadrature_passed=qscore<plan['finite_quadrature_difference_tolerance'],
                    caloric_chart_entropy_change_erg_K=float(dm.astype(ld)@S),
                    native_endpoint_difference_within_32ulp_cells=int(unresolved.sum()),
                    native_endpoint_maximum_difference_outside_32ulp=float(np.maximum(abs(native-U)-budget,0).max()))
                np.savez_compressed(OUT/f'endpoint-{steps}-{step}.npz',old_shift=shift,step_shift=delta,caloric_increment=U,
                    two_point_increment=lower,entropy_increment=S,eos=a,opacity=op,interior_Linf=L,energy_residual=residual)
                history.append(row);save(f'path-{steps}-progress.json',dict(classification='Counterexample candidate',rows=history))
                assert row['finite_quadrature_passed'],row
                assert row['caloric_chart_entropy_change_erg_K']>=0,row
                shift+=delta.astype(ld);total_energy+=energy;exchanged+=exchange;old_a=a.copy()
            record=dict(classification='Counterexample candidate',steps=steps,completed=True,history=history,
                global_increment_energy_relative_to_exchange=float(abs(total_energy.sum())/exchanged),maximum_logT_change=float(abs(shift).max()),
                physical_EOS_certified=False,continuous_caloric_error_certified=False,full_GR_evolution=False)
            if previous is not None:
                error=float(abs(shift-previous).max());record.update(time_refinement_logT_difference=error,
                    finite_refinement_passed=error<plan['finite_time_refinement_logT_tolerance'])
            np.savez_compressed(OUT/f'path-{steps}.npz',lnT_base=state['lnT'],lnT_shift=shift,eos=old_a,total_energy_increment=total_energy)
            save(f'path-{steps}.json',record);records.append(record);previous=shift.copy()
    save('result.json',dict(classification='Counterexample candidate',completed=True,paths=records,
        original_endpoint_failure_preserved=True,physical_EOS_certified=False,continuous_caloric_error_certified=False,
        fixed_density_metric_composition=True,full_GR_evolution=False))
