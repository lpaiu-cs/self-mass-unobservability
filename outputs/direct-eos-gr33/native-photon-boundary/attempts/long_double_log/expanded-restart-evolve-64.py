def evolve(n):
    assert read(OUT/'metric-result.json')['passed'];Model=initialize()
    prefix=np.load(OUT/f'sweep-1/photons/interval-15-{n}.npz');count=len(prefix['joint_stage_times'])
    recovered=np.load(OUT/f'prefix-recovered-{n}.npz');moments=list(recovered['photon_moments'][:count]);ports=list(recovered['radial_ports'][:count])
    oldrow=read(OUT/f'prefix-seed-{n}.json');cut=float(prefix['actual_step_edges'][-1])
    native=[v for v in oldrow['anchor_checks'] if v['time']<=cut+1e-18]
    branches=[v for v in oldrow['branch_checks'] if v['time']<=cut+1e-18]
    intervals=oldrow['intervals'][:15];run=Model.run;stage=run.__globals__['stages'];boundary=Model.boundary_ports
    def observed(m,t,h,x,g,lus):
        pair,mechanical=stage(m,t,h,x,g,lus)
        for v in pair:
            xx=v[0];moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*actual.AMP)
        np.savez_compressed(OUT/f'last-pair-{n}.npz',time=t,step=h,photons=[v[0] for v in pair],gas=[v[1] for v in pair],x_initial=x,g_initial=g)
        write(OUT/f'capture-{n}.json',dict(actual_steps=len(moments)//2,new_steps=(len(moments)-count)//2))
        return pair,mechanical
    def port(m,t,x):
        value=boundary(m,t,x);ports.append(value.copy()*actual.AMP);return value
    Model.boundary_ports=port;Model.run=bind(run,stages=observed);restart=f'interval-15-{n}'
    for j in range(16,17):
        m=Model(n);label=f'interval-{j}-{n}';row=m.run(n,label,j*n//16,restart);assert row['passed']
        native+=m.anchor_checks;branches+=m.branch_checks;intervals.append(dict(row))
        path=OUT/f'sweep-1/photons/{label}.npz';z=np.load(path)
        np.savez_compressed(OUT/f'recovered-{n}.npz',times=z['joint_stage_times'],weights=z['joint_stage_weights'],photon_moments=moments,radial_ports=ports,collision_rates=z['joint_collision_rates_scaled'],angular=z['accepted_angular_luminosity'],endpoint_occupation=z['restart_x']*m.scale*actual.AMP)
        checks=dict(newton=m.newton_iterations,stages=m.stage_log);restart=label;del m;gc.collect()
    dst=OUT/f'sweep-1/photons/return-{n}.npz';os.link(path,dst);os.link(path.with_suffix('.json'),dst.with_suffix('.json'))
    for k in ['actual_step_edges','joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled','joint_discard_rates_scaled','conserved_material_history','material_floor_discard_history_scaled','material_history','photon_history_scaled_occupation','radial_ports','collision_transfer','accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','moments','t']:
        assert np.array_equal(z[k][:len(prefix[k])],prefix[k]),('Restart prefix changed',k)
    assert np.array_equal(np.array(moments)[:count],recovered['photon_moments'][:count]) and np.array_equal(np.array(ports)[:count],recovered['radial_ports'][:count])
    anchor=np.load(saved(n));assert np.array_equal(z['actual_step_edges'],anchor['actual_step_edges']) and np.array_equal(z['joint_stage_times'],anchor['joint_stage_times'])
    audit,_,_,_=geometry.prior.run.verify(dst,dst);expected=np.sum(z['joint_stage_weights'][:,None,None]*ports,axis=0,dtype=LD)
    error=float(np.max(abs(expected-z['radial_ports'][-1])/np.maximum(abs(z['radial_ports'][-1]),LD('1e-290'))));assert error<1e-12,error
    row.update(audit=audit,intervals=intervals,anchor_checks=native,branch_checks=branches,
        actual_stage_photon_moments_captured=True,captured_radial_port_relative=error,
        maximum_true_stage=max(v[-1]['relative'] for v in checks['newton']),maximum_true_physical_stage=max(max(v[-1]['moments']) for v in checks['newton']),
        same_saved_stage_equation=True,actual_return_time_evolved=True,full_declared_period=True,
        self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',**checks))
