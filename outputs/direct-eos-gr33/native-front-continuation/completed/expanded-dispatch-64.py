def pilot():
    assert read(OUT/'check-result.json')['passed'];m=Model(64)
    row=m.run(64,'pilot-64',4,'restored-64');file=paths(1)[0]/'pilot-64.npz';p=dict(np.load(file))
    p.update(joint_stage_times=np.array(m.stage_t),joint_stage_weights=np.array(m.stage_h),
        joint_stage_conserved_scaled=np.array(m.stage_states),joint_native_rates_scaled=np.array(m.stage_native),
        joint_collision_rates_scaled=np.array(m.stage_collision),joint_discard_rates_scaled=np.array(m.stage_discard))
    np.savez_compressed(file,**p)
    accumulated=np.sum((p['joint_native_rates_scaled']+p['joint_collision_rates_scaled'])*p['joint_stage_weights'][:,None,None],axis=(0,1),dtype=LD)
    final=p['delta_material']/AMP*m.units+p['material_floor_discard_scaled']
    balance=(abs(np.sum(final,axis=0)-accumulated)/np.maximum(np.sum(abs(final),axis=0),LD('1e-290'))).astype(float).tolist()
    port=joint.previous.run.packets(file)[2]
    row.update(same_solution_material_balance=balance,angular_port_relative=port,
        maximum_true_stage=max(a[-1]['relative'] for a in m.newton_iterations),
        maximum_true_physical_stage=max(max(a[-1]['moments']) for a in m.newton_iterations),
        maximum_newton_iterations=max(map(len,m.newton_iterations)),full_incident_input_applied=True,
        self_GR_return_closed=False,time_comparison_completed=False,full_horizon_completed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    row['passed']=bool(row['passed'] and max(balance)<1e-8 and port<1e-12)
    write(file.with_suffix('.json'),row);write(OUT/'pilot-result.json',row)
    write(OUT/'stage-checks.json',dict(classification='Counterexample candidate',newton=m.newton_iterations,stages=m.stage_log))
    print(json.dumps(row),flush=True);assert row['passed'],row
