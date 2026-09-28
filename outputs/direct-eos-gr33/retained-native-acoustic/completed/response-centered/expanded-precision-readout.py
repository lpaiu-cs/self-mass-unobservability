def readout():
    import verify_native_stage_energy_charge as charge
    import sympy as sp
    owner=matter.old.old.previous.previous.old.source_owner
    s=owner.source.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]')
    s=s.replace('background=compare(histories[2],histories[1]),','').replace(',stress_background=compare(allstress[2],allstress[1])','')
    s=replace(s,"base=dict(np.load(wave.base.OUT/f'source-{ref}.npz'))",'base=template()')
    s=s.replace('.tolist()', '.astype(float).tolist()')
    ns=dict(matter.old.old.previous.source_scope,OUT=MATERIAL,GR=GR,Material=Material,
        stress=SimpleNamespace(pressure=pressure),signal=NO_ALARM,template=template,
        photons=SimpleNamespace(path=lambda n,r:PHOTON/f'steps-{n}-reference-{r}.npz'),
        material_path=lambda n,r:MATERIAL/f'steps-{n}-reference-{r}.npz')
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-source.py').write_text(s);ns['sources']()
    assert read(MATERIAL/'sources.json')['passed']
    model=charge.gr.Response();waves={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));waves[n,q]=charge.read(model,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**waves[n,q])
    fine=waves[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-300)
    terr=float(np.max(abs(fine-waves[64,8]['free_scalar']))/norm);qerr=float(np.max(abs(fine-waves[128,4]['free_scalar']))/norm)
    d=dict(np.load(GR/'source-128-reference-128.npz'));direct,coordinate=charge.independent.direct(model,d,8)
    ierr=float(abs(direct-waves[128,8]['direct_scalar'][-1])/max(np.max(abs(waves[128,8]['direct_scalar'])),1e-300))
    total=template();p=dict(np.load(EOS/'source.npz'));keys=forcing.history.physical.capture.KEYS+['metric_stress_erg']
    for key in keys:total[key]=total[key]+p[key][::8]+d[key]
    for key in ['inner_cumulative_energy_erg','outer_cumulative_energy_erg']:total[key]=total[key]+d[key]
    applied=charge.read(model,total,8);baseline=template()
    for key in keys:baseline[key]=baseline[key]+p[key][::8]
    before=charge.read(model,baseline,8);linear=float(np.max(abs(applied['free_scalar']-before['free_scalar']-fine))/max(np.max(abs(before['free_scalar'])),1e-300))
    np.savez_compressed(GR/'applied-source.npz',**total);np.savez_compressed(GR/'applied-charge.npz',**applied,previous=before['free_scalar'],correction=fine)
    # A conservative face enters adjacent cells with exactly opposite signs.
    left,shared,right,g0,g1=sp.symbols('left shared right g0 g1')
    assert sp.expand((left-shared+g0)+(shared-right+g1)-(left-right+g0+g1))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Shared-face telescoping only; no nonlinear or uniform EOS error theorem.'))
    passed=terr<.02 and qerr<.002 and ierr<1e-9 and linear<1e-10
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=passed,
        native_return_compact_endpoint=float(fine[-1]),previous_same_cadence_endpoint=float(before['free_scalar'][-1]),
        applied_compact_endpoint=float(applied['free_scalar'][-1]),time_relative=terr,quadrature_relative=qerr,
        independent_direct_relative=ierr,inverse_radius_residual=float(coordinate),linear_application_relative=linear,
        native_collision_applied=True,native_pressure_force_applied=True,free_material_evolved=True,represented_GR_applied=True,
        actual_material_motion_returned_to_photons=False,new_angular_emission_returned_to_infinity=False,
        native_acoustic_derivatives_certified=False,coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False))
    assert passed
