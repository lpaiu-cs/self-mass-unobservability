def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(235)
    spec=json.loads((OUT/'plan.json').read_text())
    for p,h in spec['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    assert json.loads((OUT/'thermo-control.json').read_text())['passed']
    pilot=Flow(224).run('pilot-224');forecast=1.3*pilot['seconds']*(2+16+64)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',measured_pilot_seconds=pilot['seconds'],forecast_seconds=forecast,
        assumption='Explicit CFL cost proportional to cells squared for the same fixed horizon;30percent margin. Finer-grid wave speeds and costs are not yet measured.'))
    assert forecast<235,'Production forecast over budget'
    flat=Flow(224,flat=True).run('flat-224');reference=json.loads((prior.OUT/'result.json').read_text())['gas_outside_original_cut_g']
    flat_error=abs(flat['gas_outside_original_radius_g']/reference-1)
    write(OUT/'flat-control.json',dict(classification='Counterexample candidate',native_exact_fan_mass_relative=flat_error,passed=flat_error<.03))
    assert flat_error<.03,'Independent native flat fan control failed'
    rows=[Flow(n).run(f'cells-{n}') for n in [896,1792]]
    errors={k:abs(rows[0][k]/rows[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    passed=all(r['baryon_ledger_relative']<1e-10 and r['isentropic_Killing_energy_response_relative']<.02 for r in rows) and max(errors.values())<.02
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=bool(passed),refinement=errors,seconds=time.monotonic()-start,
        gas_replaced_not_added=True,actual_radial_nonlinear_gas_evolved=True,current_elastic_scattering_work_paired=True,
        complete_photon_transport=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)
