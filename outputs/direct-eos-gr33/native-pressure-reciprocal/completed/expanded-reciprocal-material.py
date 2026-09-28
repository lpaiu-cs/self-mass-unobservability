def material(pilot):
    assert read(OUT/'packet-check.json')['passed'];initialize();folder=paths(1)[1];rows=[];histories=[];forecasts=[];started=time.monotonic()
    if not pilot:assert read(folder/'pilot.json')['eligible']
    for n in [64,128]:
        label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        mark=time.monotonic();m=c.Material(128,n)
        row=m.run(n,label,n//16 if pilot else None,None if pilot else f'pilot-{n}')
        wall=time.monotonic()-mark
        with np.load(folder/f'{label}.npz') as d:
            z=d['delta_scaled'];end=d['t'][-1]
            rates=[m.rhs(end,z,p)[0] for p in [.5,1.,2.]]
            probe=[(np.sum(abs(r-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).astype(float).tolist() for r in [rates[0],rates[2]]]
            with np.load(paths(1)[0]/f'steps-{n}-reference-128.npz') as p:
                clock=p['t'][p['t']<=end+1e-18]
                ids=[int(np.argmin(abs(d['t']-t))) for t in clock];assert np.max(abs(d['t'][ids]-clock))<1e-18
                histories.append(d['history_scaled'][ids].copy())
            if not pilot:
                with np.load(folder/f'pilot-{n}.npz') as old:
                    for key in ['t','history_scaled','ledgers_scaled','norms_scaled','discards_scaled']:
                        assert np.array_equal(d[key][:len(old[key])],old[key]),('Accepted prefix changed',key)
        row.update(worker_wall_seconds=wall,physical_branch_ratio=m.physical_branch_ratio,
            maximum_owner_error=max(v['owner_error'] for v in m.cache.values()),endpoint_probe_half_nominal_double=probe,constitutive_log_step=1e-5,gross_flux_difference_used=False,forward_probe_indicator=None)
        row['passed']=bool(row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8 and np.max(probe)<.002)
        write(folder/f'{label}.json',row);rows.append(row);assert row['passed'],row
        if pilot:
            old=read(Path('native-mixed-return162-work/sweep-1/material')/f'steps-{n}-reference-128.json')
            per_step=max(row['seconds']/row['substeps'],old['seconds']/old['substeps'])
            setup=max(old['worker_wall_seconds']-old['seconds'],wall-row['seconds'])
            forecasts.append(max(0,old['substeps']-row['substeps'])*per_step+setup+5)
        del m;gc.collect()
    errors=c.relative(*histories);result=dict(classification='Counterexample candidate',passed=max(errors)<.02,
        rows=rows,time_comparison=errors,seconds=time.monotonic()-started,
        photon_paths_recomputed=False,final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    if pilot:result.update(upper_remaining_seconds=2*sum(forecasts),eligible=result['passed'] and 2*sum(forecasts)<CAPS['production'])
    write(folder/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed'],result
