"""Read-only tracing of the failed nonlinear thermal endpoint arithmetic."""
import json,shutil,sys
import numpy as np
import gr_nonlinear_thermal as thermal

g=thermal.g;OLD=thermal.OUT;OUT=g.OUT/'gr-thermal-endpoint-replay'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failed=json.loads((OLD/'path-1-progress.json').read_text())['rows'][0]
    assert failed['global_energy_relative_to_exchange']>1e-8
    assert not (OLD/'result.json').exists()
    oldcode=g.ROOT/'verification/gr_nonlinear_thermal.py'
    shutil.copy2(oldcode,OLD/'failed-gr_nonlinear_thermal.py')
    log=g.ROOT/'outputs/gr-nonlinear-thermal33-run.log';shutil.copy2(log,OLD/'failed-run.log')
    manifest=dict(classification='Counterexample candidate',passed=False,failed_gate=failed,
        endpoint_arrays_not_captured_in_original_run=True,
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OLD.iterdir() if p.is_file()})
    (OLD/'failure-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    plan=json.loads((OLD/'plan.json').read_text())
    plan['endpoint_replay']=dict(classification='Counterexample candidate',checkpoint='9cf504d',
        original_failed_manifest_sha256=g.c.sha(OLD/'failure-manifest.json'),
        diagnostic_code_sha256=g.c.sha(g.ROOT/'verification/gr_thermal_endpoint_replay.py'),
        intervention='No physical, algorithmic or tolerance change. Observe the original Python run frame after constructing its first step record and before its failed assertion. Save the previously uncaptured arrays, and retain the expected failed gate.',
        expected_failed_step=failed)
    save('plan.json',plan)
    initial=np.array([float(2**60),4.]);increment=np.array([1.,-1.])
    endpoint=initial+increment;reconstructed=endpoint-initial
    assert increment.sum()==0 and reconstructed.sum()==-1
    save('offset-control.json',dict(classification='Proven',passed=True,
        binary64_initial=initial.tolist(),exact_increment=increment.tolist(),
        endpoint_reconstructed_increment=reconstructed.tolist(),
        statement='Two exact opposite increments conserve their sum, but storing large-offset binary64 endpoints and subtracting their initial values loses the first increment. A small local relative residual cannot certify a small global energy defect relative to the transported increment. This manufactured arithmetic control is separate from attribution of the actual failed run.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text());replay=plan['endpoint_replay']
    assert g.c.sha(OLD/'failure-manifest.json')==replay['original_failed_manifest_sha256']
    assert g.c.sha(g.ROOT/'verification/gr_thermal_endpoint_replay.py')==replay['diagnostic_code_sha256']
    for rel,digest in json.loads((OLD/'failure-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    captured=[];filename=str(g.ROOT/'verification/gr_nonlinear_thermal.py')
    def trace(frame,event,arg):
        if event=='line' and frame.f_code.co_filename==filename and frame.f_code.co_name=='run':
            variables=frame.f_locals;row=variables.get('row',{})
            if isinstance(row,dict) and 'global_energy_relative_to_exchange' in row and not captured:
                assert variables['steps']==1 and variables['step']==0
                keys=['t','old','a','mass','dm','N','dt','flux','divergence','residual','scale','change']
                np.savez_compressed(OUT/'failed-endpoint.npz',**{k:np.asarray(variables[k]).copy() for k in keys},
                    initial_lnT=variables['state']['lnT'])
                save('captured-step.json',row);captured.append(True);sys.settrace(None)
                return None
        return trace
    thermal.OUT=OUT;sys.settrace(trace)
    try:thermal.run()
    except AssertionError as error:
        assert captured and isinstance(error.args[0],dict),repr(error)
        assert error.args[0]==replay['expected_failed_step'],error.args[0]
    else:raise AssertionError('Expected original global-energy gate to fail unchanged')
    finally:sys.settrace(None)
    analyze()


def analyze():
    a=dict(np.load(OUT/'failed-endpoint.npz'));ld=np.longdouble
    observed=a['mass'].astype(ld)*(a['a'][:,2].astype(ld)-a['old'][:,2].astype(ld))
    expected=ld(a['dt'])*a['divergence'].astype(ld);residual=observed-expected
    temp_step=expected/(a['mass'].astype(ld)*a['a'][:,10].astype(ld))
    unresolved=abs(temp_step)<np.spacing(abs(a['initial_lnT']))/2
    unchanged=a['t']==a['initial_lnT'];ulp=a['mass'].astype(ld)*np.spacing(abs(a['old'][:,2])).astype(ld)
    exchange=abs(expected).sum();value=abs(residual.sum())/exchange
    row=json.loads((OUT/'captured-step.json').read_text())
    assert abs(float(value)-row['global_energy_relative_to_exchange'])<1e-12
    result=dict(classification='Counterexample candidate',failure_reproduced=True,
        original_global_energy_gate_passed=False,global_energy_relative_to_exchange=float(value),
        cells_with_unchanged_stored_lnT=int(unchanged.sum()),
        cells_with_required_increment_below_half_lnT_ulp=int(unresolved.sum()),
        fraction_absolute_energy_residual_in_temperature_unresolved_cells=float(abs(residual[unresolved]).sum()/abs(residual).sum()),
        signed_energy_residual_temperature_unresolved_relative_to_exchange=float(residual[unresolved].sum()/exchange),
        aggregate_one_ulp_energy_endpoint_budget_relative_to_exchange=float(ulp.sum()/exchange),
        maximum_energy_residual_in_units_of_endpoint_ulp=float((abs(residual)/ulp).max()),
        physical_EOS_certified=False,full_GR_evolution=False)
    save('diagnosis.json',result);print('THERMAL ENDPOINT DIAGNOSIS',result,flush=True)
    paths=[p for p in OUT.iterdir() if p.is_file()]
    save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}))


if __name__=='__main__':globals()[sys.argv[1]]()
