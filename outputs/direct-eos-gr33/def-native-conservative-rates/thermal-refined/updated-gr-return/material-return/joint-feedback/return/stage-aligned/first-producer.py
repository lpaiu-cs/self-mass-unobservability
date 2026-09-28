"""Apply the accepted source-knot repair to the large material/GR return."""
from pathlib import Path
from types import FunctionType
import inspect,json,sys
import numpy as np
import sympy as sp
import verify_native_updated_joint_feedback as old
import def_native_stage_aligned_material as fixed

OUT=old.OUT/'stage-aligned';GR=OUT/'gr';matter=old.matter;prior=old.prior
write=old.write;sha=old.sha


class Material(old.Material):
    # Reuse the corrected methods, but retain the large physical metric/transfer
    # initialization. Inheriting fixed.Material would select the tiny increment.
    closing_stage=False
    interval_fields=fixed.Material.interval_fields
    directional_rhs=fixed.Material.directional_rhs
    fields=fixed.Material.fields
    rhs=fixed.Material.rhs
    end_rhs=fixed.Material.end_rhs
    run=FunctionType(fixed.Material.run.__code__,dict(fixed.Material.run.__globals__,OUT=OUT),
                     argdefs=fixed.Material.run.__defaults__)


dispatch=FunctionType(old.dispatch.__code__,dict(old.dispatch.__globals__,OUT=OUT,__file__=__file__))
worker=FunctionType(old.worker.__code__,dict(old.worker.__globals__,OUT=OUT,Material=Material))
pilot_source=prior.replace(inspect.getsource(matter.pilot),"previous.OUT/'material-audit.json'","old.OUT/'production.json'")
scope=dict(vars(matter),OUT=OUT,Material=Material,dispatch=dispatch,__file__=__file__,old=old)
exec(compile(pilot_source,__file__,'exec'),scope);pilot=scope['pilot']
production=FunctionType(matter.production.__code__,scope)
source_scope=dict(old.sources_owner.__globals__,OUT=OUT,GR=GR,Material=Material,
                  material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
sources_owner=FunctionType(old.sources_owner.__code__,source_scope)


class GRResponse(old.GRResponse):
    run=FunctionType(old.GRResponse.run.__code__,dict(old.GRResponse.run.__globals__,OUT=GR))


fields_owner=FunctionType(old.fields_owner.__code__,dict(old.fields_owner.__globals__,OUT=OUT,GR=GR,
                            GRResponse=GRResponse,__file__=__file__))
check_jump=FunctionType(fixed.check_jump.__code__,dict(vars(fixed),OUT=OUT,base=old,Material=Material))


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir()
    assert json.loads((old.OUT/'audit.json').read_text())['passed']
    assert json.loads((fixed.OUT/'audit.json').read_text())['passed']
    paths=[Path(__file__),Path(old.__file__),Path(fixed.__file__),Path(matter.__file__),
           Path(matter.branch.base.__file__),Path(old.source_owner.__file__),
           old.OUT/'audit.json',old.OUT/'production.json',fixed.OUT/'normalized-forcing-jump.json']
    for n,r in matter.PATHS:paths.extend([old.photon_path(n,r),old.OUT/f'steps-{n}-reference-{r}.npz'])
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='93f2be745',
        claim='Apply the accepted SSP closing-stage left-interval and zero-weight probe normalization fixes to the large Phase134 material return; propagate the changed pressure/trace through compact GR and measure the remaining material/photon mismatch.',
        decision='If the old residual changes, replace the large material input to the next coupled step with this corrected trajectory. Do not iterate the negligible Phase135 GR increment to avoid this dominant error.',
        reuse='Reuse all saved physical corrected-EOS/background/atmosphere, Phase134 photons and primary metric histories. Same531 cells,64/128 paths,SSP2,CFL,3.434431ms and native owners. No new EOS calls, meshes or horizons.',
        scope='This repairs the main material response and measures its charge change. Finite waveform norms and charge differences do not certify a fixed point, all EOS derivatives, exterior/floor dynamics, nonlinear GR or final charge.',
        budgets=dict(check_seconds=15,pilot_seconds=45,production_wall_seconds=350,CPU_processes=3,threads_each=1,total_virtual_GiB=6,source_seconds=75,GR_seconds=120,production_attempts=1),
        forecast='Use completed Phase134 late-CFL raw-call counts and newly measured concurrent two-step prefixes. Require2x maximum forecast including15s setup/export within350s. Later changed trajectories remain extrapolated; hard caps stop the run.',
        gates=dict(owner=1e-8,conservation=1e-8,directional=.002,branch=.01,small_state=1e-6,time=.02,background=.02,pressure=.002,GR_quadrature=.002,GR_independent=1e-9,source_identity=1e-12),
        stop='Stop on any physics gate, forecast or wall-cap failure. No automatic refinement, extra paths, longer horizon, relaxed gate or additional production attempt.',
        bindings={str(p):sha(p) for p in paths}))
    a,b,h=sp.symbols('a b h');assert sp.simplify(h*(a+b)/2-h*a-h*(b-a)/2)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='For a step ending at a source knot, using the right instead of left derivative adds h*(b-a)/2. This identity does not bound the full EOS or coupled error.'))
    import signal
    signal.signal(signal.SIGALRM,matter.branch.base.flow.old.optical.timeout);signal.alarm(15)
    check_jump('forcing-jump.json');signal.alarm(0)


def sources():
    matter.configure();assert json.loads((OUT/'production.json').read_text())['passed']
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=75,
        claim='Apply corrected large conserved material histories to the unchanged complete pressure/trace/nonrest GR source map and original comparison gates.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.source_owner.__file__),OUT/'plan.json',OUT/'production.json',old.joint.OUT/'audit.json']}))
    sources_owner()


def fields():
    matter.configure();fields_owner()


def audit():
    a=json.loads((OUT/'production.json').read_text());b=json.loads((OUT/'sources.json').read_text());g=json.loads((GR/'result.json').read_text())
    assert a['passed'] and b['passed'] and g['passed'];checked=0;source_error=0.;changes=[]
    for plan in [OUT/'plan.json',OUT/'execution-plan.json',OUT/'source-plan.json',GR/'plan.json']:
        for p,h in json.loads(plan.read_text())['bindings'].items():assert sha(p)==h,p;checked+=1
    for n,r in matter.PATHS:
        d=np.load(GR/f'source-{n}-reference-{r}.npz');stress=np.load(OUT/f'stress-{n}-reference-{r}.npz')['material']
        rest=d['baryon_g'].astype(np.longdouble)*np.longdouble(d['cx'])*np.longdouble(prior.C)**2
        total=rest+d['gas_nonrest_energy_erg']
        checks=[total-stress[:,0],d['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
                d['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),
                d['metric_stress_erg']-(total+d['photon_energy_erg']-stress[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(x)) for x in checks)/max(np.max(abs(stress)),1e-300)))
        previous=np.load(old.OUT/f'steps-{n}-reference-{r}.npz');current=np.load(OUT/f'steps-{n}-reference-{r}.npz')
        assert np.array_equal(previous['t'],current['t'])
        x=previous['history_scaled'];y=current['history_scaled'];norm=np.maximum(np.max(np.sum(abs(y),axis=2),axis=0),1e-300)
        changes.append(dict(steps=n,reference=r,times=len(current['t']),relative=(np.max(np.sum(abs(x-y),axis=2),axis=0)/norm).tolist()))
    assert source_error<1e-12
    before=json.loads((old.OUT/'audit.json').read_text());charge=g['paths'][0]['endpoint_compact_with_metric']
    primary=float(np.load(prior.run.physical.GR/'wave-128-g8.npz')['free_scalar'][-1])
    residual=b['paths'][1]['energy_H_waveform_residual'];prior_residual=before['energy_H_waveform_residual']
    delta=charge-before['additional_compact_charge']
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,source_identity_relative=source_error,
        corrected_large_material_return_evolved=True,actual_corrected_sources_applied_to_compact_GR=True,
        additional_compact_charge=charge,charge_change_from_original_large_return=delta,
        relative_charge_change_from_original_large_return=delta/before['additional_compact_charge'],
        additional_over_primary_charge=charge/primary,charge_change_over_primary=delta/primary,
        material_history_change=changes,conserved_order=['baryon','radial_momentum_c','reference_energy','neutral_H'],
        canonical_common_times=17,previous_energy_H_waveform_residual=prior_residual,energy_H_waveform_residual=residual,
        residual_ratio_to_previous=(np.array(residual)/prior_residual).tolist(),
        material_time_max=max(b['comparisons']['time']),material_background_max=max(b['comparisons']['background']),
        stress_time_max=max(b['comparisons']['stress_time']),stress_background_max=max(b['comparisons']['stress_background']),
        corrected_motion_returned_to_photons=False,corrected_GR_returned_to_transport=False,
        uniform_contraction_bound=False,coupled_fixed_point_verified=False,full_EOS_history_error_enclosed=False,
        full_exterior_scalar=False,discarded_material_transport_closed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
