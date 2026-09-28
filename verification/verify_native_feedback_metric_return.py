"""Return the compensated GR-driven photons through free matter to charge."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,sys
import numpy as np
import sympy as sp
import def_native_feedback_metric_return as run
import verify_native_updated_joint_feedback as previous

matter=previous.matter;prior=run.prior;OUT=run.OUT/'return';GR=OUT/'gr'
write=run.write;sha=run.sha


def photon_path(n,r):return run.RESPONSE/f'steps-{n}-reference-{r}.npz'


class Material(previous.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps)
        g=dict(np.load(run.LAPSE/f'metric-{steps}-reference-{reference}-g8.npz'))
        assert np.array_equal(g['t'],self.t) and np.array_equal(g['radius_E'],self.rE);self.metric=g
        p=np.load(photon_path(steps,reference));ids=[int(np.argmin(abs(p['t']-t))) for t in self.t]
        assert np.max(abs(p['t'][ids]-self.t))<1e-18
        c=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/matter.AMP
    run=FunctionType(matter.branch.base.Material.run.__code__,dict(vars(matter.branch.base),OUT=OUT),argdefs=matter.branch.base.Material.run.__defaults__)


dispatch=FunctionType(matter.parallel.dispatch.__code__,dict(vars(matter.parallel),OUT=OUT,__file__=__file__))
worker=FunctionType(matter.worker.__code__,dict(vars(matter),OUT=OUT,Material=Material))
pilot_source=prior.replace(inspect.getsource(matter.pilot),"previous.OUT/'material-audit.json'","previous.OUT/'production.json'")
scope=dict(vars(matter),OUT=OUT,Material=Material,dispatch=dispatch,__file__=__file__,previous=previous)
exec(compile(pilot_source,__file__,'exec'),scope);pilot=scope['pilot']
production=FunctionType(matter.production.__code__,scope)

# Tiny physical increments must be assessed relative to their own norm.
source=prior.replace(previous.source_owner.source,
    'max(np.max(np.sum(abs(source[:,[3,1]]),axis=2)),1.)',
    'max(np.max(np.sum(abs(source[:,[3,1]]),axis=2)),1e-300)')
source=prior.replace(source,
    'np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1.)',
    'np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1e-300)')
source_scope=dict(previous.source_owner.scope,OUT=OUT,GR=GR,Material=Material,
    photons=SimpleNamespace(path=photon_path),material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
exec(compile(source,__file__,'exec'),source_scope)


class GRResponse(matter.previous.wave.Response):
    run=FunctionType(matter.previous.wave.base.Response.run.__code__,dict(vars(matter.previous.wave.base),OUT=GR))


fields_owner=FunctionType(previous.fields_owner.__code__,dict(previous.field_scope,OUT=OUT,GR=GR,GRResponse=GRResponse,__file__=__file__))


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir()
    assert json.loads((run.RESPONSE/'result.json').read_text())['passed']
    p,q=sp.symbols('p q',nonnegative=True)
    for old,delta,answer in [(p,q-p,q-p),(p,-q-p,-p),(-p,q+p,q),(-p,p-q,0)]:
        assert sp.simplify(sp.Max(old+delta,0)-sp.Max(old,0)-answer)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Nonnegative p,q parameterize all four old/new sign cases, including zero. Symbolic Max identities establish the compensated signed-frequency positive-part formula; they do not bound EOS or continuum errors.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='60877e975',
        claim='Apply only the compensated new GR geometry and its actual paired photon transfers to free mass/momentum/energy/H, then export its pressure/trace and solve its compact GR contribution.',
        budgets=dict(pilot_seconds=45,production_seconds=350,CPU_processes=3,threads_each=1,total_virtual_GiB=6,source_seconds=75,GR_seconds=120),
        reuse='Same updated physical background and adaptive signed material derivative owner. Its existing perturbation-size normalization handles the smaller direction without adding it to the large completed response. Reuse Phase134 full late-CFL raw-call counts for the forecast.',
        gates=dict(time=.02,background=.02,conservation=1e-8,directional=.002,pressure=.002,GR_quadrature=.002,independent=1e-9),
        limits='A finite additional GR feedback sweep, not an error enclosure. Do not use the measured gain as a uniform operator norm or disregard the remaining Phase134 waveform mismatch.',
        bindings={str(f):sha(f) for f in [Path(__file__),Path(run.__file__),Path(previous.__file__),run.OUT/'plan.json',run.LAPSE/'result.json',run.RESPONSE/'result.json',previous.OUT/'production.json']}))


def sources():
    matter.configure();assert json.loads((OUT/'production.json').read_text())['passed']
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=75,
        claim='Export actual incremental conservative sources with the same canonical geometric subtraction, actual photon moments/ports and all radial/tangential stress fields. Normalize pressure and waveform diagnostics by the tiny response itself.',
        bindings={str(f):sha(f) for f in [Path(__file__),Path(previous.source_owner.__file__),OUT/'production.json',run.RESPONSE/'result.json']}))
    (OUT/'expanded-source.py').write_text(source);source_scope['sources']()


def fields():
    matter.configure();fields_owner()


def audit():
    files=[run.LAPSE/'result.json',run.RESPONSE/'result.json',OUT/'production.json',OUT/'sources.json',GR/'result.json']
    results=[json.loads(p.read_text()) for p in files];assert all(p['passed'] for p in results)
    checked=0
    for path in [run.OUT/'plan.json',run.RESPONSE/'execution-plan.json',OUT/'plan.json',OUT/'execution-plan.json',OUT/'source-plan.json',GR/'plan.json']:
        for p,h in json.loads(path.read_text())['bindings'].items():assert sha(p)==h,p;checked+=1
    angular_errors=[];source_error=0.;gamma=1-1/np.sqrt(2)
    for n,r in matter.PATHS:
        d=np.load(photon_path(n,r));h=d['t'][-1]/n;times=d['accepted_angular_times'];L=d['accepted_angular_luminosity']
        assert np.max(abs(times-np.ravel(h*(np.arange(n)[:,None]+np.array([gamma,1.])))))<1e-18
        flux=L@((np.arange(4)*2+1)/32);integral=h*np.sum(flux.reshape(n,2)*[1-gamma,gamma])
        e=float(abs(integral-d['radial_ports'][-1,1,1])/max(h*np.sum(abs(flux)),1e-300));assert e<1e-12;angular_errors.append(e)
        s=np.load(GR/f'source-{n}-reference-{r}.npz');stress=np.load(OUT/f'stress-{n}-reference-{r}.npz')['material']
        rest=s['baryon_g'].astype(np.longdouble)*np.longdouble(s['cx'])*np.longdouble(prior.C)**2
        total=rest+s['gas_nonrest_energy_erg']
        errors=[total-stress[:,0],s['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),s['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),s['pressure_volume_erg']-stress[:,3],s['metric_stress_erg']-(total+s['photon_energy_erg']-stress[:,1]-s['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(stress)),1e-300)))
    assert source_error<1e-12
    before=json.loads((previous.OUT/'audit.json').read_text());charge=results[-1]['paths'][0]['endpoint_compact_with_metric']
    primary=np.load(prior.run.physical.GR/'wave-128-g8.npz')['free_scalar'][-1]
    old=np.load(previous.GR/'fields-128-reference-128-g8.npz');new=np.load(GR/'fields-128-reference-128-g8.npz')
    gains={k:float(np.max(abs(new[k]))/max(np.max(abs(old[k])),1e-300)) for k in ['U','delta_lambda']}
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,source_identity_relative=source_error,
        accepted_angular_port_relative=angular_errors,actual_returned_GR_applied_to_photons_and_free_material=True,
        actual_incremental_material_photon_source_applied_to_GR=True,new_compact_charge_increment=charge,
        increment_over_previous_compact_correction=charge/before['additional_compact_charge'],increment_over_primary_charge=float(charge/primary),
        measured_field_increment_ratios=gains,retained_Phase134_energy_H_waveform_residual=before['energy_H_waveform_residual'],
        new_incremental_energy_H_waveform_residual=results[-2]['paths'][1]['energy_H_waveform_residual'],
        incremental_material_motion_returned_to_photons=False,uniform_contraction_bound=False,coupled_fixed_point_verified=False,
        full_EOS_history_error_enclosed=False,discarded_material_transport_closed=False,full_exterior_scalar=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
