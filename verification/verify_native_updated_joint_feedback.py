"""Return the updated simultaneous sweep to matter/GR and measure mismatch."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,sys
import numpy as np
import def_native_updated_joint_feedback as joint
import def_native_updated_material_return as matter
import verify_native_updated_material_return as source_owner

OUT=joint.OUT/'return';GR=OUT/'gr';write=joint.write;sha=joint.sha;prior=matter.prior


def photon_path(n,r):return joint.OUT/f'steps-{n}-reference-{r}.npz'


class Material(matter.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);p=np.load(photon_path(steps,reference))
        ids=[int(np.argmin(abs(p['t']-t))) for t in self.t];assert np.max(abs(p['t'][ids]-self.t))<1e-18
        collision=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(collision[:,:,0]),p['moments'][ids,3]/self.a,collision[:,:,0],collision[:,:,1]],axis=1)/matter.AMP
    run=FunctionType(matter.branch.base.Material.run.__code__,dict(vars(matter.branch.base),OUT=OUT),argdefs=matter.branch.base.Material.run.__defaults__)


dispatch=FunctionType(matter.parallel.dispatch.__code__,dict(vars(matter.parallel),OUT=OUT,__file__=__file__))
worker=FunctionType(matter.worker.__code__,dict(vars(matter),OUT=OUT,Material=Material))
pilot_source=inspect.getsource(matter.pilot).replace("previous.OUT/'material-audit.json'","matter.OUT/'production.json'")
scope=dict(vars(matter),OUT=OUT,Material=Material,dispatch=dispatch,__file__=__file__,matter=matter)
exec(compile(pilot_source,__file__,'exec'),scope);pilot=scope['pilot']
production=FunctionType(matter.production.__code__,scope)
source_scope=dict(source_owner.scope,OUT=OUT,GR=GR,Material=Material,photons=SimpleNamespace(path=photon_path),material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
sources_owner=FunctionType(source_owner.scope['sources'].__code__,source_scope)


class GRResponse(matter.previous.wave.Response):
    run=FunctionType(matter.previous.wave.base.Response.run.__code__,dict(vars(matter.previous.wave.base),OUT=GR))


field_scope=dict(source_owner.namespace,OUT=OUT,GR=GR,GRResponse=GRResponse,__file__=__file__)
fields_owner=FunctionType(source_owner.namespace['fields'].__code__,field_scope)


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir();assert json.loads((joint.OUT/'audit.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the actual updated simultaneous photon/material collision histories back to free material and compact GR; measure the remaining common energy/H waveform mismatch and the change in all four conserved material histories.',
        reuse='Same updated531-cell background,EOS banks,GR/lapse,shared fluxes and three64/128 paths. Retain the2*T deep kinetic stress and established donor/precision checks. No physical background or completed photon replay.',
        budgets=dict(pilot_seconds=45,production_wall_seconds=350,CPU_processes=3,threads_each=1,total_virtual_memory_GiB=6,source_seconds=75,GR_seconds=120),
        forecast='Use completed Phase133 raw-call counts including late CFL work with newly measured concurrent prefix seconds/raw-call;2x maximum plus existing overhead must fit350s.',
        gates=dict(owner=1e-8,conservation=1e-8,directional=.002,time=.02,background=.02,pressure=.002,GR_quadrature=.002,GR_independent=1e-9),
        interpretation='A measured small waveform mismatch is not a universal contraction constant. Preserve incomplete fixed-point status and all exterior/floor/EOS conditions; do not repeatedly iterate solely to acquire another passing record.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(joint.__file__),Path(matter.__file__),Path(source_owner.__file__),joint.OUT/'audit.json',matter.OUT/'audit.json']}))


def sources():
    matter.configure();assert json.loads((OUT/'production.json').read_text())['passed']
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=75,
        claim='Use direct conserved native primitive variations and the same canonical-source subtraction to export the returned material and actual new photon stress; retain every GR source field consistently.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(source_owner.__file__),OUT/'production.json',joint.OUT/'audit.json']}))
    sources_owner()


def fields():
    matter.configure();fields_owner()


def audit():
    a=json.loads((OUT/'production.json').read_text());b=json.loads((OUT/'sources.json').read_text());g=json.loads((GR/'result.json').read_text())
    assert a['passed'] and b['passed'] and g['passed'];checked=0
    plans=[joint.OUT/'plan.json',joint.OUT/'angular-repair-plan.json',joint.OUT/'execution-plan.json',OUT/'plan.json',OUT/'execution-plan.json',OUT/'source-plan.json',GR/'plan.json']
    for plan in plans:
        for p,h in json.loads(plan.read_text())['bindings'].items():
            target=joint.OUT/'first-producer.py' if plan==joint.OUT/'plan.json' and Path(p)==Path(joint.__file__) else Path(p)
            assert sha(target)==h,p;checked+=1
    changes=[]
    for n,r in matter.PATHS:
        old=np.load(matter.OUT/f'steps-{n}-reference-{r}.npz');new=np.load(OUT/f'steps-{n}-reference-{r}.npz')
        assert np.array_equal(old['t'],new['t'])
        x=old['history_scaled'];y=new['history_scaled'];norm=np.maximum(np.max(np.sum(abs(y),axis=2),axis=0),1.)
        changes.append(dict(steps=n,reference=r,relative=(np.max(np.sum(abs(x-y),axis=2),axis=0)/norm).tolist()))
    before=json.loads((matter.OUT/'audit.json').read_text());charge=g['paths'][0]['endpoint_compact_with_metric']
    primary=np.load(prior.run.physical.GR/'wave-128-g8.npz')['free_scalar'][-1]
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,
        actual_updated_material_motion_returned_to_photons=True,actual_new_photon_transfer_returned_to_material=True,
        actual_returned_sources_applied_to_compact_GR=True,additional_compact_charge=charge,
        additional_over_primary_charge=float(charge/primary),compact_charge_change_from_first_material_return=charge-before['additional_compact_charge'],
        material_history_change=changes,conserved_order=['baryon','radial_momentum_c','reference_energy','neutral_H'],
        prior_energy_H_waveform_residual=before['material_energy_H_waveform_residual'],energy_H_waveform_residual=b['paths'][1]['energy_H_waveform_residual'],
        updated_feedback_GR_reapplied_to_transport=False,coupled_fixed_point_verified=False,full_EOS_history_error_enclosed=False,
        full_exterior_scalar=False,discarded_material_transport_closed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
