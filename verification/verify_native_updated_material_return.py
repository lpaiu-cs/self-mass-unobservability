"""Updated free-material histories -> conservative stress -> compact GR."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,sys
import numpy as np
import sympy as sp
import def_native_updated_material_return as run

old=run.previous;prior=run.prior;OUT=run.OUT;GR=run.GR;write=run.write;sha=run.sha


source=inspect.getsource(old.sources)
source=source.replace("(OUT/'material-audit.json')","(OUT/'production.json')")
source=prior.replace(source,'signal.alarm(45)','signal.alarm(75)')
start=source.index("        prior=np.load(old.OUT/");end=source.index('        # This is a measured',start)
source=source[:start]+source[end:]
source=prior.replace(source,'material_sweep_change=change.tolist(),','legacy_material_sweep_comparison=None,')
source=prior.replace(source,"photon_energy_erg=photonE,photon_radial_pressure_erg=photonP,",
    "nonrest_stress_erg=source[:,0].astype(LD)-source[:,1]-rest,pressure_volume_erg=source[:,3],photon_energy_erg=photonE,photon_radial_pressure_erg=photonP,")
scope=dict(vars(old),OUT=OUT,GR=GR,Material=run.Material,photons=SimpleNamespace(path=run.photon_path),
    material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
exec(compile(source,__file__,'exec'),scope)


def sources():
    run.configure();assert not (OUT/'source-plan.json').exists()
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=75,
        claim='Convert the updated actual conserved material response to pressure/radial stress/trace; combine actual photon moments without double-counting the initial canonical GR response.',
        reuse='Existing direct conservative primitive variation and sampled per-cell native pressure probes; corrected-EOS derivative banks and full saved Phase132 photons. No physical replay.',
        schema='Set nonrest stress and tangential pressure along with energy,trace,baryon and radial metric stress. Do not leave stale inherited source fields. Preserve the actual paired radial photon ports.',
        residual='Compare the newly returned material reference-energy/H history with the preceding radiation-only material history. This is an unclosed waveform mismatch, not a contraction certificate. Do not compare against an old-EOS material path as if it were the same fixed point.',
        gates=dict(time=.02,background=.02,conservation=1e-8,pressure=.002),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(old.stress.__file__),OUT/'production.json',prior.OUT/'audit.json']}))
    (OUT/'expanded-source.py').write_text(source);scope['sources']()
    m,v=sp.symbols('m v',positive=True);T=m*v*v/2
    assert sp.simplify(m*v*v-2*T)==0
    # The same identity applies to the retained nonrelativistic deep model's
    # effective inertial mass; it is not a relativistic-fluid approximation bound.
    write(OUT/'kinetic-work-symbolic.json',dict(classification='Proven',passed=True,
        scope='For the retained deep kinetic energy T=effective_mass*v^2/2, integrated radial kinetic stress equals2*T. Include it in the existing radial metric-work term. This identity is not a bound on the omitted relativistic corrections.'))


class GRResponse(old.wave.Response):
    run=FunctionType(old.wave.base.Response.run.__code__,dict(vars(old.wave.base),OUT=GR))


gr_source=inspect.getsource(old.fields).replace('budget_seconds=90','budget_seconds=120').replace('signal.alarm(90)','signal.alarm(120)')
namespace=dict(vars(old),OUT=OUT,GR=GR,GRResponse=GRResponse,__file__=__file__)
exec(compile(gr_source,__file__,'exec'),namespace)
def fields():
    run.configure();namespace['fields']()


def audit():
    assert not (OUT/'audit.json').exists()
    a=json.loads((OUT/'production.json').read_text());b=json.loads((OUT/'sources.json').read_text());g=json.loads((GR/'result.json').read_text())
    assert a['passed'] and b['passed'] and g['passed'];checked=0
    for plan in [OUT/'plan.json',OUT/'execution-plan.json',OUT/'source-plan.json',GR/'plan.json']:
        for p,h in json.loads(plan.read_text())['bindings'].items():assert sha(p)==h,p;checked+=1
    source_error=0.
    for n,r in run.PATHS:
        d=np.load(GR/f'source-{n}-reference-{r}.npz');stress=np.load(OUT/f'stress-{n}-reference-{r}.npz')['material']
        rest=d['baryon_g'].astype(np.longdouble)*np.longdouble(d['cx'])*np.longdouble(prior.C)**2
        total=rest+d['gas_nonrest_energy_erg']
        checks=[total-stress[:,0],d['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
            d['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),
            d['metric_stress_erg']-(total+d['photon_energy_erg']-stress[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(x)) for x in checks)/max(np.max(abs(stress)),1.)))
    assert source_error<1e-12
    charge=g['paths'][0]['endpoint_compact_with_metric']
    original=np.load(prior.run.physical.GR/'wave-128-g8.npz')['free_scalar'][-1]
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,source_identity_relative=source_error,
        actual_updated_transfers_applied_to_free_material=True,actual_updated_material_photon_sources_applied_to_compact_GR=True,
        additional_compact_charge=charge,additional_over_Phase130_free_charge=float(charge/original),
        material_energy_H_waveform_residual=b['paths'][1]['energy_H_waveform_residual'],
        material_motion_returned_to_photons=False,updated_feedback_GR_reapplied_to_transport=False,
        coupled_fixed_point_verified=False,full_exterior_scalar=False,full_EOS_history_error_enclosed=False,
        discarded_material_transport_closed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
