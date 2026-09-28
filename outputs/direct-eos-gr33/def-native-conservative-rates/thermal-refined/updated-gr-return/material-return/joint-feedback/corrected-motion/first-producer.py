"""Return the corrected large material history to photons/energy/H."""
from pathlib import Path
from types import FunctionType
import json,sys
import numpy as np
import sympy as sp
import def_native_corrected_joint_return as matter
import def_native_updated_joint_feedback as old
import verify_native_corrected_joint_return as execution

OUT=old.OUT/'corrected-motion';write=old.write;sha=old.sha;PATHS=old.PATHS;prior=old.prior


class Response(old.Response):
    __init__=FunctionType(old.Response.__init__.__code__,dict(vars(old),matter=matter))
    run=FunctionType(old.Response.run.__code__,dict(old.Response.run.__globals__,OUT=OUT),
                     argdefs=old.Response.run.__defaults__)


worker=FunctionType(old.worker.__code__,dict(vars(old),OUT=OUT,Response=Response))
dispatch=FunctionType(execution.dispatch.__code__,dict(vars(execution),OUT=OUT,worker=worker))
production_owner=FunctionType(old.parallel.production.__code__,dict(vars(old.parallel),OUT=OUT,dispatch=dispatch))
angular_owner=FunctionType(old.audit.__code__,dict(vars(old),OUT=OUT))


def prepare():
    assert not OUT.exists();OUT.mkdir();accepted=json.loads((matter.OUT/'audit.json').read_text());assert accepted['passed']
    files=[Path(__file__),Path(old.__file__),Path(old.base.__file__),Path(old.stage.__file__),Path(matter.__file__),
           Path(matter.fixed.__file__),Path(execution.__file__),matter.OUT/'audit.json',old.OUT/'result.json',old.OUT/'corrected-pilot.json']
    for n,r in PATHS:files.extend([matter.OUT/f'steps-{n}-reference-{r}.npz',matter.old.photon_path(n,r)])
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='245e3d611',
        claim='Apply the corrected large material baryon/momentum/inventory and noncollisional energy/H transport to actual simultaneous photon/energy/H equations, then return their new paired transfers to material and GR.',
        decision='Measure whether correcting the material input reduces the unresolved0.089750percent energy and0.143735percent H mismatch to the temporal-comparison scale. A small charge correction alone is not acceptance of coupling.',
        reuse='Same corrected physical EOS/background/atmosphere,531 cells,8 angles,152 frequencies,17 knots,64/128 clocks,3.434431ms,conservative primitive map and stage-time SDIRK. Reuse all primary GR/lapse and native banks. Change the prescribed material history only.',
        numerical_scope='One corrected waveform return; not uniform contraction, continuous-time residual control, full exterior/floor transport, nonlinear GR or final charge. Preserve the Phase135 tiny GR response separately; do not silently claim it is included in this large-motion sweep.',
        budgets=dict(pilot_wall_seconds=120,production_wall_seconds=650,CPU_processes=3,threads_each=1,total_virtual_GiB=9,new_native_states=0,production_attempts=1),
        forecast='Measure concurrent4/8/4 prefixes using the already loaded owners with Linux fork. For each path use both17-point/remaining-step extrapolation and the completed Phase134 late operator/stepping costs scaled by its matching old/new prefix rates. Require2x the larger forecast within650s. Changed late Krylov cost and shared I/O remain unmeasured; hard caps apply.',
        gates=dict(conservation=1e-8,linear=1e-12,velocity_jet=1e-4,mapping=1e-10,time=.02,background=.02,angular_port=1e-12),
        stop='Stop on a failed prefix, forecast, path, gate or cap; preserve prefixes. No automatic additional waveform iteration, refinement, longer horizon, new support, gate relaxation or physical replay.',
        bindings={str(p):sha(p) for p in files}))
    E,B,k,q=sp.symbols('E B k q')
    assert sp.expand((E-k*B)+k*B-E)==0 and sp.expand(q-q)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='The fixed-reference nonrest energy coordinate is exactly invertible; opposite paired photon/material collision transfers cancel. These identities do not certify EOS derivatives or continuum/coupling error.'))


def angular(label,n,count):
    d=np.load(OUT/f'{label}.npz');t=d['accepted_angular_times'];L=d['accepted_angular_luminosity'];gamma=1-1/np.sqrt(2)
    h=prior.run.flow.old.END/n;expected=(h*(np.arange(count)[:,None]+[gamma,1.])).ravel()
    assert len(L)==len(t)==2*count and np.max(abs(t-expected))<1e-18
    pairs=(L@((np.arange(4)*2+1)/32)).reshape(count,2)
    error=float(abs(h*np.sum(pairs*[1-gamma,gamma])-d['radial_ports'][-1,1,1])/max(h*np.sum(abs(pairs)),1.))
    assert error<1e-12;return error


def pilot():
    assert not (OUT/'pilot.json').exists()
    specs=[(64,128,'pilot-64-128',4,None),(128,128,'pilot-128-128',8,None),(128,64,'pilot-128-64',4,None)]
    rows,seconds=dispatch(specs,120);previous=json.loads((old.OUT/'corrected-pilot.json').read_text())['rows']
    full=json.loads((old.OUT/'result.json').read_text())['paths'];estimates=[];details=[]
    for row in rows:
        pair=(row['steps'],row['reference']);p=next(x for x in previous if (x['steps'],x['reference'])==pair)
        f=next(x for x in full if (x['steps'],x['reference'])==pair)
        point=row['operator_point_seconds']/row['operator_points'];step=row['stepping_seconds']/row['new_steps']
        extrapolated=17*point+(row['steps']-row['completed_steps'])*step+20
        measured_late=f['operator_point_seconds']*point/(p['operator_point_seconds']/p['operator_points'])+f['stepping_seconds']*step/(p['stepping_seconds']/p['new_steps'])+20
        estimates.append(max(extrapolated,measured_late));details.append(dict(prefix_extrapolation=extrapolated,scaled_completed_late_cost=measured_late))
    a=np.load(OUT/'pilot-64-128.npz')['moments'][-1,[0,1,2,3,5,6]];b=np.load(OUT/'pilot-128-128.npz')['moments'][-1,[0,1,2,3,5,6]]
    errors=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
    ports=[angular(label,n,count) for n,r,label,count,_ in specs]
    p=dict(classification='Counterexample candidate',rows=rows,equal_horizon=errors,angular_port_relative=ports,
           forecast_each_seconds=estimates,forecast_details=details,upper_seconds=2*max(estimates),seconds=seconds,
           eligible=all(r['passed'] for r in rows) and max(errors)<.02 and 2*max(estimates)<650)
    write(OUT/'pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=650,
        paths=[[n,r,f'pilot-{n}-{r}'] for n,r in PATHS],
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),Path(matter.__file__),Path(execution.__file__),OUT/'plan.json',OUT/'pilot.json']}))


def production():
    production_owner();angular_owner()
    p=json.loads((OUT/'audit.json').read_text());p.update(corrected_material_motion_applied=True,
        new_transfers_returned_to_free_material=False,new_material_sources_returned_to_GR=False)
    write(OUT/'audit.json',p)


if __name__=='__main__':globals()[sys.argv[1]]()
