"""Saved-state checks only; never start a new trajectory or call native EOS."""
from pathlib import Path
import inspect
import json
import time
import numpy as np
import def_resolved_scalar_pulse as s
import def_resolved_scalar_pulse_run as run
import def_companion_benchmark as benchmark


def audit():
    started=time.monotonic();plan=s.bindings()
    assert s.symbolic()==plan['symbolic']
    declared=json.loads((s.OUT/'companion-benchmark.json').read_text())
    assert declared.pop('source_sha256')==s.e.digest(Path(benchmark.__file__))
    assert declared==benchmark.define()
    result=json.loads((run.OUT/'result.json').read_text())
    for name,key in [('execution-plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((run.OUT/name).read_text())[key].items():
            assert s.e.digest(s.ROOT/rel)==digest,rel
    source=inspect.getsource(s.old.two.cones).replace('def cones(', 'def frame_cones(').replace('e.TAU',"z['tau_cond']")
    namespace=dict(vars(s.old.two));exec(compile(source,__file__,'exec'),namespace)
    cone=namespace['frame_cones'];rows=[];count=0
    initial=np.load(s.OUT/'initial.npz');V=initial['volume'];totalB=initial['baryons'].sum()
    for path in result['paths']:
        label,n=path['label'],path['steps'];folder=run.OUT/f'{label}-{n}'
        R=s.ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm']);h=4*R/s.e.C/n
        peak=s.ld(0) if label=='undriven' else s.ld('.0002')
        star,_=s.initialize(s.old.imported.CachedOnly(),0 if label=='decoupled' else -4,peak,h)
        previous_phi=np.zeros(star.n,dtype=s.ld)
        maxima=dict(cone_speed=0.,cone_imaginary=0.,baryon_inventory=0.,species_inventory=0.,native_norm=0.,
            direct_trace_over_energy_tolerance=0.,without_trace_native_norm=0.,baryon_weighted_scalar_rms=0.,innermost_scalar=0.)
        cones_pass=True
        for step in range(1,n+1):
            with np.load(folder/f'step-{step:03d}.npz') as z:
                c=cone(z);cones_pass &= c['sampled_cone_inside_light_cone']
                direct=z['gstar']*(z['psi']-previous_phi)/star.heat0
                values=dict(cone_speed=c['maximum_local_rest_characteristic_speed_over_c'],
                    cone_imaginary=c['maximum_characteristic_imaginary_part'],
                    baryon_inventory=float(abs(np.sum(V*z['dU'][:,0]))/totalB),
                    species_inventory=float(abs(np.sum(V[:,None]*z['dBX'],axis=0)).max()/totalB),
                    native_norm=float(np.max(abs(z['residual'])/s.ATOL)),
                    direct_trace_over_energy_tolerance=float(abs(direct).max()/s.ATOL[1]),
                    without_trace_native_norm=float(max(np.max(abs(z['residual'])/s.ATOL),abs(z['residual'][:,1]-direct).max()/s.ATOL[1])),
                    baryon_weighted_scalar_rms=float(np.sqrt(np.sum(initial['baryons']*z['psi']**2)/totalB)),
                    innermost_scalar=float(abs(z['psi'][0])))
                for k,v in values.items():maxima[k]=max(maxima[k],v)
                previous_phi=z['psi'].copy()
            count+=1
        # Reconstruct the actual 31 endpoint equations using only exact cached
        # native rows; CachedOnly raises if any unrecorded EOS input is needed.
        with np.load(folder/f'step-{n-1:03d}.npz') as cp:p={k:cp[k].copy() for k in cp.files}
        with np.load(folder/f'step-{n:03d}.npz') as cp:z={k:cp[k].copy() for k in cp.files}
        star.previous=p;star.field_guess=[z[k].copy() for k in ['psi','Pi','Phi']]
        star.previous_boundary=star.amplitude=s.ld(0)
        y=star.base+z['delta'];la=z['logA']
        star.material_cache.update({s.e.material_key(row):raw for row,raw in
            zip(zip(y[:,0]-3*la,y[:,1]-la,y[:,5:]),z['raw'])})
        previous=(p['delta'],p)
        value,replayed=s.residual(star,z['delta'],previous,previous,h,(s.ld(1),s.ld(-1),s.ld(0)))
        replay_norm=float(np.max(abs(value)/s.ATOL))
        replay_difference=float(np.max(abs(value-z['residual'])/s.ATOL))
        passed=bool(cones_pass and maxima['baryon_inventory']<1e-9 and maxima['species_inventory']<1e-9
            and maxima['native_norm']<=1 and replay_norm<=1 and replay_difference<1e-5)
        rows.append(dict(label=label,steps=n,passed=passed,maxima=maxima,
            endpoint_replay_norm=replay_norm,endpoint_replay_difference=replay_difference,
            endpoint_scalar_residual=replayed['scalar_residual']))
    recomputed=run.summarize(result['paths'])
    assert recomputed['readouts']==result['readouts'] and recomputed['passed']==result['passed']
    return dict(classification='Counterexample candidate',passed=all(x['passed'] for x in rows),
        time_comparison_passed=result['passed'],saved_states=count,paths=rows,seconds=time.monotonic()-started,
        source_sha256=s.e.digest(Path(__file__)),scope='All saved native residual arrays, inventories and sampled heat cones; nine cached-native endpoint equation replays. No new trajectory or native call, spatial or continuous error certificate.')


if __name__=='__main__':
    target=s.OUT/'saved-audit.json';assert not target.exists()
    result=audit();s.e.write(target,result);print(json.dumps(result));assert result['passed']
