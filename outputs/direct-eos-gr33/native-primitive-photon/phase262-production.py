"""Apply the verified returned metric to the existing exterior photon cohorts."""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import propagate_returned_moments as adapter
b=adapter.base; p=b.p; out=adapter.OUT
start=time.monotonic(); error=None; action=sys.argv[1]
resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2); p.incident.native.deadline(14400)
try:
    for f,h in b.read(out/'execution-plan.json')['bindings'].items():assert b.sha(f)==h,f
    assert b.read(out/'pilot.json')['eligible']
    if action=='collect':
        rows={k:b.read(out/k/'result.json')['rows'] for k in p.SETTINGS}
        assert all(len(v)==16 for v in rows.values())
        arrays={k:np.asarray([r['stress_weak_moments_erg'] for r in v],p.LD) for k,v in rows.items()}
        norm=np.maximum(np.max(abs(arrays['fine']),axis=0),p.LD('1e-290'))
        errors={k:np.asarray(np.max(abs(v-arrays['fine']),axis=0)/norm,float).tolist() for k,v in arrays.items() if k!='fine'}
        passed=all(np.max(v)<.002 for v in errors.values())
        result=dict(classification='Counterexample candidate',passed=passed,controls=errors,
            terminal=rows['fine'][-1],returned_metric_geometric_photon_response_computed=True,
            same_applied_outgoing_continuation=True,returned_source_time_error='Existing applied metric time gate retained; new exterior geometry has no separate64/128trajectory control in this execution.',
            reciprocal_scalar_source_complete=False,physical_GR_boundary_applied=False,physical_final_charge_solved=False,
            final_charge_conclusion='unadjudicated',full_goal_complete=False)
        np.savez_compressed(out/'stress-history.npz',geometry_fine=arrays['fine'])
        b.write(out/'result.json',result); assert passed,errors
    else:
        a,t,g=p.SETTINGS[action]; b.Moments=adapter.Moments; p.Metric=b.Metric
        p.initialize(); m=p.Photons(a,g); folder=out/action; folder.mkdir(); rows=[]
        for j,now in enumerate(m.clock[1:],1):
            z,row=m.propagate(float(now),t)
            reference,port=m.reference('low',128,float(now),t)
            z['low_reference_stress_weak_moments_erg']=reference
            z['low_physical_stress_weak_moments_erg']=reference+z['stress_weak_moments_erg']
            r=z['radius_cm']; mu=z['direction']; dr=z['delta_radius_cm']; dm=z['delta_direction']
            _,N,B,A,_=m.d.bg.metric(r/m.d.model.m.R)
            local=z['delta_log_H']-z['metric_nu']-z['metric_lambda']; E=z['background_packet_energy_erg']
            fac=p.LD(p.G)/p.LD(p.C)**4
            row.update(reference_port=port,component='returned_metric_geometry',
                photon_J_source_cm=float(fac*np.sum(E*local,dtype=p.LD)),
                photon_lapse_particular=float(-fac*np.sum(E/(r*A)*((1+mu*mu)*(local-dr/(r*B))+2*mu*dm),dtype=p.LD)),
                physical_final_charge_solved=False)
            np.savez_compressed(folder/f'snapshot-{j}.npz',**z); b.write(folder/f'snapshot-{j}.json',row); rows.append(row)
            b.write(folder/'progress.json',dict(completed=j,total=16,latest=row))
        b.write(folder/'result.json',dict(classification='Counterexample candidate',rows=rows))
except BaseException as exc:error=repr(exc); raise
finally:
    b.write(out/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=b.sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
