"""Compare actual finite-radius stress histories without promoting them to charge."""
from pathlib import Path
import json, resource, time
import numpy as np
import propagate_exterior_vacuum as adapter
p=adapter.base

out=p.OUT;start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(600)
try:
    rows={k:p.read(out/k/'result.json')['rows'] for k in p.SETTINGS}
    assert all(len(v)==16 for v in rows.values())
    for k in p.SETTINGS:assert p.read(out/f'{k}-receipt.json')['error'] is None
    geometry={k:np.array([r['stress_weak_moments_erg'] for r in v],dtype=p.LD) for k,v in rows.items()}
    physical={k:np.array([r['high_physical_stress_weak_moments_erg'] for r in v],dtype=p.LD) for k,v in rows.items()}
    controls={}
    for family,values in [('geometry',geometry),('physical',physical)]:
        norm=np.maximum(np.max(abs(values['fine']),axis=0),p.LD('1e-290'))
        controls[family]={k:np.asarray(np.max(abs(values[k]-values['fine']),axis=0)/norm,float).tolist() for k in ['angular','temporal','geometry']}
    p.initialize();m=p.Photons(8,8);time_rows=[];time_path=[]
    for j,now in enumerate(m.clock[1:],1):
        val,row=m.reference('high',64,float(now),8)
        val=val+np.load(out/f'fine/snapshot-{j}.npz')['stress_weak_moments_erg']
        time_path.append(val);time_rows.append(row)
    time_path=np.asarray(time_path,p.LD)
    norm=np.maximum(np.max(abs(physical['fine']),axis=0),p.LD('1e-290'))
    time_control=np.asarray(np.max(abs(time_path-physical['fine']),axis=0)/norm,float)
    passed=all(np.max(v)<.002 for family in controls.values() for v in family.values()) and np.max(time_control)<.02
    result=dict(classification='Counterexample candidate',passed=bool(passed),controls=controls,time_control=time_control.tolist(),
        actual_same_solution_emission_applied=True,incident_metric_photon_worldlines_propagated=True,
        physical_launch_energy_and_propagation_work_separated=True,finite_radius_stress_source_exported=True,
        terminal=rows['fine'][-1],reference_time_ports=time_rows,
        source_representation='Distributional packet worldlines and first variations; smooth weak energy/flux/pressure moments are controls, not a radial density reconstruction.',
        returned_metric_geometric_response_complete=False,reciprocal_scalar_work_complete=False,
        physical_GR_boundary_applied=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        remaining=['Apply the returned metric to background photon paths, with its own consistent exterior continuation.',
                   'Insert these physical stress distributions and reciprocal scalar stress into the same mass/lapse/scalar boundary and evolve its feedback.',
                   'Close background mass normalization, EOS/uniform errors, nonlinear/static-EFT/observational scope.'],
        bindings={str(q):p.sha(q) for q in [Path(__file__),Path(p.__file__),Path(adapter.__file__),out/'plan.json',out/'pilot.json']+[out/k/'result.json' for k in p.SETTINGS]})
    np.savez_compressed(out/'stress-history.npz',t=m.clock[1:],physical_fine=physical['fine'],physical_coarse=time_path,geometry_fine=geometry['fine'])
    p.write(out/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['bindings','reference_time_ports']}),flush=True)
    assert passed,controls
except BaseException as exc:error=repr(exc);raise
finally:p.write(out/'collect-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=p.sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
