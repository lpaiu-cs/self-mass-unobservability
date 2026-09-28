"""Localize the actual time-gate failure using saved same-solution terms."""
import time
started=time.monotonic()
import numpy as np
import complete_native_fluid_time as run
pieces=[]
for n in [64,128]:
    with np.load(run.OUT/f'sweep-1/photons/pilot-{n}.npz') as p:
        b=p['conserved_material_history'][-1,0].astype(run.LD)/run.joint.AMP
        native=np.sum(p['joint_stage_weights'][:,None].astype(run.LD)*p['joint_native_rates_scaled'][:,:,2],axis=0,dtype=run.LD)
        discard=p['material_floor_discard_scaled'][:,2].astype(run.LD)
        assert not np.any(p['joint_collision_rates_scaled'][:,:,2])
        assert np.sum(abs(b+discard-native))/max(np.sum(abs(b)),run.LD('1e-290'))<1e-8
        pieces.append((b,native,discard,p['radius_E'].copy()))
coarse,fine=pieces;difference=coarse[0]-fine[0];norm=np.sum(abs(difference));order=np.argsort(abs(difference))[::-1]
native=coarse[1]-fine[1];discard=coarse[2]-fine[2]
assert np.sum(abs(difference-native+discard))/norm<1e-8
rows=[dict(cell=int(i),radius_cm=float(fine[3][i]),share=float(abs(difference[i])/norm),
    delta_baryon_scaled=float(difference[i]),native_integral_difference=float(native[i]),floor_difference=float(discard[i])) for i in order[:12]]
result=dict(classification='Counterexample candidate',passed=True,trajectory_replayed=False,rows=rows,
    original_time_gate_passed=False,relative=float(norm/np.sum(abs(fine[0]))),
    floor_difference_L1_over_total_difference=float(np.sum(abs(discard))/norm),
    native_difference_L1_over_total_difference=float(np.sum(abs(native))/norm),
    seconds=time.monotonic()-started,source_sha256=run.sha(__file__),
    scope='Exact decomposition of the saved endpoint baryon difference into its own native-flux integral and floor removal. It localizes terms but does not prove the root of their time errors.')
run.write(run.OUT/'baryon-localization.json',result);print(result)
