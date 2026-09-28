"""Audit saved motion/transfer identities and the actual finite return mismatch."""
from pathlib import Path
import json,signal,time
import numpy as np
import def_native_finite_motion_feedback as run


def main():
    used=json.loads((run.RETURN/'sources.json').read_text())['seconds'];budget=min(40,int(90-used));assert budget>0
    run.write(run.OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Independently reconstruct photon moments and the once-only noncollision transport identity, verify all saved material balances, and measure residual against time-comparison scale.',
        budget_seconds=budget,source_seconds_used=used,source_plus_audit_cap=90,
        reuse='Only completed saved arrays; no trajectory, coefficient bank, native inverse or extra waveform iteration.',
        limits='Finite comparison and exact stored bookkeeping, not continuous EOS, coupling or nonlinear error enclosure.',
        bindings={str(p):run.sha(p) for p in [Path(__file__),Path(run.__file__),run.OUT/'result.json',run.RETURN/'production.json',run.RETURN/'sources.json',run.GR/'result.json']}))
    started=time.monotonic();signal.alarm(budget);run.configure();m=run.Response()
    LD=np.longdouble;paths=[];states=[];photon=[];old_residual=[];prefix=[]
    background=m.material.point(16);Q=background['Q'];units=np.maximum(abs(Q),1.);units[1]=np.maximum(Q[0]*run.C**2,1.)
    for n in [64,128]:
        p=dict(np.load(run.photon_path(n,128)));d=dict(np.load(run.RETURN/f'steps-{n}-reference-128.npz'))
        before=dict(np.load(run.matter.OUT/f'steps-{n}-reference-128.npz'))
        oldp=dict(np.load(run.prior.OUT/f'steps-{n}-reference-128.npz'))
        ids=np.array([int(np.argmin(abs(d['t']-t))) for t in p['t']]);assert np.max(abs(d['t'][ids]-p['t']))<1e-18
        assert len(d['t'])==n+1 and len(p['t'])==17
        z=d['history_scaled'][ids].astype(LD)*LD(run.AMP);oldz=before['history_scaled'][ids].astype(LD)*LD(run.AMP)
        # Separate old noncollision transport from the NEW paired collisions.
        expected=oldz[:,[2,3]]-oldp['collision_transfer'].transpose(0,2,1)
        actual=p['moments'][:,[1,2]].astype(LD)-p['collision_transfer'].transpose(0,2,1)
        scale=np.maximum(np.max(np.sum(abs(p['moments'][:,[1,2]]),axis=2),axis=0),1.)
        mech=np.max(np.sum(abs(actual-expected),axis=2),axis=0)/scale
        weights=m.Eweight/m.scale
        raw=p['photon_history_scaled_occupation'];energy=np.sum(raw*weights,axis=(2,3))
        radial=np.sum(raw*weights*m.model.bulk.mu2[None,None,:,None],axis=(2,3))
        moment=max(float(np.max(np.sum(abs(v-p['moments'][:,i]),axis=1))/max(np.max(np.sum(abs(p['moments'][:,i]),axis=1)),1.)) for v,i in [(energy,0),(radial,5)])
        balance=float(np.max(abs(d['history_scaled'].astype(LD).sum(2)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        residual=np.max(np.sum(abs(p['moments'][:,[1,2]]-z[:,[2,3]]),axis=2),axis=0)
        normalized=residual/np.maximum(np.max(np.sum(abs(z[:,[2,3]]),axis=2),axis=0),1.)
        olderr=np.max(np.sum(abs(oldp['moments'][:,[1,2]]-oldz[:,[2,3]]),axis=2),axis=0)
        old_residual.append(olderr);states.append(z);photon.append(p['moments'][:,[1,2]])
        pilot=dict(np.load(run.RETURN/f'pilot-{n}.npz'));assert np.array_equal(d['history_scaled'][:3],pilot['history_scaled'])
        prefix.append(True)
        relative=np.asarray(abs(z[-1])/units,float);relative[:,~background['active']]=0.
        field,cell=np.unravel_index(np.argmax(relative),relative.shape)
        peak=dict(component=['baryon','momentum','reference_energy','neutral_H'][field],cell=int(cell),relative=float(relative[field,cell]),
                  physical_change=float(z[-1,field,cell]),background=float(Q[field,cell]),deep=bool(cell<m.nb))
        paths.append(dict(steps=n,mechanical_once_relative=np.asarray(mech,float).tolist(),photon_moment_relative=moment,endpoint_peak_relative_state=peak,
            full_material_balance=balance,energy_H_residual=np.asarray(normalized,float).tolist(),
            absolute_residual=np.asarray(residual,float).tolist(),previous_absolute_residual=np.asarray(olderr,float).tolist()))
    time_material=np.max(np.sum(abs(states[0][:,[2,3]]-states[1][:,[2,3]]),axis=2),axis=0)
    time_photon=np.max(np.sum(abs(photon[0]-photon[1]),axis=2),axis=0)
    numerical=np.maximum(time_material,time_photon);res=np.array(paths[1]['absolute_residual'])
    ratio=np.asarray(res/np.maximum(numerical,LD(1e-300)),float)
    change=np.asarray(res/np.maximum(old_residual[1],LD(1e-300)),float)
    checked=0
    for path,h in json.loads((run.OUT/'plan.json').read_text())['bindings'].items():assert run.sha(path)==h,path;checked+=1
    gr=json.loads((run.GR/'result.json').read_text());assert gr['passed'] and gr['seconds']<90
    assert max(v for p in paths for v in p['mechanical_once_relative'])<1e-8
    assert max(p['photon_moment_relative'] for p in paths)<1e-12 and max(p['full_material_balance'] for p in paths)<1e-8
    result=dict(classification='Counterexample candidate',passed=True,paths=paths,plan_bindings=checked,material_pilot_prefix_exact=prefix,
        absolute_material_time_difference=np.asarray(time_material,float).tolist(),absolute_photon_time_difference=np.asarray(time_photon,float).tolist(),
        residual_over_temporal_comparison=ratio.tolist(),residual_change_ratio=change.tolist(),
        residual_below_temporal_comparison_scale=bool(max(ratio)<1),
        actual_native_forcing_and_finite_motion_in_photons=True,new_transfers_in_finite_material_and_GR=True,
        full_nonlinear_radiation=False,uniform_EOS_derivative_bound=False,uniform_coupled_error_bound=False,
        exterior_floor_feedback_closed=False,coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False,
        seconds=time.monotonic()-started)
    run.write(run.OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    signal.signal(signal.SIGALRM,run.native.forcing.history.flow.old.optical.timeout)
    main()
