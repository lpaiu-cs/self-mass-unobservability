"""Recover stored GR verdict and apply native pressure+collision sources together."""
from pathlib import Path
import json,signal,sys,time
import numpy as np
import def_native_compensated_finite_return as run
import verify_native_stage_energy_charge as charge

OUT=run.OUT;GR=run.GR;write=run.write;sha=run.sha;LD=np.longdouble


def prepare():
    assert not (GR/'result.json').exists();assert json.loads((OUT/'sources.json').read_text())['passed']
    write(GR/'serialization-failure.json',dict(error='numpy boolean in passed could not be JSON serialized after all three GR waves and independent direct integration completed',
        physical_paths_replayed=False,GR_wave_replay=False,charged_GR_seconds=90))
    write(GR/'recovery-plan.json',dict(classification='Counterexample candidate',
        claim='Use the completed three charge arrays, recover the scalar JSON verdict, independently audit source identities, then apply the sum of original+native pressure+actual collision/material sources to the same GR operator.',
        reuse='All completed photons,finite material,pressure and GR arrays. No new trajectories,EOS calls or repeated three-wave propagation.',
        budget_seconds=60,reallocation='Reallocate60s of the77.399529066s unused original source budget. Charge failed GR finalization its full90s cap. Original source+GR aggregate180s is unchanged.',
        gates=dict(time=.02,quadrature=.002,independent=1e-9,source_identity=1e-12,combined_application=1e-10,materiality=.02),
        limits='Free compact first-variation GR on saved geometry, not the full exterior or nonlinear scalar solution. Fixed native forcing and finite free material; new motion has not been returned into photons.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),OUT/'production.json',OUT/'sources.json',GR/'source-128-reference-128.npz',GR/'wave-128-g8.npz',GR/'wave-64-g8.npz',GR/'wave-128-g4.npz']}))


def main():
    started=time.monotonic();signal.alarm(60);run.configure();history=run.old.old.repaired.forcing.history
    fine=dict(np.load(GR/'wave-128-g8.npz'));coarse=dict(np.load(GR/'wave-64-g8.npz'));quad=dict(np.load(GR/'wave-128-g4.npz'))
    norm=max(np.max(abs(fine['free_scalar'])),1e-300)
    terr=float(np.max(abs(fine['free_scalar']-coarse['free_scalar']))/norm);qerr=float(np.max(abs(fine['free_scalar']-quad['free_scalar']))/norm)
    m=run.old.old.previous.GRResponse();d=dict(np.load(GR/'source-128-reference-128.npz'))
    direct,coordinate=charge.independent.direct(m,d,8);ierr=float(abs(direct/fine['direct_scalar'][-1]-1))
    source_error=0.;balances=[]
    for n in [64,128]:
        a=dict(np.load(GR/f'source-{n}-reference-128.npz'));s=dict(np.load(OUT/f'stress-{n}-reference-128.npz'))
        stress=s['material'];rest=a['baryon_g'].astype(LD)*LD(a['cx'])*LD(run.C)**2
        error=[rest+a['gas_nonrest_energy_erg']-stress[:,0],rest+a['nonrest_trace_erg']-stress[:,2],
               rest+a['nonrest_stress_erg']-(stress[:,0]-stress[:,1]),a['pressure_volume_erg']-stress[:,3],
               a['metric_stress_erg']-(stress[:,0]+a['photon_energy_erg']-stress[:,1]-a['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(x)) for x in error)/max(np.max(abs(stress)),1.)))
        z=dict(np.load(OUT/f'steps-{n}-reference-128.npz'))
        balance=float(np.max(abs(z['history_scaled'].astype(LD).sum(2)+z['discards_scaled']-z['ledgers_scaled'])/np.maximum(z['norms_scaled'],1.)))
        balances.append(balance)
    prior=json.loads((history.OUT/'result.json').read_text());base=dict(np.load(history.physical.EV/'source-128.npz'));eos=dict(np.load(history.OUT/'source.npz'))
    t=base['t'];j=np.clip(np.searchsorted(d['t'],t,side='right')-1,0,len(d['t'])-2);alpha=((t-d['t'][j])/(d['t'][j+1]-d['t'][j])).astype(LD)
    keys=history.physical.capture.KEYS+['metric_stress_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg']
    for key in keys:
        v=d[key].astype(LD);f=alpha if v.ndim==1 else alpha[:,None]
        base[key]=base[key].astype(LD)+eos[key].astype(LD)+(1-f)*v[j]+f*v[j+1]
    applied=charge.read(m,base,8);previous=dict(np.load(history.OUT/'applied-charge.npz'))
    ids=np.array([int(np.argmin(abs(applied['t']-t))) for t in fine['t']]);assert np.max(abs(applied['t'][ids]-fine['t']))<1e-18
    oldnorm=max(np.max(abs(previous['free_scalar'])),1e-300)
    linear=float(np.max(abs(applied['free_scalar'][ids]-previous['free_scalar'][ids]-fine['free_scalar']))/oldnorm)
    materiality=float(norm/oldnorm)
    np.savez_compressed(GR/'combined-source.npz',**base);np.savez_compressed(GR/'combined-charge.npz',**applied)
    passed=terr<.02 and qerr<.002 and ierr<1e-9 and source_error<1e-12 and max(balances)<1e-8 and linear<1e-10 and materiality<.02
    result=dict(classification='Counterexample candidate',passed=bool(passed),endpoint_collision_charge=float(fine['free_scalar'][-1]),
        previous_native_pressure_endpoint=prior['native_source_endpoint'],combined_native_pressure_collision_endpoint=float(applied['free_scalar'][-1]),
        maximum_collision_over_previous_signal=materiality,time_relative=terr,quadrature_relative=qerr,independent_direct_relative=ierr,
        retarded_coordinate_residual=float(coordinate),source_identity_relative=source_error,independent_material_balances=balances,
        combined_direct_application_relative=linear,combined_endpoint_positive=bool(applied['free_scalar'][-1]>0),
        native_collision_propagated_through_finite_material_and_GR=True,pressure_check_kind='16eps arithmetic resolution diagnostic on actual finite primitive pressure; not a rigorous EOS enclosure',
        compact_potential_included=False,material_motion_returned_to_photons=False,uniform_EOS_derivative_bound=False,
        full_floor_exterior_feedback_enclosed=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-started)
    write(GR/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':
    signal.signal(signal.SIGALRM,run.old.old.repaired.forcing.history.flow.old.optical.timeout)
    try:prepare();main()
    except Exception as exc:write(GR/'recovery-failure.json',dict(error=repr(exc)));raise
