"""Read the actually retained material into the same initial GR operator."""
from pathlib import Path
from types import FunctionType
import json,sys,time,signal
import numpy as np
import sympy as sp
import evolve_native_retained_tail as run
import verify_native_stage_energy_charge as original

OUT=run.OUT;EV=run.EV;GR=run.GR;read=run.read;write=run.write;sha=run.sha


def main():
    assert read(EV/'result.json')['passed'] and not (GR/'result.json').exists()
    spent=read(OUT/'resumed-pilot.json')['aggregate_spent_seconds']+read(EV/'result.json')['seconds']
    cap=min(90.,1300-spent)
    assert cap>=40,'Insufficient aggregate budget for GR readout'
    write(GR/'plan.json',dict(classification='Counterexample candidate',hard_seconds=cap,
        claim='Apply the completed lower-floor actual photon/material source histories to the identical initial characteristic GR operator and compare to the stored thermal-refined64/128 background.',
        physics='Only the original numerical deletion support changes. Previous native pressure/collision corrections on the older background are not silently inherited on the changed trajectory. This is a same-operator compact first-variation readout, not a new final null-infinity/fully nonlinear result.',
        gates=dict(energy=1e-8,baryon=1e-10,time=.02,quadrature=.002,independent=1e-9),
        budget='Use at most90s of the actual unspent revised1300s aggregate. No further evolution, EOS state or resolution.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in [Path(__file__),Path(run.__file__),
            EV/'result.json',EV/'source-64.npz',EV/'source-128.npz',run.physical.GR/'wave-128-g8.npz',run.physical.GR/'wave-64-g8.npz']}))
    start=time.monotonic();run.support.owner.deadline(start,cap)
    model=original.gr.Response();paths={}
    for n,order in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(EV/f'source-{n}.npz'));wave=original.read(model,d,order);paths[n,order]=wave
        np.savez_compressed(GR/f'wave-{n}-g{order}.npz',**wave)
        run.support.owner.deadline(start,cap)
    fine=paths[128,8];norm=max(np.max(abs(fine['free_scalar'])),1e-300)
    time_error=float(np.max(abs(fine['free_scalar'][::2]-paths[64,8]['free_scalar']))/norm)
    quadrature=float(np.max(abs(fine['free_scalar']-paths[128,4]['free_scalar']))/norm)
    d=dict(np.load(EV/'source-128.npz'));direct,residual=original.independent.direct(model,d,8)
    agreement=abs(direct/fine['direct_scalar'][-1]-1)
    comparisons=[];changes={};inventories=[]
    for n in [64,128]:
        old=np.load(run.physical.GR/f'wave-{n}-g8.npz');now=paths[n,8];delta=now['free_scalar']-old['free_scalar'];changes[n]=delta
        z=np.load(EV/f'checkpoint-{n}.npz');assert int(z['completed'])==n
        history=np.load(EV/f'history-{n}.npz');previous=np.load(run.physical.EV/f'history-{n}.npz')
        current_floor=float((history['discard'][-1,2]+model.model.m.a0*model.model.cx*history['discard'][-1,0])*model.model.gas_scale)
        prior_floor=float((previous['discard'][-1,2]+model.model.m.a0*model.model.cx*previous['discard'][-1,0])*model.model.gas_scale)
        tail=read(EV/f'tail-{n}.json')['rows'];assert len(tail)==n+1
        comparisons.append(dict(steps=n,scalar_endpoint=float(now['free_scalar'][-1]),prior_endpoint=float(old['free_scalar'][-1]),
            endpoint_change=float(delta[-1]),change_relative_to_prior=float(delta[-1]/old['free_scalar'][-1]),
            discarded_Killing_energy_erg=current_floor,prior_discarded_Killing_energy_erg=prior_floor,
            remaining_discard_ratio=current_floor/prior_floor,
            maximum_retained_tail_mass_g=max(v['mass_g'] for v in tail),
            maximum_retained_tail_pressure_volume_erg=max(v['pressure_volume_erg'] for v in tail),
            maximum_retained_tail_speed_c=max(v['maximum_speed_c'] for v in tail),
            maximum_retained_tail_cells=max(v['cells'] for v in tail),
            minimum_retained_tail_temperature_K=min(v['minimum_temperature_K'] for v in tail if v['cells'])))
        src=np.load(EV/f'source-{n}.npz');defect=abs(np.sum(src['baryon_g'][-1],dtype=np.longdouble)+z['discard'][0]*model.model.gas_scale/run.C**2)
        f=model.model.flow;m=model.model.m;initial=np.sum(f.initial[0]*m.vol,dtype=np.longdouble)*model.model.gas_scale/run.C**2
        inventories.append(dict(steps=n,baryon_relative_to_initial=float(defect/initial)))
        assert len(np.load(EV/f'accepted-ports-{n}.npz')['angular_luminosity'])==n
    delta_norm=max(np.max(abs(changes[128])),1e-300)
    delta_time=float(np.max(abs(changes[128][::2]-changes[64]))/delta_norm)
    old_quad=np.load(run.physical.GR/'wave-128-g4.npz')['free_scalar']
    delta_quad=float(np.max(abs((paths[128,4]['free_scalar']-old_quad)-changes[128]))/delta_norm)
    # Linear response superposition is an identity of this fixed operator;
    # it gives no contraction/continuum certificate for the coupled model.
    A,x,y=sp.symbols('A x y');assert sp.expand(A*(x+y)-A*x-A*y)==0
    run.support.owner.deadline(start,cap)
    passed=time_error<.02 and quadrature<.002 and agreement<1e-9 and max(r['baryon_relative_to_initial'] for r in inventories)<1e-10
    result=dict(classification='Counterexample candidate',passed=bool(passed),comparisons=comparisons,baryon=inventories,
        time_relative=time_error,quadrature_relative=quadrature,independent_direct_relative=float(agreement),inverse_residual_cm=float(residual),
        tail_change_time_relative=delta_time,tail_change_quadrature_relative=delta_quad,
        tail_change_resolved_on_original_gates=delta_time<.02 and delta_quad<.002,
        retained_native_material_full_horizon=True,actual_current_sources_applied_to_compact_GR=True,
        original_failed_pilot_passed=False,uniform_EOS_derivative_bound=False,full_floor_feedback_enclosed=False,
        coupled_fixed_point_verified=False,new_emission_at_null_infinity=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
        symbolic_linearity_check=True,seconds=time.monotonic()-start,aggregate_scientific_seconds=spent+time.monotonic()-start)
    write(GR/'result.json',result);signal.setitimer(signal.ITIMER_REAL,0.);print(json.dumps(result),flush=True)


def coarse_gr():
    assert read(EV/'result-64.json')['passed'] and not (GR/'coarse-result.json').exists()
    spent=read(OUT/'resumed-pilot.json')['aggregate_spent_seconds']+read(OUT/'production-receipt.json')['seconds']
    cap=min(25.,1300-spent);assert cap>0
    write(GR/'coarse-plan.json',dict(classification='Counterexample candidate',hard_seconds=cap,
        claim='Apply the completed64-step retained-material source to the same compact first-variation GR operator;128-step evolution and time acceptance remain incomplete.',
        budget='Only the actual unspent portion of the revised1300s aggregate, at most25s. No new fluid path or EOS grid.',
        bindings={str(p):sha(p) for p in [Path(__file__),EV/'source-64.npz',EV/'result-64.json',OUT/'production-receipt.json']}))
    start=time.monotonic();run.support.owner.deadline(start,cap)
    try:
        model=original.gr.Response();run.support.owner.deadline(start,cap);d=dict(np.load(EV/'source-64.npz'))
        wave=original.read(model,d,8);np.savez_compressed(GR/'coarse-wave-g8.npz',**wave)
        q=original.read(model,d,4);run.support.owner.deadline(start,cap)
        previous=np.load(run.physical.GR/'wave-64-g8.npz');delta=wave['free_scalar']-previous['free_scalar'];norm=max(abs(wave['free_scalar']).max(),1e-300)
        direct,residual=original.independent.direct(model,d,8);agreement=float(abs(direct/wave['direct_scalar'][-1]-1))
        quad=float(max(abs(wave['free_scalar']-q['free_scalar']))/norm)
        history=np.load(EV/'history-64.npz');before=np.load(run.physical.EV/'history-64.npz');m=model.model
        def loss(z):return float((z['discard'][-1,2]+m.m.a0*m.cx*z['discard'][-1,0])*m.gas_scale)
        result=dict(classification='Counterexample candidate',passed=quad<.002 and agreement<1e-9,
            coarse_endpoint=float(wave['free_scalar'][-1]),prior_same_clock_endpoint=float(previous['free_scalar'][-1]),
            endpoint_change=float(delta[-1]),relative_change=float(delta[-1]/previous['free_scalar'][-1]),
            quadrature_relative=quad,independent_direct_relative=agreement,inverse_residual_cm=float(residual),
            discarded_Killing_energy_erg=loss(history),prior_discarded_Killing_energy_erg=loss(before),
            remaining_discard_ratio=loss(history)/loss(before),seconds=time.monotonic()-start,
            actual_completed_coarse_source_applied_to_GR=True,coarse_fine_time_acceptance=False,
            uniform_EOS_derivative_bound=False,full_floor_feedback_enclosed=False,new_emission_at_null_infinity=False,
            final_charge_solved=False,full_goal_complete=False)
        run.support.owner.deadline(start,cap);write(GR/'coarse-result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(GR/'coarse-failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


if __name__=='__main__':globals()[sys.argv[1] if len(sys.argv)>1 else 'main']()
