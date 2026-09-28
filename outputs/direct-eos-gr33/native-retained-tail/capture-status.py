"""Preserve accepted prefix data and classify the stopped execution."""
from pathlib import Path
import io,json,time
import numpy as np
import sympy as sp

ROOT=Path('outputs/direct-eos-gr33');OUT=ROOT/'native-retained-tail';S=OUT/'supported-temperature';EV=S/'evolution'
read=lambda p:json.loads(p.read_text())
load=lambda p:np.load(io.BytesIO(p.read_bytes()))
start=time.monotonic();coarse=load(EV/'checkpoint-64.npz');fine=load(EV/'checkpoint-128.npz');hist=load(EV/'history-128.npz')
assert int(coarse['completed'])==64 and int(fine['completed'])==16 and len(hist['t'])==33
assert read(S/'production-receipt.json')['returncode']==124
assert read(EV/'result-64.json')['passed'] and read(EV/'capture-64.json')['passed']
z=load(EV/'state-128-16.npz')
keys=['U','I','u','theta','eta','Pi','h','j']
assert all(np.array_equal(z[k],fine[k]) for k in keys)
history=read(EV/'tail-128.json')['rows'];assert len(history)==33
assert len(load(EV/'accepted-ports-128.npz')['angular_luminosity'])==32
assert abs(float(hist['t'][16])-float(hist['t'][-1])/2)<1e-18
# Source export is the original capture identity, from existing accepted
# moments. This creates no new fluid state or extrapolated history.
d=dict(load(EV/'source-64.npz'));a=hist['moments']
for i,key in enumerate(['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','nonrest_stress_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg']):d[key]=a[:,i]
C=2.99792458e10
d.update(t=hist['t'],inner_cumulative_energy_erg=hist['ports'][:,0],outer_cumulative_energy_erg=hist['ports'][:,1],
    metric_stress_erg=a[:,0]*np.longdouble(d['cx'])*np.longdouble(C)**2+a[:,3]+a[:,5]-a[:,6])
stream=io.BytesIO();np.savez_compressed(stream,**d)
temp=EV/'source-prefix-128-32.saving';temp.write_bytes(stream.getbuffer());temp.replace(EV/'source-prefix-128-32.npz')
old=load(ROOT/'def-native-conservative-rates/thermal-refined/evolution/history-64.npz')
current=load(EV/'history-64.npz');ratio=float(current['discard'][-1,0]/old['discard'][-1,0])
bank=load(S/'support/repaired-bank.npz');scale=4*np.pi*float(d['RJ'])**2*float(bank['rho0'])
tails=read(EV/'tail-64.json')['rows'];assert len(tails)==65
D,v,r,tau,p,cx=sp.symbols('D v r tau p cx',real=True)
internal=(tau+p-cx*D*v*v/(r*(1+r)))*r*r-p
assert sp.simplify((-cx*D*v*v/(1+r)+internal-3*p-(tau-v*v*(cx*D+tau+p)-3*p)).subs(v*v,1-r*r))==0
result=dict(classification='Counterexample candidate',overall_phase_passed=False,
    contribution='Loophole progress: native low-density matter is retained in actual coupled photon/material evolution for the complete coarse horizon; full two-clock/GR conclusion remains open.',
    native_bank_passed=True,native_controls_passed=True,full_coarse_steps=64,full_fine_steps=128,
    coarse_completed=64,fine_recorded_physical_step=32,fine_complete_restart_step=16,
    fine_complete_restart_bitwise_against_canonical=True,fine_prefix_source_exported=32,
    truncated_checkpoint_recovery_passed=read(S/'checkpoint-recovery.json')['passed'],
    retained_tail_endpoint_g=tails[-1]['mass_g'],maximum_retained_tail_g=max(z['mass_g'] for z in tails),
    minimum_retained_tail_temperature_K=min(z['minimum_temperature_K'] for z in tails if z['cells']),
    coarse_discarded_baryon_g=float(current['discard'][-1,0]*scale),
    old_coarse_discarded_baryon_g=float(old['discard'][-1,0]*scale),coarse_remaining_discard_ratio=ratio,
    coarse_conserved_Killing_balance=read(EV/'result-64.json')['conserved_Killing_balance'],
    physical_density_floor_cgs=float(bank['rho0']*np.exp(bank['x'][0])),
    production=read(S/'production-receipt.json'),GR_attempt=read(S/'coarse_gr-receipt.json'),
    source_complete_on_both_clocks=False,new_GR_charge_available=False,uniform_EOS_derivative_bound=False,
    full_floor_feedback_enclosed=False,coupled_fixed_point_verified=False,nonlinear_GR=False,
    final_charge_solved=False,full_goal_complete=False,
    original_trace_identity_symbolically_rechecked=True,post_run_storage_repair_test=read(S/'checkpoint-control.json'),
    scientific_action_timer_sum=read(S/'resumed-pilot.json')['aggregate_spent_seconds']+read(S/'production-receipt.json')['seconds']+read(S/'coarse_gr-receipt.json')['seconds'],
    bookkeeping_seconds=time.monotonic()-start)
(OUT/'status.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
