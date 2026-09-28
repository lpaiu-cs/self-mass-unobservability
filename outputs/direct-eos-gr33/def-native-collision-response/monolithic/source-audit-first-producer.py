"""Read saved moments independently; export the additional GR stress source."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_monolithic_response as run

OUT=run.OUT;old=run.old;write=run.write;sha=run.sha


def main():
    assert not (OUT/'source-audit.json').exists();start=time.monotonic()
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',
        claim='Independently reconstruct saved energy/H invariants, local native thermochemical increments and the additional gas/photon GR stress moments.',
        budget=dict(seconds=25,new_native_bank_calls=0,new_evolution_steps=0),
        gates=dict(energy=1e-8,species=1e-8,pressure_reconstruction=1e-10,small_log_increment=1e-6),
        subtraction='Photon moments subtract the canonical volume/momentum contribution already present in Phase122. Collision-generated material moments are additional to its canonical gas response.',
        limits='No material advection, momentum response, exterior/deep scalar or GR re-evolution is inferred from this export. No final charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(old.__file__),OUT/'result.json',OUT/'steps-128-reference-128.npz']}))
    signal.signal(signal.SIGALRM,old.flow.old.optical.timeout);signal.alarm(25)
    m=old.Response(128);z=np.load(OUT/'steps-128-reference-128.npz');indices=[int(np.argmin(abs(z['t']-t))) for t in m.t]
    assert np.max(abs(z['t'][indices]-m.t))<1e-18
    rows=z['moments'][indices];x=z['delta_packet_scaled_occupation'];gas=z['delta_material'];units=np.load(old.OUT/'bank-128/units.npz')
    energy=np.sum(x*m.weights*m.E)+np.sum(gas[:,0]*units['energy'])
    number=np.sum(x*m.weights)-np.sum(gas[:,1]*units['neutral'])
    norms=[max(np.sum(abs(x)*m.weights*m.E),np.sum(abs(gas[:,0])*units['energy']),1.),max(np.sum(abs(x)*m.weights),np.sum(abs(gas[:,1])*units['neutral']),1.)]
    defects=[float(abs(energy-z['ledger'][1])/norms[0]),float(abs(number-z['ledger'][0])/norms[1])]
    dt=[];dy=[];pr=[];trace=[];pressure_error=0.;unsupported=0.
    for k in range(17):
        c=dict(np.load(old.OUT/f'bank-128/point-{k}.npz'));active=c['dt_energy']>0
        yy=np.divide(rows[k,2],c['neutral'],out=np.zeros(m.n),where=c['neutral']>0)
        tt=np.divide(rows[k,1]-c['dy_energy']*yy,c['dt_energy'],out=np.zeros(m.n),where=active)
        pp=(c['dt_p']*tt+c['dy_p']*yy)*m.volume
        pressure_error=max(pressure_error,float(np.sum(abs(pp-rows[k,6]))/max(np.sum(abs(rows[k,6])),1.)))
        unsupported=max(unsupported,float(np.sum(abs(rows[k,1,~active]))/max(np.sum(abs(rows[k,1])),1.)))
        energyproper=rows[k,1]/m.a
        radial=pp.copy();v=c['beta'][m.nb:]
        radial[m.nb:]=v*v*energyproper[m.nb:]+(1+v*v)*pp[m.nb:]
        dt.append(tt);dy.append(yy);pr.append(radial);trace.append(energyproper-radial-2*pp)
    E0=np.sum(m.I*m.weights[None]*m.E,axis=(2,3));P0=np.sum(m.I*m.weights[None]*m.E*m.model.bulk.mu2[None,None,:,None],axis=(2,3))
    mu4=(m.edges_mu[1:]**5-m.edges_mu[:-1]**5)/(5*np.diff(m.edges_mu))
    R4=np.sum(m.I*m.weights[None]*m.E*mu4[None,None,:,None],axis=(2,3))
    u=m.g['delta_u'];lam=m.g['delta_lambda'];canonicalE=-u*E0-lam*P0;canonicalP=-u*P0+lam*(R4-2*P0)
    additionalE=(rows[:,0]-canonicalE)/m.a;additionalP=(rows[:,5]-canonicalP)/m.a
    np.savez_compressed(OUT/'additional-stress-128.npz',t=m.t,radius_E=m.r,volume=m.volume,reference_lapse=m.a,
        additional_gas_energy_erg=rows[:,1]/m.a,additional_gas_radial_pressure_erg=pr,additional_gas_trace_erg=trace,
        additional_photon_energy_erg=additionalE,additional_photon_radial_pressure_erg=additionalP,
        canonical_photon_energy_ref_erg=canonicalE,canonical_photon_radial_pressure_ref_erg=canonicalP,
        cumulative_material_momentum_impulse_ref_erg=rows[:,3],delta_logT_collision=dt,delta_log_neutral_collision=dy)
    result=dict(classification='Counterexample candidate',passed=max(defects)<1e-8 and pressure_error<1e-10 and unsupported==0 and max(np.max(abs(dt)),np.max(abs(dy)))<1e-6,
        independent_endpoint_energy_relative=defects[0],independent_endpoint_photon_minus_neutral_relative=defects[1],pressure_reconstruction_relative=pressure_error,inactive_material_response_fraction=unsupported,
        max_delta_logT=float(np.max(abs(dt))),max_delta_log_neutral=float(np.max(abs(dy))),
        maximum_material_trace_L1_erg=float(np.max(np.sum(abs(trace),axis=1))),maximum_additional_photon_energy_L1_erg=float(np.max(np.sum(abs(additionalE),axis=1))),
        maximum_additional_photon_radial_pressure_L1_erg=float(np.max(np.sum(abs(additionalP),axis=1))),
        endpoint_material_trace_erg=float(np.sum(trace[-1])),seconds=time.monotonic()-start,
        additional_material_motion_evolved=False,additional_GR_source_exported=True,additional_GR_source_applied=False,final_charge_solved=False)
    write(OUT/'source-audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
