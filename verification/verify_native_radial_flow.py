"""Read-only audit of the evolved gas and its inner/photon exchange ports."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_radial_thermo as task

OUT=task.OUT


def main():
    assert not (OUT/'flow-audit.json').exists();start=time.monotonic();signal.alarm(30)
    task.task.write(OUT/'flow-audit-plan.json',dict(classification='Counterexample candidate',
        claim='Independently check the native initial gas state, integrated baryon balance and Lorentz-transformed photon work; export the exact accumulated inner transfer for the next bulk coupling.',
        gates=dict(initial_EOS=.002,baryon=1e-10,independent_scattering_work=1e-8),
        budget_seconds=30,native_calls=0,new_evolution_steps=0,
        bindings={str(p.relative_to(task.task.old.ROOT)):task.task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'cells-1792.npz',OUT/'cells-1792.json',OUT/'result.json']}))
    m=task.Flow(1792);d=np.load(OUT/'cells-1792.npz');row=json.loads((OUT/'cells-1792.json').read_text())
    U=d['U'];initial=d['initial'];history=d['history'];tail=d['baryon_discard'];C=task.task.C
    rho,v,sigma=m.primitive(U);p,u,gamma,T,kap=m.eos(rho,sigma)
    scale=4*np.pi*m.RJ**2*m.eos.rho0
    residual=np.sum((U[0]-initial[0])*m.vol)-history[-1,3]+tail[0]
    baryon=abs(residual)/np.sum(initial[0]*m.vol)
    inp=m.eos.d
    pp,uu,gg,tt,kk=m.eos(inp['initial_raw'][:,0]/m.eos.rho0,inp['initial_sigma'])
    init_errors=[float(max(abs(pp*m.eos.rho0*C*C/inp['initial_raw'][:,1]-1))),
        float(max(abs(uu*C*C/inp['initial_raw'][:,2]-1))),float(max(abs(tt/inp['initial_T']-1)))]
    # Independent angular quadrature of the Lorentz boost, followed by a
    # forward triangular sweep of the luminosity-work debit. The production
    # source uses analytic moments and two simultaneous fixed-point sweeps.
    gx,gw=np.polynomial.legendre.leggauss(16);area=(m.r/m.RJ)**2
    mu0=np.sqrt(np.maximum(0,1-(m.RJ/m.r)**2*(m.a/m.a0)**2));mu=mu0[:,None]+(1-mu0[:,None])*(gx+1)/2
    weights=gw[None,:]*(1-mu0[:,None])/2
    baseF=m.F0/(area*(m.a/m.a0)**2);W=1/np.sqrt(1-v*v);debit=0.;work=[]
    for i in range(len(rho)):
        F=baseF[i]+debit/(m.a[i]**2*area[i]*C)
        intensity=2*F/(1-mu0[i]**2)
        transformed=np.sum(weights[i]*intensity*W[i]**2*(mu[i]-v[i])*(1-v[i]*mu[i]))
        force=m.eos.rho0*rho[i]*kap[i]*transformed
        power=-m.vol[i]*m.a[i]**2*C*W[i]*v[i]*force
        work.append(power);debit+=power
    _,ports,_=m.rhs(U,float(history[-1,0]));work_error=abs(-debit/ports[2]-1)
    # Export cumulative, already accepted RK ledgers, rather than differentiating
    # sparse snapshots or silently adding the shell to the old bulk gas.
    assert np.max(abs(d['snapshots'][:,:,-1]))==0,'Gas reached the outer boundary'
    times=np.r_[0.,history[:,0]];baryon_port=np.r_[0.,history[:,3]]*scale
    energy_port=np.r_[0.,history[:,4]-history[:,5]]*scale*C*C
    radiation_work=np.r_[0.,history[:,5]]*scale*C*C
    np.savez_compressed(OUT/'bulk-radiation-ports.npz',time_seconds=times,
        baryon_into_layer_g=baryon_port,baryon_into_bulk_g=-baryon_port,
        Killing_nonrest_energy_into_layer_erg=energy_port,Killing_nonrest_energy_into_bulk_erg=-energy_port,
        scattering_work_into_gas_erg=radiation_work,scattering_work_into_photons_erg=-radiation_work,
        reference_redshifted_rest_energy_per_g=m.a0*m.eos.cx*C*C,
        cumulative_unresolved_tail_baryon_g=float(tail[0]*scale))
    e0=m.energy(initial);energy_residual=float(np.sum((m.energy(U)-e0)*m.vol)-history[-1,4])*scale*C*C
    source_error_scale=abs(energy_residual/row['integrated_trace_energy_erg'])
    passed=max(init_errors)<.002 and baryon<1e-10 and work_error<1e-8
    result=dict(classification='Counterexample candidate',passed=bool(passed),initial_EOS_relative=init_errors,
        baryon_ledger_relative=float(baryon),independent_angular_and_serial_scattering_work_relative=float(work_error),
        instantaneous_gas_scattering_power_erg_s=float(-debit*scale*C*C),
        integrated_gas_scattering_work_erg=float(radiation_work[-1]),
        inner_baryon_into_layer_g=float(baryon_port[-1]),inner_Killing_nonrest_energy_into_layer_erg=float(energy_port[-1]),
        discarded_dilute_baryon_g=float(tail[0]*scale),discarded_material_is_a_ledger_not_physical_vacuum=True,
        actual_energy_balance_residual_erg=energy_residual,energy_residual_over_trace_change=source_error_scale,
        source_accuracy_scope='The accepted energy gate is normalized by kinetic-plus-pressure response. The residual relative to the smaller trace change is reported separately. No defect is added as physical heat, and no pairwise tolerance is called a rigorous charge enclosure.',
        inner_port_applied_back_to_bulk=False,full_retarded_scattering=False,absorption_chemical_kinetics=False,
        full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-start)
    task.task.write(OUT/'flow-audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':main()
