"""Export collision occupation and additional, not double-counted, GR sources."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_compensated_transport as run

OUT=run.OUT;write=run.write;sha=run.sha;flow=run.flow


def main():
    assert not (OUT/'collision-input-audit.json').exists()
    write(OUT/'collision-input-plan.json',dict(classification='Counterexample candidate',
        claim='Convert the conserved packet increment to the occupation seen by actual collisions, and remove the canonical metric-volume photon moments already included in Phase122 before exporting additional GR forcing.',
        identity='delta_F=delta_N/reference_phase_weight-(3*alpha*delta_phi+delta_lambda)*F0. Fixed-canonical integrated energy changes by -u*E0-lambda*Pr0; integrated radial pressure changes by -u*Pr0+lambda*(R4_0-2*Pr0).',
        limitations='This prepares the actual collision and next GR inputs; it is not their evolution. Additional fluid density/velocity/temperature/H responses and their collision derivatives remain.',
        budget=dict(seconds=15,new_native_calls=0,new_fluid_steps=0,new_transport_steps=0),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),OUT/'steps-128-reference-128.npz',run.metric.OUT/'corrected/metric-128-g8.npz',flow.OUT/'coupled-128.npz']}))
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(15);start=time.monotonic();m=run.Transport(128)
    d=np.load(OUT/'steps-128-reference-128.npz');g=m.g;dv=3*g['delta_u']+g['delta_lambda']
    x=d['delta_count_scaled_occupation'];occupation=x-dv[-1,:,None,None]*m.I[-1]
    residual=float(np.max(abs(occupation+dv[-1,:,None,None]*m.I[-1]-x))/max(np.max(abs(x)),1e-300));assert residual<1e-12
    mu4=(m.edges_mu[1:]**5-m.edges_mu[:-1]**5)/(5*np.diff(m.edges_mu))
    E0=np.sum(m.I*m.weights*m.E,axis=(2,3));P0=np.sum(m.I*m.weights*m.E*m.model.bulk.mu2[None,None,:,None],axis=(2,3))
    R4=np.sum(m.I*m.weights*m.E*mu4[None,None,:,None],axis=(2,3))
    canonicalE=-g['delta_u']*E0-g['delta_lambda']*P0
    canonicalP=-g['delta_u']*P0+g['delta_lambda']*(R4-2*P0)
    a=np.r_[m.model.bulk.d['a'],m.model.m.a]
    extraE=(d['moments'][:,1]-canonicalE)/a;extraP=(d['moments'][:,2]-canonicalP)/a
    np.savez_compressed(OUT/'collision-input-128.npz',delta_occupation=occupation,delta_log_proper_volume=dv[-1],
        reference_t=g['t'][-1],reference_coupled_sha256=sha(flow.OUT/'coupled-128.npz'))
    np.savez_compressed(OUT/'additional-photon-sources-128.npz',t=g['t'],radius_E=m.r,additional_energy_erg=extraE,
        additional_radial_pressure_erg=extraP,canonical_reference_energy_erg=canonicalE,canonical_reference_pressure_erg=canonicalP,
        full_collision_or_hydrodynamic_feedback=False)
    row=dict(classification='Counterexample candidate',passed=True,occupation_reconstruction_relative=residual,
        maximum_additional_photon_energy_L1_erg=float(np.max(np.sum(abs(extraE),axis=1))),
        maximum_additional_radial_pressure_L1_erg=float(np.max(np.sum(abs(extraP),axis=1))),
        maximum_canonical_reference_energy_L1_erg=float(np.max(np.sum(abs(canonicalE),axis=1))),
        seconds=time.monotonic()-start,collision_input_exported=True,canonical_double_count_removed=True,
        collision_response_evolved=False,full_GR_feedback=False,final_charge_solved=False)
    write(OUT/'collision-input-audit.json',row);print(json.dumps(row),flush=True);signal.alarm(0)


if __name__=='__main__':main()
