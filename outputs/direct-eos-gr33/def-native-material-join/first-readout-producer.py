"""Counterexample candidate: reuse actual-source retarded charge and native audit."""
from pathlib import Path
import inspect
import json
import signal
import sys
import textwrap
import time
import numpy as np
import def_native_material_join as run
import verify_native_interior_feedback as previous

OUT=run.OUT;write=run.write;sha=run.sha
replace=run.replace

source=textwrap.dedent(inspect.getsource(previous.source))
source=replace(source,"'h','j','t']", "'h','j','Pi','t']")
source=replace(source,"model.h=data['h'][k];", "model.Pi=data['Pi'][k];model.h=data['h'][k];")
audit=textwrap.dedent(inspect.getsource(previous.audit))
audit=replace(audit,"m.h=z['h'];", "m.Pi=z['Pi'];m.h=z['h'];")
audit=audit.replace('run.chem','run.previous.chem')
namespace=dict(vars(previous),run=run,OUT=OUT)
exec(compile(source+'\n'+textwrap.dedent(inspect.getsource(previous.readout))+'\n'+audit,__file__,'exec'),namespace)


def prepare():
    assert not (OUT/'readout-plan.json').exists()
    plan=json.loads((run.previous.OUT/'readout-plan.json').read_text())
    plan.update(claim='Use actual shared-material64/128 histories to judge the conditional direct charge after releasing the copied ghost. Cell momenta come from this run, not the old staggered face interpolation.',
        limits='Acoustic interior mechanics and first-order density/inventory; background-reconstructed native face at34km center-to-face gap. Fixed metric and unresolved spatial/frequency/interior-angle errors. Not final physical charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(run.previous.__file__),Path(previous.__file__)]})
    write(OUT/'readout-plan.json',plan)
    (OUT/'readout-producer.py').write_text(Path(__file__).read_text())
    (OUT/'expanded-source.py').write_text(source)


def readout():namespace['readout']()


def check():
    # Native endpoint and instantaneous number checks retain their old limits.
    namespace['audit']()
    start=time.monotonic();signal.signal(signal.SIGALRM,run.old.optical.timeout);signal.alarm(15)
    m=run.Coupled();z=np.load(OUT/'coupled-128.npz');m.h=z['h'];m.Pi=z['Pi'];m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(run.old.END)
    f=m.flow;b=m.bulk;u=z['u'];theta=z['theta'];eta=z['eta'];state=m.material_state(u,eta)
    m.recover_material(state,run.old.END,theta);V=f.primitive(z['U']);L,R=f.reconstruct(V,run.old.END)
    flux=m.join_flux(R[:,0],V[3,0]);physical=m.mflux.copy();deep=m.material_rhs(u,theta,eta)
    # Independently convert the atmosphere's finite-volume boundary units.
    scale=4*np.pi*m.m.RJ**2*f.eos.rho0
    converted=flux*np.array([run.C,run.C**2,run.C**3,run.C*f.eos.nH])*scale
    unit=float(np.max(abs(converted-physical)/np.maximum(abs(physical),1.)))
    # Telescoping mass/species/energy flux sums independently read the actual RHS.
    mass=float(abs(-np.diff(deep[0]).sum()+physical[0])/max(abs(physical[0]),1.))
    energy=float(abs(deep[2].sum()+physical[2])/max(abs(physical[2]),1.))
    neutral=float(abs(deep[3].sum()+physical[3])/max(abs(physical[3]),1.))
    passed=max(unit,mass,energy,neutral)<1e-8
    row=dict(classification='Counterexample candidate',passed=bool(passed),shared_face_conversion_relative=unit,instantaneous_global_mass_relative=mass,
        instantaneous_global_energy_relative=energy,instantaneous_global_neutral_relative=neutral,
        integrated_shared_mass_g=float(z['scalar_join_mass']),integrated_shared_momentum_g_cm_s=float(z['scalar_join_momentum']),
        integrated_shared_Killing_energy_erg=float(z['scalar_join_energy']),integrated_shared_neutral=float(z['scalar_join_neutral']),
        face_velocity_over_c=float(f.join_state[1]),face_temperature_K=float(np.exp(f.join_state[2])),
        face_density_cgs=float(f.join_state[0]*f.eos.rho0),face_neutral_fraction=float(f.join_state[3]),
        actual_shared_mass_momentum_energy_species_evolved=True,spatial_interface_reconstruction_certified=False,
        full_neutral_trajectory_ledger=False,full_GR_feedback=False,final_charge_solved=False,seconds=time.monotonic()-start)
    write(OUT/'join-audit.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert passed


if __name__=='__main__':globals()[sys.argv[1]]()
