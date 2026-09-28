"""Apply the completed energy-conserving gas histories to the prior readout."""
from pathlib import Path
from types import SimpleNamespace,FunctionType
import json
import signal
import time
import numpy as np
import def_native_energy_flow as evolution
import def_native_release_charge as prior

OUT=evolution.OUT/'charge'


def adapter(n):
    m=evolution.Flow(n)
    # Freeze the final native support. No native call or table extension is
    # permitted while reading the already completed histories.
    m.eos=evolution.temperature.parent.EOS(evolution.OUT/'runtime-columns.npz')
    m.eos.fan=SimpleNamespace(calls=0)
    base=m.base;base.eos=m.eos;base.primitive=m.primitive
    return base


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic();signal.alarm(30)
    paths=[Path(__file__),Path(evolution.__file__),Path(prior.__file__),evolution.OUT/'runtime-columns.npz']
    paths += [evolution.OUT/f'cells-{n}.npz' for n in [896,1792]]
    prior.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the completed conservative-energy histories to exactly the same direct Green and inner acoustic readout. Compare both registered gas grids and the preserved prior result.',
        seconds=30,native_calls=0,new_fluid_steps=0,grid_gate=.02,
        unchanged_scope='Direct prescribed-source scalar response. Full metric/matter/chemical/radiative closure is still required.',
        bindings={str(p.relative_to(prior.old.ROOT)):prior.old.photons.digest(p) for p in paths}))
    env=dict(vars(prior));env['OUT']=OUT;env['prior']=SimpleNamespace(Flow=adapter,OUT=evolution.OUT)
    readout=FunctionType(prior.readout.__code__,env,argdefs=prior.readout.__defaults__)
    old=np.load(prior.OUT/'cells-1792-linear-g12.npz');times=old['u_seconds'];raw=np.load(prior.OUT/'native-bulk.npz')['raw']
    fine,fr=readout(1792,'linear',times,raw);coarse,cr=readout(896,'linear',times,raw)
    error=float(max(abs(fine-coarse))/max(abs(fine)))
    change=float(max(abs(fine-old['normalized_charge']))/max(abs(fine)))
    result=dict(classification='Counterexample candidate',passed=error<.02,grid_wave_relative=error,change_from_entropy_evolution_relative=change,
        endpoint_direct_normalized_charge=fr['endpoint_direct_normalized_charge'],endpoint_inner_acoustic_energy_mismatch_erg=fr['endpoint_bulk_energy_mismatch_erg'],
        seconds=time.monotonic()-start,native_calls=0,new_fluid_steps=0,energy_conserving_input_applied=True,
        final_charge_solved=False,full_GR_scalar_feedback=False,full_goal_complete=False)
    prior.write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
