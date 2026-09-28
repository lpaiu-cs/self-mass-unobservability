"""Correct only the signed scalar boundary term; reuse every photon ray."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
import def_native_dynamic_lapse as run

OUT=run.OUT/'corrected';write=run.write;sha=run.sha


def repair():
    assert not OUT.exists();OUT.mkdir()
    # Independent weak-scalar vacuum test: f=exp[-(r-1)]/r, Phi=K/r^2.
    # Direct integration of -delta_nu_prime fixes the boundary sign.
    K=.01
    direct=quad(lambda r:K*np.exp(1-r)/r**2,1,np.inf,epsabs=1e-14)[0]
    tail=quad(lambda r:2*K*np.exp(1-r)/r**3,1,np.inf,epsabs=1e-14)[0]
    correct=K-tail;wrong=-K-tail
    assert abs(correct-direct)<1e-13 and abs(wrong-direct)>.01
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        correction='The scalar lapse boundary is +Phi0*U0. The derivative identity was correct but the minus sign of the integral from r0 to infinity was applied incorrectly. Preserve the old producer, all raw results and symbolic text; do not relabel their numerical gate as validation of this sign.',
        action='Add2*Phi0*U0 uniformly to lapse values and the same correction to lapse time derivatives. Radial lapse derivatives, all actual packet rays, mass constraints and scalar fields are unchanged and reused.',
        original_symbolic_angular_failure='Earlier Hamiltonian verification also caught the opposite sign of the time-dependent angular drift; that failure occurred before ray runs and was already repaired.',
        budget=dict(seconds=10,new_ray_integrations=0,new_fluid_steps=0,new_native_calls=0),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),run.OUT/'lapse-producer-before-boundary-sign.py',run.OUT/'result.json']}))
    model=run.Lapse();Phi0=model.bg.fields(np.array([model.r0]))['Phi'][0];rows=[]
    for steps,order in [(128,8),(128,4),(64,8)]:
        path=run.OUT/f'metric-{steps}-g{order}.npz';d=dict(np.load(path));U=np.load(run.prior.new.OUT/f'fields-{steps}-g8.npz')['U'][:,-1]
        correction=2*Phi0*U
        for key in ['delta_nu','delta_nu_faces','delta_log_lapse','delta_log_speed']:d[key]+=correction[:,None]
        d['delta_nu_interval_rate']+=np.diff(correction)[:,None]/np.diff(d['t'])[:,None]
        d['outer_scalar_lapse']+=correction
        np.savez_compressed(OUT/path.name,**d)
        rows.append(dict(steps=steps,order=order,source_sha256=sha(path),maximum_boundary_sign_correction=float(np.max(abs(correction))),
            endpoint_outer_scalar_lapse=float(d['outer_scalar_lapse'][-1]),maximum_delta_nu=float(np.max(abs(d['delta_nu'])))))
    write(OUT/'audit.json',dict(classification='Counterexample candidate',passed=True,independent_boundary_sign_absolute_error=abs(correct-direct),
        old_wrong_sign_absolute_error=abs(wrong-direct),paths=rows,all_ray_results_reused=True,full_spatial_feedback=False,final_charge_solved=False))
    write(OUT/'symbolic.json',run.symbolic())
    print(json.dumps(rows),flush=True)


if __name__=='__main__':repair()
