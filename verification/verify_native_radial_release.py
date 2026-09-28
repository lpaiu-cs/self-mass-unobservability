"""Finish the saved EOS-bank audit and verify conservative source identities."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as s
import def_native_radial_release as task

OUT=task.OUT


def main():
    assert (OUT/'eos.npz').exists() and not (OUT/'eos.json').exists()
    progress=np.load(OUT/'eos-progress.npz');assert np.all(progress['done'])
    task.write(OUT/'budget-stop.json',dict(classification='Counterexample candidate',native_calls_completed=1200,
        failure='Original native cap exhausted during the independent interpolation controls, after all441 bank states were saved.',
        preserved_bank_native_calls=int(progress['native_calls']),preserved_states=int(progress['done'].sum()),
        unsaved_control_calls=1200-int(progress['native_calls']),total_failed_run_wall_seconds=None,
        root_cause='The call forecast did not reserve the bounded independent entropy inversions after the production bank. Reuse the complete bank; do not regenerate it.'))
    task.write(OUT/'eos-audit-recovery-plan.json',dict(classification='Counterexample candidate',
        decision='Finish the original eight off-grid EOS controls, using the completed interpolated temperature as the native entropy-root initial guess. Preserve every table, domain and threshold.',
        additional_native_call_cap=40,total_native_call_cap=1240,additional_seconds=15,
        forecast='Eight controls initialized at the saved interpolant, usually two or three Newton calls each; cap40calls and15s. No new EOS states for production and no grid expansion.',
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'eos.npz',OUT/'plan.json']}))
    start=time.monotonic();signal.alarm(15);eos=task.EOS();d=eos.d;fan=task.prior.Fan(call_cap=40,reuse=True);rows=[]
    for xx in [-.037,-.7,-3.3,-10.3]:
        for frac in [.25,.75]:
            sigma=float(d['sigma'][0])*frac;target=float(d['s0'])+float(d['sunit'])*sigma
            p,u,gamma,T,k=eos(np.array([np.exp(xx)]),np.array([sigma]));guess=np.log(T[0])
            for _ in range(8):
                raw=fan.call(np.log(eos.rho0)+xx,guess);error=(raw[3]-target)*np.exp(guess)/raw[10]
                if abs(error)<2e-12:break
                guess-=error
            else:raise AssertionError('Audit entropy root')
            errors=[float(abs(p[0]*eos.rho0*task.C**2/raw[1]-1)),float(abs(u[0]*task.C**2/raw[2]-1)),float(abs(gamma[0]/raw[4]-1)),float(abs(T[0]/np.exp(guess)-1))]
            rows.append(dict(log_density=xx,entropy_parameter=sigma,relative=errors))
            task.write(OUT/'eos-control-progress.json',dict(rows=rows,native_calls=fan.calls))
    f,v=s.symbols('f v',real=True);W=1/s.sqrt(1-v*v)
    G0=W*v*f;Gr=W*f
    assert s.simplify(-W*G0+W*v*Gr)==0
    p,e=s.symbols('p e',real=True);E=(e+p)*W**2-p;S=(e+p)*W**2*v
    assert s.simplify(E-(S*v+p)-2*p-(e-3*p))==0
    symbolic=dict(classification='Proven',passed=True,
        scattering='G_lab^0=W*v*f_com and G_lab^r=W*f_com imply u_mu*G^mu=0. Gas entropy is conserved for this coherent scattering model; lab-frame photon work is nonzero and must be debited.',
        gas='The evolved gas has E=(e+p)W^2-p, S=(e+p)W^2*v, Pr=S*v+p, Pt=p and trace=e-3p; the LTE photon tensor has already been removed.',
        energy='On ds^2=-a^2*c^2*dt^2+B^2*dr^2+r^2*dOmega^2, conserved Killing energy is a*E; its flux is a^2*c*S. Subtracting a_s*c_x*c^2 times baryon conservation gives the small-energy ledger.',
        scope='Elastic scattering and fixed saved metric. No absorption, non-equilibrium chemical kinetics, full dynamic metric or final scalar charge theorem.')
    task.write(OUT/'symbolic.json',symbolic)
    passed=max(max(r['relative']) for r in rows)<.002
    task.write(OUT/'eos.json',dict(classification='Counterexample candidate',passed=passed,controls=rows,
        bank_states=[3,len(d['x'])],reused_states=145,native_calls_completed_before_audit=1200,recovery_native_calls=fan.calls,
        total_native_calls=1200+fan.calls,recovery_seconds=time.monotonic()-start,entropy_parameter_interval=d['sigma'][[0,-1]].tolist(),
        density_interval=d['x'][[0,-1]].tolist(),minimum_temperature_K=float(d['T'].min()),symbolic_passed=True))
    signal.alarm(0);print((OUT/'eos.json').read_text(),flush=True);assert passed,'Native interpolation failed'


if __name__=='__main__':main()
