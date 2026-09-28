"""Adjudicate Request 13 evidence without promoting convergence to a proof."""
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/research-remediation'


def main():
    # Central difference truncation coefficient for a cubic is exactly 1/6.
    x,h,a,b,c,d=sp.symbols('x h a b c d',real=True)
    f=a*x**3+b*x*x+c*x+d
    error=sp.expand((f.subs(x,x+h)-f.subs(x,x-h))/(2*h)-sp.diff(f,x))
    assert sp.simplify(error-h*h*sp.diff(f,x,3)/6)==0
    ti,t0,di,d0,D,spin,spindot,turn=sp.symbols('ti t0 di d0 D spin spindot turn',real=True)
    tref=sp.symbols('tref',real=True)
    T=ti*D-tref; T0=t0*D-tref
    original=spin/D*(T-di)+spindot/(2*D*D)*(T-di)**2-spin/D*(T0-d0)-spindot/(2*D*D)*(T0-d0)**2-turn
    stable=spin*(ti-t0)-turn+spindot/2*(ti-t0)*(ti+t0-2*tref/D)
    stable-=(spin/D+spindot*T/D**2)*di-(spin/D+spindot*T0/D**2)*d0
    stable+=spindot/(2*D**2)*(di**2-d0**2)
    assert sp.simplify(original-stable)==0
    base=np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    sw=1/base['errs']
    rows=[]
    for label in ['control','integration','mesh','inversion','stable']:
        for j in [21,22,27]:
            full=OUT/f'{label}-jac-{j:02d}-1.0.npz'
            half=OUT/f'{label}-jac-{j:02d}-0.5.npz'
            if not full.exists() or not half.exists(): continue
            u=np.load(full)['dcol']; v=np.load(half)['dcol']
            ref=np.load(OUT/f'control-jac-{j:02d}-0.5.npz')['dcol']
            rows.append(dict(label=label,column=j,
                h_vs_half_relative=float(np.linalg.norm(sw*(u-v))/np.linalg.norm(sw*v)),
                half_vs_control_relative=float(np.linalg.norm(sw*(v-ref))/np.linalg.norm(sw*ref))))
    assert json.loads((OUT/'control-preflight.json').read_text())['archived_max_difference_us']==0
    stellar=json.loads((OUT/'stellar-matching.json').read_text())
    native=[r for r in stellar['stars'] if r['interpolation']=='lal']
    final=native[-1]
    assert abs(final['gravitational_mass_solar']-1.4378144085)<1e-6
    independent=stellar['independent_LAL_refinement'][-1]
    assert abs(independent['mass_relative_difference'])<1e-6
    assert abs(independent['radius_relative_difference'])<1e-6
    assert abs(native[0]['static']['susceptibility_m']/final['static']['susceptibility_m']-1)<2e-6
    poles=np.array([[p['omega_real_per_s'],p['omega_imag_per_s']] for p in stellar['poles']])
    assert all(p['success'] and p['dimensionless_mismatch']<1e-7 for p in stellar['poles'])
    pole_change=float(np.max(np.linalg.norm(poles-poles[0],axis=1))/np.linalg.norm(poles[0]))
    assert pole_change<1e-6
    matches=[]
    for path in sorted(OUT.glob('nonlinear-*.json')):
        result=json.loads(path.read_text())
        history=result['history']; dev=np.array([r['deviance'] for r in history])
        assert np.all(np.diff(dev)<=1e-8)
        assert all(0<=r['eccentricity_extra']<1 for r in history)
        matches.append(dict(label=result['label'],iterations=len(history),deviance=float(dev[-1]),
            eccentricity_extra=history[-1]['eccentricity_extra'],stationarity_certified=result['stationarity_certified'],
            termination='three failed local line searches' if all(not r['accepted'] for r in history[-3:]) else 'iteration budget'))
    gradients=[]
    checkpoint=OUT/'nonlinear-zero.npz'
    for j in range(28):
        path=OUT/f'gradient-zero-{j:02d}.npz'
        z=np.load(path)
        assert str(z['checkpoint_sha256'])==hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        assert np.isfinite(z['dcol']).all()
        gradients.append(j)
    result=dict(status='Imported from prior work',derivative_diagnostics=rows,
        stellar_numerical_checks_pass=True,scalar_pole_contour_relative_change=pole_change,
        nonlinear_fits=matches,fresh_gradient_columns=gradients,
        rigorous_derivative_certificate=dict(pass_gate=False,function_error_upper_bound=None,
            third_derivative_upper_bound=None,reason='No continuum ODE/interpolation/roundoff enclosure is supplied. Measured differences cannot fill these fields.'),
        full_physical_inference_pass=False)
    (OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
