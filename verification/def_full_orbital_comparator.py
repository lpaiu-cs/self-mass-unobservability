"""Remove known interior wave propagation from a static-response comparison.

The comparator retains the measured static stellar impedance plus the exact
test-scalar propagation on the same frozen interior metric. This separates a
frequency-independent material response from geometrical wave storage.
"""
from pathlib import Path
import json
import time
import signal
import numpy as np
from scipy.integrate import solve_ivp
import def_full_orbital_exterior as full
import def_orbital_charge_audit as audit

OUT=full.OUT;prior=full.prior;write=full.write


def propagation(omega,saved,tol,step,flat=False):
    rs=float(saved['r'][-1]);radius=saved['r']/rs;mass=saved['m']/rs
    lapse=saved['N'];fc=float(lapse[0])
    def F(r):
        if flat:return 1.
        m=mass[1]*(r/radius[1])**3 if r<radius[1] else np.interp(r,radius,mass)
        return np.interp(r,radius,lapse)*np.sqrt(1-2*m/r) if r else fc
    # psi=1+omega^2*u, r^2*F*psi'=omega^2*w. Scaling avoids
    # subtracting a nearly constant incident field before differentiating.
    def rhs(r,y):
        f=F(r)
        return [y[1]/(r*r*f) if r else 0.,-r*r/f*(1+omega*omega*y[0])]
    sol=solve_ivp(rhs,[0,1],[0.,0.],method='DOP853',rtol=tol,atol=tol*1e-3,max_step=step)
    assert sol.success
    u,w=sol.y[:,-1]
    return omega*omega*w/(F(1)*(1+omega*omega*u)),sol.nfev


def main():
    assert not (OUT/'comparator-plan.json').exists()
    write(OUT/'comparator-plan.json',dict(classification='Counterexample candidate',
        reassessment='The full frozen-Robin contrast is dominated by a term of order omega^2/3, also present in a flat empty ball. It must not be labeled a microscopic relaxation signal.',
        comparator='At physical stellar surface use Z_static_star+Z_test_wave(omega). Test wave solves (r^2 F psi_prime)_prime+omega^2*r^2/F*psi=0 on the same frozen metric with regular centre. Z_test_wave(0)=0. Exterior scalar+metric propagation remains identical. This is one explicitly specified static-material comparator, not all static EFTs.',
        controls='Flat exact omega*cot(omega)-1; original p2/p4/outer/quadrature data and two geometric IVP tolerances. Test material susceptibility stays frequency independent; metric is fixed.',
        gates=dict(flat_relative=1e-9,geometry_relative_on_retained_signal=.002,spatial=.02,outer=.002,quadrature=.002),
        budget=dict(seconds=120,new_EOS_calls=0,new_FEM_solves=0,new_time_steps=0),
        decision='Stop on a failed contrast gate; no extra harmonics or longer integration. Keep the frozen-Robin result as the diagnostic that motivated this separate declared comparator.',
        bindings={str(p):prior.go.task.digest(p) for p in [Path(__file__),OUT/'result.json',OUT/'plan.json',prior.BENCH,prior.patch.base.surface.OUT/'background.npz']}))
    signal.alarm(120);start=time.monotonic();audit.mp.mp.dps=70
    original=json.loads((OUT/'result.json').read_text());saved=np.load(prior.patch.base.surface.OUT/'background.npz')
    mu,flux,M,K,bg=full.background();F=float(bg.y[1,-1]);pilot=json.loads((prior.OUT/'pilot.json').read_text());ratio=pilot['R_m']/pilot['ADM_geom_m']
    base=json.loads((prior.OUT/'result.json').read_text());omegas=[r['omega_R_over_c'] for r in base['cases']['p4']]
    geometric={};flat=[]
    for setting,tol,step in [('standard',2e-11,1/128),('tight',2e-12,1/256)]:
        geometric[setting]=[0.]
        for w in omegas[1:]:geometric[setting].append(propagation(w,saved,tol,step)[0])
    for w in omegas[1:]:
        exact=float(audit.mp.mpf(str(w))/audit.mp.tan(audit.mp.mpf(str(w)))-1)
        value,_=propagation(w,saved,2e-12,1/256,True);flat.append(abs(value/exact-1))
    assert max(flat)<1e-9
    cases={};amplitudes=[audit.mp.mpf(str(r['drive_amplitude'])) for r in base['rows']]
    for label,oldrows in original['rows'].items():
        for setting in ['standard','tight']:
            rows=[];z0=oldrows[0]['surface_impedance']
            for n,row in enumerate(oldrows):
                w=omegas[n];wave=full.ext.outgoing(mu,flux,w,tol=2e-13);h=wave['h'];Zo=wave['impedance']
                z=row['surface_impedance'];Zc=z0+geometric[setting][n];dz=z-Zc;D=Zc-Zo
                drive=np.exp(-1j*w)/(F*h);tail=-drive*dz/(D*(D+dz)*h);gain=-ratio*tail
                rows.append(dict(harmonic=n,charge=[float(gain.real),float(gain.imag)],
                    residual_surface_impedance=dz,test_wave_impedance=geometric[setting][n],
                    delta_alpha_radiative_over_phi0=float(abs(gain)*float(amplitudes[n])/.001)))
            cases[label+'-'+setting]=rows
    compare=[];ref=cases['p4-tight']
    for n in [1,2,3]:
        a=complex(*ref[n]['charge']);item=dict(harmonic=n)
        for label in ['p2-tight','outer-tight','quadrature-tight','p4-standard']:
            item[label]=abs(complex(*cases[label][n]['charge'])-a)/abs(a)
        compare.append(item)
    passed=all(r['p2-tight']<.02 and r['outer-tight']<.002 and r['quadrature-tight']<.002 and r['p4-standard']<.002 for r in compare)
    fits={label:[audit.projection(rows,d,amplitudes,audit.mp.mpf('.001')) for d in [0,1,2,3,4,5]] for label,rows in cases.items()}
    error4=max(float(audit.mp.norm(audit.mp.matrix(rows[4]['residual'])-audit.mp.matrix(fits['p4-tight'][4]['residual']))) for label,rows in fits.items() if label!='p4-tight')
    result=dict(classification='Counterexample candidate',passed=passed,cases=cases,comparisons=compare,
        flat_exact_relative=max(flat),fits=fits,degree4_control_difference=error4,
        sum_first_three_amplitudes=sum(r['delta_alpha_radiative_over_phi0'] for r in ref[1:]),
        interpretation='Full adiabatic stellar response minus declared static-material plus same-metric propagation comparator. Includes matter/metric/scalar response, not only free-minus-clamped motion. Conditional radiation readout, not thermally matched binary sensitivity or observed signal.',
        full_objective_complete=False,seconds=time.monotonic()-start)
    write(OUT/'comparator-result.json',result)
    print('COMPARATOR',json.dumps({k:v for k,v in result.items() if k not in ['cases','fits']}),flush=True)
    print('FINE',ref,'FITS',[(r['degree'],r['residual_norm']) for r in fits['p4-tight']],flush=True)


if __name__=='__main__':main()
