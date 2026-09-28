"""Exact obstruction to closing dynamic transport with static mean data.

No stellar evolution is run. The examples establish input nonuniqueness,
not physical opacity models for this star or a no-go for dynamic charge.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import def_reactive_cauchy as coupled

h=coupled.h
OUT=coupled.OUT.parent/'def-thermal-closure-identifiability'


def symbolic():
    # Two strictly positive four-group spectra, same two independently
    # prescribed static moments. These equal group weights are explicit.
    a=[sp.Integer(1),sp.Integer(1),sp.Rational(1,2),sp.Integer(2)]
    b=[sp.Rational(3,2),sp.Rational(2,3),(7+sp.sqrt(13))/6,(7-sp.sqrt(13))/6]
    assert all(x.is_positive for x in a+b)
    mean=lambda q:sp.simplify(sum(q)/4)
    inverse=lambda q:sp.simplify(sum(1/x for x in q)/4)
    assert mean(a)==mean(b)==sp.Rational(9,8)
    assert inverse(a)==inverse(b)==sp.Rational(9,8)
    s=sp.symbols('s')
    # Laplace-domain angular dipole relaxation, with c*rho scaled to one.
    # Its DC conductivity is the same inverse-opacity mean for both sets.
    transfer=lambda q:sp.factor(sp.simplify(sum(1/(x+s) for x in q)/4))
    A,B=transfer(a),transfer(b);difference=sp.factor(A-B)
    assert difference!=0 and sp.simplify(difference.subs(s,0))==0
    assert sp.simplify(difference.subs(s,sp.I))!=0
    C,K,tau,k,z=sp.symbols('C K tau k z',positive=True)
    # Entropic Cattaneo current with fixed EOS heat capacity and DC K.
    generator=sp.Matrix([[0,-sp.I*k/C],[-sp.I*k*K/tau,-1/tau]])
    entropy_metric=sp.diag(C,tau/K)
    dissipation=sp.simplify(entropy_metric*generator+sp.conjugate(generator.T)*entropy_metric)
    assert dissipation==sp.diag(0,-2/K)
    polynomial=sp.factor((z*sp.eye(2)-generator).det())
    assert sp.simplify(polynomial-(z*z+z/tau+K*k*k/(C*tau)))==0
    # A mean-free-time is an independent constitutive assumption: positivity,
    # entropy production and subluminal characteristic speeds give no upper
    # bound on tau or unique transient response from C,K.
    return dict(classification='Proven',passed=True,
        static_double_mean=dict(first=list(map(str,a)),second=list(map(str,b)),
            Planck_equal_weight_mean=str(mean(a)),Rosseland_equal_weight_mean=str(1/inverse(a)),
            first_dynamic_transfer=str(A),second_dynamic_transfer=str(B),difference=str(difference),
            real_frequency_one_difference=str(sp.simplify(difference.subs(s,sp.I))),
            domain='Positive four-group relaxation model with explicitly equal normalized Planck and Rosseland weights. Not a fit to this star or an assertion of its physical group weights.'),
        conductive_family=dict(equations='C*T_t+q_x=0; tau*q_t+q=-K*T_x',
            quadratic_energy='(C*T^2+tau*q^2/K)/2',dissipation='-q^2/K, apart from boundary flux',
            mode_polynomial=str(polynomial),characteristic_speed='sqrt(K/(C*tau))',
            causal_requirement='tau>=K/(C*c^2); every larger tau obeys this requirement and has the same static EOS and conductivity.',
            driven_transfer='T/source=(1+s*tau)/(C*s*(1+s*tau)+K*k^2)',
            limitation='A mathematical constitutive family, not a statement that arbitrary tau is realized by this stellar electron plasma.'),
        consequence='Even exact equilibrium EOS, positive DC conductivity, Planck and Rosseland means, passivity and a speed ceiling do not in general select a unique dynamic heat closure. Actual frequency/angle collision information or an explicitly declared, independently justified constitutive model is needed.',
        not_proved=['No physical stellar dynamic-charge no-go','No numerical lower bound on the opacity uncertainty of this star','No failure of the accepted Phase52 source-frozen partial response'])


def main():
    assert not OUT.exists();OUT.mkdir()
    sources=[Path(__file__),Path(coupled.__file__),coupled.thermal.OUT/'coefficients.npz',
             h.ROOT/'verification/gr_two_carrier_evolution.py',
             h.ROOT/'outputs/direct-eos-gr33/gr-radiative-boundary/spectral-nonuniqueness.json']
    h.write(OUT/'plan.json',dict(classification='Proven',
        claim='Determine whether the saved static thermal inputs and generic admissibility conditions uniquely supply the missing dynamic radiation/conduction closure before extending stellar evolution.',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in sources},
        budget=dict(hard_seconds=30,new_native_calls=0,new_stellar_time_steps=0),
        decision='If nonunique, do not turn the existing guessed relaxation time or Rosseland mean into a physically certified full-star calculation. Specify exactly which missing inputs each evolution block needs.'))
    proof=symbolic();h.write(OUT/'symbolic.json',proof)
    d=np.load(coupled.thermal.OUT/'coefficients.npz');rho=d['raw'][:,0];capacity=rho*d['thermo'][:,3]/np.exp(d['lnT'])
    # This is a characteristic bound for the declared scalar Cattaneo block,
    # not a relativistic coupled-fluid causality certificate or a calibration.
    lower=d['K']/(capacity[:,None]*(h.gr.C*100)**2)
    assert np.min(lower)>0 and np.all(np.isfinite(lower))
    np.savez_compressed(OUT/'conditional-causal-floor.npz',heat_capacity_volume=capacity,
                        tau_lower_seconds=lower,radius_cm=d['radius_cm'])
    result=dict(classification='Counterexample candidate',same_saved_background_cells=len(rho),
        conditional_photon_tau_floor_range_seconds=[float(lower[:,0].min()),float(lower[:,0].max())],
        conditional_conductive_tau_floor_range_seconds=[float(lower[:,1].min()),float(lower[:,1].max())],
        physical_relaxation_times_identified=False,full_radiation_hydrodynamics_closed=False,
        physical_dynamic_charge_identified=False,full_dynamic_charge_solved=False,
        decision='No further full physical photon/conduction evolution is warranted by these static coefficients alone. An uncalibrated admissible closure can produce a numerical path but cannot discharge the physical completion requirement.')
    h.write(OUT/'result.json',result);print(json.dumps(dict(symbolic_passed=True,**result)),flush=True)


if __name__=='__main__':main()
