"""Nonzero-heat characteristic audit and conditional initial GR tangent.

The Cattaneo constitutive model and its two uncalibrated heat times are
explicit. An initial derivative is not a finite GR evolution trajectory.
"""
import json,sys
from fractions import Fraction as F
import numpy as np
import sympy as sp
import gr_heat_characteristics as char

g=char.g;OUT=g.OUT/'gr-heat-initial-tangent'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();char.verify()
    paths=[g.ROOT/'verification/gr_heat_initial_tangent.py',char.OUT/'manifest.json',char.OUT/'coefficients.npz',
        g.OUT/'gr-heat-initial-constraints/initial-data.npz',g.OUT/'gr-heat-fluid-frame/manifest.json',
        g.OUT/'gr-heat-fluid-frame/rate-coefficients.npz',g.OUT/'gr-heat-relaxation-boundary/thresholds.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e6aa3f1',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        trial_multipliers=[2,3],root_bracket_halfwidth_powers=[45,40,35,30,25,20],
        method='Use floating eigenvalues only to propose disjoint rational sign brackets. Four opposite-sign intervals strictly inside (-1,1) certify all roots of the exact frozen quartic. No unbracketed approximate eigenvalue is accepted.',
        initial_constitutive_condition='The discrete initial Fourier target Q_F is set equal to the stored initial proper Q=heat_flux/c. This matches the existing face-to-cell heat current; it is an explicit collocation assumption, not an independent continuum gradient or atmosphere validation.',
        physical_heat_time_calibrated=False,full_GR_evolution=False))
    b,r,d,e,j,lam,v=sp.symbols('b r d e j lam v',real=True)
    M=sp.Matrix([[1,0,0,0],[0,b,2*j,0],[0,0,1,1],[0,0,1,lam]])
    N=sp.Matrix([[0,0,1,0],[0,0,e,1],[r,d,2*j,0],[0,1,0,0]])
    C=lam*(b*r+d*e)+1-e-d
    poly=b*(lam-1)*v**4+2*j*(lam*(d-b)-1)*v**3-C*v*v+2*j*v+r
    assert sp.simplify((N-v*M).det()-poly)==0
    assert sp.simplify(M.det()-b*(lam-1))==0
    tau,h,w,Q,QF,Lapse,R=sp.symbols('tau h w Q QF Lapse R',real=True)
    Qdot=(Lapse*(QF-Q)-h*R)/(tau-h)
    vdot=(tau*R-Lapse*(QF-Q))/(w*(tau-h))
    assert sp.simplify(Qdot+w*vdot-R)==0
    assert sp.simplify(tau*Qdot+w*h*vdot-Lapse*(QF-Q))==0
    # At initial material velocity zero the coordinate radial heat component
    # is Q/a. The metric connection cancels its basis derivative.
    a,adot,Qt,qd=sp.symbols('a adot Qt qd',real=True)
    assert sp.simplify(a*(qd/a-Qt*adot/a**2+adot/a*(Qt/a))-qd)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        constitutive_model='tau*Delta^mu_nu*u^alpha*nabla_alpha q^nu+q^mu=-K*(Delta^mu_nu*nabla^nu T+T*a^mu), with consistent c factors restored in the physical-second initial equations. No viscosity, baryon diffusion or extra Israel-Stewart expansion terms are included.',
        nonzero_heat_principal_part='With j=Q/w and initial material v=0, energy has the added 2*j*v_dot and momentum the added 2*c*j*v_x. The projected Cattaneo time/spatial principal coefficients are unchanged at this instant. Coefficients and background are frozen only for this local characteristic test.',
        polynomial='b*(lambda-1)*v^4+2*j*(lambda*(d-b)-1)*v^3-[lambda*(b*r+d*e)+1-e-d]*v^2+2*j*v+r; v is characteristic speed/c.',
        root_certificate='A nonzero real quartic with opposite endpoint signs in four disjoint intervals contained in (-1,1) has exactly four distinct real subluminal roots, by the intermediate value theorem and its degree. This certifies the frozen principal matrix, not lower-order finite-k stability on a non-equilibrium background.',
        covariant_initial_reduction='At v=0, the projected derivative of radial Q/a is Q_dot/(a*N), because the time derivative of the spatial basis cancels the connection. The material acceleration includes (ln N)_r/a + v_dot/(N*c). Thus tau*Q_dot+w*h*v_dot=N*(Q_F-Q), h=K*T/(w*c^2), where tau is a proper time and dots use the same coordinate seconds.',
        coupled_initial_rates='Let R=2*c*N*A*Q. Momentum gives Q_dot+w*v_dot=R. Hence Q_dot=[N*(Q_F-Q)-h*R]/(tau-h), v_dot=[tau*R-N*(Q_F-Q)]/[w*(tau-h)].',
        initial_target_specialization='For Q_F=Q, Q_dot=-h*R/(tau-h) and v_dot=tau*R/[w*(tau-h)]. Insert this v_dot into the earlier corrected material energy and temperature rates. It is generally nonzero despite initial v=0.',
        limitations='This is a specified initial tangent and pointwise cone audit. No finite time step, spatial truncation certificate, lapse update, nonlinear stability, physical transport calibration, reaction/convection/atmosphere or complete GR solution is claimed.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    fields=dict(np.load(char.OUT/'coefficients.npz'))['coefficients']
    initial=dict(np.load(g.OUT/'gr-heat-initial-constraints/initial-data.npz'))
    frame=dict(np.load(g.OUT/'gr-heat-fluid-frame/rate-coefficients.npz'))
    w=dict(np.load(g.OUT/'gr-heat-relaxation-boundary/thresholds.npz'))['enthalpy_density']
    c=g.c.gr.C*100;Q=initial['proper_heat_flux_cgs']/c;j=Q/w
    # This initial momentum forcing uses the same saved pointwise geometry.
    R=2*c*initial['lapse']*initial['radial_extrinsic_curvature_per_cm']*Q
    bracket_data=[];rate_data=[];maximum_residual=F(0);cone_failures=[]
    for index,values in enumerate(fields):
        B,P,D,E,H=[F(float(x)) for x in values];J=F(float(j[index]))
        sound=P+D*E/B
        lamstar=max(F(1),(B+1-E-D-P)/(B*(1-sound)),(2*B+1-E-D)/(B*(2-sound)))
        for multiplier in plan['trial_multipliers']:
            lam=multiplier*lamstar
            coeff=[B*(lam-1),2*J*(lam*(D-B)-1),-(lam*(B*P+D*E)+1-E-D),2*J,P]
            def poly(x):
                value=F(0)
                for a in coeff:value=value*x+a
                return value
            Mf=np.array([[1,0,0,0],[0,float(B),2*float(J),0],[0,0,1,1],[0,0,1,float(lam)]])
            Nf=np.array([[0,0,1,0],[0,0,float(E),1],[float(P),float(D),2*float(J),0],[0,1,0,0]])
            approximate=np.linalg.eigvals(np.linalg.solve(Mf,Nf))
            assert np.all(approximate.imag==0),(index,multiplier,approximate)
            centers=sorted(F(float(x)) for x in approximate.real);brackets=None
            for power in plan['root_bracket_halfwidth_powers']:
                half=F(1,2**power);trial=[(x-half,x+half) for x in centers]
                if (trial[0][0]>-1 and trial[-1][1]<1 and
                    all(trial[i][1]<trial[i+1][0] for i in range(3)) and
                    all(poly(lo)*poly(hi)<0 for lo,hi in trial)):
                    brackets=trial;break
            if brackets is None:cone_failures.append(dict(cell=index,multiplier=multiplier));continue
            # Save exact dyadic endpoints as strings; float display would not
            # preserve the rational interval membership certificate.
            bracket_data.append(dict(cell=index,multiplier=multiplier,
                bounds=[[str(lo),str(hi)] for lo,hi in brackets]))
            TAU=H*lam;RR=F(float(R[index]));WW=F(float(w[index]))
            qdot=-H*RR/(TAU-H);vdot=TAU*RR/(WW*(TAU-H))
            residual=max(abs(qdot+WW*vdot-RR),abs(TAU*qdot+WW*H*vdot))
            maximum_residual=max(maximum_residual,residual)
            correction=F(float(frame['coefficient_vdot_in_logT_dot'][index]))*vdot
            rate_data.append([index,multiplier,float(TAU),float(qdot),float(vdot),float(correction),
                float(F(float(frame['conditional_zero_vdot_logT_dot'][index]))+correction)])
    save('root-brackets.json',dict(classification='Proven',rows=bracket_data))
    rates=np.array(rate_data);np.savez_compressed(OUT/'initial-rates.npz',values=rates,Q=Q,j=j,momentum_forcing=R)
    result=dict(classification='Proven',completed=True,cells=len(fields),
        exact_frozen_cones_passed=not cone_failures,cone_failures=cone_failures,
        certified_cones=len(bracket_data),exact_initial_coupled_rate_residual=str(maximum_residual),
        maximum_abs_initial_vdot_per_coordinate_second=float(abs(rates[:,4]).max()),
        maximum_abs_added_logT_rate_per_coordinate_second=float(abs(rates[:,5]).max()),
        maximum_abs_Qdot_erg_cm3_per_coordinate_second=float(abs(rates[:,3]).max()),
        scope='Exact arithmetic on the declared frozen coefficient and forcing numbers. Rates and maxima saved as floating displays. Q_F=Q is an explicit discrete target assumption; no continuous gradient error is certified.',
        actual_time_trajectory_computed=False,physical_heat_time_calibrated=False,
        nonzero_heat_finite_wavenumber_stability_certified=False,physical_EOS_certified=False,full_GR_evolution=False)
    save('result.json',result)
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('NONZERO HEAT INITIAL TANGENT',result,flush=True)
    assert not cone_failures and maximum_residual==0,result
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    result=json.loads((OUT/'result.json').read_text())
    assert result['exact_frozen_cones_passed'] and result['exact_initial_coupled_rate_residual']=='0'
    print('PASS NONZERO HEAT INITIAL TANGENT SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
