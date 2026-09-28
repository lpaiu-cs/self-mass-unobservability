"""Conditional boundary for promoting Fourier heat flow to a dynamical GR fluid.

This reproduces a restricted Eckart instability mechanism, not a new physical
relaxation calibration or a statement about every first-order fluid theory.
"""
import json,sys,urllib.request
import numpy as np
import sympy as sp
import gr_microphysics as micro

g=micro.g;OUT=g.OUT/'gr-heat-relaxation-boundary'
URL='https://ccom.ucsd.edu/~lindblom/Publications/24_PhysRevD.31.725.pdf'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir()
    with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'hiscock-lindblom-1985.pdf').write_bytes(response.read())
    assert (OUT/'hiscock-lindblom-1985.pdf').read_bytes().startswith(b'%PDF')
    paths=[g.ROOT/'verification/gr_heat_relaxation_boundary.py',g.OUT/'gr-heat-fluid-frame/manifest.json',
        g.OUT/'initial-state-17-4.npz',g.OUT/'gr-microphysics/auxiliaries.npz',
        g.OUT/'gr-opacity/evaluation.npz',OUT/'hiscock-lindblom-1985.pdf']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='a4fb34c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source=dict(classification='Imported from prior work',url=URL,
            title='Generic instabilities in first-order dissipative relativistic fluid theories',
            authors=['William A. Hiscock','Lee Lindblom'],year=1985,doi='10.1103/PhysRevD.31.725',
            use='The primary article identifies instabilities of dissipative relativistic fluid theories. The restricted homogeneous equations below are independently derived and checked; no new-discovery claim is made.'),
        scope='Linear homogeneous local inertial perturbations in the baryon frame of an isotropic zero-flux equilibrium, constant positive enthalpy/temperature/conductivity; no viscosity or particle diffusion. The necessary homogeneous condition is not a sufficient finite-wave-number, nonlinear or causal stability theorem.',
        physical_relaxation_calibrated=False,full_GR_evolution=False))
    w,K,T,c,tau,s=sp.symbols('w K T c tau s',positive=True)
    # For physical velocity V and heat flux f, at zero wave number:
    # w*V_dot + f_dot=0, tau*f_dot+f=-(K*T/c^2)*V_dot.
    matrix=sp.Matrix([[w*s,s],[K*T*s/c**2,1+tau*s]])
    determinant=sp.factor(matrix.det());tmin=K*T/(w*c*c)
    assert sp.simplify(determinant-w*s*(1+(tau-tmin)*s))==0
    growth=sp.simplify(-1/(tau-tmin))
    assert sp.simplify(growth.subs(tau,0)-1/tmin)==0
    assert sp.simplify(growth.subs(tau,2*tmin)+1/tmin)==0
    assert sp.simplify(growth.subs(tau,3*tmin)+1/(2*tmin))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        homogeneous_equations='w*dV/dt+df/dt=0; tau*df/dt+f=-(K*T/c^2)*dV/dt. Here w=epsilon+P, K is thermal conductivity, and f is physical heat flux. The acceleration term is the local baryon-frame Fourier/Eckart term.',
        nonzero_pole='s=-1/(tau-tau_min), tau_min=K*T/(w*c^2). Tau=0 produces a positive pole 1/tau_min. Tau=tau_min makes the evolution matrix singular; positive decay of this homogeneous mode requires tau>tau_min.',
        nonuniqueness='Tau=2*tau_min and tau=3*tau_min have the same equilibrium EOS and stationary conductivity but different stable homogeneous poles. EOS plus stationary opacity/conductivity do not determine the missing transient coefficient.',
        limitations='Only this specified homogeneous subsystem is analyzed. The inequality does not prove full spatial causality/stability, select a physical tau, or identify a stellar/macroscopic relaxation time with the local heat-flux time. Other constitutive frames/theories are not excluded.'))
    state,aux=micro.inputs();a=aux['eos'];rho=np.exp(state['lnd']);temperature=np.exp(state['lnT'])
    cvT=a[:,10];c0=g.c.gr.C*100;CX=(state['X']/g.c.A)@g.c.W
    enthalpy=rho*(CX*c0*c0+a[:,2])+a[:,1]
    opacity=dict(np.load(g.OUT/'gr-opacity/evaluation.npz'))['values'][:,0]
    conductivity=16*5.670400e-5*temperature**3/(3*rho*opacity)
    minimum=conductivity*temperature/(enthalpy*c0*c0)
    assert np.all(minimum>0) and np.all(np.isfinite(minimum))
    # A separate fixed-fluid Cattaneo heat equation has speed^2=K/(C*tau).
    # Its elementary subluminal condition is a different, also conditional,
    # test. It cannot replace the coupled-fluid characteristic analysis.
    heat_capacity_volume=rho*cvT/temperature
    fixed_fluid_bound=conductivity/(heat_capacity_volume*c0*c0)
    result=dict(classification='Counterexample candidate',completed=True,cells=len(rho),
        effective_conductivity_cgs_range=[float(conductivity.min()),float(conductivity.max())],
        homogeneous_necessary_tau_min_seconds_range=[float(minimum.min()),float(minimum.max())],
        fixed_fluid_Cattaneo_subluminal_tau_bound_seconds_range=[float(fixed_fluid_bound.min()),float(fixed_fluid_bound.max())],
        these_are_actual_tau_measurements=False,physical_relaxation_calibrated=False,
        nonlinear_or_finite_wavenumber_stability_certified=False,full_GR_evolution=False)
    np.savez_compressed(OUT/'thresholds.npz',effective_conductivity=conductivity,
        enthalpy_density=enthalpy,homogeneous_tau_min_seconds=minimum,
        fixed_fluid_Cattaneo_tau_bound_seconds=fixed_fluid_bound)
    save('result.json',result)
    save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('HEAT RELAXATION BOUNDARY',result,flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS HOMOGENEOUS HEAT RELAXATION BOUNDARY SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
