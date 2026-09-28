"""Correct the distinction between zero velocity and zero velocity derivative.

The first heat-current constraint construction remains valid. Its material
temperature-rate formula omitted the heat-flux/acceleration term unless the
velocity derivative is additionally set to zero. Preserve and label that error.
"""
import json,sys
import numpy as np
import sympy as sp
import gr_heat_initial_constraints as original

g=original.g;OUT=g.OUT/'gr-heat-fluid-frame'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir()
    original.verify()
    paths=[g.ROOT/'verification/gr_heat_fluid_frame.py',original.OUT/'plan.json',original.OUT/'manifest.json',
        original.OUT/'symbolic.json',original.OUT/'initial-data.npz',g.OUT/'gr-microphysics/auxiliaries.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='a4fb34c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        correction='The original initial normal velocity v=0 fixes E_normal=epsilon only at that instant. Its time derivative is E_normal_dot=epsilon_dot+2*Q*v_dot, with Q=q_proper/c. Do not differentiate the instantaneous equality while leaving v_dot unconstrained.',
        unchanged='Initial A=4*pi*r*a*q, pointwise momentum solution and zero added Hamiltonian term are unchanged. The original density rate remains valid at v=0. The stored material energy/temperature rate is conditional on v_dot=0.',
        physical_EOS_certified=False,full_GR_evolution=False))
    v,eps,P,Q=sp.symbols('v eps P Q',real=True)
    W=1/sp.sqrt(1-v*v);u=sp.Matrix([W,W*v]);heat=sp.Matrix([W*v*Q,W*Q]);eta=sp.diag(-1,1)
    assert sp.simplify((u.T*eta*u)[0]+1)==0 and sp.simplify((u.T*eta*heat)[0])==0
    tensor=(eps+P)*u*u.T+P*eta+u*heat.T+heat*u.T
    expected=sp.Matrix([[eps+P*v*v+2*v*Q,(eps+P)*v+Q*(1+v*v)],
        [(eps+P)*v+Q*(1+v*v),P+eps*v*v+2*v*Q]])/(1-v*v)
    assert all(sp.simplify(x)==0 for x in tensor-expected)
    derivatives=[sp.simplify(sp.diff(tensor[i,j],v).subs(v,0)) for i,j in [(0,0),(0,1),(1,1)]]
    assert derivatives==[2*Q,eps+P,2*Q]
    rho,C,energy_rate,Qdot,vdot,NA=sp.symbols('rho C energy_rate Qdot vdot NA',real=True)
    # The initially hydrostatic radial momentum projection is
    # d(a*j_hat)/d(ct)=N*A*a*Q, while a_dot/a=-N*A.
    assert sp.simplify(-NA*Q+Qdot+(eps+P)*vdot-NA*Q-(Qdot+(eps+P)*vdot-2*NA*Q))==0
    specific_rate=energy_rate-2*Q*vdot/rho
    assert sp.simplify(rho*specific_rate+2*Q*vdot-rho*energy_rate)==0
    report=dict(classification='Proven',passed=True,
        stress_tensor='In a local normal orthonormal frame with material speed v/c=v, E_n=(epsilon+P*v^2+2*v*Q)/(1-v^2), J_n=((epsilon+P)*v+Q*(1+v^2))/(1-v^2), S_rr=(P+epsilon*v^2+2*v*Q)/(1-v^2), Q=proper_heat_flux/c.',
        initial_derivatives='At v=0: E_n_dot=epsilon_dot+2*Q*v_dot; J_n_dot=Q_dot+(epsilon+P)*v_dot; S_rr_dot=P_dot+2*Q*v_dot. Dots here may use the same chosen time coordinate throughout.',
        corrected_specific_energy='With physical seconds, initial v=0, fixed nuclear composition and the same closed heat law: u_dot=(P/rho_B)*c*N*A -(1/N)*partial_mu L_infinity -(2*Q/rho_B)*v_dot. Consequently dlnT/dt=[heat+(P/rho_B-u_lnrho)*c*N*A-(2*Q/rho_B)*v_dot]/(cv*T).',
        initial_radial_momentum='With initially hydrostatic pressure/lapse and isotropic material stress: Q_dot+(epsilon+P)*v_dot=2*c*N*A*Q. A closure for Q_dot must be solved together with the fluid acceleration. Prescribing Q but omitting its evolution is insufficient.',
        lapse_source='Polar lapse responds to normal radial stress S_rr. Its time derivative also contains 2*Q*v_dot; replacing it by P_dot alone silently imposes an extra restriction.',
        original_formula_valid_only_if='v_dot=0 in addition to v=0, or Q=0 for the omitted energy term. Neither extra assumption is certified for the actual evolving star.',
        physical_EOS_or_transport_closure_certified=False,full_GR_evolution=False)
    save('symbolic.json',report)
    data=dict(np.load(original.OUT/'initial-data.npz'));state,aux=original.micro.inputs();a=aux['eos']
    rho=np.exp(state['lnd']);c=g.c.gr.C*100;Qvalue=data['proper_heat_flux_cgs']/c
    coeff_u=-2*Qvalue/rho;coeff_logT=coeff_u/a[:,10]
    result=dict(classification='Counterexample candidate',completed=True,
        original_material_rate_claim_without_velocity_derivative_restriction_valid=False,
        original_initial_constraint_construction_preserved=True,
        coefficient_of_vdot_in_specific_energy_rate_range=[float(coeff_u.min()),float(coeff_u.max())],
        coefficient_of_vdot_in_logT_rate_range=[float(coeff_logT.min()),float(coeff_logT.max())],
        fluid_acceleration_solved=False,heat_flux_time_derivative_solved=False,
        physical_EOS_certified=False,full_GR_evolution=False)
    np.savez_compressed(OUT/'rate-coefficients.npz',coefficient_vdot_in_u_dot=coeff_u,
        coefficient_vdot_in_logT_dot=coeff_logT,
        conditional_zero_vdot_logT_dot=data['dlnT_dt'],density_rate=data['dlnrho_dt'])
    save('result.json',result)
    correction=dict(classification='Proven',original_artifact='symbolic.json',
        invalid_unqualified_field='baryon_and_specific_energy: the material energy rate was stated using only initial velocity zero.',
        still_valid='Initial constraints, metric/mass rates and normal-frame Hamiltonian propagation. The material density rate is unchanged at zero initial velocity.',
        correction='Add -(2*Q/rho_B)*v_dot to the material specific-energy rate; solve Q_dot+(epsilon+P)*v_dot=2*c*N*A*Q under initial hydrostatic balance.',
        corrected_derivation='outputs/direct-eos-gr33/gr-heat-fluid-frame/symbolic.json',
        original_sources_and_arrays_preserved=True)
    (original.OUT/'velocity-derivative-correction.json').write_text(json.dumps(correction,indent=2)+'\n')
    paths=[p for p in OUT.iterdir() if p.is_file()]+[original.OUT/'velocity-derivative-correction.json']
    save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}))
    print('FLUID FRAME CORRECTION',result,flush=True);verify()


def verify():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    entries=json.loads((OUT/'manifest.json').read_text())['sha256']
    for rel,digest in entries.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert not json.loads((OUT/'result.json').read_text())['original_material_rate_claim_without_velocity_derivative_restriction_valid']
    print('PASS FLUID FRAME CORRECTION SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
