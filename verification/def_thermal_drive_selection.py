"""Proven: DEF source and thermal-readout power counting on a smooth branch.

This checks the actual matter-source coupling, not a fitted relaxation time.
The resolvent is assumed bounded as phi0 tends to zero; critical/scalarized
branches and scalar-independent tidal/heating drives are outside the result.
"""
import json
from pathlib import Path

import sympy as sp


def frame_mapping_check():
    # Matter is minimally coupled to the Jordan metric A^2*g_E. Baryon
    # density and temperature therefore map as n_E=A^3*n_J, T_E=A*T_J.
    A, n, T = sp.symbols('A n T', positive=True)
    free_energy = sp.Function('f')
    F = A**4*free_energy(n/A**3, T/A)  # free energy per Einstein volume
    energy = F-T*sp.diff(F, T)
    pressure = n*sp.diff(F, n)-F
    assert sp.simplify(A*sp.diff(F, A)-(energy-3*pressure)) == 0
    assert sp.simplify(F.subs(A, 1)-free_energy(n, T)) == 0
    # Photon gas is traceless; a cold conserved rest-mass gas is not.
    mass, radiation = sp.symbols('mass radiation', positive=True)
    cold = A**4*mass*n/A**3
    photons = -A**4*radiation*(T/A)**4/3
    assert sp.simplify(A*sp.diff(cold, A)-cold) == 0
    assert sp.diff(photons, A) == 0
    return dict(classification='Proven', passed=True,
        definitions='A=exp(beta*phi^2/2), n_E=A^3*n_J, T_E=A*T_J, f_E(n_E,T_E,phi)=A^4*f_J(n_E/A^3,T_E/A,composition).',
        identity='At fixed Einstein-frame baryon density, temperature and composition: d f_E/d ln A = epsilon_E-3P_E; d f_E/d phi = beta*phi*(epsilon_E-3P_E).',
        primitive_mapping='P_E=A^4*P_J; epsilon_E=A^4*epsilon_J; specific total energy per conserved baryon mass scales as A. The nuclear rest-energy contribution must scale along with thermal energy.',
        controls=['Symbolic arbitrary differentiable free energy.', 'A=1 exact identity.',
            'Radiation free energy has zero scalar trace source.', 'Cold rest-mass gas has a nonzero trace source.'],
        implementation_boundary='EOS frame map only. Radiation/conduction constitutive transport, metric/scalar dynamics, boundary matching and conservative time discretization still require consistent coupling. Do not merely add a scalar body force while holding rest energy fixed.')


def check():
    beta, trace, phi0, epsilon = sp.symbols('beta trace phi0 epsilon')
    psi, delta_psi, grad_psi, grad_delta = sp.symbols('psi delta_psi grad_psi grad_delta')
    # Covector form of the Einstein-frame matter balance:
    # nabla_mu T^mu_nu = beta*phi*T*d_nu(phi).
    phi = phi0*(psi+epsilon*delta_psi)
    gradient = phi0*(grad_psi+epsilon*grad_delta)
    source = beta*trace*phi*gradient
    linear = sp.diff(source, epsilon).subs(epsilon, 0)
    expected = beta*trace*phi0**2*(psi*grad_delta+delta_psi*grad_psi)
    assert sp.expand(linear-expected) == 0
    assert source.subs(phi0, 0) == 0
    # Scalar-independent GR stress corrections start at phi0^2. Their
    # contribution to this already quadratic source starts at phi0^4.
    stress_correction = sp.symbols('stress_correction')
    corrected = source.subs(trace, trace+phi0**2*stress_correction)
    assert sp.Poly(sp.expand(corrected-source), phi0).monoms() == [(4,)]
    # Keep the common orbital drive explicit: delta_phi_ext=phi0*delta_u.
    resolvent, f, h, companion, delta_u = sp.symbols('resolvent f h companion delta_u')
    fluid_change = phi0**2*resolvent*f*delta_u
    thermal_charge_change = phi0*h*fluid_change
    force_change = phi0*companion*thermal_charge_change
    assert sp.expand(force_change-phi0**4*companion*h*resolvent*f*delta_u) == 0
    independent_heat = sp.symbols('independent_heat')
    independent_force = phi0*companion*(phi0*h*independent_heat)
    assert sp.Poly(independent_force, phi0).monoms() == [(2,)]
    # A bounded resolvent is substantive: approaching a pole can compensate
    # the prefactor. Do not convert power counting into a uniform small bound.
    s, rate = sp.symbols('s rate')
    scalar_response = phi0**2/(s+rate)
    assert sp.cancel(scalar_response.subs(s, 0).subs(rate, phi0**2)) == 1
    return dict(classification='Proven', passed=True,
        theory='Massless DEF with A(phi)=exp(beta*phi^2/2), smooth unscalarized branch; not scalarized or critical.',
        source_covector='nabla_mu T^mu_nu = beta*phi*T*d_nu(phi)',
        companion_drive='phi=phi0*(psi0+epsilon*delta_psi); delta_phi_ext=phi0*delta_u. Companion feedback must remain regular as phi0 tends to zero.',
        linear_source=str(expected),
        thermal_force_order='With bounded fluid/metric and scalar solution operators, deltaU=O(phi0^2*delta_u), deltaq_thermal=O(phi0^3*delta_u), and delta(a_scalar/g)_thermal=O(phi0^4*delta_u). The ordinary direct scalar force starts at phi0^2.',
        independent_heating_order='If deltaU is supplied independently of phi0, the corresponding scalar-force modulation starts at phi0^2 instead. Do not identify a free thermal transient with companion-scalar forcing.',
        pole_counterexample='phi0^2/(s+lambda) is order one at s=0, lambda=phi0^2. Without a uniform spectral/resolvent condition no absolute suppression bound follows.',
        integration_requirements=['Transform Jordan-frame EOS and baryon variables consistently to the Einstein frame.',
            'Include scalar stress-energy in metric constraints and matter-scalar energy exchange at the same order.',
            'Use compatible scalar initial/boundary data, not a GR-only state labeled a finite-scalar equilibrium.',
            'Compare driven and undriven paths on the same background. A changing background requires a two-time response kernel; no time-invariant transfer function is assumed.',
            'Remove only preregistered, physically justified instantaneous/derivative nuisance directions.'],
        limits='No numerical residue, actual orbital thermal mode, physical EOS error bound or observed signal is established.')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--frame-check', action='store_true')
    frame = parser.parse_args().frame_check
    result = frame_mapping_check() if frame else check()
    name = 'gr-driven-frame-map.json' if frame else 'gr-driven-source-selection.json'
    output = Path(__file__).resolve().parents[1]/'outputs/direct-eos-gr33'/name
    assert not output.exists(), 'Preserve previous result.'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
