"""Conjectural physical benchmark, separated from the fast pulse control.

Declared model parameters, not measurements or observational constraints.
Newtonian wide-binary drive and leading weak-self-gravity scalar charge only.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import def_resolved_scalar_pulse as s


def define():
    x,a,ecc=sp.symbols('x a eccentricity',positive=True)
    phi,beta,u=sp.symbols('phi beta delta_u')
    assert sp.expand((beta*phi*(phi*u))-beta*phi**2*u)==0
    assert sp.simplify(1/(a*(1-ecc))-1/a-ecc/(a*(1-ecc)))==0
    data=np.load(s.OUT/'initial.npz');M=s.e.GRAV*data['energy'].sum();R=data['rf'][-1]
    ph0=s.ld('.001');b=s.ld(-4);orbital_a=10*R;ecc=s.ld('.1')
    alpha_companion=b*ph0
    period=2*np.pi*np.sqrt(orbital_a**3/(2*M*s.e.C**2))
    field_excursion=abs(alpha_companion)*M/orbital_a*ecc/(1-ecc)
    du=field_excursion/ph0
    target=s.ld('1e-9');error=s.ld('1e-10')
    return dict(classification='Conjectural',kind='Declared equal-mass weak-field binary benchmark; not a measured binary or an executed orbit.',
        parameters=dict(beta=float(b),background_phi=float(ph0),companion_mass_over_star=1,
            semimajor_axis_over_star_radius=10,eccentricity=float(ecc)),
        inputs=dict(material_mass_geom_cm=float(M),radius_cm=float(R),initial_sha256=s.e.digest(s.OUT/'initial.npz')),
        leading_drive=dict(alpha_companion=float(alpha_companion),period_seconds=float(period),
            omega_R_over_c=float(2*np.pi*R/(period*s.e.C)),
            maximum_scalar_excursion_from_inverse_semimajor_reference=float(field_excursion),
            maximum_delta_u=float(du),maximum_linear_delta_logA=float(abs(b)*ph0*field_excursion)),
        drive_definition='phi_ext(t)=-alpha_B M_B/r_AB(t), r=a(1-e cos E), E-e sin E=Omega t. Include the DC field in the stationary background; actual source includes nonzero phi0. alpha_B=beta*phi0 is a leading weak-body approximation, not a matched charge certificate.',
        primary_readout='delta(alpha_A/phi0) from a scalar exterior matched to infinity, after paired undriven background subtraction; report the direct field response and material-mediated response separately. Force readout is alpha_B delta(alpha_A). Temperature is an intermediate state only.',
        comparator='Same-inventory quasistatic equilibrium response with exterior matching; optional free real derivatives through order four must be explicitly labeled as a restricted nuisance class. Three or more orbital harmonics, with withheld-frequency prediction. No novelty claim against arbitrary unbounded derivative order.',
        numerical_target=dict(charge_over_phi0_residual=float(target),combined_error_budget=float(error),
            target_is_not_observational_sensitivity=True,
            normalized_thermal_response_needed_for_target=float(target/(ph0**2*du))),
        execution_rule='First solve normalized matter/charge sensitivities or bound them. Test background evolution over the intended orbital interval before a stationary transfer-function approximation. Do not run an orbit merely to subtract nearly equal native states. If the certified residual upper bound is below the declared 1e-9 computational target, close this benchmark without a long run; this is not a universal no-go.',
        fast_pulse_relation=dict(pulse_duration_seconds=float(2*R/s.e.C),pulse_amplitude=2e-4,
            pulse_is_companion_matched=False,period_over_pulse_duration=float(period/(2*R/s.e.C)),
            pulse_amplitude_over_leading_orbital_excursion=float(s.ld('.0002')/field_excursion)),
        unclosed=['Matched companion charge and nonzero-background stellar initial data',
            'Normalized charge response and physical orbital interval background evolution',
            'Physical heat-transport coefficients, EOS and atmosphere uncertainty',
            'Static/derivative subtraction and observational likelihood'],
        source='Khalil et al. 2022, arXiv:2206.13233v2, effective point-particle scalar charge and Newtonian binary Hamiltonian; parameters here are newly declared, not their neutron-star results.')


if __name__=='__main__':
    target=s.OUT/'companion-benchmark.json';assert not target.exists()
    result=define();result['source_sha256']=s.e.digest(Path(__file__))
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
