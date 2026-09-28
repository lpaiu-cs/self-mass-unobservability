"""Independent integral and density-coordinate checks of the matched envelope."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.integrate import simpson
import def_native_radiative_envelope as run


def main():
    out=run.OUT
    assert not (out/'audit.json').exists()
    run.write(out/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Reconstruct luminosity from the integrated Tolman temperature gradient, hydrostatic pressure drop and baryon mass from saved radial profiles. Check five states using native density coordinates rather than the producer pressure coordinates.',
        gates=dict(integral_relative=.001,native_pressure_relative=1e-7,native_adiabatic_relative=1e-7),
        limits='Only integral/numerical consistency of the declared grey envelope. Convective stability and exact incoming/outgoing transport are separate physical gates.',
        budget=dict(hard_seconds=60,maximum_native_calls=5),
        bindings={str(p):run.digest(p) for p in [Path(__file__),Path(run.__file__),out/'plan.json',out/'result.json',out/'fine.npz',out/'thin.npz']}))
    signal.alarm(60);start=time.monotonic()
    result=json.loads((out/'result.json').read_text())
    assert result['passed']
    p=run.Envelope();d=np.load(out/'fine.npz');r=d['r'];L=float(d['Linfinity'])
    theta=d['A']*d['N']*d['T']
    integral=simpson(d['optical']*d['N']**2/r**2,x=r)
    reconstructed=4*np.pi*p.arad*p.c/3*(theta[0]**4-theta[-1]**4)/integral
    hydro=simpson((d['energy']+d['Ptotal'])*d['g'],x=r)
    baryon=simpson(4*np.pi*r*r*d['A']**3*d['rho']/np.sqrt(d['b']),x=r)/p.total_baryon
    scores=dict(luminosity=abs(reconstructed/L-1),hydrostatic=abs(hydro/(d['Ptotal'][0]-d['Ptotal'][-1])-1),
        baryon=abs(baryon/d['states'][6,-1]-1))
    assert max(scores.values())<.001,scores
    controls=[]
    for i in [0,40,150,280,400]:
        a=p.eos(2,float(np.log(d['rho'][i])),float(np.log(d['T'][i])),p.X)
        ad=(a[1]/a[0]-a[9])/a[10]/a[4]
        row=dict(index=i,pressure_relative=abs(a[1]/d['Ptotal'][i]-1),
            nabla_ad_relative=abs(ad/d['nabla_ad'][i]-1))
        assert row['pressure_relative']<1e-7 and row['nabla_ad_relative']<1e-7,row
        controls.append(row)
    tau=d['states'][7]-d['states'][7,-1]
    proper_temperature_scale=(d['A']/np.sqrt(d['b']))/abs(d['logT_prime']+d['g'])
    Kn=1/(d['rho']*d['opacity']*proper_temperature_scale)
    excess=d['nabla']-d['nabla_ad']
    thick=(tau>10)&(Kn<.01)
    unstable=thick&(excess>0)
    pr=d['Prad'];force=d['optical']*d['F']/p.c
    pressure_identity=d['gas_pressure_prime']+4*pr*d['logT_prime']+(d['energy']+d['Ptotal'])*d['g']
    pressure_score=float(max(abs(pressure_identity)/((d['energy']+d['Ptotal'])*d['g'])))
    assert pressure_score<1e-12
    grids=p.opacity.data
    ranges={suffix:[float(max(grids[k][0] for k in grids if k.startswith('kap_z_tables-') and k.endswith(suffix))),
                    float(min(grids[k][-1] for k in grids if k.startswith('kap_z_tables-') and k.endswith(suffix)))] for suffix in ['-R','-T']}
    opaque_domain=bool(min(d['opacity_logR'])>=ranges['-R'][0] and max(d['opacity_logR'])<=ranges['-R'][1]
        and min(np.log10(d['T']))>=ranges['-T'][0] and max(np.log10(d['T']))<=ranges['-T'][1])
    original=p.d
    # Compare only overlapping pressures, avoiding the arbitrary geometric
    # identification of the old finite-pressure point with the new photosphere.
    oldP=original['raw'][:p.i+1,1]
    newT=np.interp(np.log(oldP),np.log(d['Ptotal'])[::-1],d['T'][::-1])
    oldT=np.exp(original['lnT'][:p.i+1])
    mass_fraction=d['states'][6,-1]
    original_removed=(original['dm'][:p.i].sum()+original['dm'][p.i]/2)/p.total_baryon
    value=dict(classification='Counterexample candidate',passed=True,
        integral_relative_errors={k:float(v) for k,v in scores.items()},native_density_coordinate_controls=controls,
        gas_plus_radiation_pressure_identity=pressure_score,opacity_unclipped_in_common_table_domain=opaque_domain,
        common_table_ranges=ranges,maximum_heat_stress_over_energy=float(max(d['F']/(p.c*d['energy']))),
        maximum_nabla_minus_adiabatic=float(max(excess)),opaque_adiabatic_parcel_test_samples=int(thick.sum()),
        opaque_Schwarzschild_unstable_samples=int(unstable.sum()),
        opaque_maximum_nabla_minus_adiabatic=float(max(excess[thick])) if thick.any() else None,
        Schwarzschild_scope='Homogeneous-composition adiabatic parcel criterion, reported only where tau>10 and transport Kn<0.01; not an evolved convection model or non-LTE stability proof.',
        old_outer_temperature_K=float(oldT[0]),new_temperature_at_old_outer_pressure_K=float(newT[0]),
        maximum_native_overlap_temperature_relative=float(max(abs(newT/oldT-1))),
        reconstructed_envelope_baryon_fraction=float(mass_fraction),original_native_prefix_baryon_fraction=float(original_removed),
        difference_against_native_prefix_only=float(mass_fraction-original_removed),
        inventory_scope='Native cells0..131 plus half cell132 only; the previous added isentropic atmosphere is excluded from this reference. This is not a fixed-inventory replacement comparison.',
        unchanged_interior_luminosity_matched=bool(abs(result['required_minus_original_L_relative'])<.002),
        same_inventory_whole_star=False,full_time_radial_Einstein_equation_solved=False,spectral_atmosphere=False,full_goal_complete=False,
        seconds=time.monotonic()-start,total_compute_seconds=result['total_compute_seconds']+time.monotonic()-start)
    assert value['total_compute_seconds']<600
    run.write(out/'audit.json',value);print(json.dumps(value),flush=True)


if __name__=='__main__':main()
