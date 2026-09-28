"""Use the full original carrier inventory and signed archived luminosities."""
from pathlib import Path
from types import FunctionType
import json
import shutil
import time
import numpy as np
import def_free_surface_thermal as old
import def_thermal_native_carrier as failed

OUT=old.OUT


def inputs():
    bg,raw,state=old.inputs();data,_=old.h.inputs()
    state['L']=data['L'].copy()
    assert np.all(np.isfinite(state['L'])) and np.all(state['L']!=0)
    return bg,raw,state


ns=dict(vars(old),inputs=inputs)
source=FunctionType(old.source.__code__,ns)


def main():
    assert not (OUT/'inventory-source-plan.json').exists()
    old.h.write(OUT/'inventory-source-plan.json',dict(classification='Counterexample candidate',
        bindings={str(p.relative_to(old.h.ROOT)):old.h.digest(p) for p in [Path(__file__),Path(old.__file__),Path(failed.__file__),old.h.OLD/'molecular-state-17-8.npz']},
        preserved_failures=['32-cell carrier failed set_vars after native evaluations; original log retained.',
            'The attempted positive-L carrier guard rejected 435 negative archived luminosities before any new native run. Signed luminosities are valid carrier data; no sign change or absolute-value replacement.'],
        intervention='Keep all 5735 original carrier baryon cells and signed archived luminosities. Only requested rho,T,X and actual new EOS auxiliaries define source values. Carrier luminosities do not define physical transport.',
        budget=dict(full_batches=2,hard_timeout_seconds=60,CPU_threads=1,expected_seconds=[18,40],native_EOS_calls=0),
        cost_basis='The identical full-grid tracer previously took 8.50 seconds. The 32-cell carrier was not a valid speed pilot; use the required full-grid baseline as first measured batch and reuse it.',
        stop='Stop on source/profile/input failure, or first batch above 25 seconds; no larger run.'))
    rows=np.arange(5735);native,value,first=source('inventory-actual',rows,[0,0])
    old.h.write(OUT/'inventory-source-pilot.json',dict(classification='Counterexample candidate',seconds=first,cells=len(rows),within_budget=first<25))
    assert first<25
    other,control,second=source('inventory-control',rows,[1,-1])
    equal={k:bool(np.array_equal(value[k],control[k])) for k in ['dxdt','heat','neutrino']};assert all(equal.values()),equal
    _,profile=old.g.c.mesa(OUT/'inventory-actual-profile.data.gz');thermal=profile['non_nuc_neu']
    coeff=np.load(OUT/'coefficients.npz');R=value['dxdt'];dm=coeff['dm'];A=coeff['A'];N=coeff['N']
    rest=(old.g.c.W/old.g.c.A-1)*(old.h.gr.C*100)**2
    heating=-R@rest-value['neutrino']-thermal
    baryon=float(np.max(abs(R.sum(1))/np.maximum(abs(R).sum(1),1e-100)));assert baryon<1e-10,baryon
    np.savez_compressed(OUT/'sources.npz',dxdt=R,neutrino=value['neutrino'],thermal_neutrino=thermal,
        rest_to_internal_heating=heating,heat_Q_reference=value['heat'])
    result=dict(classification='Counterexample candidate',cells=len(R),rate_values_bitwise_independent_of_derivative_arguments=equal,
        baryon_source_relative=baryon,proper_heating_range_erg_g_s=[float(heating.min()),float(heating.max())],
        redshifted_net_power_erg_s=float(dm@(A*A*N*N*heating)),maximum_species_rate=float(abs(R).max()),
        seconds=first+second,internal_energy_composition_derivative_included=False,returned_Jacobian_used=False,
        finite_thermal_reactive_step=False,full_dynamic_charge_solved=False)
    old.h.write(OUT/'sources.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
