"""Stored-result audit; no native calls or additional response/evolution solves."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import def_reactive_structure_balanced as balanced

old=balanced.old
h=old.h
OUT=old.chemical.paired.old.thermal.OUT


def main():
    assert not (OUT/'completion.json').exists()
    # The homogeneous first-order mass-defect block with J(0)=0 has the exact
    # solution J=0. A componentwise relative residual of its roundoff-only
    # numerical values is one; eliminate that known solution, then test every
    # original row. Do not loosen the residual denominator or gate.
    source=balanced.source
    anchor='    scale=np.asarray(abs(matrix).sum(1)).ravel()'
    assert source.count(anchor)==1
    source=source[:source.index(anchor)]+'    return matrix,rhs,grid\n'
    source=source.replace('def solve(', 'def assemble(',1)
    ns=dict(balanced.ns);exec(compile(source,__file__,'exec'),ns)
    fn,proof=old.symbolic();saved=np.load(old.surface.OUT/'background.npz')
    matrix,rhs,grid=ns['assemble'](1,1,fn,old.exterior_rows(saved),drive=1.)
    y=np.load(balanced.OUT/'adiabatic-1-1-zero.npz')['response'].copy()
    original_J=float(np.max(abs(y[:,4])));y[:,4]=0
    residual=float(np.max(abs(matrix@y.ravel()-rhs)/(np.asarray(abs(matrix)@abs(y.ravel())).ravel()+abs(rhs)+1e-100)))
    assert residual<1e-9,residual
    np.savez_compressed(balanced.OUT/'adiabatic-exact-zero-defect.npz',grid=grid,response=y)
    r,m,p,e,v,A=sp.symbols('r m p e v A',positive=True);b=1-2*m/r
    M=4*sp.pi*r*r*A**4*e+r*r*b*v*v/2
    nu=m/(r*r*b)+4*sp.pi*r*A**4*p/b+r*v*v/2
    lam=(M/r-m/r**2)/b
    assert sp.simplify(nu+lam-(4*sp.pi*r*A**4*(e+p)/b+r*v*v))==0
    bindings=[]
    for file in sorted(OUT.rglob('*.json')):
        record=json.loads(file.read_text())
        if not isinstance(record,dict):continue
        for rel,sha in record.get('bindings',{}).items():
            assert h.digest(h.ROOT/rel)==sha,(str(file),rel)
            bindings.append((str(file.relative_to(h.ROOT)),rel))
    result=json.loads((balanced.OUT/'result.json').read_text())
    chemical=json.loads((old.chemical.OUT/'result.json').read_text())
    native=json.loads((OUT/'sources.json').read_text())
    assert chemical['forcing_gate_passed'] and native['baryon_source_relative']<1e-10
    assert result['spatial_relative_difference']['surface_radial_fraction_rate']<.02
    assert result['mass_first_law_relative_difference']<.02
    assert result['adiabatic_displacement_relative_difference']<.001
    row=dict(classification='Counterexample candidate',progress_class='loophole progress',
        same_background_native_reaction_values=True,chemical_energy_direction_connected=True,
        chemical_direction_cells=chemical['cells'],chemical_direction_gate_passed=chemical['forcing_gate_passed'],
        quasi_static_free_surface_mass_and_displacement_gates_passed=True,
        linear_residual_max_after_exact_homogeneous_elimination=max(residual,max(r['linear_residual'] for r in result['rows'])),
        adiabatic_control=dict(original_roundoff_J_max=original_J,exact_J_zero_residual=residual,
            explanation='J has zero forcing and zero central boundary. Eliminate its exact zero solution and recheck all original matrix rows with unchanged residual normalization; retain the original residual=1 failure.'),
        integrating_factor_symbolic_passed=True,verified_plan_bindings=len(bindings),
        reactive_charge_grid_relative=result['spatial_relative_difference']['normalized_charge_rate'],
        reactive_charge_spatial_gate_passed=False,physical_time_evolution=False,
        thermal_transport_evolved=False,physical_radiative_atmosphere=False,full_dynamic_charge_solved=False,
        scope='Actual same-background native reaction/chemical-energy and quasi-static free-surface structure connection. No finite physical trajectory or certified reactive charge.',
        remaining='Evolve thermal/reactive matter with inertia, heat flux and physical surface closure on this same background. A well-balanced charge readout is needed before using the tiny reactive charge; do not automatically refine the mesh or start orbit-long integration.')
    h.write(OUT/'completion.json',row)
    print(json.dumps(row),flush=True)


if __name__=='__main__':main()
