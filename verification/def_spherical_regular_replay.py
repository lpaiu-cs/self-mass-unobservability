"""Counterexample candidate: saved finite constraints and trace-resolution audit.

No EOS calls or trajectory rerun. This is not an independent native residual
evaluation: those residuals come from the completed, source-bound stages.
"""
from pathlib import Path
import json
import numpy as np
import def_spherical_regular as s


def check():
    initial=np.load(s.old.imported.initial.OUT/'initial.npz')
    rows=[]
    for cells in [16,5735]:
        folder=s.OUT/f'cells-{cells}'
        for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
            for rel,sha in json.loads((folder/name).read_text())[key].items():
                assert s.e.digest(s.ROOT/rel)==sha,rel
        V=initial['volume'][:cells];r=initial['r'][:cells];rf=initial['rf'][:cells+1]
        heat0=np.exp(initial['base'][:cells,0])*initial['aux'][:cells,10]
        for index in range(2):
            d=np.load(folder/f'path-{index}/endpoint.npz')
            dm=d['dmf'];term=s.e.GRAV*V*d['dEt'];error=np.diff(dm)-term
            scale=np.maximum(abs(dm[:-1])+abs(dm[1:])+abs(term),s.ld('1e-100'))
            row_error=float(np.max(abs(error)/scale));assert row_error<2e-15
            a0=d['a']-d['da'];b0=1/a0**2
            f=(r**3-rf[:-1]**3)/(rf[1:]**3-rf[:-1]**3)
            midmass=(1-f)*dm[:-1]+f*dm[1:]
            rhs=4/(np.sqrt(b0)+np.sqrt(d['b']))**2*midmass/r
            secant_error=float(np.max(abs(d['da']/(a0+d['da']/2)-rhs)/np.maximum(abs(rhs),s.ld('1e-100'))))
            assert secant_error<2e-15
            direct=d['gstar']*d['psi']
            removed=d['residual'].copy();removed[:,1]-=direct/heat0
            rows.append(dict(cells=cells,path=index,radial_row_scaled_error=row_error,
                absolute_radial_constraint_error_cm=float(abs(error).max()),metric_secant_relative_error=secant_error,
                weighted_direct_trace_work_erg=float(np.sum(d['H']*V*direct)),
                maximum_direct_term_over_native_energy_tolerance=float(np.max(abs(direct)/heat0)/s.ATOL[1]),
                native_norm_with_only_material_trace_term_removed=float(np.max(abs(removed)/s.ATOL)),
                wave_travel_over_outer_cell_width=float(s.e.C*s.e.TAU/128/(rf[-1]-rf[-2]))))
    return dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Frozen file bindings and radial constraints of four saved endpoints. No fresh native EOS calls.',
        trace_feedback_numerically_isolated=False,
        limitation='Deleting the material trace term at the saved state still passes the native tolerance. The full first step checks assembly and finite mass budget, not resolved reciprocal matter response.',
        no_added_evolution=True)


if __name__=='__main__':
    target=s.OUT/'saved-check.json';assert not target.exists()
    result=check();result['source_sha256']=s.e.digest(Path(__file__))
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
