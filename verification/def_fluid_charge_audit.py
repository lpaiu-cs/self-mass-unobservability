"""No native calls: trace positive control and conserved-response decomposition."""
from pathlib import Path
import json
import time
import numpy as np
import def_fluid_charge as f


def audit():
    started=time.monotonic();out=f.OUT
    plan=json.loads((out/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert f.e.digest(f.s.ROOT/rel)==sha,rel
    for rel,sha in json.loads((out/'manifest.json').read_text())['sha256'].items():
        assert f.e.digest(f.s.ROOT/rel)==sha,rel
    result=json.loads((out/'result.json').read_text());grad=dict(np.load(out/'gradients.npz'))
    star,_=f.s.initialize(f.s.old.imported.CachedOnly(),-4,f.ld('.001'),f.ld(1))
    base=dict(np.load(out/'baseline.npz'));g=f.geometry(star,np.zeros(star.n,dtype=f.ld))
    B=base['a']*base['rho']*base['W']*g.volume
    C0,_=f.readout(g,base['E'],base['R'],base['trace'])
    rows=[]
    for saved in result['endpoints']:
        label,n=saved['label'],saved['steps'];z=dict(np.load(out/f'{label}-{n}.npz'))
        source=f.production.OUT/f'{label}-{n}'/f'step-{n:03d}.npz'
        assert f.e.digest(source)==saved['source_sha256']
        eta=z['a']*z['rho']*z['W']*g.volume/B-1
        theta=z['delta'][:,1]-base['delta'][:,1]
        velocity=z['v']-base['v']
        terms={key:float(grad[key]@value) for key,value in [('baryon',eta),('temperature',theta),('velocity',velocity)]}
        C,_=f.readout(g,z['E'],z['R'],z['trace']);actual=float(C-C0)
        rows.append(dict(label=label,steps=n,actual_increment=actual,linear_terms=terms,
            remaining_increment=actual-sum(terms.values()),
            relative_total_baryon_change=float(abs(np.sum(B*eta))/np.sum(B)),
            maximum_composition_change=float(abs(z['delta'][:,5:]-base['delta'][:,5:]).max())))
    # Real baseline heat flux is too small to resolve its linear velocity term
    # in the native controls. Use a declared nonzero-Q kinematic control; this
    # is algebra only, not a physically evolved heat-flow background.
    eps=base['E'];P=base['P'];heat=f.ld('.01')*(eps+P)
    direction=np.sin(np.pi*star.r/star.rf[-1])
    def moving(t):
        v=t*direction
        E=(eps+P*v*v+2*heat*v)/(1-v*v)
        R=(eps*v*v+P+2*heat*v)/(1-v*v)
        trace=-eps+3*P
        value,_=f.readout(g,E,R,trace)
        wrong,_=f.readout(g,E,R,-E+3*P)
        return value,wrong
    expected=float(moving(1e-24j)[0].imag/1e-24)
    values=[float((moving(f.ld(h))[0]-moving(-f.ld(h))[0])/(2*f.ld(h))) for h in ['.001','.0003']]
    assert abs(expected)>1e-3 and max(abs(v-expected) for v in values)<1e-4*abs(expected)
    wrong=float(moving(1e-24j)[1].imag/1e-24)
    assert abs(wrong-expected)>1e-3
    contrasts=[]
    for other in ['undriven','decoupled']:
        a=next(r for r in rows if r['label']=='driven' and r['steps']==192)
        b=next(r for r in rows if r['label']==other and r['steps']==192)
        terms={k:a['linear_terms'][k]-b['linear_terms'][k] for k in a['linear_terms']}
        actual=a['actual_increment']-b['actual_increment']
        remainder=actual-sum(terms.values())
        contrasts.append(dict(name='driven_minus_'+other,actual=actual,linear_terms=terms,
            remainder=remainder,remainder_over_actual=abs(remainder/actual)))
    return dict(classification='Counterexample candidate',passed=True,
        frozen_bindings_verified=True,new_native_calls=0,new_evolution_steps=0,
        kinematic_heat_positive_control=dict(expected=expected,central_differences=values,
            incorrect_rest_trace_derivative=wrong,scope='Manufactured Q=0.01(eps+P), fixed primitive density/temperature. Verifies Lorentz stress and trace handling, not a physical heat-flow solution.'),
        baseline_velocity_native_control_unresolved=True,
        baseline_velocity_derivative_l1=result['gradient_l1']['velocity'],
        native_derivative_absolute_floor=plan['gates']['derivative_absolute'],
        endpoints=rows,finest_contrasts=contrasts,seconds=time.monotonic()-started,
        interpretation='Linear baryon and thermal decomposition is at the original background; the remainder retains composition, nonlinear and cross contributions. Radius and Eulerian baryon terms are alternative perturbation coordinates, not additive descriptions of the same material displacement.')


if __name__=='__main__':
    target=f.OUT/'audit.json';assert not target.exists()
    result=audit();result['source_sha256']=f.e.digest(Path(__file__))
    f.e.write(target,result);print(json.dumps(result,indent=2))
