"""Replay nonzero-background native endpoint equations without new EOS calls."""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import json
import time
import numpy as np
import def_nonstationary_interior as n

s,e,ld=n.s,n.e,n.ld


def audit():
    start=time.monotonic();plan=n.bindings();assert n.symbolic()==plan['symbolic']
    result=json.loads((n.OUT/'result.json').read_text());rows=[]
    V=np.load(s.OUT/'initial.npz')['volume'];R=ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm'])
    for path in result['paths']:
        label,steps=path['label'],path['steps'];folder=n.OUT/f'{label}-{steps}';last=steps//2
        h=2*R/e.C/steps;star,_=s.initialize(s.old.imported.CachedOnly(),-4,ld('.001'),h)
        star.equilibrium_pressure_flux=np.zeros(star.n+1,dtype=ld)
        star.equilibrium_gravity=np.zeros(star.n,dtype=ld)
        star.equilibrium_pressure=np.zeros(star.n,dtype=ld)
        initial=dict(np.load(folder/'initial.npz'));totalB=np.sum(initial['B']*V)
        maxima=dict(native_norm=0.,baryon=0.,isotope=0.,cone_speed=0.,cone_imaginary=0.)
        for j in range(1,last+1):
            z=dict(np.load(folder/f'step-{j:03d}.npz'))
            cone=FunctionType(s.old.two.cones.__code__,dict(s.old.two.cones.__globals__,
                e=SimpleNamespace(C=e.C,TAU=z['tau_cond'])))(z)
            values=dict(native_norm=float(np.max(abs(z['residual'])/s.ATOL)),
                baryon=float(abs(np.sum((z['B']-initial['B'])*V)/totalB)),
                isotope=float(abs(np.sum((z['dBX']-initial['dBX'])*V[:,None],axis=0)/totalB).max()),
                cone_speed=cone['maximum_local_rest_characteristic_speed_over_c'],
                cone_imaginary=cone['maximum_characteristic_imaginary_part'])
            for key,v in values.items():maxima[key]=max(maxima[key],v)
        previous=dict(np.load(folder/f'step-{last-1:03d}.npz'))
        star.previous=previous;star.field_guess=[z[k].copy() for k in ['psi','Pi','Phi']]
        boundary=initial['psi'][-1]+initial['Phi'][-1]*star.distance[-1]
        star.amplitude=boundary
        sign={'minus':-1,'undriven':0,'plus':1}[label]
        t=ld(last-1)*h;tc=R/e.C
        star.previous_boundary=boundary+sign*ld(str(plan['drive_amplitude']))*np.sin(ld(str(np.pi))*t/tc)**8
        y=star.base+z['delta'];la=z['logA']
        star.material_cache.update({e.material_key(row):raw for row,raw in
            zip(zip(y[:,0]-3*la,y[:,1]-la,y[:,5:]),z['raw'])})
        p=(previous['delta'],previous)
        value,replayed=s.residual(star,z['delta'],p,p,h,(ld(1),ld(-1),ld(0)))
        replay_norm=float(np.max(abs(value)/s.ATOL))
        difference=float(np.max(abs(value-z['residual'])/s.ATOL))
        # The stored initial at=0 was an unused initialization placeholder.
        # Its physical value follows from initial energy flux; do not call it
        # a stationary-metric constraint or a replayed initial derivative.
        initial_flux=s.m.fluxes(star,initial)[1]
        initial_at=e.GRAV*(-e.C*4*np.pi*star.rf**2*initial_flux)
        initial_at=initial['a']**3*(initial_at[:-1]+star.fraction*np.diff(initial_at))/star.r
        row=dict(label=label,steps=steps,saved_states=last,maxima=maxima,
            endpoint_replay_norm=replay_norm,endpoint_replay_difference=difference,
            initial_physical_max_abs_at_per_second=float(abs(initial_at).max()),
            initial_at_placeholder_unused=True,
            passed=bool(maxima['native_norm']<=1 and maxima['baryon']<1e-9 and maxima['isotope']<1e-9
                and maxima['cone_speed']<1 and maxima['cone_imaginary']<1e-10 and replay_norm<=1 and difference<1e-5))
        rows.append(row)
    return dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),
        time_comparison_passed=result['passed'],paths=rows,seconds=time.monotonic()-start,
        new_native_calls=0,new_evolution_steps=0,source_sha256=e.digest(Path(__file__)))


if __name__=='__main__':
    target=n.OUT/'saved-audit.json';assert not target.exists()
    result=audit();e.write(target,result);print(json.dumps(result));assert result['passed']
