"""Saved native equation replay and independent exponential work controls."""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import json
import time
import numpy as np
import def_exponential_coupled_run as run

x,s,e,ld=run.x,run.s,run.e,run.ld


def matrix_controls():
    Q=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]]).astype(ld)
    rows=[]
    for diagonal in [[2,0,3,0],[2,0,-3,0]]:
        M=np.diag(np.array(diagonal,dtype=ld))
        for h in [ld('.1'),ld(str(np.pi))/np.sqrt(ld(6))]:
            aug=np.zeros((8,8),dtype=ld);aug[:4,:4]=h*Q@M;aug[:4,4:]=h*Q
            exp=x.exponential(aug);A,B=exp[:4,:4],exp[:4,4:]
            y=np.array([.2,-.1,.4,.3],dtype=ld);offset=np.array([.3,.2,-.2,.1],dtype=ld)
            new=A@y+B@offset;mid=(new+y)/2
            work=(M@mid+offset)@(new-y)
            energy=(new@M@new-y@M@y)/2+offset@(new-y)
            scale=max(float(abs(new@M@new)+abs(y@M@y)+abs(offset@(new-y))),1.)
            errors=[float(abs(A.T@M@A-M).max()),float(abs(A.T@M@B-(np.eye(4)-A.T)).max()),
                float(abs(B.T@M@B+B+B.T).max()),float(abs(work)/scale),float(abs(energy)/scale)]
            assert max(errors)<2e-15,errors
            rows.append(dict(diagonal=diagonal,h=float(h),errors=errors))
    return dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Singular positive/indefinite reference Hamiltonians, affine discrete gradient, and a half-period where the alternative A+I inverse is singular. No material or continuum certificate.')


def audit():
    start=time.monotonic();plan=run.bindings();result=json.loads((run.OUT/'result.json').read_text());rows=[]
    R=ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm'])
    for path in result['paths']:
        label,steps=path['label'],path['steps'];folder=run.OUT/f'{label}-{steps}';last=steps//2;h=2*R/e.C/steps
        star,_=s.initialize(s.old.imported.CachedOnly(),-4,ld('.001'),h);star.__class__=x.ExponentialStar
        star.equilibrium_pressure_flux=np.zeros(star.n+1,dtype=ld)
        star.equilibrium_gravity=np.zeros(star.n,dtype=ld);star.equilibrium_pressure=np.zeros(star.n,dtype=ld)
        initial=dict(np.load(folder/'initial.npz'));star.previous=initial
        boundary=initial['psi'][-1]+initial['Phi'][-1]*star.distance[-1]
        amplitude=ld(str(plan['amplitude']))*{'minus':-1,'undriven':0,'plus':1}[label]
        x.setup(star,initial,boundary,amplitude)
        totalB=np.sum(initial['B']*star.volume);before=initial
        maxima=dict(native_norm=0.,baryon=0.,isotope=0.,mass_identity=0.,cone_speed=0.,cone_imaginary=0.,boundary_error=0.)
        for j in range(1,last+1):
            z=dict(np.load(folder/f'step-{j:03d}.npz'))
            cone=FunctionType(s.old.two.cones.__code__,dict(s.old.two.cones.__globals__,e=SimpleNamespace(C=e.C,TAU=z['tau_cond'])))(z)
            dy=z['canonical']-before['canonical'];mid=(z['canonical']+before['canonical'])/2
            work=-R/e.GRAV*(star.clock_matrix@mid)@dy
            dm=(z['dmf'][-1]-before['dmf'][-1])/e.GRAV
            res=np.sum(z['H']*star.volume*star.heat0*z['residual'][:,1])
            scale=np.sum(star.volume*(abs(z['Ephi'])+abs(before['Ephi'])+abs(z['dEm']-before['dEm'])))+abs(work)
            boundary_actual=star.boundary_basis@z['canonical'][star.n:star.n+9]
            expected=boundary+amplitude*np.sin(ld(str(np.pi))*2*ld(j)/steps)**8
            values=dict(native_norm=float(np.max(abs(z['residual'])/s.ATOL)),
                baryon=float(abs(np.sum((z['B']-initial['B'])*star.volume)/totalB)),
                isotope=float(abs(np.sum((z['dBX']-initial['dBX'])*star.volume[:,None],axis=0)/totalB).max()),
                mass_identity=float(abs(dm-work-res)/max(scale,ld('1e-100'))),
                cone_speed=cone['maximum_local_rest_characteristic_speed_over_c'],
                cone_imaginary=cone['maximum_characteristic_imaginary_part'],boundary_error=float(abs(boundary_actual-expected)))
            for k,v in values.items():maxima[k]=max(maxima[k],v)
            before=z
        p=dict(np.load(folder/f'step-{last-1:03d}.npz'));star.previous=p
        star.field_guess=[z[k].copy() for k in ['psi','Pi','Phi']]
        star.previous_boundary=star.boundary_basis@p['canonical'][star.n:star.n+9]
        star.amplitude=star.boundary_basis@z['canonical'][star.n:star.n+9]
        y=star.base+z['delta'];la=z['logA']
        star.material_cache.update({e.material_key(row):raw for row,raw in zip(zip(y[:,0]-3*la,y[:,1]-la,y[:,5:]),z['raw'])})
        previous=(p['delta'],p)
        value,replayed=s.residual(star,z['delta'],previous,previous,h,(ld(1),ld(-1),ld(0)))
        norm=float(np.max(abs(value)/s.ATOL));difference=float(np.max(abs(value-z['residual'])/s.ATOL))
        passed=bool(norm<=1 and difference<1e-5 and maxima['native_norm']<=1 and maxima['baryon']<1e-9
            and maxima['isotope']<1e-9 and maxima['mass_identity']<2e-15 and maxima['cone_speed']<1
            and maxima['cone_imaginary']<1e-10 and maxima['boundary_error']<2e-19)
        rows.append(dict(label=label,steps=steps,passed=passed,maxima=maxima,endpoint_replay_norm=norm,
            endpoint_replay_difference=difference,endpoint_scalar_residual=replayed['scalar_residual']))
    compared=run.compare(result['paths']);assert compared['comparisons']==result['comparisons'] and compared['passed']==result['passed']
    return dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),
        time_comparison_passed=result['passed'],rows=rows,matrix_controls=matrix_controls(),
        new_native_calls=0,new_evolution_steps=0,seconds=time.monotonic()-start,source_sha256=e.digest(Path(__file__)))


if __name__=='__main__':
    target=run.OUT/'saved-audit.json';assert not target.exists()
    result=audit();e.write(target,result);print(json.dumps(result));assert result['passed']
