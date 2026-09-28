"""Independent invariant and transfer audit of stored native forcing."""
from pathlib import Path
import json,sys,time
import numpy as np
import sympy as sp
import repair_native_collision_defect as run


def main():
    start=time.monotonic();run.configure();m=run.response.Response();N=m.Nweight/m.scale;W=m.Eweight/m.scale
    rows=[];changes=[];old=[];moments=[];times=[];bank_difference=0.
    for k in range(17):
        z=dict(np.load(run.OUT/f'point-{k}.npz'));p=z['photon'];b=z['bound'];e=z['escape'];d=z['number_projection'];before=p-d
        norm=np.maximum(np.sum((abs(p)+abs(b))*N,axis=(1,2))+abs(e[0]),1.)
        number=float(np.max(abs(np.sum((p-b)*N,axis=(1,2))+e[0])/norm))
        prior_number=float(np.max(abs(np.sum((before-b)*N,axis=(1,2))+e[0])/norm))
        energy_scale=np.maximum(np.sum(abs(p)*W,axis=(1,2)),1.)
        ep=float(np.max(abs(np.sum(d*W,axis=(1,2)))/energy_scale))
        mp=float(np.max(abs(np.sum(d*W*m.mu[None,:,None],axis=(1,2)))/energy_scale))
        size=float(np.max(np.sum(abs(d)*W,axis=(1,2))/energy_scale))
        recovered=np.array([np.sum(p*W,axis=(1,2))+e[1],np.sum(b*N,axis=(1,2)),
            -(np.sum(p*W*m.mu[None,:,None],axis=(1,2))+e[2])/m.a])
        export=float(np.max(np.sum(abs(recovered-z['defect_moments']),axis=1)/np.maximum(np.sum(abs(recovered),axis=1),1.)))
        assert number<1e-12 and max(ep,mp,export)<1e-12 and size<1e-8
        rows.append(dict(k=k,number_before_projection=prior_number,number=number,energy_projection=ep,momentum_projection=mp,projection_L1=size,export_relative=export))
        times.append(float(z['t']));old.append(z['original_moments']);moments.append(recovered)
        new=dict(np.load(run.response.OUT/f'bank-128/point-{k}.npz'));previous=dict(np.load(run.response.updated.run.COLL/f'bank-128/point-{k}.npz'))
        for key in new:
            bank_difference=max(bank_difference,float(np.max(abs(new[key]-previous[key]))/max(np.max(abs(previous[key])),1.)))
    old=np.array(old);moments=np.array(moments)
    ratios={name:(np.trapezoid(np.sum(abs(moments[:,:,part]),axis=2),times,axis=0)/np.maximum(np.trapezoid(np.sum(abs(old[:,:,part]),axis=2),times,axis=0),1.)).tolist()
        for name,part in [('deep',slice(0,19)),('atmosphere',slice(19,None)),('all',slice(None))]}
    assert max(x for row in ratios.values() for x in row)<.02
    R,L,H,mu=sp.symbols('R L H mu');lo=-R*H/(H-L);hi=R*L/(H-L)
    assert sp.simplify(lo+hi+R)==0 and sp.simplify(lo*L+hi*H)==0 and sp.simplify(mu*(lo*L+hi*H))==0
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,integrated_net_relative=ratios,
        projection_identity=dict(classification='Proven',scope='Two-node number correction has zero energy and same-angle radial momentum in exact arithmetic.'),
        regenerated_bank_relative_difference=bank_difference,actual_velocity_owner=m.model.velocity.__func__.__module__,
        old_missing_j_diagnosis_retracted=True,bookkeeping_correction='Producer field export_number_before_projection accidentally records the second correction residual; use independent number_before_projection here.',
        seconds=time.monotonic()-start,full_uniform_error_bound=False)
    run.write(run.OUT/'audit.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


def material_gate():
    import def_native_collision_material_charge as physical
    physical.configure();m=physical.Material(128,64);z=dict(np.load(physical.OUT/'pilot-64.npz'))
    t=float(z['time']);j,f,_,_=m.fields(t);Q=(1-f)*m.point(j)['Q']+f*m.point(j+1)['Q']
    active=m.active(t);units=np.maximum(abs(Q),1.);units[1]=np.maximum(Q[0]*physical.previous.prior.C**2,1.)
    relative=abs(z['delta_scaled'])*physical.AMP/units;relative[:,~active]=0
    component,cell=np.unravel_index(np.argmax(relative),relative.shape)
    pilot=json.loads((physical.OUT/'pilot-64.json').read_text());assert not pilot['passed']
    result=dict(classification='Counterexample candidate',linear_material_admissible=False,
        completed_steps=int(z['completed']),t=t,component=['baryon','momentum_c','reference_energy','neutral_H'][component],
        cell=int(cell),region='deep' if cell<m.nb else 'atmosphere',radius_E=float(m.rE[cell]),
        relative=float(relative[component,cell]),gate=1e-6,background=float(Q[component,cell]),
        actual_delta=float(z['delta_scaled'][component,cell]*physical.AMP),component_maxima=np.max(relative,axis=1).tolist(),
        conserved_baryon_relative_max=float(np.max(relative[0])),neutral_inventory_positive=bool(np.all(Q[3,active]+z['delta_scaled'][3,active]*physical.AMP>0)),
        numerical_conservation_passed=max(pilot['balance_relative'])<1e-8,directional_check_passed=pilot['directional_relative']<.002,
        donor_branch_passed=pilot['physical_branch_ratio']<.01,material_production_started=False,GR_readout=False,
        next='A finite-amplitude material/constitutive response or an explicit nonlinear remainder bound is required before applying this correction to GR. Do not reduce forcing or relabel the failed linear small-state gate as passed.',
        frozen_plan=run.sha(physical.OUT/'plan.json'),native_forcing_changed=False)
    physical.write(physical.OUT/'pilot-failure.json',dict(error='Original small-state gate failed',charged_pilot_seconds=60,completed_paths=['pilot-64'],physical_replay=False))
    physical.write(physical.OUT/'gate-result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':(material_gate if sys.argv[1:]==['material'] else main)()
