"""Read the saved rejected pair; test the native energy-coordinate round trip."""
from pathlib import Path
from types import FunctionType
import json,os,time,numpy as np
import continue_precise_native as prior
OUT=Path('native-energy-coordinate201-work');OLD=prior.OUT
assert not OUT.exists();OUT.mkdir()
for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
    dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst)
prior.joint.previous.original.inf.incident.native.deadline(120)
start=time.monotonic();init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT));init(False)
m=prior.owner.Model(64);engine=prior.joint.previous.engine;LD=prior.LD
z=dict(np.load(OLD/'last-accepted-64.npz'));f=dict(np.load(OLD/'failed-linear-64.npz'))
t,h=z['next_time'][()],z['next_step'][()];pairs=[m.unpack(row)[1] for row in f['solution'].reshape(2,-1)]
old=np.array([m.native(t+c*h,g) for c,g in zip(prior.joint.C,pairs)])
mark='et=z.astype(LD).copy();et[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*et[0]'
assert engine.source.count(mark)==1
source=engine.source.replace(mark,'et=z.astype(LD).copy();et[2]=m.native_etilde')
ns=dict(engine.tangent.__globals__);exec(compile(source,__file__,'exec'),ns)
new=[];energy=[]
for c,g in zip(prior.joint.C,pairs):
    exact=g[:,0]*m.eu;q=m.conserved(g)
    reconstructed=q[2]-(m.material.a.astype(LD)-m.material.model.m.a0)*m.material.model.cx*prior.joint.previous.engine.C**2*q[0]
    m.material.native_etilde=exact
    new.append(m.native(t+c*h,g,tangent=ns['tangent']))
    error=reconstructed-exact
    energy.append(dict(cell261_direct=str(exact[261]),cell261_roundtrip_error=str(error[261]),
        cell261_relative=float(abs(error[261])/max(abs(exact[261]),LD('1e-290')))))
new=np.array(new);scale=np.linalg.norm(f['rhs']);A=prior.joint.A
gas=np.array(pairs);initial=z['g']
measure=lambda x:float(np.linalg.norm(x)/scale)
oldB=gas[:,:,2]-initial[:,2]-h*(A@old[:,:,2]);newB=gas[:,:,2]-initial[:,2]-h*(A@new[:,:,2])
result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
    energy=energy,old_actual_B=measure(oldB),direct_energy_actual_B=measure(newB),
    B_stage_change=measure(h*A@(new[:,:,2]-old[:,:,2])),original_gate=1e-12,
    new_physical_steps=0,actual_stage_accepted=False,final_charge_conclusion='unadjudicated')
np.savez_compressed(OUT/'coordinate-comparison.npz',old_native=old,direct_energy_native=new,old_B=oldB,new_B=newB)
prior.write(OUT/'result.json',result);print(json.dumps(result,indent=2))
