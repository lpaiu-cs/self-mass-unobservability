"""Read the stored rejected photon proposal against the original joint equation."""
from pathlib import Path
from types import FunctionType
import json,resource,time
import numpy as np
import recover_joint_photon_stages as rec

out=Path('native-exact-stage205-work');receipt=out/'proposal-audit-receipt.json'
assert not receipt.exists();start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));rec.joint.previous.original.inf.incident.native.deadline(180)
error=None
try:
    initialize=FunctionType(rec.base.prior.initialize.__code__,dict(rec.base.prior.initialize.__globals__,OUT=out));initialize()
    m=rec.base.owner.Model(128);p=dict(np.load(out/'rejected-128.npz'));archive=dict(np.load(rec.saved(128)))
    step=int(np.load(out/'accepted-128.npz')['step']);t,h=p['time'][()],p['step_size'][()]
    q=archive['joint_stage_conserved_scaled'][2*step-1]
    g0=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]);g0[~m.material.active(t)]=0
    gas=p['gas'];photons=p['photon_stage_solution'];times=p['stage_times'];v=m.pack(p['x_initial'],g0)
    cs=[m.local(now) for now in times];ss=[m.source(now) for now in times];maps=[m.jacobian(now,g) for now,g in zip(times,gas)]
    src=[];rates=[];cerrors=[];nerrors=[];proof=[]
    for j,(xx,g,c,s,(J,base)) in enumerate(zip(photons,gas,cs,ss,maps)):
        affine=base-(J@g.ravel()).reshape(m.n,4)
        src.append(m.pack(s[0]/(m.scale*rec.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+affine))
        p,q,*_=m.collision(c,xx,g,True);native=m.native(times[j],g)
        rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+p+s[0]/(m.scale*rec.AMP),q+native))
        for values,key,output in [(q,'joint_collision_rates_scaled',cerrors),(native,'joint_native_rates_scaled',nerrors)]:
            old=archive[key][2*step+j];actual=values*m.units
            output.append((np.sum(abs(actual-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),rec.LD('1e-290'))).astype(float).tolist())
        proof.append(m.pack(xx,g))
    rhs=(v+ h*(rec.A@np.array(src))).ravel();sol=np.array(proof).ravel()
    defect=(np.array(proof)-v-h*(rec.A@np.array(rates))).ravel()
    scales=rec.joint.scales(m,rhs,sol);relative=float(np.linalg.norm(defect)/np.linalg.norm(rhs));physical=(rec.joint.physical_norm(m,defect)/scales).astype(float)
    row=dict(classification='Counterexample candidate',original_joint_stage_passed=relative<1e-12 and max(physical)<1e-13,
        relative=relative,physical_relative=physical.tolist(),original_gates=dict(stage=1e-12,physical=1e-13),
        collision_rate_identity=cerrors,native_rate_identity=nerrors,step=step,
        exact_archival_rate_identity_passed=False,old_failure_preserved=True,new_physical_steps=0,
        final_charge_conclusion='unadjudicated')
    rec.write(out/'proposal-audit.json',row);print(json.dumps(row),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:rec.write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=rec.sha(__file__)))
