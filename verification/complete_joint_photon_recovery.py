"""Counterexample candidate: complete reconstruction under original equations.

205's extra archival collision-rate identity remains FAILED. Reuse its five
strict-passing steps, check the actual joint equation on every remaining pair,
then check the original saved endpoint and same-solution integrated ledgers.
"""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import recover_joint_photon_stages as base
import recover_exact_joint_stages as prior

OUT=Path('native-equation-recovery208-work');OLD=prior.OUT
read,write,sha=base.read,base.write,base.sha
LD,AMP,A=base.LD,base.AMP,base.A
CAPS=dict(prepare=180,fine=1800,audit=180)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'proposal-audit.json')['original_joint_stage_passed']
    assert 'collision_relative' in read(OLD/'fine-receipt.json')['error']
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['accepted-128.npz','rejected-128.npz','expanded-recovery.py','recovered-64.npz','recovered-64.json','fine-receipt.json','proposal-audit.json']]+[base.saved(128)]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='f749d11ef',
        claim='Complete the missing fine photon history using the original joint-equation, physical-moment, saved-endpoint and integrated same-solution ledger gates, explicitly distinct from205strict rate-identity claim.',
        reason='205rejected pair independently passes the original full coupled equation1.11e-15 and physical moments3.84e-15 while a small H-rate relative identity is1.61e-12. The archive-identity failure remains; do not fit moments or accept it as an exact-rate reconstruction.',
        method='Resume all five strict-passing steps bit-identically. Use the retained sixth photon pair as a proposal and re-evaluate its original equation. Solve only the remaining two conditional photon blocks. Audit the full original joint equation at every remaining step and all integrated reconstructed collision/packet/port balances at the original stored endpoint.',
        scope='Original205T/32and128clock,8steps; no extra native-fluid evolution or new physical grid. Coarse remains205strict-passing reconstruction. Prefix5 retains its own strict identity evidence, not a claim of new independent full-state audits on those old steps.',
        gates=dict(conditional=1e-14,conditional_physical=1e-13,original_stage=1e-12,original_physical=1e-13,endpoint=1e-12,packet=1e-12,radial_port=1e-12,material_ledger=1e-8),
        strict_archival_rate_identity_admitted=False,budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='205conditional steps took17..21s; one whole-joint proposal audit25.84s including constructor. Two new photon solves plus three audits expected80..180s; allow30minutes. Accepted five-step prefix and retained pair are reused.',
        stop='Any unchanged original equation/endpoint/ledger/port gate or generous wall cap. Preserve205strict failure and any new failure. No full-period reconstruction until this original endpoint is tested.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))


def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
    q=z['joint_stage_conserved_scaled'][2*step-1]
    initial=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]);initial[~m.material.active(t)]=0
    v=m.pack(x,initial);src=[];rates=[];states=[];identity=[]
    for j,(xx,g,now,c,s) in enumerate(zip(photons,gas,times,cs,ss)):
        J,b=m.jacobian(now,g);affine=b-(J@g.ravel()).reshape(m.n,4)
        src.append(m.pack(s[0]/(m.scale*AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+affine))
        ph,q,*_=m.collision(c,xx,g,True);native=m.native(now,g)
        rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+ph+s[0]/(m.scale*AMP),q+native));states.append(m.pack(xx,g))
        old=z['joint_native_rates_scaled'][2*step+j];identity.append(float(np.max(abs(native*m.units-old))))
    rhs=(v+h*(A@np.array(src))).ravel();sol=np.array(states).ravel();defect=(np.array(states)-v-h*(A@np.array(rates))).ravel()
    rel=float(np.linalg.norm(defect)/np.linalg.norm(rhs));physical=(base.joint.physical_norm(m,defect)/base.joint.scales(m,rhs,sol)).astype(float).tolist()
    row=dict(step=step,relative=rel,physical=physical,native_identity_absolute=identity,passed=bool(rel<1e-12 and max(physical)<1e-13 and max(identity)==0))
    write(OUT/f'original-equation-{step}.json',row);assert row['passed'],row
    return row


def run():
    source=(OLD/'expanded-recovery.py').read_text()
    old="zero=np.zeros((m.n,4),LD);zeroJ=sparse.csr_matrix((4*m.n,4*m.n));logs=[];moments=[];collisions=[];packets=[];timings=[];ports=[]"
    new=old+"\n    checkpoint=dict(np.load(OLD/'accepted-128.npz'));begin=int(checkpoint['step']);assert begin==5\n    x=checkpoint['x'].copy();logs=json.loads(str(checkpoint['logs']));moments=list(checkpoint['moments']);collisions=list(checkpoint['collisions']);packets=list(checkpoint['packets']);ports=list(checkpoint['ports'])\n    for row in logs:assert np.max(row['collision_relative'])<1e-12 and row['packet_relative']<1e-12\n    seed=dict(np.load(OLD/'rejected-128.npz'));assert np.array_equal(seed['x_initial'],x)"
    changes=[("if n==128:assert read(OLD/'recovered-64.json')['passed']","assert n==128 and read(OLD/'recovered-64.json')['passed']"),(old,new),('for step in range(stop):','for step in range(begin,stop):'),
        ('sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]',"sol=seed['photon_stage_solution'].copy().ravel() if step==begin else np.tile(x[None],(2,1,1,1)).ravel();calls=[]"),
        ('if np.max(errors)>=1e-12 or packet_error>=1e-12:','if packet_error>=1e-12:'),
        ('x=pairs[-1].copy();timings.append(time.monotonic()-started)',"row['original_joint_equation']=original_equation(m,z,step,t,h,x,gas,pairs,times,cs,ss)\n        row['strict_archival_rate_identity_passed']=bool(np.max(errors)<1e-12)\n        x=pairs[-1].copy();timings.append(time.monotonic()-started)"),
        ("same_material_history_unchanged=True,new_physical_steps=0", "same_material_history_unchanged=True,reused_strict_prefix_steps=begin,reused_rejected_proposal=True,strict_archival_rate_identity_passed=False,original_equations_audited_steps=list(range(begin,stop)),new_physical_steps=0")]
    for old,new in changes:assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
    ns=dict(base.run.__globals__,OUT=OUT,OLD=OLD,original_equation=original_equation);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-recovery.py').write_text(source);ns['run'](128)


def audit():
    result=read(OUT/'recovered-128.json');assert result['passed']
    z=dict(np.load(base.saved(128)));p=dict(np.load(OUT/'recovered-128.npz'));count=len(p['times']);assert count==16
    q=z['conserved_material_history'][1]/AMP;k=z['energy_offset_reference'][1]/AMP
    actual=np.column_stack([q[2]-k,q[3],q[0],q[1]])+z['material_floor_discard_history_scaled'][1]
    expected=np.sum(p['weights'][:,None,None]*(z['joint_native_rates_scaled'][:count]+p['collision_rates']),axis=0,dtype=LD)
    error=np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))
    ports=np.sum(p['weights'][:,None,None]*p['radial_ports'],axis=0,dtype=LD);old=z['radial_ports'][1]
    port=float(np.max(abs(ports-old)/np.maximum(abs(old),LD('1e-290'))))
    packet=float(np.sum(abs(p['angular']-z['accepted_angular_luminosity'][:count]))/max(np.sum(abs(z['accepted_angular_luminosity'][:count])),LD('1e-290')))
    row=dict(classification='Counterexample candidate',original_endpoint_and_ledger_passed=bool(max(error)<1e-8 and port<1e-12 and packet<1e-12),material_ledger=error.astype(float).tolist(),radial_port_relative=port,packet_relative=packet,
        reconstructed_steps=8,reused_strict_steps=5,original_equation_audited_steps=[5,6,7],endpoint_relative=result['endpoint_relative'],strict_archival_rate_identity_passed=False,original205failure_preserved=True,full_period_recovery=False,GR_source_time_admission_passed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',row);print(json.dumps(row),flush=True);assert row['original_endpoint_and_ledger_passed'],row


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        (run if action=='fine' else globals()[action])()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
