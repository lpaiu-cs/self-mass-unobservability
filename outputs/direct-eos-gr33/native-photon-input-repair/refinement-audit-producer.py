"""Read the retained extra corrections without accepting their failed tighter gate."""
from pathlib import Path
from types import FunctionType
import json,os,resource,time
import numpy as np
import restart_photon_recovery_from_archive as owner

prior=owner.prior;OLD=owner.OUT;OUT=OLD/'refinement-audit'
read,write,sha=owner.read,owner.write,owner.sha
assert not OUT.exists();OUT.mkdir();start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.base.joint.previous.original.inf.incident.native.deadline(300)
for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-128']:(OUT/folder).mkdir(parents=True)
for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:os.link(p,OUT/p.relative_to(OLD))
write(OUT/'plan.json',dict(classification='Conjectural',
    claim='Measure whether the retained extra photon corrections change the failed physical port; do not accept or hide the failed1e-18target.',
    method='Reuse220accepted15steps, initial state and actual gas; evaluate the saved final correction with the unchanged original whole equation, native identity, endpoint, material ledger, angular and radial ports. No new linear solve or physical step.',
    seconds=300,CPU_threads=1,virtual_GiB=6,bindings={str(p):sha(p) for p in [Path(__file__),OLD/'refinement/rejected-linear-128.npz',OLD/'current-pair.npz',OLD/'accepted-128.npz',OLD/'refine-receipt.json',OLD/'clock-128/snapshot-02.json',prior.prior.saved(128)]},final_charge_conclusion='unadjudicated'))
try:
    FunctionType(prior.base.base.prior.initialize.__code__,dict(prior.base.base.prior.initialize.__globals__,OUT=OUT))()
    z=dict(np.load(prior.prior.saved(128)));m=prior.interval_model(z,128,15)
    p=dict(np.load(OLD/'current-pair.npz'));new=dict(np.load(OLD/'refinement/rejected-linear-128.npz'));saved=dict(np.load(OLD/'accepted-128.npz'))
    assert int(p['step'])==int(new['step'])==int(saved['step'])==15 and np.array_equal(p['x_initial'],new['x_initial']) and np.array_equal(new['x_initial'],saved['x'])
    pairs=new['solution'].reshape(p['solution'].shape);gas=p['gas'];times=z['joint_stage_times'][30:32]
    cs=[m.local(t) for t in times];ss=[m.source(t) for t in times]
    ports=list(saved['ports']);collisions=list(saved['collisions']);packets=list(saved['packets'])
    for xx,g,c,t in zip(pairs,gas,cs,times):
        collisions.append(m.collision(c,xx,g,True)[1]*m.units)
        ports.append(m.boundary_ports(float(t),xx)*owner.AMP);packets.append(m.angular[-1])
    ns=dict(prior.prior.prior.original_equation.__globals__,OUT=OUT/'clock-128',restored_gas=prior.restored_gas)
    exec(compile((OLD/'expanded-original-128.py').read_text(),__file__,'exec'),ns)
    whole=ns['original_equation'](m,z,15,p['time'][()],p['step_size'][()],p['x_initial'],gas,pairs,times,cs,ss)
    snap=FunctionType(prior.prior.snapshot.__code__,dict(prior.prior.snapshot.__globals__,OUT=OUT))
    try:row=snap(128,m,z,15,pairs[-1],collisions,ports,packets)
    except AssertionError:row=read(OUT/'clock-128/snapshot-02.json')
    oldrow=read(OLD/'clock-128/snapshot-02.json')
    weights=z['joint_stage_weights'][:32]
    actual=np.sum(weights[:,None,None]*np.array(ports),axis=0,dtype=owner.LD)
    original=np.sum(weights[:,None,None]*np.r_[saved['ports'],[m.boundary_ports(float(t),xx)*owner.AMP for t,xx in zip(times,p['solution'])]],axis=0,dtype=owner.LD)
    result=dict(classification='Counterexample candidate',original_equation=whole,snapshot=row,
        old_radial_port_relative=oldrow['radial_port_relative'],outer_number_port_change=float(actual[1,0]-original[1,0]),
        extra_linear_correction_solution_relative=float(np.linalg.norm(new['solution']-p['solution'].ravel())/np.linalg.norm(p['solution'])),
        conditional_tighter_gate_failed=True,original_port_gate_passed=row['passed'],new_physical_steps=0,new_Krylov_iterations=0,final_charge_conclusion='unadjudicated')
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
