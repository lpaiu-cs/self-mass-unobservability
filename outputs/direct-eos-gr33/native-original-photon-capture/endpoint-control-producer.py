"""Compare actual saved gas/photon states at the two ends of the failed interval."""
from pathlib import Path
from types import FunctionType
import gc,json,os,resource,sys,time
import numpy as np
import restart_photon_recovery_from_archive as owner

prior=owner.prior;LD,AMP=owner.LD,owner.AMP
OUT=Path('native-original-state221-work');OLD=owner.OUT
read,write,sha=owner.read,owner.write,owner.sha
assert not OUT.exists();OUT.mkdir();start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.base.joint.previous.original.inf.incident.native.deadline(300)
files=[]
for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True)
for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
    os.link(p,OUT/p.relative_to(OLD));files.append(p)
checkpoints=[owner.INPUT,Path('native-front-continuation184-work/sweep-1/photons/pilot-128.npz')]
files += checkpoints+[prior.prior.saved(128),owner.SEED/'recovered-128.npz',OLD/'current-pair.npz',OLD/'refinement/rejected-linear-128.npz']
files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
write(OUT/'plan.json',dict(classification='Conjectural',
    claim='Separate the effects of conserved-to-normalized gas inversion, photon-state differences and boundary evaluation at actual stored endpoints of the failed recovery interval.',
    method='Use original restart_guide (the saved pre-floor final gas stage) and restart_x atT/32andT/16. Require their conserved state, stage time and collision history to match the actual archive. Evaluate a two-by-two gas/photon substitution without any solve or fitted coordinate. Independently evaluate the linear outer photon-number difference before subtracting large totals.',
    decision='Only a demonstrated coordinate/operator defect can justify another reconstruction repair. Preserve the port failure; do not infer every-stage or final-charge certification from these two endpoints.',
    cap_seconds=300,CPU_threads=1,virtual_GiB=6,new_physical_steps=0,new_linear_solves=0,
    forecast='Two original per-interval model constructions14..16s each plus saved-point contractions; allow5minutes. Stop on mismatched original state/rates.',
    bindings={str(p):sha(p) for p in dict.fromkeys(files+[Path(__file__)])},final_charge_conclusion='unadjudicated'))
try:
    FunctionType(prior.base.base.prior.initialize.__code__,dict(prior.base.base.prior.initialize.__globals__,OUT=OUT))()
    z=dict(np.load(prior.prior.saved(128)));rows=[]
    for step,path in zip([7,15],checkpoints):
        original=dict(np.load(path));m=prior.interval_model(z,128,step);q=z['joint_stage_conserved_scaled'][2*step+1];t=z['joint_stage_times'][2*step+1]
        assert original['t'][-1]==z['actual_step_edges'][step+1] and np.array_equal(original['joint_stage_conserved_scaled'][-1],q)
        g=original['restart_guide'];restored=prior.restored_gas(m,q);assert np.array_equal(m.conserved(g),q)
        x=original['restart_x']
        if step==7:recovered=np.load(owner.SEED/'recovered-128.npz')['endpoint_occupation']/(m.scale*AMP)
        else:recovered=np.load(OLD/'refinement/rejected-linear-128.npz')['solution'].reshape((2,*x.shape))[-1]
        c=m.local(t);m.source(t);values={}
        for xn,xx in [('original',x),('recovered',recovered)]:
            for gn,gg in [('original',g),('inverted',restored)]:values[xn,gn]=m.collision(c,xx,gg,True)[1]*m.units
        archived=z['joint_collision_rates_scaled'][2*step+1];anchor=values['original','original']
        identity=bool(np.array_equal(anchor,archived))
        differences={f'{a}_{b}':np.sum(abs(v-anchor),axis=0,dtype=LD).astype(float).tolist() for (a,b),v in values.items()}
        number_scale=max(np.sum(abs(archived[:,1]),dtype=LD),LD('1e-290'))
        bt=z['accepted_angular_times'][2*step+1]
        before=m.boundary_ports(bt,x)*AMP;after=m.boundary_ports(bt,recovered)*AMP
        delta=(recovered-x)*m.scale;mask=m.mu>0
        # The same boundary functional; subtract occupation before large port sums.
        term=delta[-1,mask]*(m.w*m.mu)[mask,None]*m.num
        direct=LD(4)*LD(np.pi)*LD(prior.base.base.owner.joint.previous.original.C)*LD(m.area[-1])*np.sum(term,dtype=LD)*AMP
        changed=np.argwhere(g!=restored)
        row=dict(classification='Counterexample candidate',step=step+1,time=float(t),original_conserved_identity=True,
            original_collision_bit_identity=identity,original_collision_max_abs=float(np.max(abs(anchor-archived))),
            changed_gas_coordinates=changed.tolist(),maximum_gas_change=float(np.max(abs(g-restored))),
            collision_L1_changes=differences,collision_number_L1_scale=float(number_scale),
            boundary_port_original=before.astype(float).tolist(),boundary_port_recovered=after.astype(float).tolist(),
            outer_number_difference_by_subtraction=float(after[1,0]-before[1,0]),outer_number_difference_direct=float(direct),
            photon_N_L1_relative=float(np.sum(abs(recovered-x)*m.Nweight,dtype=LD)/max(np.sum(abs(x)*m.Nweight,dtype=LD),LD('1e-290'))),
            photon_E_L1_relative=float(np.sum(abs(recovered-x)*m.Eweight,dtype=LD)/max(np.sum(abs(x)*m.Eweight,dtype=LD),LD('1e-290'))))
        write(OUT/f'endpoint-{step+1}.json',row);rows.append(row)
        assert identity,('Original endpoint collision identity',step+1,row)
        del m,original;gc.collect()
    write(OUT/'result.json',dict(classification='Counterexample candidate',rows=rows,new_physical_steps=0,new_linear_solves=0,final_charge_conclusion='unadjudicated'));print(json.dumps(rows),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
