"""Recover six missing photon pairs; reuse the actual final two captures.

Counterexample candidate. No native-fluid reintegration or paired-time claim.
The original whole-equation and exact native-history gates remain mandatory.
"""
from pathlib import Path
from types import FunctionType
import gc,inspect,json,os,resource,sys,time
import numpy as np
import continue_captured_photon_history as recovery
import continue_native_flux_precision as flux

OUT=Path('native-coarse-photon240-work')
FULL=Path('native-true-momentum238-work')
CAPTURE=Path('native-thermal-continuation235-work')
prior=recovery.prior;archive=prior.prior;base=prior.base
read,write,sha,LD=prior.read,prior.write,prior.sha,prior.LD
AMP,A=base.AMP,base.A
CAPS=dict(prepare=300,recover=3600,assemble=300)
saved=lambda n:FULL/f'sweep-1/photons/complete-{n}.npz'
seed_path=lambda n:recovery.OUT/f'accepted-{n}.npz'


def bind(fn,**values):
    return FunctionType(fn.__code__,dict(fn.__globals__,**values),argdefs=fn.__defaults__)


def prepare():
    assert read(FULL/'path-64.json')['passed']
    assert read(FULL/'prefix-result.json')['passed']
    assert read(FULL/'controller-status.json')['state']=='completed'
    assert read(recovery.OUT/'result.json')['same_solution_accepted_history_recovered']
    assert not OUT.exists();OUT.mkdir();files=[]
    for part in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64']:
        (OUT/part).mkdir(parents=True)
    for p in list((FULL/'sweep-0').rglob('*.npz'))+[FULL/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(FULL);os.link(p,dst);files += [p,dst]
    z=np.load(saved(64));old=np.load(archive.saved(64));seed=np.load(seed_path(64))
    assert int(seed['step'])==111 and len(z['actual_step_edges'])==120
    for k in ['joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled']:
        assert np.array_equal(z[k][:222],old[k]),k
    for k in ['moments','collisions','ports','packets']:assert len(seed[k])==222
    checkpoint=np.load(CAPTURE/'last-accepted-64.npz')
    assert len(checkpoint['actual_edges'])==118
    assert np.array_equal(checkpoint['actual_edges'],z['actual_step_edges'][:118])
    captures=[CAPTURE/f'captured-64-{i:03d}.npz' for i in [234,235]]
    captures += [FULL/f'captured-64-{i:03d}.npz' for i in [236,237]]
    for i,p in enumerate(captures,234):
        c=np.load(p)
        assert c['time']==z['joint_stage_times'][i] and c['weight']==z['joint_stage_weights'][i]
        assert np.array_equal(c['collision_rates'],z['joint_collision_rates_scaled'][i])
    files += captures+[saved(64),archive.saved(64),seed_path(64),CAPTURE/'last-accepted-64.npz',
        recovery.OUT/'expanded-recovery-64.py',recovery.OUT/'recovered-64.json',recovery.OUT/'recovered-64.npz',
        FULL/'path-64.json',FULL/'prefix-result.json',FULL/'coarse-receipt.json',CAPTURE/'coarse-receipt.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Supply the complete original119-step coarse photon/material/port history to the same-solution GR consumer, without replaying accepted material evolution.',
        reuse='Reuse223steps1..111 exactly and actual235/238captures of118/119. Recover only112..117 by conditioning the original photon Radau block on the saved accepted gas. Compare the117photon endpoint with its actual saved checkpoint before joining the captured tail.',
        arithmetic='Steps112..114 use the original native owner.115..117 use the original20260-digit B arithmetic, including covariance, RHS and true B defect. Restore each archived conserved coordinate exactly. Require exact original native-rate identity and the unchanged whole coupled equation; no fitted rate or archived-RHS substitution.',
        gates=read(archive.OUT/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=2,virtual_GiB=8,
        forecast='223measured conditional steps were about20..45seconds before later branches. Six missing pairs suggest2..5minutes plus setup; later costs unmeasured. Allow1hour,12linear corrections. No new physical step, fine path, resolution or period.',
        stop='First original conditional/whole-vector/native-identity/endpoint/ledger/angular/radial gate or cap; preserve failed pairs and accepted photon checkpoint. Do not relax gates or replay the material path.',
        decision='Only the resulting same coarse history is assembled.239fine paired-time acceptance remains required before full-period GR admission.236actual GR return and239running sources/plans stay unchanged.',
        new_material_steps=0,missing_photon_steps=6,reused_prefix_steps=111,captured_tail_steps=2,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    from fractions import Fraction as F
    for k in range(3):assert F(3,4)*F(1,3)**k+F(1,4)==F(1,k+1)
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Original Radau quadrature moments0..2 only. No physical or uniform-error certificate.'))
    write(OUT/'reuse-check.json',dict(classification='Counterexample candidate',passed=True,
        original111stage_arrays_exact=True,captured_tail_times_weights_collision_exact=True,new_material_steps=0))


def interval_model(z,n,step):
    m=prior.interval_model(z,n,step)
    with flux.precision.mp.workdps(60):m.precise_tangent=flux.precision.build(flux.joint,flux.owner)
    m.precise_values={}
    return m


def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
    initial=prior.restored_gas(m,z['joint_stage_conserved_scaled'][2*step-1])
    initial[~m.material.active(t)]=0
    v=m.pack(x,initial);src=[];rates=[];states=[];identity=[];maps=[]
    try:
        for j,(xx,g,now,c,s) in enumerate(zip(photons,gas,times,cs,ss)):
            J,b=m.jacobian(now,g);maps.append((J,b));affine=b-(J@g.ravel()).reshape(m.n,4)
            src.append(m.pack(s[0]/(m.scale*AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+affine))
            ph,q,*_=m.collision(c,xx,g,True)
            native=m.native(now,g,details=True)
            if step>=114:native=flux.precise_native(m,now,g,native)
            rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+ph+s[0]/(m.scale*AMP),q+native[0]))
            states.append(m.pack(xx,g))
            identity.append(float(np.max(abs(native[0]*m.units-z['joint_native_rates_scaled'][2*step+j]))))
        rhs=(v+h*(A@np.array(src))).ravel();sol=np.array(states).ravel()
        defect=(np.array(states)-v-h*(A@np.array(rates))).ravel()
        if step>=114:
            cache=m.precise_values
            rhs=flux.precise_rhs(m,t,h,v,maps,gas,rhs);m.precise_values=cache
            defect=flux.precise_defect(m,t,h,v,sol,defect)
        rel=float(np.linalg.norm(defect)/np.linalg.norm(rhs))
        physical=(base.joint.physical_norm(m,defect)/base.joint.scales(m,rhs,sol)).astype(float).tolist()
        row=dict(step=step,relative=rel,physical=physical,native_identity_absolute=identity,
            original202_B_arithmetic=step>=114,passed=rel<1e-12 and max(physical)<1e-13 and max(identity)==0)
        write(OUT/f'clock-64/original-equation-{step}.json',row)
        assert row['passed'],row
        return row
    except BaseException:
        np.savez_compressed(OUT/'rejected-original-64.npz',step=step,time=t,step_size=h,x_initial=x,
            photon_stage_solution=photons,gas=gas,stage_times=times)
        raise


def recover():
    source=(recovery.OUT/'expanded-recovery-64.py').read_text()
    changes=[("target=float(z['t'][-1])","target=float(z['actual_step_edges'][117])"),
        ("expected=z['photon_history_scaled_occupation'][-1];actual=x*m.scale*AMP",
         "expected=np.load(CAPTURE/'last-accepted-64.npz')['x']*m.scale*AMP;actual=x*m.scale*AMP"),
        ('snapshot(n,m,z,step,pairs[-1],collisions,ports,packets,initial_x=x,photon_pair=pairs,gas=gas,t=t,h=h)',
         'snapshot(n,m,z,step,pairs[-1],collisions,ports,packets)')]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    ns=dict(base.run.__globals__,OUT=OUT,CAPTURE=CAPTURE,saved=saved,seed_folder=lambda n:recovery.OUT,
        seed=bind(prior.seed,seed_path=seed_path),seed_path=seed_path,RESUME=True,
        checkpoint=bind(archive.checkpoint,OUT=OUT),snapshot=bind(archive.snapshot,OUT=OUT),
        restored_gas=prior.restored_gas,original_equation=original_equation,
        interval_model=interval_model,interval_starts=prior.interval_starts,gc=gc)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-recovery-64.py').write_text(source)
    ns['run'](64)


def assemble():
    row=read(OUT/'recovered-64.json');assert row['passed'] and row['recovered_steps']==117
    z=np.load(saved(64));d=dict(np.load(OUT/'recovered-64.npz'))
    initial=np.load(recovery.OUT/'recovered-64.npz')
    for k in ['times','weights','photon_moments','collision_rates','angular','radial_ports']:
        assert np.array_equal(d[k][:222],initial[k]),k
    captures=[np.load(folder/f'captured-64-{i:03d}.npz') for folder,ids in [(CAPTURE,[234,235]),(FULL,[236,237])] for i in ids]
    for key,source in [('photon_moments','photon_moments'),('collision_rates','collision_rates'),('radial_ports','radial_ports')]:
        d[key]=np.concatenate([d[key],np.array([p[source] for p in captures])])
    d['angular']=np.concatenate([d['angular'],z['accepted_angular_luminosity'][234:]])
    d['times']=z['joint_stage_times'];d['weights']=z['joint_stage_weights']
    d['endpoint_occupation']=z['photon_history_scaled_occupation'][-1]
    rel=lambda x,y:float(np.sum(abs(x-y))/max(np.sum(abs(y)),LD('1e-290')))
    port=np.sum(d['weights'][:,None,None]*d['radial_ports'],axis=0,dtype=LD)
    radial=float(np.max(abs(port-z['radial_ports'][-1])/np.maximum(abs(z['radial_ports'][-1]),LD('1e-290'))))
    angular=rel(d['angular'],z['accepted_angular_luminosity'])
    q=z['conserved_material_history'][-1]/AMP;k=z['energy_offset_reference'][-1]/AMP
    actual=np.column_stack([q[2]-k,q[3],q[0],q[1]])+z['material_floor_discard_history_scaled'][-1]
    target=np.sum(d['weights'][:,None,None]*(z['joint_native_rates_scaled']+d['collision_rates']),axis=0,dtype=LD)
    ledger=(np.sum(abs(actual-target),axis=0)/np.maximum(np.sum(abs(actual)+abs(target),axis=0),LD('1e-290'))).astype(float).tolist()
    result=dict(classification='Counterexample candidate',passed=radial<1e-12 and angular<1e-12 and max(ledger)<1e-8,
        clock=64,recovered_steps=119,horizon_seconds=float(z['actual_step_edges'][-1]),reused_prefix_steps=111,
        new_conditional_steps=6,actual_captured_tail_steps=2,endpoint117_relative=row['endpoint_relative'],
        final_endpoint_from_same_actual_coupled_checkpoint=True,radial_port_relative=radial,angular_relative=angular,
        material_ledger=ledger,original_equation_audited_steps=list(range(111,117)),
        new_material_steps=0,paired_time_admitted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'assembly-check.json',result);assert result['passed'],result
    for ext in ['npz','json']:os.rename(OUT/f'recovered-64.{ext}',OUT/f'conditional-recovered-64.{ext}')
    np.savez_compressed(OUT/'recovered-64.npz',**d);write(OUT/'recovered-64.json',result)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    os.sched_setaffinity(0,{2});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    base.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
