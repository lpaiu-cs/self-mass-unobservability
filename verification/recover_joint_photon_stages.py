"""Recover missing photon stages conditional on the SAME saved material stages.

Counterexample candidate: solve only the linear photon block of the original
Radau equations, then compare its collisions, packets and saved endpoint.
There is no fluid replay, changed time grid or new physical solution.
"""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator,splu,gmres
import complete_full_incident_horizon as base

OUT=Path('native-stage-recovery204-work');INPUT=base.OUT
read,write,sha=base.read,base.write,base.sha
LD,AMP=base.LD,base.AMP;joint=base.owner.joint
A,B,C=joint.A,joint.B,joint.C
CAPS=dict(prepare=120,coarse=1200,coarse_retry=1200,fine=1800,audit=60)


def saved(n):return INPUT/f'sweep-1/photons/interval-02-{n}.npz'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(INPUT/'comparison-02.json')['passed']
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:dst=OUT/p.relative_to(INPUT);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [saved(n) for n in [64,128]]+[INPUT/'comparison-02.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='97935eaf2',
        claim='Recover the accepted joint solution photon stages missing from its archive, to supply time-resolved SAME-solution GR sources without repeating the stiff native-fluid solve.',
        evidence='The archive retains every material Radau state/rate and actual angular packet but only canonical/pilot photon snapshots. The current GR consumer uses17snapshots, so it does not control source interpolation between them.',
        method='Condition the original two-stage linear photon block on the stored actual four-component gas stages. Reuse original transport, collision, incident source, exact step edges and weights. Start at the saved zero photon state; compare the first saved nonzero endpoint, every collision moment and angular packet. No native-fluid Jacobian or Newton step is recomputed.',
        scope='First saved nonzero pilot endpoint only, using existing64/128clocks. No extra clock, new physical evolution or final charge. Both paths must preserve the original joint-history identities before any longer source recovery is admitted.',
        gates=dict(linear=1e-14,physical_linear=1e-13,endpoint=1e-12,collision=1e-12,packet=1e-12),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='The existing model constructor costs about15seconds. Conditional photon-only stage solves are unmeasured and omit the expensive native Jacobian/Newton system. Allow20/30minutes for the bounded existing prefix; measure coarse first, then use its observed cost. No full-period recovery admitted by this plan.',
        stop='Any original equation, stored collision/packet/endpoint identity,12linear refinements or wall cap. Preserve failures, do not widen gates or replay the whole accepted fluid solution.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    x,g,L,K,f,h,a=sp.symbols('x g L K f h a',commutative=False)
    assert sp.expand(x-h*a*(L*x+K*g+f)-(x-h*a*L*x-h*a*(K*g+f)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Conditioning a linear photon block on saved gas moves its known gas term to the RHS; uniqueness and archival accuracy remain numerical checks.'))


def run(n):
    if n==128:assert read(OUT/'recovered-64.json')['passed']
    initialize=FunctionType(base.prior.initialize.__code__,dict(base.prior.initialize.__globals__,OUT=OUT));initialize()
    m=base.owner.Model(n);z=dict(np.load(saved(n)))
    target=float(z['t'][1]);edges=z['actual_step_edges'];stop=int(np.argmin(abs(edges-target)))
    assert stop>0 and abs(edges[stop]-target)<1e-18
    x=np.asarray(z['photon_history_scaled_occupation'][0]/(m.scale*AMP),LD);assert not np.any(x)
    zero=np.zeros((m.n,4),LD);zeroJ=sparse.csr_matrix((4*m.n,4*m.n));logs=[];moments=[];collisions=[];packets=[];timings=[]
    def stream(v):return (m.A@v.reshape(m.n*m.q,m.nf)).reshape(v.shape)
    def relative(v,reference):return float(np.sum(abs(v-reference),dtype=LD)/max(np.sum(abs(reference),dtype=LD),LD('1e-290')))
    for step in range(stop):
        started=time.monotonic();t,h=edges[step],edges[step+1]-edges[step];times=t+C*h
        ids=np.arange(2*step,2*step+2);assert np.max(abs(z['joint_stage_times'][ids]-times))<1e-18
        gas=[]
        for q in z['joint_stage_conserved_scaled'][ids]:gas.append(np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]))
        cs=[m.local(now) for now in times];ss=[m.source(now) for now in times]
        force=np.array([m.collision(c,np.zeros_like(x),g,True)[0]+s[0]/(m.scale*AMP) for c,g,s in zip(cs,gas,ss)])
        rhs=(x[None]+h*np.einsum('ij,jnqf->inqf',A,force)).ravel()
        shape=(2,*x.shape)
        def mat(value):
            v=value.reshape(shape);rates=np.array([stream(vv)+m.collision(c,vv,zero)[0] for c,vv in zip(cs,v)])
            return (v-h*np.einsum('ij,jnqf->inqf',A,rates)).ravel()
        inv=[];lus=[]
        for j,c in enumerate(cs):
            # Reuse the existing photon rank-two angular inverse, with gas
            # feedback disabled only in this photon-only preconditioner.
            bare=dict(c,**{k:np.zeros_like(c[k]) for k in ['B','Bb','Be']})
            inv.append(m.inverse_pair(bare,h*A[j,j],zeroJ))
            lus.append(splu(sparse.eye(m.n*m.q,format='csc')-h*A[j,j]*m.A))
        def pre(value):
            rows=[]
            for j,v in enumerate(value.reshape(shape)):
                streamed=lus[j].solve(np.asarray(v,float).reshape(m.n*m.q,m.nf)).reshape(x.shape)
                rows.append(m.unpack(inv[j](streamed,np.zeros((m.n,4))))[0])
            return np.asarray(rows,float).ravel()
        op=LinearOperator((len(rhs),)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
        sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]
        for attempt in range(12):
            residual=rhs-op.matvec(sol);err=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),LD('1e-290')))
            r=residual.reshape(shape);s=sol.reshape(shape);b=rhs.reshape(shape)
            physical=[float(np.sum(abs(r)*w)/max(np.sum(abs(b)*w),np.sum(abs(s)*w),1.)) for w in [m.Nweight,m.Eweight]]
            if err<1e-14 and max(physical)<1e-13:break
            assert attempt<11,('Conditional photon residual',err,physical)
            history=[];start=time.monotonic()
            delta,info=gmres(op,np.asarray(residual,float),M=P,rtol=1e-12,atol=0.,restart=40,maxiter=10,callback=history.append,callback_type='pr_norm')
            sol+=np.asarray(delta,LD);calls.append(dict(info=int(info),iterations=len(history),seconds=time.monotonic()-start))
            write(OUT/f'linear-{n}.json',dict(step=step,calls=calls))
        pairs=sol.reshape(shape);errors=[]
        for j,(xx,g,c) in enumerate(zip(pairs,gas,cs)):
            _,q,*_=m.collision(c,xx,g,True);actual=q*m.units;stored=z['joint_collision_rates_scaled'][ids[j]]
            errors.append([relative(actual[:,k],stored[:,k]) for k in range(4)])
            collisions.append(actual)
            moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*AMP)
            m.boundary_ports(float(times[j]),xx);packets.append(m.angular[-1])
        packet_error=relative(np.array(packets[-2:]),z['accepted_angular_luminosity'][ids])
        row=dict(step=step,linear_relative=err,physical_relative=physical,collision_relative=errors,packet_relative=packet_error,calls=calls)
        logs.append(row);write(OUT/f'progress-{n}.json',dict(completed=step+1,total=stop,last=row))
        assert np.max(errors)<1e-12 and packet_error<1e-12,row
        x=pairs[-1].copy();timings.append(time.monotonic()-started)
    expected=z['photon_history_scaled_occupation'][1];actual=x*m.scale*AMP
    endpoint=relative(actual,expected)
    result=dict(classification='Counterexample candidate',passed=endpoint<1e-12,clock=n,recovered_steps=stop,
        horizon_seconds=target,endpoint_relative=endpoint,rows=logs,step_seconds=timings,
        same_material_history_unchanged=True,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/f'recovered-{n}.npz',times=z['joint_stage_times'][:2*stop],weights=z['joint_stage_weights'][:2*stop],
        photon_moments=np.array(moments),collision_rates=np.array(collisions),angular=np.array(packets),endpoint_occupation=actual)
    write(OUT/f'recovered-{n}.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);assert result['passed'],result


def audit():
    rows=[read(OUT/f'recovered-{n}.json') for n in [64,128]];assert all(r['passed'] for r in rows)
    assert rows[0]['horizon_seconds']==rows[1]['horizon_seconds']
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,full_period_recovery=False,
        GR_source_interpolation_resolved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    started=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                actual=OUT/'prepared-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(actual)==h,p
            assert 'Cannot cast array data' in read(OUT/'coarse-receipt.json')['error']
            assert not (OUT/'progress-64.json').exists() or action!='coarse_retry'
            write(OUT/'dispatch-amendment.json',dict(classification='Counterexample candidate',
                change='Cast only approximate preconditioner input/output to the binary64 sparse-LU contract. The true photon equation and extended residual retain their original precision. Initial dispatch failed before any recovered stage.',
                original_failure_seconds=read(OUT/'coarse-receipt.json')['seconds'],
                retry_cap_seconds=CAPS['coarse_retry'],prepared_source_sha256=sha(OUT/'prepared-producer.py'),executed_source_sha256=sha(__file__)))
        if action in ['coarse','coarse_retry','fine']:run(128 if action=='fine' else 64)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-started,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
