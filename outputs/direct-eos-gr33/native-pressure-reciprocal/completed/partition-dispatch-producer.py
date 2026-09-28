"""Apply the accepted same-solution material motion back to actual photons.

Counterexample candidate. One new reciprocal sweep, unchanged Radau equation,
original gates and fixed local clock. Previous complete paths are immutable.
"""
from pathlib import Path
import gc,inspect,json,os,resource,shutil,sys,time
import numpy as np
import continue_native_material_flux_tangent as material
import continue_native_pressure_response as photon

engine=material.prior;run=engine.run;front=run.front;original=run.original
OUT=Path('native-pressure-reciprocal175-work');BEFORE=run.OUT
read,write,sha=run.read,run.write,run.sha;LD=np.longdouble;AMP=run.AMP
CAPS=dict(prepare=30,check=60,pilot=220,run=4800,material_pilot=60,material_production=450,residual=30,block=300)
TOTAL=sum(CAPS.values())
CAPS['recheck']=50  # Remaining part of the original60s check, not a new budget.


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def initialize():
    run.OUT=OUT;engine.face.BASE_PATHS=paths;engine.paths=paths
    engine.initialize()
    Parent=run.c.Response
    class Response(Parent):
        def __init__(self,n):
            super().__init__(n)
            # Retain the small noncollisional remainder before binary64 maps.
            with np.load(paths(0)[1]/f'steps-{n}-reference-128.npz') as d,np.load(paths(0)[0]/f'steps-{n}-reference-128.npz') as p:
                ids=[np.argmin(abs(d['t']-t)) for t in self.t]
                motion=d['history_scaled'][ids].astype(LD)
                offset=(self.a.astype(LD)-self.model.m.a0)*self.model.cx*LD(original.C)**2*motion[:,0]
                motion[:,2]-=offset
                collision=p['collision_transfer'].astype(LD).transpose(0,2,1)/LD(AMP)
            self.motion=np.asarray(motion,float);self.energy_offset=np.asarray(offset,float)
            self.mechanical=np.asarray(motion[:,[2,3]]-collision,float)
    run.c.Response=Response


def status(v):
    p=OUT/'status.json';temp=p.with_suffix('.tmp');write(temp,v);temp.replace(p)


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[Path(__file__),Path(material.__file__),Path(engine.__file__),
        Path(engine.face.__file__),Path(run.__file__),Path(photon.__file__),Path(front.__file__),
        Path(front.radau.__file__),BEFORE/'result.json',BEFORE/'full-material-audit.json',
        BEFORE/'sweep-1/material-flux-tangent/production.json']
    assert read(files[-1])['passed'] and read(BEFORE/'full-material-audit.json')['passed']
    for s in [0,1]:
        for p in paths(s):p.mkdir(parents=True)
    for name in ['normalization.json','photon-conservation-plan.json']:
        shutil.copyfile(BEFORE/name,OUT/name);files.append(BEFORE/name)
    for n in [64,128]:
        for src,dest in [(BEFORE/'sweep-1/photons',paths(0)[0]),(BEFORE/'sweep-1/material-flux-tangent',paths(0)[1])]:
            p=src/f'steps-{n}-reference-128.npz';q=dest/p.name
            shutil.copyfile(p,q);assert sha(p)==sha(q);files.append(p)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='184582b28',
        previous_turn='Progress: same173photon history was applied to actual free matter; atmospheric derivative bottleneck repaired and full64/128material paths passed original gates.',
        claim='Return that actual B,S,inventory displacement and noncollisional Etilde/H to the simultaneous corrected photon/thermal/H equations. Return the resulting collisions to the same free material operator and test the reciprocal finite block before GR.',
        decision='Whether the corrected same-solution input reaches the original0.2percent block criterion. Photon/material passes alone do not decide final charge. A failed reciprocal block stops this registered sweep; no automatic third sweep.',
        equation='Existing operator L on the same retained background; affine drive includes actual B,S,xi(B),dM/dt with M=(Etilde,H)-collision transfer. Etilde=Eref-(a-a_surface)*cx*c^2*B. No manual work addition to charge or mass.',
        reuse='Copy completed173photons and174material byte-for-byte as the preceding iterate. No background/EOS/grid/waveform replay or changed physical amplitude. New paths are necessary because their actual free-material input changes.',
        time='Same64/128base clocks, input-defined single local bisection, direct physical Radau stages and actual quadrature weights. Equal-horizon4/8macro-step prefixes; resume only accepted prefixes. Compare every canonical interval.',
        gates=dict(time=.02,block=.002,paired=.002,directional=.002,conservation=1e-8,branch=.01,stage=1e-12,physical_stage_moment=1e-13,port=1e-12,mapping=1e-10),
        budgets=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3,max_new_reciprocal_sweeps=1,
        forecast='Photon173completed in2630.61action seconds. Measure the new first interval against its same-clock accepted predecessor, then scale the measured old interval-by-interval continuation costs by the largest observed cost ratio and1.75margin. Must fit4800s, rechecked after each completed pair. Changed late Krylov costs remain assumptions. New material uses physical SSP-step cost and the450s original cap.',
        stop='Original gate failure, nondecreasing reciprocal residual, memory or action/total cap. No automatic extra grid, clock, pulse, bisection, sweep, or weakened criterion.',
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    status(dict(state='prepared',final_charge_conclusion='unadjudicated'))


def check(retry=False):
    if retry:
        failed=read(OUT/'check-receipt.json');assert failed['error'] and failed['seconds']+CAPS['recheck']<60
        write(OUT/'partition-precision-plan.json',dict(classification='Counterexample candidate',
            failure=failed['error'],repair='Compute Etilde and subtract actual collision transfer in longdouble before casting the small M to binary64. Apply in the common active Response constructor for both clocks and later block checks; same physical arrays and equation.',
            original_check_cap=60,remaining_check_cap=CAPS['recheck'],no_evolution_repeated=True,
            source_sha256=sha(__file__),original_source_sha256=sha(OUT/'initial-producer.py')))
    initialize();m=run.c.Response(128);rows=[]
    assert np.any(m.motion[:,0]) and np.any(m.motion[:,1]) and np.any(m.mechanical)
    data=np.load(paths(0)[1]/'steps-128-reference-128.npz');p=np.load(paths(0)[0]/'steps-128-reference-128.npz')
    ids=[np.argmin(abs(data['t']-t)) for t in m.t];z=data['history_scaled'][ids].astype(LD)
    offset=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*LD(original.C)**2*z[:,0];z[:,2]-=offset
    assert np.array_equal(m.motion,np.asarray(z,float))
    partition=float(np.max(abs(m.mechanical-(z[:,[2,3]]-p['collision_transfer'].transpose(0,2,1)/AMP)))/max(np.max(abs(m.mechanical)),1e-290))
    assert partition<1e-8,partition
    current={key:getattr(m,key) for key in ['motion','mechanical','xi','energy_offset']}
    for t in [m.t[-1]/96,m.t[-1]/16]:
        actual=m.local(t)
        try:
            for key,v in current.items():setattr(m,key,np.zeros_like(v))
            zero=m.local(t)
        finally:
            for key,v in current.items():setattr(m,key,v)
        for key in ['loss','sc','B','Bb','Be','pressure_map']:assert np.array_equal(actual[key],zero[key]),key
        for key in ['S','esc']:
            diff=actual[key]-zero[key];assert not (np.any(diff.data) if hasattr(diff,'nnz') else np.any(diff))
        assert any(np.any(actual[key]) for key in ['q','qb','qe','mechanical'])
        s,l,e=m.source(t);expected,ports,error=original.correction(m,t)
        assert np.array_equal(s,expected) and np.array_equal(l,ports) and max(e,error)<1e-12
        rows.append(dict(t=float(t),actual_material_input=True,operator_unchanged=True,external_source_unchanged=True))
    write(OUT/'input-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        partition_relative=partition,mapping=m.mapping_error,velocity_jet=m.velocity_jet_error,
        final_charge_conclusion='unadjudicated'))
    import sympy as sp
    E,B,a,a0,c,C,h=sp.symbols('E B a a0 c C h')
    Et=E-(a-a0)*c*B
    assert sp.expand(Et+(a-a0)*c*B-E)==0
    assert sp.expand((Et+h)-(C+h)-(Et-C))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Conserved energy-coordinate inversion and cancellation of a shared collision transfer from M=Etilde-C. No full physical or nonlinear error bound.'))


def historical(k,n):
    if k==1:return read(front.OUT/f'sweep-1/photons/pilot-{n}.json')
    return next(v for v in read(photon.OUT/f'comparison-{k:02d}.json')['rows'] if v['steps']==n)


def remaining_cost(first,ratio,overhead):
    return 1.75*sum(max(v.get('continuation_wall_seconds',v['seconds']),ratio*v['stepping_seconds']+overhead)
                    for k in range(first,17) for n in [64,128] for v in [historical(k,n)])


def pilot():
    assert read(OUT/'input-check.json')['passed'];initialize();rows=[];hist=[];ratio=1.;overhead=0.;start=time.monotonic()
    for n in [64,128]:
        mark=time.monotonic();m=run.c.Response(n);r=m.run(n,f'pilot-{n}',n//16);wall=time.monotonic()-mark
        err,moments,count=photon.verify_packet(paths(1)[0]/f'pilot-{n}.npz',paths(1)[0]/f'pilot-{n}.npz')
        assert r['passed'] and count==({64:6,128:12}[n]),r
        r.update(actual_angular_quadrature_relative=err,worker_wall_seconds=wall);rows.append(r);hist.append(moments)
        ratio=max(ratio,r['stepping_seconds']/historical(1,n)['stepping_seconds']);overhead=max(overhead,wall-r['stepping_seconds'])
        write(paths(1)[0]/f'pilot-{n}.json',r);del m;gc.collect()
    errors=run.c.relative(*hist);upper=remaining_cost(2,ratio,overhead)
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,rows=rows,time_comparison=errors,
        upper_remaining_seconds=upper,eligible=max(errors)<.02 and upper<CAPS['run'],cost_ratio=ratio,overhead=overhead,
        seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated')
    write(paths(1)[0]/'pilot.json',result);status(dict(state='pilot_completed' if result['eligible'] else 'pilot_rejected',result=result))
    print(json.dumps(result),flush=True);assert result['eligible'],result


def evolve():
    pilot_result=read(paths(1)[0]/'pilot.json');assert pilot_result['eligible'];initialize()
    lock=OUT/'run.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
    identity=dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],started_unix=time.time(),source_sha256=sha(__file__))
    labels={n:f'pilot-{n}' for n in [64,128]};completed={64:6,128:12};ratio=pilot_result['cost_ratio'];overhead=pilot_result['overhead'];start=time.monotonic()
    try:
        for k in range(2,17):
            rows=[];hist=[]
            for n in [64,128]:
                status(dict(identity,state='running',active_base_clock=n,target_interval=k,completed_intervals=k-1,total_intervals=16,actual_completed_steps=completed,elapsed_seconds=time.monotonic()-start,cap_seconds=CAPS['run'],final_charge_conclusion='unadjudicated'))
                mark=time.monotonic();m=run.c.Response(n);label=f'steps-{n}-reference-128' if k==16 else f'interval-{k:02d}-{n}'
                row=m.run(n,label,k*(n//16),restart=labels[n]);wall=time.monotonic()-mark;assert row['passed'],row
                err,moments,count=photon.verify_packet(paths(1)[0]/f'{label}.npz',paths(1)[0]/f'pilot-{n}.npz')
                row.update(actual_angular_quadrature_relative=err,accepted_prefix_preserved=True,continuation_wall_seconds=wall)
                write(paths(1)[0]/f'{label}.json',row);rows.append(row);hist.append(moments);completed[n]=count;labels[n]=label
                ratio=max(ratio,row['stepping_seconds']/historical(k,n)['stepping_seconds']);overhead=max(overhead,wall-row['stepping_seconds'])
                del m;gc.collect()
            errors=run.c.relative(*hist);elapsed=time.monotonic()-start;upper=remaining_cost(k+1,ratio,overhead)
            result=dict(classification='Counterexample candidate',passed=max(errors)<.02,interval=k,time_comparison=errors,rows=rows,
                        elapsed_seconds=elapsed,upper_remaining_seconds=upper,actual_completed_steps=completed.copy())
            write(OUT/f'comparison-{k:02d}.json',result);print(json.dumps({a:b for a,b in result.items() if a!='rows'}),flush=True)
            assert result['passed'],result
            assert elapsed+upper<CAPS['run'],('Cost admission',elapsed,upper)
        result.update(full_photon_horizon_completed=True,actual_previous_free_material_input_applied=True,
            new_free_material_return_completed=False,reciprocal_block_accepted=False,GR_return_completed=False,
            final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
        write(paths(1)[0]/'result.json',result);write(OUT/'photon-result.json',result);status(dict(identity,state='completed',result=result))
    except BaseException as exc:
        failure=dict(error=repr(exc),last_labels=labels,actual_completed_steps=completed,elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated')
        write(OUT/'failure.json',failure);status(dict(identity,state='failed',result=failure));raise


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                frozen=OUT/'initial-producer.py' if Path(p)==Path(__file__) else p
                assert sha(frozen)==h,p
            if action!='recheck':assert sha(__file__)==read(OUT/'partition-precision-plan.json')['source_sha256']
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        if action=='run':evolve()
        elif action=='recheck':check(True)
        elif action in ['prepare','check','pilot']:globals()[action]()
        else:raise AssertionError('Downstream actions require completed photons and a measured material admission')
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
