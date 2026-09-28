"""Counterexample candidate: resolve slow-cell pulse onset in the coupled solve.

Use one fixed, input-defined local bisection. No method/source change or
automatic resolution ladder. All ledgers use the actual substep times.
"""
from pathlib import Path
import gc,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
import apply_native_direct_radau as prior

OUT=Path('native-pressure-front172-work');OLD=prior.OUT
radau=prior.radau;original=prior.original
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=30,check=60,pilot=220,finish=30);TOTAL=420
BASE_INITIALIZE=radau.initialize


def split_flags(m,steps):
    # Fixed across the pair: slow radial transport in the native deep core.
    reference_h=m.t[-1]/64;rate=(-m.A.diagonal()).reshape(m.n,m.q)
    cells=np.flatnonzero(np.max(reference_h*rate[:m.nb],axis=1)<1.)
    arrivals=-np.asarray(m.redshift_driver.xc[cells],float)/original.C
    starts=np.r_[arrivals,arrivals+m.redshift_driver.D]
    width=2*reference_h;h=m.t[-1]/steps
    flags=np.array([any(t< (k+1)*h and t+width>k*h and 0<=t<m.t[-1] for t in starts) for k in range(steps)])
    return flags,cells,starts


def initialize():
    radau.OUT=OUT;BASE_INITIALIZE();model=original.c.Response
    runner=(OUT/'sweep-1/expanded-radau-run.py').read_text()
    def change(a,b):
        nonlocal runner
        assert runner.count(a)==1,(a,runner.count(a));runner=runner.replace(a,b)
    change("h=self.t[-1]/steps;gamma=1/4", "h=self.t[-1]/steps;base_h=h;gamma=1/4;flags,front_cells,front_times=split_flags(self,steps)")
    change("lus=[splu(sparse.eye(self.n*self.q,format='csc')-v*h*self.A) for v in [5/12,1/4]]",
        "lu_by_parts={parts:[splu(sparse.eye(self.n*self.q,format='csc')-v*base_h/parts*self.A) for v in [5/12,1/4]] for parts in [1,2]}")
    change('times=[];begin=0;moment=0.', 'times=[];begin=0;moment=0.;stage_weights=[];actual_edges=[0.]')
    anchor="self.angular=list(z['accepted_angular_luminosity']);"
    change(anchor,anchor+"stage_weights=list(z['accepted_angular_quadrature_weights']);actual_edges=list(z['actual_step_edges']);")
    a=runner.index('        t=k*h\n');b=runner.index('        if (k+1)%',a)
    body=runner[a:b].replace('        t=k*h\n','        t=k*base_h+sub*h\n')
    body=body.replace('self.lift((k+1)*h)','self.lift(t+h)')
    runner=runner[:a]+'''        parts=2 if flags[k] else 1;h=base_h/parts;lus=lu_by_parts[parts]
        for sub in range(parts):
'''+''.join('    '+line+'\n' for line in body.splitlines())+'''            stage_weights.extend(h*RK_B);actual_edges.append(t+h)
'''+runner[b:]
    runner=runner.replace('count*h','count*base_h').replace('begin*h','begin*base_h').replace('record((k+1)*h)','record((k+1)*base_h)')
    change('accepted_angular_luminosity=self.angular)',
        'accepted_angular_luminosity=self.angular,accepted_angular_quadrature_weights=stage_weights,actual_step_edges=actual_edges)')
    change("accepted_angular_quadrature_weights=h*np.tile(RK_B,count)",
        "accepted_angular_quadrature_weights=stage_weights,actual_step_edges=actual_edges,split_macro_steps=flags,front_cells=front_cells,front_times=front_times")
    change("    row['passed']=", "    row.update(actual_completed_steps=len(stage_weights)//2,actual_new_steps=len(stage_weights)//2-sum(1+int(v) for v in flags[:begin]),split_completed_macro_steps=int(sum(flags[:count])),base_steps=steps)\n    row['passed']=")
    namespace=dict(model.run.__globals__,split_flags=split_flags)
    exec(compile(runner,__file__,'exec'),namespace);model.run=namespace['run']
    (OUT/'sweep-1/expanded-front-run.py').write_text(runner)


def toy(n,local):
    T=read(OUT/'inspection.json')['prefix_end_seconds'];a=read(OUT/'inspection.json')['deep_cell_arrivals'][15]['seconds']/T
    x=p=0.;count=0
    for k in range(n):
        parts=2 if local and k/n<a+.5 and (k+1)/n>a else 1;h=1/n/parts
        for j in range(parts):
            t=k/n+j*h;f=np.maximum(t+h*radau.RK_C-a,0.)**3
            stages=x+h*(radau.RK_A@f);p+=h*(radau.RK_B@stages);x=stages[-1];count+=1
    return dict(steps=count,x=float(x),pressure=float(p))


def prepare():
    assert OUT.exists() and not (OUT/'plan.json').exists() and not read(OLD/'result.json')['passed']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    files=[Path(__file__),Path(prior.__file__),Path(radau.__file__),Path(original.__file__),OLD/'result.json',OUT/'inspection-plan.json',OUT/'inspection.json',OUT/'inspection-receipt.json']
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    inspection=Path('pressure-mode172-inspection.py');shutil.copyfile(inspection,OUT/'inspection-producer.py');files.append(OUT/'inspection-producer.py')
    t=sp.symbols('t',positive=True);x=t**4/4;p=t**5/20
    assert sp.diff(x,t)==t**3 and sp.diff(p,t)==x
    rows=[toy(n,local) for local in [False,True] for n in [4,8]]
    errors=[abs(rows[i]['pressure']-rows[i+1]['pressure'])/abs(rows[i+1]['pressure']) for i in [0,2]]
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='Zero-initial xprime=t^3, Pprime=x gives x=t^4/4 and P=t^5/20; leading compact-pulse derivative is cubic.',
        scope='Leading source cascade only; not an error bound or proof of dominance for the full coupled problem.'))
    write(OUT/'onset-toy.json',dict(classification='Conjectural',rows=rows,time_differences=errors,
        scope='Same measured arrival and leading source power, no fitted amplitude. The toy motivates one local bisection, not a physical acceptance verdict.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='4f469fdfe71e6eeb2fee53dd02cd716664284963',
        claim='Resolve the localized pressure-rise time error in the SAME direct-source coupled photon/gas equations using one input-defined local step bisection.',
        evidence='171pressure mismatch is97.93percent concentrated in cell15. Its gas-block eigenvalue*h is about0.32, radial diagonal*h below0.22, and pressure-map change below2e-9. It receives the primary front after58.39percent of the prefix. Slow-photon-to-gas response requires two time integrations of a leading cubic source. The saved toy predicts9.32percent versus1.41percent paired pressure differences before/after local bisection; no guarantee for the full model.',
        rule='For all existing native deep cells with max radial loss diagonal*(T/64)<1, use their fixed primary onset and trailing-front times. Bisect once every base interval overlapping [front,front+2*T/64]. Same physical rule and window on both base64/128 clocks. No fitted cell list, new equation, method, pulse or background. Record actual step counts honestly.',
        decision='Only all six unchanged2percent time gates plus original physical checks permit considering full coupled continuation and final-charge readout. Final-charge stability remains unadjudicated until the same corrected photon/free-material/GR solution and boundary history are completed.',
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12),
        budget=dict(actions=CAPS,inspection_cap_seconds=60,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3),
        forecast='171actual action103.33s, with68.72s of stepping for12steps. Expected18steps here gives about130-165s with setup/IO;220s hard pilot cap. Full-path late costs and counts are not admitted by this forecast. One pair only.',
        representation='Canonical17state outputs remain unchanged. Integrate collisions, gas, impulse, spectral ghosts/work and actual boundary packets at each actual substep. Save explicit angular weights and step edges; consumers must use them, not a uniform-clock formula.',
        stop='Stop on any gate or cap. No second bisection, extra pair, full path or automatic grid/horizon/method expansion.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def check():
    initialize();m=original.c.Response(128);rows=[]
    for n in [64,128]:
        flags,cells,fronts=split_flags(m,n);h=m.t[-1]/n
        edges=np.r_[0.,np.concatenate([k*h+np.arange(1,(2 if flags[k] else 1)+1)*h/(2 if flags[k] else 1) for k in range(n)])]
        assert np.all(np.diff(edges)>0) and abs(edges[-1]-m.t[-1])<1e-18
        weights=(np.diff(edges)[:,None]*radau.RK_B).ravel()
        assert abs(sum(weights)-m.t[-1])<1e-18
        assert all(np.min(abs(edges-t))<1e-18 for t in m.t)
        rows.append(dict(base_steps=n,full_planned_actual_steps=len(edges)-1,prefix_actual_steps=int(sum(1+flags[:n//16])),
            slow_deep_cells=cells.tolist(),split_macro_steps=np.flatnonzero(flags).tolist()))
    assert [r['prefix_actual_steps'] for r in rows]==[6,12],rows
    write(OUT/'grid-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows))
    print(json.dumps(rows),flush=True)


def pilot():
    assert read(OUT/'grid-check.json')['passed'];initialize();folder=OUT/'sweep-1/photons';rows=[];then=time.monotonic()
    for n in [64,128]:
        m=original.c.Response(n);row=m.run(n,f'pilot-{n}',n//16)
        with np.load(folder/f'pilot-{n}.npz') as p:
            edges=p['actual_step_edges'];width=np.diff(edges);expected=(edges[:-1,None]+width[:,None]*radau.RK_C).ravel()
            tt=p['accepted_angular_times'];w=p['accepted_angular_quadrature_weights']
            assert len(tt)==2*row['actual_completed_steps'] and np.max(abs(tt-expected))<1e-18
            assert np.max(abs(w-(width[:,None]*radau.RK_B).ravel()))<1e-18
            flux=p['accepted_angular_luminosity']@(np.arange(1,8,2)/32)
            err=float(abs(w@flux-p['radial_ports'][-1,1,1])/max(w@abs(flux),1e-290));assert err<1e-12
        row['actual_angular_quadrature_relative']=err;row['time_integrator']='RadauIIA2 with one prescribed local bisection'
        write(folder/f'pilot-{n}.json',row);rows.append(row);assert row['passed'],row
        del m;gc.collect()
    data=[np.load(folder/f'pilot-{n}.npz')['moments'][:,[0,1,2,3,5,6]] for n in [64,128]]
    errors=np.max(np.sum(abs(data[0]-data[1]),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(data[1]),axis=2),axis=0),1e-290)
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02),rows=rows,time_comparison=errors.astype(float).tolist(),
        seconds=time.monotonic()-then,full_horizon_completed=False,physical_final_charge_solved=False,full_goal_complete=False)
    write(folder/'pilot.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def finish():
    radau.OUT=OUT;radau.finish();result=read(OUT/'result.json');result.pop('export_repair',None)
    result['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    result['actual_time_grid']=read(OUT/'grid-check.json')['rows'];write(OUT/'result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
