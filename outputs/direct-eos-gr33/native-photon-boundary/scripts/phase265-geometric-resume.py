"""Finish the 575-time geometric lapse after stopping four slow workers for cost.

Workers 6-9 ran the same module on the same inputs at about 47 s per launch cell
against about 15 s for the other ten, which projects past their registered 3 h
cap. They are stopped between items and their completed items are kept. Their
remaining items run through a shared queue: first on idle CPUs, then on the
CPU of each original worker as it finishes. Packets, quadrature, projection and
gates are the module's own; only the scheduling changes. The module's collect
cannot run (four workers have no receipt), so the same checks are applied here
to the merged items, and the stop evidence is kept.
"""
from pathlib import Path
import fcntl,hashlib,json,os,signal,subprocess,sys,time
import numpy as np

OUT=Path('native-geometric-clock265-work');STATUS=Path('native-geometric-clock265-status.json')
QUEUE=OUT/'queue.json';LOCK=OUT/'queue.lock';STOPPED=[6,7,8,9]
module=Path('verification/propagate_geometric_boundary_clock.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='/home/lpaiu/work/nutimo_pilot/request13_deps:verification')


def alive(pid,ticks):
    try:
        stat=Path(f'/proc/{pid}/stat').read_text();fields=stat[stat.rindex(')')+2:].split()
        return fields[19]==ticks and fields[0]!='Z'
    except (FileNotFoundError,ProcessLookupError):return False


def stop():
    """Stop the controller, then the four slow workers, keeping identity evidence."""
    assert not (OUT/'cost-stop.json').exists()
    old=read(STATUS);assert old['state']=='running' and old['action']=='work'
    boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    procs={}
    for d in Path('/proc').iterdir():
        if d.name.isdigit():
            try:cmd=(d/'cmdline').read_bytes().replace(b'\0',b' ').decode()
            except (FileNotFoundError,ProcessLookupError,PermissionError):continue
            if '.phase265-geometric-controller.py' in cmd or str(module) in cmd:
                procs[int(d.name)]=dict(cmd=cmd.strip(),ticks=(d/'stat').read_text().split()[21])
    controller=[p for p,v in procs.items() if '.phase265-geometric-controller.py' in v['cmd']];assert len(controller)==1
    workers={f'work-{k}':old['workers'][f'work-{k}'] for k in range(14)}
    for k in range(14):
        pid=workers[f'work-{k}'];assert pid in procs and procs[pid]['cmd'].endswith(f'work-{k}'),k
    record=dict(classification='Counterexample candidate',boot_id=boot,controller_pid=controller[0],controller_ticks=procs[controller[0]]['ticks'],
        workers={a:dict(pid=p,ticks=procs[p]['ticks']) for a,p in workers.items()},stopped=[f'work-{k}' for k in STOPPED],original_status=old,
        progress={f'work-{k}':read(OUT/f'work-{k}-progress.json') for k in range(14)},
        reason='Same module and inputs: workers6-9 measured about47s per launch cell, the other ten about15s. At that rate their remaining items take about4.3h, beyond the registered3h cap. Intentional cost stop between items, not a physical or numerical gate failure.',
        time_unix=time.time())
    write(OUT/'cost-stop-intent.json',record)
    os.kill(controller[0],signal.SIGTERM)
    for _ in range(100):
        if not alive(controller[0],procs[controller[0]]['ticks']):break
        time.sleep(.1)
    assert not alive(controller[0],procs[controller[0]]['ticks'])
    for k in STOPPED:
        pid=workers[f'work-{k}'];ticks=procs[pid]['ticks']
        # Stop between items: wait until the last progress write is older than1s and the item count is stable.
        before=read(OUT/f'work-{k}-progress.json')['completed']
        os.kill(pid,signal.SIGSTOP);time.sleep(1.5)
        after=read(OUT/f'work-{k}-progress.json')['completed']
        os.kill(pid,signal.SIGKILL)
        for _ in range(100):
            if not alive(pid,ticks):break
            time.sleep(.1)
        record['progress'][f'work-{k}']['completed_at_stop']=after;record['progress'][f'work-{k}']['completed_before_stop']=before
    record['stopped_unix']=time.time();write(OUT/'cost-stop.json',record)


def completed(k):
    """Saved items of worker k that also have their progress row; anything unreadable is recomputed."""
    try:
        z=np.load(OUT/f'work-{k}.npz');index=[int(v) for v in z['index']];values=z['values'];launch=z['launch']
        rows=read(OUT/f'work-{k}-progress.json')['rows']
    except Exception:return [],np.zeros((0,5)),np.zeros(0)
    n=min(len(index),len(rows))
    if [r['index'] for r in rows[:n]]!=index[:n]:return [],np.zeros((0,5)),np.zeros(0)
    return index[:n],values[:n],launch[:n]


def prepare_queue():
    assert read(OUT/'cost-stop.json') and not QUEUE.exists()
    plan=read(OUT/'execution-plan.json');items=[]
    for k in STOPPED:
        done=set(completed(k)[0]);items+=[v for v in plan['assignments'][k] if v['index'] not in done]
    items=sorted(items,key=lambda v:-v['cells'])
    write(QUEUE,dict(items=items,claimed={}))


def claim(tag):
    with LOCK.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        q=read(QUEUE);pending=[v for v in q['items'] if str(v['index']) not in q['claimed']]
        if not pending:return None
        item=pending[0];q['claimed'][str(item['index'])]=tag;write(QUEUE,q);return item


def worker(tag):
    sys.path.insert(0,'verification');import resource
    import propagate_geometric_boundary_clock as g
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);g.p.incident.native.deadline(g.CAPS['work'])
    for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
    plan=read(OUT/'execution-plan.json')
    for q,h in plan['bindings'].items():assert sha(q)==h,q
    start=time.monotonic();error=None;done=0
    try:
        g.p.initialize();a,order,geo=g.FINE;m=g.Photons(a,geo);assert np.array_equal(m.clock,np.array(plan['clock']))
        while (item:=claim(tag)) is not None:
            value,row=g.evaluate(m,item['time'],order)
            np.savez_compressed(OUT/f"resume-{item['index']}-tmp.npz",value=value,launch=np.array(row['physical_launch_energy_increment_erg']))
            os.replace(OUT/f"resume-{item['index']}-tmp.npz",OUT/f"resume-{item['index']}.npz")
            write(OUT/f"resume-{item['index']}.json",dict(item,worker=tag,packets=row['packets'],seconds=row['seconds'],
                energy_identity_relative=row['energy_identity_relative'],angular_invariant_relative=row['angular_invariant_relative']))
            done+=1
    except BaseException as exc:error=repr(exc);raise
    finally:write(OUT/f'queue-{tag}-receipt.json',dict(seconds=time.monotonic()-start,items=done,error=error,source_sha256=sha(__file__),
        module_sha256=sha(module),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


def collect():
    sys.path.insert(0,'verification');import propagate_geometric_boundary_clock as g
    plan=read(OUT/'execution-plan.json');stop=read(OUT/'cost-stop.json');t=g.metric_clock();n=len(t);LD=g.LD
    values=np.zeros((n,5),LD);launch=np.zeros(n);owner=[None]*n;rows=[]
    evaluation=np.array(t,float);knot=np.full(n,-1)
    for group in plan['assignments']:
        for item in group:
            evaluation[item['index']]=item['time']
            if item['knot'] is not None:knot[item['index']]=item['knot']
    def put(i,value,lval,who):
        assert owner[i] is None,(i,owner[i],who);owner[i]=who;values[i]=value;launch[i]=lval
    for k in range(14):
        index,vals,lau=completed(k)
        if k in STOPPED:assert not (OUT/f'work-{k}-receipt.json').exists()
        else:
            progress=read(OUT/f'work-{k}-progress.json');assert read(OUT/f'work-{k}-receipt.json')['error'] is None
            assert progress['completed']==progress['total']==len(plan['assignments'][k])==len(index)
        for i,v,l in zip(index,vals,lau):put(i,v,l,f'work-{k}')
        rows+=read(OUT/f'work-{k}-progress.json')['rows'][:len(index)]
    q=read(QUEUE)
    for item in q['items']:
        i=item['index'];z=np.load(OUT/f'resume-{i}.npz');r=read(OUT/f'resume-{i}.json')
        assert r['index']==i and r['time']==item['time'];put(i,z['value'],float(z['launch']),r['worker']);rows.append(r)
    for f in OUT.glob('queue-*-receipt.json'):assert read(f)['error'] is None,f
    assert all(owner[i] is not None for i in range(1,n)) and owner[0] is None
    saved=np.load(g.SOURCE/'photon-boundary/source.npz');exact=[]
    for j in range(1,17):
        i=int(np.flatnonzero(knot==j)[0]);v=values[i]
        exact.append(bool(v[0]==saved['photon_J_source_cm'][j] and v[1]==saved['photon_lapse_source'][j] and np.array_equal(v[2:],saved['lapse_energy_radius_angle_parts'][j-1])))
    pilot=np.load(OUT/'pilot.npz');repeat=bool(all(np.array_equal(values[i],v) for i,v in zip(pilot['index'],pilot['values'])))
    identity=max(r['energy_identity_relative'] for r in rows);invariant=max(r['angular_invariant_relative'] for r in rows)
    lapse=np.asarray(values[:,1],float);applied=np.asarray(np.load(g.ACTUAL/'metric/metric-128-g8.npz')['delta_nu_faces'][:,-1],float)
    parts=np.max(abs(np.asarray(values[:,2:],float)),axis=0)
    passed=all(exact) and repeat and identity<.002 and invariant<1e-10
    np.savez_compressed(OUT/'boundary-575.npz',t=t,evaluation_t=evaluation,knot=knot,
        photon_geometric_mass_cm=values[:,0],photon_geometric_lapse=values[:,1],lapse_energy_radius_angle_parts=values[:,2:],
        physical_launch_energy_increment_erg=launch)
    counts={}
    for who in owner[1:]:counts[who]=counts.get(who,0)+1
    result=dict(classification='Counterexample candidate',passed=passed,knot_values_exact=exact,pilot_repeat_exact=repeat,
        output_times=n,new_times=int(np.count_nonzero(knot[1:]<0)),maximum_energy_identity_relative=identity,
        maximum_angular_invariant_relative=invariant,
        maximum_geometric_lapse_over_applied_outer_lapse=float(np.max(abs(lapse))/np.max(abs(applied))),
        maximum_lapse_energy_radius_angle_parts=parts.tolist(),
        terminal_photon_geometric_mass_cm=float(values[-1,0]),terminal_photon_geometric_lapse=float(values[-1,1]),
        inherited_phase261_quadrature_controls=read(g.SOURCE/'photon-boundary/result.json')['controls'],
        scheduling=dict(items_by_process=counts,cost_stop=str(OUT/'cost-stop.json'),queue=str(QUEUE),collector=str(Path(__file__)),
            collector_sha256=sha(__file__),module_collect_not_run='Workers6-9have no receipt; the module collect requires all fourteen.'),
        boundary_applied_to_matter=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed,result


def control():
    stop();prepare_queue();start=time.monotonic();children={};source=sha(__file__)
    old=read(OUT/'cost-stop.json')['workers']
    running={a:v for a,v in old.items() if int(a.split('-')[1]) not in STOPPED}
    cpus={f'work-{k}':cpu for k,cpu in zip(range(14),read(OUT/'execution-plan.json')['CPU_affinities'])}
    free=[0,15]+[cpus[f'work-{k}'] for k in STOPPED];serial=[0]
    def spawn(cpu):
        assert sha(__file__)==source;tag=f'q{serial[0]}-cpu{cpu}';serial[0]+=1
        with open(f'.phase265-geometric-{tag}.stdout.log','w') as o,open(f'.phase265-geometric-{tag}.stderr.log','w') as e:
            children[tag]=(subprocess.Popen(['taskset','-c',str(cpu),sys.executable,__file__,'worker',tag],stdout=o,stderr=e,env=env),cpu)
    def pending():
        q=read(QUEUE);return len(q['items'])-len(q['claimed'])
    for cpu in free:spawn(cpu)
    status=lambda state,**kw:write(STATUS,dict(state=state,action='resume',seconds=time.monotonic()-start,
        original_running=list(running),queue_workers={k:v[1] for k,v in children.items()},pending_items=pending(),**kw))
    status('running')
    while True:
        for a,v in list(running.items()):
            if (OUT/f'{a}-receipt.json').exists() or not alive(v['pid'],v['ticks']):
                del running[a]
                if pending()>0:spawn(cpus[a])
        for tag,(child,cpu) in list(children.items()):
            code=child.poll()
            if code is not None and code!=0:
                for c,_ in children.values():
                    if c.poll() is None:c.terminate()
                status('failed',error=f'{tag} returned {code}');raise SystemExit(1)
        if not running and all(c.poll() is not None for c,_ in children.values()):break
        status('running');time.sleep(10)
    assert pending()==0
    status('running',phase='collect')
    code=subprocess.run([sys.executable,__file__,'collect'],env=env).returncode
    if code!=0:status('failed',error='collect');raise SystemExit(code)
    status('completed',result=read(OUT/'result.json'))


if __name__=='__main__':
    action=sys.argv[1]
    if action=='control':control()
    elif action=='worker':worker(sys.argv[2])
    elif action=='collect':collect()
    else:raise SystemExit(action)
