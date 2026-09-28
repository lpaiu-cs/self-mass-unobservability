"""Finish the 575-time geometric lapse through a shared queue after all workers ended.

The first attempt (.phase265-geometric-resume.py) stopped the original
controller in order to stop four slow workers (about 47 s per launch cell
against 15 s). Stopping that WSL-launched controller also ended all fourteen
workers of its session, so none wrote a receipt. Every item saved before that
(npz values plus its progress row) is kept. The remaining items run through a
shared queue, so that slower CPUs take fewer items. The module's packets,
quadrature, projection and gates are unchanged; only the scheduling changes.
The module's collect needs all fourteen receipts, so the same checks are
applied here to the merged items, with the stop record preserved.
"""
from pathlib import Path
import fcntl,hashlib,json,os,subprocess,sys,time
import numpy as np

OUT=Path('native-geometric-clock265-work');STATUS=Path('native-geometric-clock265-status.json')
QUEUE=OUT/'queue.json';LOCK=OUT/'queue.lock';ORIGINAL=list(range(14));CPUS=list(range(1,15))
module=Path('verification/propagate_geometric_boundary_clock.py')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='/home/lpaiu/work/nutimo_pilot/request13_deps:verification')


def completed(k):
    """Saved items of original worker k that also have their progress row."""
    try:
        z=np.load(OUT/f'work-{k}.npz');index=[int(v) for v in z['index']];values=z['values'];launch=z['launch']
        rows=read(OUT/f'work-{k}-progress.json')['rows']
    except Exception:return [],np.zeros((0,5)),np.zeros(0)
    n=min(len(index),len(rows))
    if [r['index'] for r in rows[:n]]!=index[:n]:return [],np.zeros((0,5)),np.zeros(0)
    return index[:n],values[:n],launch[:n]


def prepare():
    assert not QUEUE.exists() and not (OUT/'termination.json').exists()
    assert (OUT/'cost-stop-intent.json').exists() and not (OUT/'cost-stop.json').exists()
    assert not any((OUT/f'work-{k}-receipt.json').exists() for k in ORIGINAL)
    plan=read(OUT/'execution-plan.json');items=[];kept={}
    for k in ORIGINAL:
        done=completed(k)[0];kept[f'work-{k}']=len(done)
        items+=[v for v in plan['assignments'][k] if v['index'] not in set(done)]
    items=sorted(items,key=lambda v:-v['cells'])
    write(OUT/'termination.json',dict(classification='Counterexample candidate',
        event='Stopping the original WSL-launched controller (to stop workers6-9 for cost) also ended all fourteen workers of its session at about11:53:25KST. None wrote a receipt. This is a scheduling failure, not a physical or numerical gate failure.',
        intent=str(OUT/'cost-stop-intent.json'),kept_items=kept,remaining_items=len(items),remaining_cells=sum(v['cells'] for v in items),
        lesson='In this WSL setup, ending a launched controller ends its child workers. Stop individual workers only through a controller that stays alive.',
        time_unix=time.time()))
    write(QUEUE,dict(items=items,claimed={}))


def claim(tag):
    with LOCK.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        q=read(QUEUE);pending=[v for v in q['items'] if str(v['index']) not in q['claimed']]
        if not pending:return None
        item=pending[0];q['claimed'][str(item['index'])]=tag;write(QUEUE,q);return item


def worker(tag):
    import resource
    sys.path.insert(0,'verification');import propagate_geometric_boundary_clock as g
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);g.p.incident.native.deadline(4*3600)
    for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
    plan=read(OUT/'execution-plan.json')
    for q,h in plan['bindings'].items():assert sha(q)==h,q
    start=time.monotonic();error=None;done=0
    try:
        g.p.initialize();a,order,geo=g.FINE;m=g.Photons(a,geo);assert np.array_equal(m.clock,np.array(plan['clock']))
        while (item:=claim(tag)) is not None:
            value,row=g.evaluate(m,item['time'],order)
            np.savez_compressed(OUT/f"queue-{item['index']}-tmp.npz",value=value,launch=np.array(row['physical_launch_energy_increment_erg']))
            os.replace(OUT/f"queue-{item['index']}-tmp.npz",OUT/f"queue-{item['index']}.npz")
            write(OUT/f"queue-{item['index']}.json",dict(item,worker=tag,packets=row['packets'],seconds=row['seconds'],
                energy_identity_relative=row['energy_identity_relative'],angular_invariant_relative=row['angular_invariant_relative']))
            done+=1
    except BaseException as exc:error=repr(exc);raise
    finally:write(OUT/f'queue-{tag}-receipt.json',dict(seconds=time.monotonic()-start,items=done,error=error,source_sha256=sha(__file__),
        module_sha256=sha(module),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


def collect():
    sys.path.insert(0,'verification');import propagate_geometric_boundary_clock as g
    plan=read(OUT/'execution-plan.json');t=g.metric_clock();n=len(t);LD=g.LD
    values=np.zeros((n,5),LD);launch=np.zeros(n);owner=[None]*n;rows=[]
    evaluation=np.array(t,float);knot=np.full(n,-1)
    for group in plan['assignments']:
        for item in group:
            evaluation[item['index']]=item['time']
            if item['knot'] is not None:knot[item['index']]=item['knot']
    def put(i,value,lval,who):
        assert owner[i] is None,(i,owner[i],who);owner[i]=who;values[i]=value;launch[i]=lval
    for k in ORIGINAL:
        assert not (OUT/f'work-{k}-receipt.json').exists()
        index,vals,lau=completed(k)
        for i,v,l in zip(index,vals,lau):put(i,v,l,f'work-{k}')
        rows+=read(OUT/f'work-{k}-progress.json')['rows'][:len(index)]
    q=read(QUEUE)
    for item in q['items']:
        i=item['index'];z=np.load(OUT/f'queue-{i}.npz');r=read(OUT/f'queue-{i}.json')
        assert r['index']==i and r['time']==item['time'];put(i,z['value'],float(z['launch']),r['worker']);rows.append(r)
    receipts=sorted(OUT.glob('queue-q*-receipt.json'));assert receipts
    for f in receipts:assert read(f)['error'] is None,f
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
        scheduling=dict(items_by_process=counts,termination=str(OUT/'termination.json'),queue=str(QUEUE),collector=str(Path(__file__)),
            collector_sha256=sha(__file__),module_collect_not_run='The original fourteen workers have no receipt; the module collect requires them.'),
        boundary_applied_to_matter=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed,result


def control():
    prepare();start=time.monotonic();children={};source=sha(__file__)
    def pending():
        q=read(QUEUE);return len(q['items'])-len(q['claimed'])
    for cpu in CPUS:
        tag=f'q{cpu}'
        with open(f'.phase265-geometric-{tag}.stdout.log','w') as o,open(f'.phase265-geometric-{tag}.stderr.log','w') as e:
            children[tag]=subprocess.Popen(['taskset','-c',str(cpu),sys.executable,__file__,'worker',tag],stdout=o,stderr=e,env=env)
    status=lambda state,**kw:write(STATUS,dict(state=state,action='queue',seconds=time.monotonic()-start,
        workers={k:v.pid for k,v in children.items()},pending_items=pending(),**kw))
    try:
        while any(c.poll() is None for c in children.values()):
            for tag,child in children.items():
                code=child.poll()
                if code is not None and code!=0:raise AssertionError((tag,code))
            status('running');time.sleep(15)
        for tag,child in children.items():assert child.returncode==0,(tag,child.returncode)
        assert pending()==0 and sha(__file__)==source
        status('running',phase='collect')
        code=subprocess.run([sys.executable,__file__,'collect'],env=env).returncode;assert code==0,('collect',code)
        status('completed',result=read(OUT/'result.json'))
    except BaseException as exc:
        for child in children.values():
            if child.poll() is None:child.terminate()
        status('failed',error=repr(exc));raise


if __name__=='__main__':
    action=sys.argv[1]
    if action=='control':control()
    elif action=='worker':worker(sys.argv[2])
    elif action=='collect':collect()
    else:raise SystemExit(action)
