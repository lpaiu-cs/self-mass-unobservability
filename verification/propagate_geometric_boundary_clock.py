"""Evaluate the accepted incident-metric background-photon lapse at every applied GR time.

Counterexample candidate. The phase-261 vacuum photon propagation, fine packet
quadrature (angular 8, launch-time 8, geometry 8), background launch cohorts and
particular mass/lapse projection of .phase261-boundary.py are reused unchanged.
Only the output times change: the 575-point union stage clock of the applied
returned metric replaces the 17 background knots, and the last launch interval
is cut at the output time. A metric time within 1e-18 s of a knot is evaluated
at that knot, exactly as the stage driver aligns times. The result is an input
for the complete one-return outer lapse, not a solved boundary, scalar closure
or charge.
"""
from pathlib import Path
import json,os,resource,sys,time
import numpy as np
from numpy.polynomial import legendre as leg
import propagate_exterior_vacuum as adapter

p=adapter.base
SOURCE=adapter.OUT
ACTUAL=p.ACTUAL
OUT=Path('native-geometric-clock265-work')
p.OUT=adapter.OUT=OUT
read,write,sha,LD,G,C=p.read,p.write,p.sha,p.LD,p.G,p.C
FINE=p.SETTINGS['fine']
KNOT=1e-18
WORKERS=14
CAPS=dict(prepare=300,pilot=1800,work=10800,collect=900)


class Photons(p.Photons):
    def cohorts(self,now,order,cells=None):
        if np.any(self.clock==now):return super().cohorts(now,order,cells)
        assert cells is None
        edges=np.r_[self.clock[self.clock<now],now]
        gx,gw=leg.leggauss(order);ids=np.arange(len(edges)-1)
        dt=np.diff(edges)[ids];te=(edges[ids,None]+dt[:,None]*(gx+1)/2).ravel()
        tw=(dt[:,None]*gw/2).ravel();lum=np.column_stack([np.interp(te,self.clock,np.asarray(v,float)) for v in self.lum.T])
        owner=np.tile(np.arange(len(self.mu)),len(te));launch=np.repeat(te,len(self.mu))
        energy=(tw[:,None]*lum[:,self.bins]*self.mu*self.mw).ravel().astype(LD)
        return launch,owner,energy


def evaluate(m,now,order):
    """The .phase261-boundary.py projection, applied to one live snapshot."""
    z,row=m.propagate(now,order);d=m.d
    r=z['radius_cm'];mu=z['direction'];dr=z['delta_radius_cm'];dm=z['delta_direction']
    _,N,b,a,_=d.bg.metric(r/d.model.m.R);K=1/(r*a)
    local=z['delta_log_H']-z['metric_nu']-z['metric_lambda'];energy=z['background_packet_energy_erg']
    mass=p.LD(p.G)/p.LD(p.C)**4*np.sum(energy*local,dtype=p.LD)
    pieces=np.array([-p.LD(p.G)/p.LD(p.C)**4*np.sum(energy*K*v,dtype=p.LD) for v in
        [(1+mu*mu)*local,-(1+mu*mu)*dr/(r*b),2*mu*dm]])
    return np.r_[mass,pieces.sum(dtype=p.LD),pieces],row


def metric_clock():
    t=np.load(ACTUAL/'metric/metric-128-g8.npz')['t']
    for n,q in [(64,8),(128,4)]:assert np.array_equal(np.load(ACTUAL/f'metric/metric-{n}-g{q}.npz')['t'],t)
    return t


def prepare():
    assert not OUT.exists()
    assert read(SOURCE/'result.json')['passed'] and read(SOURCE/'photon-boundary/result.json')['passed']
    assert read(ACTUAL/'result.json')['passed'] and read(ACTUAL/'controller-status.json')['state']=='completed'
    OUT.mkdir();t=metric_clock()
    files=[Path(__file__),Path(adapter.__file__),Path(p.__file__),SOURCE/'plan.json',SOURCE/'result.json',
        SOURCE/'photon-boundary/result.json',SOURCE/'photon-boundary/source.npz',SOURCE/'photon-boundary/plan.json',
        ACTUAL/'result.json']+[ACTUAL/f'metric/metric-{n}-g{q}.npz' for n,q in [(64,8),(128,4),(128,8)]]
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Evaluate the accepted phase-261 incident-metric background-photon mass and particular outer lapse at every time of the applied returned metric, as the missing term of the complete one-return outer boundary.',
        method='Unchanged vacuum propagation, fine packet quadrature (angular8, launch-time8, geometry8), background luminosity knots and .phase261-boundary.py projection. Output time changes from the17background knots to the575union stage clock; the last launch interval is cut at the output time with the same8-point Gauss rule. Times within1e-18s of a knot use that knot, as the stage driver aligns times.',
        decision='Only exact reproduction of all16saved knot values and the original per-snapshot energy-identity and angular-invariant gates admit this input to the new applied metric. The launch-quadrature and angular/geometry controls of phase261 (max2.64e-9) are inherited, not re-run.',
        gates=dict(knot_reproduction='bitwise',energy_identity=.002,angular_invariant=1e-10),
        budgets=CAPS,workers=WORKERS,virtual_GiB_per_worker=16,output_times=len(t),
        forecast='Phase261fine took484s for its16knots (0.0027..0.0163s per packet). The575times need1308928packet propagations: about5.0h CPU at the mean rate. The pilot measures a late partial time and sets the parallel forecast with a2x margin before any production.',
        stop='Any knot mismatch, original propagation gate or cap. No new quadrature order, grid or clock beyond the applied metric clock.',
        bindings={str(q):sha(q) for q in files},physical_final_charge_solved=False,full_goal_complete=False))


def assignments(m,t):
    clock=m.clock;items=[]
    for i,now in enumerate(t):
        if i==0:assert now==0;continue
        j=int(np.argmin(abs(clock-now)));knot=abs(clock[j]-now)<=KNOT
        at=float(clock[j]) if knot else float(now)
        cells=j if knot else int(np.count_nonzero(clock<now))
        items.append(dict(index=i,time=at,metric_time=float(now),knot=j if knot else None,cells=cells))
    order=sorted(items,key=lambda v:-v['cells']);load=[0]*WORKERS;groups=[[] for _ in range(WORKERS)]
    for item in order:
        k=int(np.argmin(load));groups[k].append(item);load[k]+=item['cells']
    return items,groups,load


def pilot():
    assert not (OUT/'pilot.json').exists()
    for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
    p.initialize();a,order,g=FINE;m=Photons(a,g);t=metric_clock()
    items,groups,load=assignments(m,t);saved=np.load(SOURCE/'photon-boundary/source.npz')
    knots=[v for v in items if v['knot'] is not None];assert sorted(v['knot'] for v in knots)==list(range(1,17))
    rows=[];values={}
    first=next(v for v in knots if v['knot']==1)
    late=max((v for v in items if v['knot'] is None),key=lambda v:v['cells'])
    for item in [first,late]:
        start=time.monotonic();value,row=evaluate(m,item['time'],order);row['wall_seconds']=time.monotonic()-start
        rows.append(dict(item,**row));values[item['index']]=value
    v=values[first['index']];j=first['knot']
    exact=bool(v[0]==saved['photon_J_source_cm'][j] and v[1]==saved['photon_lapse_source'][j] and np.array_equal(v[2:],saved['lapse_energy_radius_angle_parts'][j-1]))
    per=max(r['wall_seconds']/r['packets'] for r in rows)
    packets=256*sum(v['cells'] for v in items);upper=2*per*max(load)*256+120
    result=dict(classification='Counterexample candidate',rows=rows,knot1_exact=exact,
        per_packet_seconds=per,total_packets=packets,forecast_upper_seconds_per_worker=upper,
        worker_cells=load,eligible=exact and upper<CAPS['work'])
    np.savez_compressed(OUT/'pilot.npz',index=np.array(list(values)),values=np.array(list(values.values()),LD))
    write(OUT/'pilot.json',result);print(json.dumps(result),flush=True)
    assert result['eligible'],result
    write(OUT/'execution-plan.json',dict(classification='Conjectural',workers=WORKERS,assignments=groups,
        clock=[float(v) for v in m.clock],forecast_upper_seconds_per_worker=upper,CPU_affinities=list(range(1,1+WORKERS)),
        bindings={str(q):sha(q) for q in [Path(__file__),OUT/'plan.json',OUT/'pilot.json',OUT/'pilot.npz']}))


def work(k):
    plan=read(OUT/'execution-plan.json')
    for q,h in plan['bindings'].items():assert sha(q)==h,q
    group=plan['assignments'][k];p.initialize();a,order,g=FINE;m=Photons(a,g)
    assert np.array_equal(m.clock,np.array(plan['clock']))
    index=[];values=[];launch=[];diagnostics=[];start=time.monotonic()
    for item in group:
        value,row=evaluate(m,item['time'],order)
        index.append(item['index']);values.append(value);launch.append(row['physical_launch_energy_increment_erg'])
        diagnostics.append(dict(index=item['index'],time=item['time'],packets=row['packets'],seconds=row['seconds'],
            energy_identity_relative=row['energy_identity_relative'],angular_invariant_relative=row['angular_invariant_relative']))
        np.savez_compressed(OUT/f'work-{k}.npz',index=np.array(index),values=np.array(values,LD),launch=np.array(launch))
        write(OUT/f'work-{k}-progress.json',dict(completed=len(index),total=len(group),seconds=time.monotonic()-start,rows=diagnostics))


def collect():
    plan=read(OUT/'execution-plan.json');t=metric_clock();n=len(t)
    values=np.zeros((n,5),LD);launch=np.zeros(n);filled=np.zeros(n,bool);rows=[]
    evaluation=np.array(t,float);knot=np.full(n,-1)
    for group in plan['assignments']:
        for item in group:
            evaluation[item['index']]=item['time']
            if item['knot'] is not None:knot[item['index']]=item['knot']
    for k in range(plan['workers']):
        assert read(OUT/f'work-{k}-receipt.json')['error'] is None
        z=np.load(OUT/f'work-{k}.npz');progress=read(OUT/f'work-{k}-progress.json')
        assert progress['completed']==progress['total']==len(plan['assignments'][k])
        assert not filled[z['index']].any();values[z['index']]=z['values'];launch[z['index']]=z['launch'];filled[z['index']]=True
        rows+=progress['rows']
    assert filled[1:].all() and not filled[0]
    saved=np.load(SOURCE/'photon-boundary/source.npz');exact=[]
    for j in range(1,17):
        i=int(np.flatnonzero(knot==j)[0]);v=values[i]
        exact.append(bool(v[0]==saved['photon_J_source_cm'][j] and v[1]==saved['photon_lapse_source'][j] and np.array_equal(v[2:],saved['lapse_energy_radius_angle_parts'][j-1])))
    pilot=np.load(OUT/'pilot.npz');repeat=bool(all(np.array_equal(values[i],v) for i,v in zip(pilot['index'],pilot['values'])))
    identity=max(r['energy_identity_relative'] for r in rows);invariant=max(r['angular_invariant_relative'] for r in rows)
    lapse=np.asarray(values[:,1],float);applied=np.asarray(np.load(ACTUAL/'metric/metric-128-g8.npz')['delta_nu_faces'][:,-1],float)
    parts=np.max(abs(np.asarray(values[:,2:],float)),axis=0)
    passed=all(exact) and repeat and identity<.002 and invariant<1e-10
    np.savez_compressed(OUT/'boundary-575.npz',t=t,evaluation_t=evaluation,knot=knot,
        photon_geometric_mass_cm=values[:,0],photon_geometric_lapse=values[:,1],lapse_energy_radius_angle_parts=values[:,2:],
        physical_launch_energy_increment_erg=launch)
    result=dict(classification='Counterexample candidate',passed=passed,knot_values_exact=exact,pilot_repeat_exact=repeat,
        output_times=n,new_times=int(np.count_nonzero(knot[1:]<0)),maximum_energy_identity_relative=identity,
        maximum_angular_invariant_relative=invariant,
        maximum_geometric_lapse_over_applied_outer_lapse=float(np.max(abs(lapse))/np.max(abs(applied))),
        maximum_lapse_energy_radius_angle_parts=parts.tolist(),
        terminal_photon_geometric_mass_cm=float(values[-1,0]),terminal_photon_geometric_lapse=float(values[-1,1]),
        inherited_phase261_quadrature_controls=read(SOURCE/'photon-boundary/result.json')['controls'],
        boundary_applied_to_matter=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed,result


if __name__=='__main__':
    action=sys.argv[1];key=action.split('-')[0];assert key in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(CAPS[key])
    receipt=OUT/f'{action}-receipt.json'
    try:
        if action!='prepare':
            assert not receipt.exists()
            for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
        if key=='work':work(int(action.split('-')[1]))
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
