"""Native conserved EOS on the corrected stored trajectory, applied to GR.

Counterexample candidate: 17 actual snapshots, not continuous EOS certification.
Reuse inverse and retarded owners; no physical trajectory is replayed.
"""
from pathlib import Path
import json, multiprocessing as mp, resource, signal, sys, time
import numpy as np
import sympy as sp
import verify_native_updated_gr_return as updated
import def_native_conservative_eos_readout as deep
import def_native_atmosphere_inverse as atmosphere

physical=updated.physical;flow=physical.flow;chem=physical.chem
charge=deep.charge;C=flow.C;write=flow.write;sha=flow.sha
OUT=updated.OUT/'conserved-history-charge'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(updated.__file__),Path(deep.__file__),Path(atmosphere.__file__),
           Path(chem.__file__),Path(chem.old.__file__),Path(charge.__file__),
           updated.run.BG/'coupled-128.npz',physical.OUT/'bank.npz',physical.EV/'source-128.npz',
           physical.GR/'result.json',deep.OUT/'production-samples.json',atmosphere.OUT/'production-samples.json',
           atmosphere.OUT/'source-128.npz',flow.INPUT/'balanced-20.npz',
           chem.old.old.cold.CACHE/'stable/gas.so',chem.old.old.cold.CACHE/'libfree_eos_native_cold_stable.so']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='7642e22e2',
        claim='Invert the SAME native EOS at all active deep/atmospheric conserved states at all17 saved canonical times of the corrected fine trajectory, and insert full pressure/radial-stress/trace differences into the actual retarded charge.',
        decision='If source replacement changes charge by2percent or constitutive pressure by0.2percent, repair the EOS owner before further physical interpretation. Otherwise quantify this sampled EOS source uncertainty and advance the still-open derivative/feedback terms.',
        reuse='No new trajectory, mesh, horizon, EOS table or waveform iteration. Exact conserved-coordinate matches may reuse native roots; near matches are not accepted as cached proof. All initial pressure offsets are retained.',
        closure='Deep anchors retain original other ionic/molecular inventories and explicit additive first-order inventory corrections. Atmosphere fixes D,S,K,H and lapse. Below-floor cells remain excluded by the original EOS owner; the correction there is NOT certified zero.',
        gates=dict(energy=1e-12,momentum=1e-12,density=1e-10,reseed=1e-10,pressure=.002,charge=.02,cadence=.002,quadrature=.002,linearity=1e-10),
        budget=dict(pilot_seconds=90,production_wall_seconds=900,production_worker_seconds=2400,readout_seconds=90,
                    workers=3,threads_per_worker=1,virtual_GiB_per_process=3,native_calls_per_worker_job=24000,total_new_native_call_cap=80000),
        forecast='Previous223 atmospheric roots68.17s and323 old deep roots within120s. Pilot actual current initial/middle/final states, forecast at twice the slowest measured root per region plus30s setup. Three process-isolated native workers; no shared Fortran state. Interpreter imports/binding IO reported outside driver caps.',
        limits='17 saved states do not certify all129 integration states or continuous time. Source replacement fixes photons, metric, conserved trajectory and the original initial GR coefficients; it is not EOS re-evolution, derivative enclosure, coupled gain or final physical charge.',
        stop='Stop on failed root, forecast, pressure/materiality gate or cap. Preserve successful states; do not automatically enlarge budget, physical paths or support.',
        bindings={str(p):sha(p) for p in files}))
    S,Q,p,dp=sp.symbols('S Q p dp');v=S/(Q+p);vn=S/(Q+p+dp)
    assert sp.factor(-S*(vn-v)-dp-(v*vn-1)*dp)==0
    assert sp.factor(-S*(vn-v)-3*dp-(v*vn-3)*dp)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        premise='Fixed D,S,tau,H and lapse; Q=cx*D+tau. Deep velocity is fixed by its separate conserved momentum owner.',
        atmosphere_trace='(v_old*v_native-3)*delta_p',atmosphere_metric_stress='(v_old*v_native-1)*delta_p',
        deep_trace='-3*delta_p',deep_metric_stress='-delta_p',initial_offsets_retained=True))


def restore(model,z,k):
    model.Pi=z['snapshot_Pi'][k].copy();model.h=z['snapshot_h'][k].copy();model.j=z['snapshot_j'][k].copy()
    model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(float(z['snapshot_t'][k]))


def load_inputs():
    updated.configure();model=updated.Model();f=model.flow;e=model.bulk.eos
    source=np.load(updated.run.BG/'coupled-128.npz')
    z={k:source[k] for k in ['snapshot_t','snapshot_U','snapshot_Pi','snapshot_h','snapshot_j','snapshot_theta','snapshot_eta','snapshot_u']}
    times=z['snapshot_t'];assert len(times)==17 and np.max(abs(times-np.linspace(0,flow.old.END,17)))<1e-17
    interior=[];exterior=[]
    for it,t in enumerate(times):
        restore(model,z,it);theta=z['snapshot_theta'][it];eta=z['snapshot_eta'][it];p,u,ut,*_=e.gas(theta,eta)
        interior.append(dict(t=float(t),x=np.log1p(e.density_shift)+np.log1p(e.x),xi=e.xi.copy(),
            lt=np.log(e.base.d['T'])+e.theta0+theta,y=e.base.d['y0']*(1+eta),
            target=z['snapshot_u'][it]+e.inventory[1]*e.xi,p=p.copy(),u=z['snapshot_u'][it],ut=ut.copy()))
        U=z['snapshot_U'][it];rho,v,lt,y=f.primitive(U);p,u,_,_,_,cv,_=f.eos.evaluate(rho,lt)
        active=np.flatnonzero(rho>=f.eos.floor)
        # The collision and EOS owners use rest density, not lab D.
        assert np.array_equal(active,np.flatnonzero(U[0]>=f.eos.floor)),('Floor classification differs',it)
        tau=(U[2]-(model.m.a-model.m.a0)*f.eos.cx*U[0])/model.m.a
        exterior.append(dict(U=U,rho=rho,v=v,lt=lt,y=y,p=p,u=u,cv=cv,tau=tau,active=active))
    return model,interior,exterior


def decorate(kind,it,j,r,model,ds,ats):
    r=dict(r,kind=kind,it=it,cell=j,t=ds[it]['t'])
    if kind=='deep':
        d=ds[it];e=model.bulk.eos
        r.update(x=float(d['x'][j]),xi=float(d['xi'][j]),y=float(d['y'][j]),target_u=float(d['target'][j]),
            saved_u=float(d['u'][j]),old_lt=float(d['lt'][j]),old_pressure=float(d['p'][j]))
        r['pressure']=r['raw'][1]-float(e.inventory[0,j]*d['xi'][j]);r['pressure_relative']=r['pressure']/r['old_pressure']-1
        r['density_relative']=abs(r['raw'][0]/(float(e.base.d['rho'][j])*np.exp(r['x']))-1)
    else:
        d=ats[it];r.update(old_pressure=float(d['p'][j]),old_velocity=float(d['v'][j]),old_lt=float(d['lt'][j]))
        r['pressure_relative']=r['p']/r['old_pressure']-1
    assert abs(r['pressure_relative'])<.002 and r['density_relative']<1e-10
    return r


def cached(model,ds,ats):
    rows=[]
    for r in json.loads((deep.OUT/'production-samples.json').read_text()):
        if r['it']!=0:continue
        j=r['cell'];d=ds[0]
        if all(float(d[k][j])==r[v] for k,v in [('x','x'),('y','y'),('target','target_u')]):
            rows.append(decorate('deep',0,j,r,model,ds,ats))
    old=dict(np.load(atmosphere.OUT/'source-128.npz'))
    for r in json.loads((atmosphere.OUT/'production-samples.json').read_text()):
        j=r['cell']
        if r['steps']==128 and np.array_equal(old['U'][:,j],ats[16]['U'][:,j]):
            rows.append(decorate('atmosphere',16,j,r,model,ds,ats))
    return rows


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.alarm(90)
    model,ds,ats=load_inputs();rows=cached(model,ds,ats);reused=len(rows);durations={'deep':[],'atmosphere':[]};checks=[]
    native=chem.old.Native(cap=1200);native.prefix=dict(native.prefix)
    for it in [0,8,16]:
        active=ats[it]['active']
        for j in [int(active[0]),int(active[len(active)//2]),int(active[-1])]:
            if any(r['kind']=='atmosphere' and r['it']==it and r['cell']==j for r in rows):continue
            tick=time.monotonic();r=atmosphere.root(native,model.flow,ats[it],j);durations['atmosphere'].append(time.monotonic()-tick)
            rows.append(decorate('atmosphere',it,j,r,model,ds,ats))
            if j==int(active[-1]):
                other=atmosphere.root(native,model.flow,ats[it],j,.001);err=abs(other['p']/r['p']-1);assert err<1e-10
                checks.append(dict(kind='atmosphere',it=it,cell=j,reseed_relative=err))
    for it,j in [(0,0),(8,9),(16,18)]:
        chem.setup(native,model.bulk.eos.base.d,j);tick=time.monotonic();r=deep.root(native,ds[it],j)
        durations['deep'].append(time.monotonic()-tick)
        rows=[v for v in rows if not(v['kind']=='deep' and v['it']==it and v['cell']==j)]
        rows.append(decorate('deep',it,j,r,model,ds,ats))
        other=deep.root(native,ds[it],j,.001);err=abs(other['raw'][1]/r['raw'][1]-1);assert err<1e-10
        checks.append(dict(kind='deep',it=it,cell=j,reseed_relative=err))
    totals={'deep':19*17,'atmosphere':sum(len(a['active']) for a in ats)}
    remaining={k:totals[k]-sum(r['kind']==k for r in rows) for k in totals}
    work={k:2*max(durations[k])*remaining[k] for k in totals}
    forecast=work['deep']+work['atmosphere']/3+30
    write(OUT/'pilot-samples.json',rows)
    result=dict(classification='Counterexample candidate',passed=True,eligible=forecast<900 and sum(work.values())<2400,
        totals=totals,remaining=remaining,exact_cached_states=reused,native_calls=native.ion.calls,root_seconds=durations,
        forecast_upper_wall_seconds=forecast,forecast_upper_worker_seconds=sum(work.values()),reseed=checks,
        active_atmosphere_per_time=[len(a['active']) for a in ats],seconds=time.monotonic()-start,
        source_sha256=sha(__file__))
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def worker(job):
    kind,group=job;model,ds,ats=INPUTS;name=f'{kind}-{group}';start=time.monotonic();rows=[]
    native=chem.old.Native(cap=24000);native.prefix=dict(native.prefix)
    try:
        if kind=='deep':
            for j in range(19):
                chem.setup(native,model.bulk.eos.base.d,j)
                for it in range(17):
                    if (kind,it,j) not in DONE:rows.append(decorate(kind,it,j,deep.root(native,ds[it],j),model,ds,ats))
                write(OUT/f'{name}-partial.json',rows)
        else:
            assert abs(np.exp(native.lr)/model.flow.eos.rho0-1)<1e-10
            for it in range(group,17,3):
                for j in map(int,ats[it]['active']):
                    if (kind,it,j) not in DONE:rows.append(decorate(kind,it,j,atmosphere.root(native,model.flow,ats[it],j),model,ds,ats))
                write(OUT/f'{name}-partial.json',rows)
        write(OUT/f'{name}-samples.json',rows)
        result=dict(job=name,states=len(rows),native_calls=native.ion.calls,seconds=time.monotonic()-start)
        write(OUT/f'{name}.json',result);return result
    except Exception as exc:
        write(OUT/f'{name}-partial.json',rows);write(OUT/f'{name}-failure.json',dict(error=repr(exc),completed=len(rows),native_calls=native.ion.calls,seconds=time.monotonic()-start));raise


def production():
    global INPUTS,DONE
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['eligible'] and pilot['source_sha256']==sha(__file__)
    assert not (OUT/'production.json').exists();start=time.monotonic();signal.alarm(900)
    INPUTS=load_inputs();rows=json.loads((OUT/'pilot-samples.json').read_text());DONE={(r['kind'],r['it'],r['cell']) for r in rows}
    # ponytail: fork isolates mutable Fortran globals; never use native threads.
    with mp.get_context('fork').Pool(3) as pool:jobs=pool.map(worker,[('deep',0),('atmosphere',0),('atmosphere',1),('atmosphere',2)])
    for job in jobs:rows.extend(json.loads((OUT/f"{job['job']}-samples.json").read_text()))
    keys={(r['kind'],r['it'],r['cell']) for r in rows};expected={('deep',it,j) for it in range(17) for j in range(19)}
    expected|={('atmosphere',it,int(j)) for it,d in enumerate(INPUTS[2]) for j in d['active']}
    assert keys==expected and len(rows)==len(keys)
    calls=sum(j['native_calls'] for j in jobs)+pilot['native_calls'];worker_seconds=sum(j['seconds'] for j in jobs)
    assert calls<80000 and worker_seconds<2400
    write(OUT/'production-samples.json',rows)
    result=dict(classification='Counterexample candidate',passed=True,states=len(rows),jobs=jobs,native_calls_including_pilot=calls,
        maximum_pressure_relative=max(abs(r['pressure_relative']) for r in rows),
        maximum_energy_relative=max(abs(r.get('energy_relative',r.get('energy_residual'))) for r in rows),
        maximum_momentum_relative=max(r.get('momentum_relative',0.) for r in rows),
        maximum_density_relative=max(r['density_relative'] for r in rows),seconds=time.monotonic()-start,
        worker_seconds=worker_seconds,all17_saved_active_states=True,continuous_EOS_bound=False,EOS_re_evolution=False)
    write(OUT/'production.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def source(model,ds,ats,rows,stride=1,region='all'):
    d=dict(np.load(physical.EV/'source-128.npz'));n=len(d['radius']);pV=np.zeros((17,n));trace=pV.copy();stress=pV.copy()
    volume=4*np.pi*model.m.RJ**2*model.m.vol;unit=model.flow.eos.rho0*C*C*volume
    for r in rows:
        it,j=r['it'],r['cell']
        if r['kind']=='deep' and region!='atmosphere':
            dp=(r['pressure']-r['old_pressure'])*model.bulk.volume[j];pV[it,j]=dp;trace[it,j]=-3*dp;stress[it,j]=-dp
        elif r['kind']=='atmosphere' and region!='deep':
            dp=(r['p']-r['old_pressure'])*unit[j];vv=r['v']*r['old_velocity'];k=19+j
            pV[it,k]=dp;trace[it,k]=(vv-3)*dp;stress[it,k]=(vv-1)*dp
    for key in physical.capture.KEYS+['metric_stress_erg']:d[key]=np.zeros_like(d[key])
    d['inner_cumulative_energy_erg']=np.zeros_like(d['t']);d['outer_cumulative_energy_erg']=np.zeros_like(d['t'])
    t=np.array([v['t'] for v in ds])
    for key,value in [('pressure_volume_erg',pV),('nonrest_trace_erg',trace),('nonrest_stress_erg',stress),('metric_stress_erg',stress)]:
        d[key]=np.array([np.interp(d['t'],t[::stride],col[::stride]) for col in value.T]).T
    return d,dict(t=t,pressure_volume_erg=pV,trace_erg=trace,metric_stress_erg=stress)


def readout():
    assert json.loads((OUT/'production.json').read_text())['passed'];assert not (OUT/'result.json').exists()
    start=time.monotonic();signal.alarm(90);model,ds,ats=load_inputs();rows=json.loads((OUT/'production-samples.json').read_text())
    m=charge.gr.Response();waves={}
    for label,stride,region,order in [('fine',1,'all',8),('nine',2,'all',8),('g4',1,'all',4),('deep',1,'deep',8)]:
        d,knots=source(model,ds,ats,rows,stride,region);wave=charge.read(m,d,order);waves[label]=wave
        np.savez_compressed(OUT/f'wave-{label}.npz',**wave)
        if label=='fine':np.savez_compressed(OUT/'source.npz',**d);np.savez_compressed(OUT/'source-knots.npz',**knots)
    d=dict(np.load(OUT/'source.npz'));total=dict(np.load(physical.EV/'source-128.npz'))
    for key in physical.capture.KEYS+['metric_stress_erg']:total[key]=total[key]+d[key]
    applied=charge.read(m,total,8);old=dict(np.load(physical.GR/'wave-128-g8.npz'));a=waves['fine']['free_scalar'];norm=float(np.max(abs(old['free_scalar'])))
    linear=float(np.max(abs(applied['free_scalar']-old['free_scalar']-a))/norm)
    cadence=float(np.max(abs(a-waves['nine']['free_scalar']))/norm);quad=float(np.max(abs(a-waves['g4']['free_scalar']))/norm)
    direct,error=charge.independent.direct(m,d,8);agreement=float(abs(direct-waves['fine']['direct_scalar'][-1])/max(np.max(abs(waves['fine']['direct_scalar'])),1e-300))
    materiality=float(np.max(abs(a))/norm);passed=linear<1e-10 and cadence<.002 and quad<.002 and agreement<1e-9 and materiality<.02
    np.savez_compressed(OUT/'applied-charge.npz',**applied,previous=old['free_scalar'],correction=a)
    result=dict(classification='Counterexample candidate',passed=passed,old_endpoint=float(old['free_scalar'][-1]),
        native_source_endpoint=float(applied['free_scalar'][-1]),correction_endpoint=float(a[-1]),
        deep_correction_endpoint=float(waves['deep']['free_scalar'][-1]),atmosphere_correction_endpoint=float(a[-1]-waves['deep']['free_scalar'][-1]),
        maximum_correction_relative_to_signal=materiality,correction_cadence_relative_to_signal=cadence,
        correction_quadrature_relative_to_signal=quad,independent_application_relative=linear,independent_direct_relative=agreement,
        inverse_radius_error=error,seconds=time.monotonic()-start,constitutive_repair_required=materiality>=.02,
        actual_saved_history_EOS_source_applied=True,initial_pressure_offsets_retained=True,
        native_correction_time_convergence_tested=False,full129_state_EOS_audit=False,uniform_EOS_derivative_bound=False,
        full_floor_feedback_enclosed=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,flow.old.optical.timeout)
    action=sys.argv[1]
    try:globals()[action]()
    except Exception as exc:
        if OUT.exists():write(OUT/f'{action}-failure.json',dict(error=repr(exc)))
        raise
