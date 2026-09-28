"""Read a tiny material return directly between its two saved finite states."""
from pathlib import Path
import json,resource,signal,sys,time
import numpy as np
import def_retained_motion_return as run


def prepare():
    r=run.read(run.OUT/'infinity/result.json');assert not r['passed']
    run.write(run.OUT/'direct-readout-plan.json',dict(classification='Counterexample candidate',
        failure=r['motion_return_controls'],claim='Resolve the pressure subtraction at the actual two material states before rejecting the physical return for its3.46percent charge time difference.',
        method='Keep Phase150 charge/source as the accepted baseline. Recover both saved finite material states with the existing extended spline/primitive owner, subtract their pressures before GR integration, and use the separately saved photon increment. No trajectory, native state, bank or clock replay.',
        scope='This repairs the represented source difference; it does not improve the absolute old baseline or certify primitive/EOS errors.',
        pilot_knots=[0,8,16],gates=dict(time=.02,quadrature=.002,independent=1e-9),
        cap_seconds=130,original_readout_cap_seconds=180,charged_first_readout_seconds=run.read(run.OUT/'readout-receipt.json')['seconds'],
        forecast='Measure the three representative source knots, reuse them, and require twice worst knot cost for14 remaining plus35s for saved GR/infinity reads inside the remaining130s repair allowance.',
        stop='Retain the original failed result. Stop without another evolution or clock on time/support/precision/budget failure.',
        source_sha256=run.sha(__file__),bindings={str(p):run.sha(p) for p in [run.OUT/'execution-plan.json',run.OUT/'infinity/result.json',run.OUT/'result.json',run.OUT/'material/sources.json']}))


def main():
    from types import SimpleNamespace
    import verify_native_stage_energy_charge as charge
    import def_native_global_scalar_closure as exterior
    start=time.monotonic();plan=run.read(run.OUT/'direct-readout-plan.json')
    assert run.sha(__file__)==plan['source_sha256']
    for p,h in plan['bindings'].items():assert run.sha(p)==h,p
    run.initialize();s=run.State();m=s.m;mat=m.material;LD=np.longdouble;C=run.C;AMP=run.AMP
    originals={n:dict(np.load(run.BEFORE/f'material/steps-{n}-reference-128.npz')) for n in [64,128]}
    latest={n:dict(np.load(run.MATERIAL/f'steps-{n}-reference-128.npz')) for n in [64,128]}
    states={n:(originals[n]['history_scaled'][::n//16],latest[n]['history_scaled'][::n//16]) for n in [64,128]}
    cache={};costs=[]
    def pressure():
        model=m.model;f=model.flow;nb=m.nb;f.eos.y=m.y
        p,u,*_=f.eos.evaluate(m.rho,m.lt);pg=p*LD(f.eos.rho0)*LD(C)**2
        deep=model.bulk.eos.gas(m.theta,m.eta)[0]
        gas=np.r_[deep,pg]*mat.V;radial=gas.copy();radial[:nb]+=2*model.kinetic()
        H=(m.rho*(LD(f.eos.cx)+u)*LD(f.eos.rho0)*LD(C)**2+pg)
        radial[nb:]+=H*m.beta*m.beta/(1-m.beta*m.beta)*mat.V[nb:]
        return np.array([gas,radial],LD)
    def point(k):
        mark=time.monotonic();s.precision(False);m.model.flow.seed=np.asarray(m.model.flow.seed,float);mat.point(k)
        s.precision(True);m.coefficients=pressure
        for n,(before,after) in states.items():
            a=s.coefficients(k,before[k],1.);b=s.coefficients(k,after[k],1.)
            cache[n,k]=(b-a,16*np.finfo(LD).eps*(abs(a)+abs(b)))
        s.precision(False);costs.append(time.monotonic()-mark)
    for k in [0,8,16]:point(k)
    forecast=2*max(costs)*14+35
    run.write(run.OUT/'direct-readout-pilot.json',dict(classification='Counterexample candidate',point_seconds=costs,upper_remaining_seconds=forecast,remaining_seconds=130-(time.monotonic()-start),eligible=forecast<130-(time.monotonic()-start)))
    assert forecast<130-(time.monotonic()-start)
    for k in range(17):
        if (128,k) not in cache:point(k)
    dest=run.OUT/'direct-motion';dest.mkdir();model=charge.gr.Response();waves={};pressure_rows=[]
    for n,(a,b) in states.items():
        z=(b.astype(LD)-a.astype(LD))*LD(AMP);p=np.load(run.PHOTON/f'steps-{n}-reference-128.npz')
        dp=np.array([cache[n,k][0] for k in range(17)]);floor=np.array([cache[n,k][1] for k in range(17)])
        e=(z[:,2]+LD(mat.rest)*z[:,0])/mat.a;rest=z[:,0]*LD(mat.model.cx)*LD(C)**2
        ph=p['moments'][:,0]/mat.a;pr=p['moments'][:,5]/mat.a
        d=run.prior.template();d.update(baryon_g=z[:,0],gas_nonrest_energy_erg=e-rest,
            nonrest_trace_erg=e-dp[:,1]-2*dp[:,0]-rest,nonrest_stress_erg=e-dp[:,1]-rest,pressure_volume_erg=dp[:,0],
            photon_energy_erg=ph,photon_radial_pressure_erg=pr,metric_stress_erg=e+ph-dp[:,1]-pr,
            inner_cumulative_energy_erg=p['radial_ports'][:,0,1],outer_cumulative_energy_erg=p['radial_ports'][:,1,1])
        np.savez_compressed(dest/f'source-{n}.npz',**d)
        for q in ([8,4] if n==128 else [8]):
            waves[n,q]=charge.read(model,d,q);np.savez_compressed(dest/f'wave-{n}-g{q}.npz',**waves[n,q])
        pressure_rows.append(dict(steps=n,arithmetic_over_pressure_increment=float(np.max(np.sum(floor,axis=2))/max(np.max(np.sum(abs(dp),axis=2)),1.))))
    fine=waves[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-300)
    direct,coord=charge.independent.direct(model,dict(np.load(dest/'source-128.npz')),8)
    independent=float(abs(direct-waves[128,8]['direct_scalar'][-1])/max(np.max(abs(waves[128,8]['direct_scalar'])),1e-300))
    controls=dict(compact_time=float(np.max(abs(fine-waves[64,8]['free_scalar']))/norm),compact_quadrature=float(np.max(abs(fine-waves[128,4]['free_scalar']))/norm),independent=independent)
    result=dict(classification='Counterexample candidate',controls=controls,pressure=pressure_rows,compact_endpoint=float(fine[-1]),seconds=time.monotonic()-start)
    # First measure, then decide; do not rewrite the original failed verdict.
    final=[]
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        now=dict(np.load(run.OUT/f'infinity/charge-{n}-a{a}-r{r}.npz'));old=dict(np.load(run.BEFORE/f'infinity/charge-{n}-a{a}-r{r}.npz'))
        alpha=-float(model.K/model.M) if hasattr(model,'K') else -float(np.load(run.GR/'applied-source.npz')['K_cm']/np.load(run.GR/'applied-source.npz')['M_cm'])
        mass=float(np.load(run.GR/'applied-source.npz')['M_cm'])
        de=exterior.G/C**4*(now['arrived_parts_erg'][:,1].astype(LD)-old['arrived_parts_erg'][:,1].astype(LD))/mass
        ds=now['normalized_exterior_parts'][:,1].astype(LD)-old['normalized_exterior_parts'][:,1].astype(LD)
        inc=(waves[n,8]['free_scalar'].astype(LD)+ds+(LD(alpha)+old['normalized'])*de)/(1-now['epsilon'])
        now.update(motion_return=inc,normalized=old['normalized'].astype(LD)+inc,compact=old['compact'].astype(LD)+waves[n,8]['free_scalar'],native_return=old['native_return'].astype(LD)+inc)
        now['scalar']=now['compact']+now['normalized_exterior'];final.append(now)
        np.savez_compressed(dest/f'charge-{n}-a{a}-r{r}.npz',**now)
    f=final[0]['motion_return'];norm=max(np.max(abs(f)),1e-300)
    result.update(motion_return_controls={key:float(np.max(abs(z['motion_return']-f))/norm) for key,z in zip(['angular','radial','time'],final[1:])},
        endpoint_normalized=float(final[0]['normalized'][-1]),motion_return_endpoint=float(f[-1]),original_failure_preserved=True,
        actual_material_and_photon_paths_replayed=False,final_charge_solved=False,full_goal_complete=False)
    result['passed']=controls['compact_time']<.02 and controls['compact_quadrature']<.002 and independent<1e-9 and result['motion_return_controls']['time']<.02 and max(result['motion_return_controls'][k] for k in ['angular','radial'])<.002
    run.write(dest/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        def timeout(*_):raise TimeoutError('130s remaining readout repair cap')
        signal.signal(signal.SIGALRM,timeout);signal.alarm(130);start=time.monotonic();cpu=time.process_time();error=None
        try:main()
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.alarm(0);run.write(run.OUT/'direct-readout-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=run.sha(__file__)))
