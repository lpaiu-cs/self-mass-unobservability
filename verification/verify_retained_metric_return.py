"""Current GR-driven photon/material increment -> compact and infinity charge."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,resource,signal,sys,time
import numpy as np
import sympy as sp
import def_retained_metric_return as run
import repair_retained_metric_material as repair
import verify_native_stage_energy_charge as charge

run.MATERIAL=repair.MATERIAL
OUT=run.OUT;GR=run.GR;read,write,sha=run.read,run.write,run.sha
AMP=run.AMP;C=run.C;LD=np.longdouble
NO_ALARM=SimpleNamespace(alarm=lambda n:0,signal=signal.signal,SIGALRM=signal.SIGALRM)


def prepare():
    assert read(run.MATERIAL/'production.json')['passed'];assert not (OUT/'readout-plan.json').exists()
    files=[Path(__file__),Path(run.__file__),Path(repair.__file__),OUT/'material-repair-execution.json',Path(charge.__file__),Path(run.aligned.base.__file__),
           Path(run.aligned.base.matter.previous.stress.__file__),Path(run.fields.previous.prior.__file__),
           run.PHOTON/'result.json',run.MATERIAL/'production.json',OUT/'response-execution-plan.json']
    for n in [64,128]:files += [run.PHOTON/f'steps-{n}-reference-128.npz',run.MATERIAL/f'steps-{n}-reference-128.npz']
    files += list((run.fields.BEFORE/'direct-motion').glob('charge-*.npz'))
    write(OUT/'readout-plan.json',dict(classification='Counterexample candidate',
        claim='Convert actual GR-driven material states and photons to the same canonical first-variation source, then apply this separate increment once to the accepted Phase151 charge.',
        source='Reuse the conservative primitive-to-pressure map and independent cell probes; remove the initial gas/photon canonical geometry already in the GR operator. Retain all energy, baryon, trace, radial and tangential stress fields.',
        infinity='Reuse exact causal packet fronts and the actual signed SDIRK angular packets. The previous compact/exterior/arrived-energy state is Phase151 direct-motion, including its accepted pressure-subtraction repair.',
        identity='delta_q=[delta_compact+delta_exterior+(alpha+q_previous)*delta_epsilon]/[1-epsilon_previous-delta_epsilon]. Store this increment separately; its controls use its own norm.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,independent_GR=1e-9,normalization=1e-12),
        budgets=dict(readout=180,audit=60),
        limits='Preserves the initial-slice GR response operator and finite packet representation. No uniform native Jacobian, continuum, exterior-scattering, fixed-point or final physical error enclosure.',
        bindings={str(p):sha(p) for p in files}))


def compact():
    owner=run.aligned.base;s=owner.source.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]')
    s=s.replace('background=compare(histories[2],histories[1]),','').replace(',stress_background=compare(allstress[2],allstress[1])','')
    s=run.fields.previous.prior.replace(s,"base=dict(np.load(wave.base.OUT/f'source-{ref}.npz'))",'base=template()')
    ns=dict(owner.source_scope,OUT=run.MATERIAL,GR=GR,Material=run.Material,signal=NO_ALARM,
        template=run.fields.previous.prior.template,photons=SimpleNamespace(path=lambda n,r:run.PHOTON/f'steps-{n}-reference-{r}.npz'),
        material_path=lambda n,r:run.MATERIAL/f'steps-{n}-reference-{r}.npz')
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-source.py').write_text(s);ns['sources']()
    assert read(run.MATERIAL/'sources.json')['passed']
    model=charge.gr.Response();waves={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));waves[n,q]=charge.read(model,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**waves[n,q])
    fine=waves[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-300)
    direct,coordinate=charge.independent.direct(model,dict(np.load(GR/'source-128-reference-128.npz')),8)
    controls=dict(time=float(np.max(abs(fine-waves[64,8]['free_scalar']))/norm),
        quadrature=float(np.max(abs(fine-waves[128,4]['free_scalar']))/norm),
        independent=float(abs(direct-waves[128,8]['direct_scalar'][-1])/max(np.max(abs(waves[128,8]['direct_scalar'])),1e-300)))
    result=dict(classification='Counterexample candidate',controls=controls,endpoint_metric_compact_increment=float(fine[-1]),
        inverse_radius_residual=float(coordinate),actual_current_metric_applied_to_photons_and_free_material=True,
        passed=controls['time']<.02 and controls['quadrature']<.002 and controls['independent']<1e-9,
        final_charge_solved=False,full_goal_complete=False)
    write(GR/'result.json',result);assert result['passed'],result


def infinity(remaining):
    prior=run.fields.previous.prior;s=inspect.getsource(prior.infinity)
    before=s.index("    original=dict(np.load(EOS/'applied-charge.npz'))")
    after=s.index('    for n,a,r in [(128,8,8)',before)
    s=s[:before]+"    prior=dict(np.load(BEFORE/'direct-motion/charge-128-a8-r8.npz'));alpha=-m.K/m.M;paths=[];rows=[]\n"+s[after:]
    s=prior.replace(s,"background=np.load(CURRENT/f'infinity/completed/retained-128-a{a}-r{r}.npz')",
        "background=np.load(BEFORE/f'direct-motion/charge-{n}-a{a}-r{r}.npz')")
    s=s.replace("background['exterior']","background['normalized_exterior']").replace("background['arrived']","background['arrived_energy_erg']")
    s=prior.replace(s,"previous=(baseline+d['normalized_exterior_parts'][:,0]+alpha*eps0)/(1-eps0)","previous=background['normalized']")
    s=prior.replace(s,"compact=baseline+wave['free_scalar']","compact=background['compact'].astype(LD)+wave['free_scalar']")
    s=s.replace('native_return','metric_return')
    ns=dict(vars(prior),OUT=OUT,PHOTON=run.PHOTON,GR=GR,BEFORE=run.fields.BEFORE,CAPS=dict(infinity=remaining))
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-infinity.py').write_text(s);ns['infinity']()
    result=read(OUT/'infinity/result.json')
    result.update(actual_current_GR_applied_to_photons_and_free_material=True,Phase151_direct_baseline_preserved=True,
        new_geometric_material_motion_returned_to_photons=False,uniform_native_derivative_certificate=False,
        uniform_contraction_bound=False,nonlinear_GR=False,full_source_error_enclosed=False)
    write(OUT/'infinity/result.json',result)
    print(json.dumps(result),flush=True)


def main():
    started=time.monotonic();repair.initialize();compact();infinity(180-(time.monotonic()-started))


def audit():
    started=time.monotonic();repair.initialize()
    files=[run.fields.FIELDS/'result.json',run.fields.METRIC/'result.json',run.PHOTON/'result.json',run.MATERIAL/'production.json',run.MATERIAL/'sources.json',GR/'result.json',OUT/'infinity/result.json']
    assert all(read(p)['passed'] for p in files)
    checks=0
    for name in ['plan.json','response-plan.json','response-execution-plan.json','readout-plan.json']:
        for p,h in read(OUT/name)['bindings'].items():
            path=OUT/'registered-response-producer.py' if name=='response-plan.json' and Path(p).name==Path(run.__file__).name else Path(p)
            assert sha(path)==h,p;checks+=1
    assert sha(run.__file__)==read(OUT/'response-execution-plan.json')['source_sha256']
    plan=read(OUT/'material-repair-execution.json');assert sha(repair.__file__)==plan['source_sha256']
    for p,h in plan['bindings'].items():assert sha(p)==h,p;checks+=1
    assert not read(repair.ORIGINAL/'steps-128-reference-128.json')['passed']
    assert sha(run.MATERIAL/'steps-64-reference-128.npz')==plan['inherited_coarse_sha256']
    assert max(read(run.MATERIAL/'steps-128-reference-128.json')['reverse_probe_8_16'])<.002
    fine=np.load(run.fields.FIELDS/'fields-128-g8.npz');field_controls={}
    for key,n,q in [('time',64,8),('quadrature',128,4)]:
        other=np.load(run.fields.FIELDS/f'fields-{n}-g{q}.npz')
        field_controls[key]={name:float(np.max(abs(fine[name]-other[name]))/max(np.max(abs(fine[name])),1e-300))
            for name in ['U','delta_mass_cm','delta_lambda']}
    assert max(field_controls['time'].values())<.02 and max(field_controls['quadrature'].values())<.002
    source_error=0.;port_error=0.;balance=0.;rows=[];gamma=1-1/np.sqrt(2)
    for n in [64,128]:
        p=np.load(run.PHOTON/f'steps-{n}-reference-128.npz');m=np.load(run.MATERIAL/f'steps-{n}-reference-128.npz')
        d=np.load(GR/f'source-{n}-reference-128.npz');s=np.load(run.MATERIAL/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-s[:,0],d['nonrest_trace_erg']+rest-(s[:,0]-s[:,1]-2*s[:,3]),
            d['nonrest_stress_erg']+rest-(s[:,0]-s[:,1]),d['pressure_volume_erg']-s[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(e)) for e in errors)/max(np.max(abs(s)),1e-300)))
        h=p['t'][-1]/n;times=h*(np.arange(n)[:,None]+[gamma,1.]).ravel()
        assert np.max(abs(times-p['accepted_angular_times']))<1e-18
        packet=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*p['accepted_angular_luminosity']
        aw=np.arange(1,8,2)/32;flux=packet@aw;acc=np.r_[0.,np.cumsum(flux,dtype=LD)]
        err=float(np.max(abs(acc[2*np.rint(p['t']/h).astype(int)]-p['radial_ports'][:,1,1]))/max(np.sum(abs(flux)),1e-300))
        port_error=max(port_error,err)
        balance=max(balance,float(np.max(abs(np.sum(m['history_scaled'],axis=2,dtype=LD)+m['discards_scaled']-m['ledgers_scaled'])/np.maximum(m['norms_scaled'],1.))))
        assert np.any(m['delta_scaled']) and np.any(p['delta_packet_scaled_occupation'])
        metric=np.load(run.fields.METRIC/f'corrected/metric-{n}-g8.npz');assert np.any(metric['delta_lambda'])
        rows.append(dict(steps=n,maximum_metric_lambda=float(np.max(abs(metric['delta_lambda']))),
            photon_energy_erg=float(np.sum(p['moments'][-1,0])),material_reference_energy_erg=float(np.sum(m['delta_scaled'][2])*AMP)))
    assert source_error<1e-12 and port_error<1e-12 and balance<1e-8
    norm_error=0.;baseline_error=0.;alpha=-float(d['K_cm']/d['M_cm'])
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        now=np.load(OUT/f'infinity/charge-{n}-a{a}-r{r}.npz');old=np.load(run.fields.BEFORE/f'direct-motion/charge-{n}-a{a}-r{r}.npz')
        dc=np.load(GR/f'wave-{n}-g8.npz')['free_scalar'];de=charge.G/C**4*now['arrived_parts_erg'][:,1]/float(d['M_cm'])
        inc=(dc+now['normalized_exterior_parts'][:,1]+(alpha+old['normalized'])*de)/(1-old['epsilon']-de)
        norm_error=max(norm_error,float(np.max(abs(inc-now['metric_return']))/max(np.max(abs(inc)),1e-300)))
        baseline_error=max(baseline_error,float(np.max(abs(now['previous']-old['normalized']))/max(np.max(abs(old['normalized'])),1e-300)))
        assert np.array_equal(now['normalized_exterior_parts'][:,0],old['normalized_exterior'])
        assert np.array_equal(now['arrived_parts_erg'][:,0],old['arrived_energy_erg'])
    assert norm_error<1e-12 and baseline_error<1e-14
    a,q,e,dc,ds,de=sp.symbols('a q e dc ds de');s=q*(1-e)-a*e
    assert sp.factor((s+dc+ds+a*(e+de))/(1-e-de)-q-(dc+ds+(a+q)*de)/(1-e-de))==0
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checks,source_identity_relative=source_error,
        angular_port_relative=port_error,material_balance_relative=balance,increment_normalization_relative=norm_error,
        previous_charge_relative=baseline_error,field_controls=field_controls,paths=rows,seconds=time.monotonic()-started,
        failed_fine_material_preserved=True,accepted_coarse_material_reused=True,
        actual_current_GR_photon_material_GR_return=True,nonzero_geometry_consumed=True,
        new_geometric_material_motion_returned_to_photons=False,coupled_fixed_point_verified=False,
        full_source_error_enclosed=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'verification.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    else:
        stage='readout' if action=='main' else action;receipt=OUT/f'{stage}-receipt.json';assert not receipt.exists()
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        def timeout(*_):raise TimeoutError('Registered readout/audit budget')
        signal.signal(signal.SIGALRM,timeout);signal.alarm(run.fields.CAPS[stage]);start=time.monotonic();cpu=time.process_time();error=None
        try:
            for p,h in read(OUT/'readout-plan.json')['bindings'].items():assert sha(p)==h,p
            globals()[action]()
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.alarm(0);write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=sha(__file__)))
