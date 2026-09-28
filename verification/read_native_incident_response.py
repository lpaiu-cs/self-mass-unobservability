"""Read the actual externally driven material-mediated compact GR return.

Counterexample candidate. This does not include the full direct incoming pulse
scattering at null infinity, exterior photon response, or a coupled fixed point.
"""
from pathlib import Path
from types import SimpleNamespace
import json,resource,signal,sys,time
import numpy as np
import solve_native_incident_material as fixed
import verify_native_stage_energy_charge as charge

run=fixed.lift
OUT=run.OUT;PHOTON=run.PHOTON;MATERIAL=fixed.MATERIAL;GR=run.base.GR
read,write,sha=run.read,run.write,run.sha;AMP=run.AMP;LD=np.longdouble
NO_ALARM=SimpleNamespace(alarm=lambda n:0,signal=signal.signal,SIGALRM=signal.SIGALRM)


def prepare():
    assert read(PHOTON/'result.json')['passed'] and read(MATERIAL/'production.json')['passed']
    paths=[Path(__file__),Path(fixed.__file__),Path(run.__file__),Path(run.base.__file__),Path(run.base.drive.__file__),
           OUT/'budget-reassessment.json',PHOTON/'result.json',MATERIAL/'production.json']
    paths += [folder/f'steps-{n}-reference-128.npz' for folder in [PHOTON,MATERIAL] for n in [64,128]]
    write(OUT/'readout-plan.json',dict(classification='Counterexample candidate',
        claim='Export the actual externally driven photon/thermal/H and free-material increments, remove the initial canonical geometry already in the scalar operator, and compute their compact retarded scalar return.',
        distinction='The primary incident scalar and its first potential return are stored separately. Here free_scalar means the material source response of the existing compact Green operator, not the complete outgoing waveform of the initially exterior scalar packet.',
        controls='Original17 canonical output times;64/128 evolution clocks;4/8 spatial quadrature; independent direct-trace readout; pressure probe and conserved source/port identities.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,independent_GR=1e-9,conservation=1e-8,identities=1e-12),
        scope='One physical incident-input sweep on the retained time-dependent background. Mechanical response has not yet returned into photons. No complete native Jacobian, nonlinear GR, stationary transfer function, orbital signal, static EFT nonabsorption or complete null-infinity charge.',
        budget_seconds=180,audit_seconds=60,new_background_steps=0,new_native_calls=0,
        stop='Keep original failures, waveform, physical amplitude, gates, grid and horizon. No automatic cadence/grid or runtime expansion.',
        bindings={str(p):sha(p) for p in paths}))


def compact():
    plan=read(OUT/'readout-plan.json')
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    fixed.initialize();owner=run.base.old.aligned.base;prior=run.base.drive.native.prior
    s=owner.source.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]')
    s=s.replace('background=compare(histories[2],histories[1]),','').replace(',stress_background=compare(allstress[2],allstress[1])','')
    # JSON diagnostics use ordinary scalars; retained conserved arrays stay LD.
    s=s.replace('.tolist()', '.astype(float).tolist()')
    s=prior.replace(s,"base=dict(np.load(wave.base.OUT/f'source-{ref}.npz'))",'base=template()')
    ns=dict(owner.source_scope,OUT=MATERIAL,GR=GR,Material=run.base.Material,signal=NO_ALARM,
        template=prior.template,photons=SimpleNamespace(path=lambda n,r:PHOTON/f'steps-{n}-reference-{r}.npz'),
        material_path=lambda n,r:MATERIAL/f'steps-{n}-reference-{r}.npz',write=fixed.write)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-incident-source.py').write_text(s);ns['sources']()
    assert read(MATERIAL/'sources.json')['passed']
    model=charge.gr.Response();waves={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));waves[n,q]=charge.read(model,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**waves[n,q])
    fine=waves[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-290)
    direct,coordinate=charge.independent.direct(model,dict(np.load(GR/'source-128-reference-128.npz')),8)
    controls=dict(time=float(np.max(abs(fine-waves[64,8]['free_scalar']))/norm),
        quadrature=float(np.max(abs(fine-waves[128,4]['free_scalar']))/norm),
        independent=float(abs(direct-waves[128,8]['direct_scalar'][-1])/max(np.max(abs(waves[128,8]['direct_scalar'])),1e-290)))
    result=dict(classification='Counterexample candidate',controls=controls,compact_return_endpoint=float(fine[-1]),
        compact_return_maximum=float(norm),incident_amplitude=run.base.drive.ETA,
        endpoint_over_incident_amplitude=float(fine[-1]/run.base.drive.ETA),inverse_radius_residual=float(coordinate),
        actual_external_input_consumed=True,actual_free_material_response_evolved=True,
        passed=controls['time']<.02 and controls['quadrature']<.002 and controls['independent']<1e-9,
        full_null_infinity_charge=False,complete_reciprocal_fixed_point=False,full_goal_complete=False)
    write(GR/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def audit():
    paths=[run.base.drive.METRIC/'result.json',OUT/'source-check.json',OUT/'lift-check.json',PHOTON/'result.json',
           MATERIAL/'production.json',MATERIAL/'sources.json',GR/'result.json']
    assert all(read(p)['passed'] for p in paths)
    fixed.initialize();maximum_source=0.;maximum_port=0.;maximum_balance=0.;rows=[];gamma=1-1/np.sqrt(2)
    for n in [64,128]:
        p=np.load(PHOTON/f'steps-{n}-reference-128.npz');m=np.load(MATERIAL/f'steps-{n}-reference-128.npz')
        d=np.load(GR/f'source-{n}-reference-128.npz');s=np.load(MATERIAL/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(run.base.C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-s[:,0],d['nonrest_trace_erg']+rest-(s[:,0]-s[:,1]-2*s[:,3]),
            d['nonrest_stress_erg']+rest-(s[:,0]-s[:,1]),d['pressure_volume_erg']-s[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        maximum_source=max(maximum_source,float(max(np.max(abs(e)) for e in errors)/max(np.max(abs(s)),1e-290)))
        h=p['t'][-1]/n;times=h*(np.arange(n)[:,None]+[gamma,1.]).ravel()
        assert np.max(abs(times-p['accepted_angular_times']))<1e-18
        packet=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*p['accepted_angular_luminosity']
        flux=packet@(np.arange(1,8,2)/32);acc=np.r_[0.,np.cumsum(flux,dtype=LD)]
        maximum_port=max(maximum_port,float(np.max(abs(acc[2*np.rint(p['t']/h).astype(int)]-p['radial_ports'][:,1,1]))/max(np.sum(abs(flux)),1e-290)))
        maximum_balance=max(maximum_balance,float(np.max(abs(np.sum(m['history_scaled'],axis=2,dtype=LD)+m['discards_scaled']-m['ledgers_scaled'])/np.maximum(m['norms_scaled'],1.))))
        pilot=np.load(PHOTON/f'pilot-{n}.npz')
        for key in ['t','moments','photon_history_scaled_occupation','material_history','radial_ports']:
            assert np.array_equal(p[key][:len(pilot[key])],pilot[key]),key
        assert np.any(p['photon_history_scaled_occupation']) and np.any(m['delta_scaled'])
        rows.append(dict(steps=n,photon_reference_energy_erg=float(np.sum(p['moments'][-1,0])),
            material_reference_energy_erg=float(np.sum(m['delta_scaled'][2])*AMP),
            outgoing_signed_energy_erg=float(p['radial_ports'][-1,1,1]),free_material_substeps=int(m['substeps'])))
    assert maximum_source<1e-12 and maximum_port<1e-12 and maximum_balance<1e-8
    write(OUT/'audit.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        source_identity_relative=maximum_source,angular_port_relative=maximum_port,material_conservation=maximum_balance,
        accepted_photon_prefixes_reused=True,original_pilot_failure_preserved=not read(OUT/'photons/pilot.json')['passed'],
        claim='Actual compact matter-mediated response to the external scalar packet, within the declared retained finite model.',
        full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','compact','audit'];resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    run.base.drive.native.deadline(180 if action=='compact' else 60)
    start=time.monotonic();cpu=time.process_time();error=None
    try:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        p=OUT/f'readout-{action}-receipt.json';assert not p.exists()
        write(p,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
