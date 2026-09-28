"""Counterexample candidate: native pressure at saved conserved deep states.

Keep the declared first-order advected-inventory closure explicit. This is a
source replacement on saved paths, not a new coupled EOS/fluid trajectory.
"""
from pathlib import Path
import json, signal, sys, time
import numpy as np
import sympy as sp
import verify_native_stage_energy_charge as charge

flow=charge.flow;chem=flow.previous.chem;C=flow.C;write=flow.write;sha=flow.sha
OUT=flow.OUT.parent/'def-native-conservative-eos-readout'
ROOT_TOL=1e-12


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(flow.__file__),Path(chem.__file__),Path(chem.old.__file__),
        Path(charge.__file__),flow.OUT/'coupled-64.npz',flow.OUT/'coupled-128.npz',
        charge.run.OUT/'source-64.npz',charge.run.OUT/'source-128.npz',charge.OUT/'result.json',
        flow.INPUT/'balanced-20.npz',chem.old.old.cold.CACHE/'stable/gas.so',
        chem.old.old.cold.CACHE/'libfree_eos_native_cold_stable.so',chem.old.LEVELS]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='91eeec8e9',
        claim='Evaluate native pressure at fixed saved density, specific internal energy and neutral H, then apply its full pressure/trace/stress difference to the actual retarded scalar charge.',
        decision='If native source replacement or represented interior force changes exceed2percent, correct the constitutive owner before another physical conclusion. Otherwise retain this as a measured, limited deep-source uncertainty.',
        reuse='Original64/128 paths,17 canonical snapshots each,19 deep cells; common initial states and6 pilot points reused. No new fluid steps, atmosphere EOS bank or physical grid.',
        closure='Other ionic/molecular inventories stay at each original native anchor. Retain the declared additive first-order advected-inventory u/P corrections. Native target u=saved_u+inventory_u*xi; pressure=native_P-inventory_P*xi. Density includes the installed initial volume shift. Do not subtract initial native pressure offset.',
        gates=dict(root_energy_relative=ROOT_TOL,root_iterations=6,root_logT_change=.02,root_reseed_pressure=1e-10,
            density_relative=1e-10,source_materiality=.02,force_materiality=.02,correction_cadence_relative_to_signal=.002,quadrature_relative_to_signal=.002),
        budget=dict(pilot_seconds=30,pilot_native_calls=500,production_seconds=120,production_native_calls=12000,readout_seconds=45,CPU_threads=1,memory_GB=3),
        forecast='Measure six representative native inverse states; allow2x per-state time for remaining states plus10s. Late state cost is unmeasured; hard call/time caps apply.',
        stop='Preserve failures. No automatic finer clocks, larger support or longer horizon. Material source error triggers constitutive repair, not a time-only fitted pressure patch.',
        limits='Saved-state comparison is not coupled EOS feedback, atmosphere certification, advected-composition exactness, uniform derivative/time error or nonlinear GR.',
        bindings={str(p):sha(p) for p in paths}))
    E,P,V,dp=sp.symbols('E P V dp')
    assert sp.expand((E-3*(P+dp))*V-(E-3*P)*V)==-3*dp*V
    assert sp.expand((E-P-dp)*V-(E-P)*V)==-dp*V
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        premise='Fixed conserved energy/density and declared acoustic momentum, isotropic material pressure.',
        trace_change='-3 delta_P V',metric_stress_change='-delta_P V',energy_change=0,
        limitations='No fluid feedback or exact GR-fluid primitive equivalence is inferred.'))


def inputs(model,steps):
    b=model.bulk;e=b.eos;z=np.load(flow.OUT/f'coupled-{steps}.npz')
    wanted=np.linspace(0,flow.old.END,17);ids=np.array([np.argmin(abs(z['snapshot_t']-t)) for t in wanted])
    assert np.max(abs(z['snapshot_t'][ids]-wanted))<1e-17
    rows=[]
    for it,k in enumerate(ids):
        model.Pi=z['snapshot_Pi'][k].copy();model.h=z['snapshot_h'][k].copy();model.j=z['snapshot_j'][k].copy()
        model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(float(z['snapshot_t'][k]))
        theta=z['snapshot_theta'][k];eta=z['snapshot_eta'][k];p,u,ut,*_=e.gas(theta,eta)
        rows.append(dict(t=float(wanted[it]),x=np.log1p(e.density_shift)+np.log1p(e.x),xi=e.xi.copy(),
            lt=np.log(e.base.d['T'])+e.theta0+theta,y=e.base.d['y0']*(1+eta),
            target=z['snapshot_u'][k]+e.inventory[1]*e.xi,p=p.copy(),u=z['snapshot_u'][k],ut=ut.copy()))
    return rows,z,ids


def root(native,inp,j,offset=0.):
    lt=float(inp['lt'][j]+offset);target=float(inp['target'][j]);start=native.ion.calls
    for iteration in range(6):
        state=native.state(float(inp['x'][j]),lt,float(inp['y'][j]));raw=state['raw']
        residual=float((raw[2]-target)/max(abs(target),1.))
        if abs(residual)<=ROOT_TOL:break
        # Interpolated ut accelerates inversion only. Acceptance uses the
        # actual native energy, never equilibrium native derivative columns.
        lt-=float((raw[2]-target)/inp['ut'][j])
        assert abs(lt-inp['lt'][j])<.02,('Native inverse support',j,lt)
    else:raise AssertionError(('Native conserved-energy inverse',j,residual))
    return dict(lt=lt,raw=raw.tolist(),energy_residual=residual,population_error=float(state['population_error']),
        fraction=state['fraction'].tolist(),affinity=float(state['affinity']),iterations=iteration+1,calls=native.ion.calls-start)


def native_samples(pilot=False):
    name='pilot' if pilot else 'production';assert not (OUT/f'{name}.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    cap=30 if pilot else 120;signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(cap)
    started=time.monotonic();model=flow.Coupled();b=model.bulk;e=b.eos
    data={steps:inputs(model,steps)[0] for steps in [64,128]}
    native=chem.old.Native(cap=500 if pilot else 12000)
    rows=[] if pilot else json.loads((OUT/'pilot-samples.json').read_text())
    done={(r['steps'],r['it'],r['cell']) for r in rows};count0=len(rows);durations=[];reseed=[]
    if not pilot:assert json.loads((OUT/'pilot.json').read_text())['eligible']
    try:
        for j in ([0,9,18] if pilot else range(b.n)):
            chem.setup(native,e.base.d,j)
            for steps in ([128] if pilot else [128,64]):
                for it in ([0,16] if pilot else range(17)):
                    key=(0 if it==0 else steps,it,j)
                    if key in done:continue
                    inp=data[steps][it];tick=time.monotonic();a=root(native,inp,j);durations.append(time.monotonic()-tick)
                    pressure=a['raw'][1]-float(e.inventory[0,j]*inp['xi'][j]);rho=float(e.base.d['rho'][j]*np.exp(inp['x'][j]))
                    density_error=abs(a['raw'][0]/rho-1);assert density_error<1e-10
                    a.update(steps=key[0],it=it,cell=j,t=inp['t'],x=float(inp['x'][j]),xi=float(inp['xi'][j]),
                        y=float(inp['y'][j]),target_u=float(inp['target'][j]),saved_u=float(inp['u'][j]),
                        old_lt=float(inp['lt'][j]),old_pressure=float(inp['p'][j]),pressure=pressure,
                        pressure_correction=pressure-float(inp['p'][j]),density_error=density_error)
                    rows.append(a);done.add(key)
                    if pilot:
                        independent=root(native,inp,j,.001)
                        err=abs(independent['raw'][1]/a['raw'][1]-1);assert err<1e-10
                        reseed.append(dict(cell=j,it=it,pressure_relative=err,independent=independent))
        write(OUT/f'{name}-samples.json',rows)
    except Exception as exc:
        write(OUT/f'{name}-partial.json',rows)
        write(OUT/f'{name}-failure.json',dict(error=repr(exc),completed=len(rows),native_calls=native.ion.calls,seconds=time.monotonic()-started));raise
    remaining=19*(17*2-1)-len(rows);forecast=2*max(durations)*remaining+10 if pilot else 0
    result=dict(classification='Counterexample candidate',passed=True,unique_states=len(rows),new_states=len(rows)-count0,
        maximum_energy_residual=max(abs(r['energy_residual']) for r in rows),maximum_density_error=max(r['density_error'] for r in rows),
        maximum_pressure_relative=max(abs(r['pressure_correction']/r['old_pressure']) for r in rows),
        maximum_logT_change=max(abs(r['lt']-r['old_lt']) for r in rows),native_calls=native.ion.calls,
        seconds=time.monotonic()-started,forecast_upper_seconds=forecast,eligible=forecast<120,reseed=reseed,
        coupled_EOS_evolution=False,uniform_EOS_error_bound=False)
    write(OUT/f'{name}.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def correction_data(steps,rows,volumes,stride=1,relative=False):
    d=dict(np.load(charge.run.OUT/f'source-{steps}.npz'));t=np.linspace(0,flow.old.END,17)
    by={(r['steps'],r['it'],r['cell']):r for r in rows}
    dp=np.array([[by[(0 if i==0 else steps,i,j)]['pressure_correction'] for j in range(19)] for i in range(17)])
    if relative:dp=dp-dp[0]
    pV=np.array([np.interp(d['t'],t[::stride],dp[::stride,j]) for j in range(19)]).T*volumes
    for key in charge.run.KEYS+['metric_stress_erg']:d[key]=np.zeros_like(d[key])
    d['inner_cumulative_energy_erg']=np.zeros_like(d['t']);d['outer_cumulative_energy_erg']=np.zeros_like(d['t'])
    d['pressure_volume_erg'][:,:19]=pV;d['nonrest_trace_erg'][:,:19]=-3*pV
    d['nonrest_stress_erg'][:,:19]=-pV;d['metric_stress_erg'][:,:19]=-pV
    return d,dp


def forces(model,rows):
    """Same accepted interior force owner; hold the shared port fixed.

    The last cell is excluded: changing the deep EOS also changes the join
    primitive/radiation, which requires a coupled solve rather than this test.
    """
    inp,z,ids=inputs(model,128);b=model.bulk;e=b.eos
    by={(r['steps'],r['it'],r['cell']):r for r in rows};old=[];delta=[]
    for it,k in enumerate(ids):
        model.Pi=z['snapshot_Pi'][k].copy();model.h=z['snapshot_h'][k].copy();model.j=z['snapshot_j'][k].copy()
        model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(inp[it]['t']);model.mflux=np.zeros(4)
        args=(z['snapshot_u'][k],z['snapshot_theta'][k],z['snapshot_eta'][k])
        before=model.material_rhs(*args)[1][:-1]
        dp=np.array([by[(0 if it==0 else 128,it,j)]['pressure_correction'] for j in range(19)])
        gas=e.gas
        def corrected(theta,eta):
            values=list(gas(theta,eta));values[0]=values[0]+dp;return values
        e.gas=corrected
        try:after=model.material_rhs(*args)[1][:-1]
        finally:e.gas=gas
        old.append(before);delta.append(after-before)
    old=np.array(old);delta=np.array(delta)
    norm=max(np.sum(abs(old),axis=1).max(),1e-300)
    np.savez_compressed(OUT/'interior-force.npz',t=[r['t'] for r in inp],old=old,delta=delta)
    return dict(maximum_L1_relative=float(np.sum(abs(delta),axis=1).max()/norm),
        maximum_cell_peak_relative=float(np.max(np.max(abs(delta),axis=0)/np.maximum(np.max(abs(old),axis=0),1e-300))),
        shared_port_held_fixed=True,last_cell_excluded=True,coupled_force_feedback=False)


def readout():
    assert not (OUT/'result.json').exists();started=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(45)
    rows=json.loads((OUT/'production-samples.json').read_text());assert len(rows)==627
    m=charge.gr.Response();waves={}
    for steps,stride,relative,order,label in [(128,1,False,8,'fine'),(64,1,False,8,'coarse'),
            (128,2,False,8,'nine'),(128,1,False,4,'g4'),(128,1,True,8,'dynamic')]:
        d,dp=correction_data(steps,rows,m.model.bulk.volume,stride,relative);wave=charge.read(m,d,order);waves[label]=wave
        np.savez_compressed(OUT/f'correction-{label}.npz',**wave,pressure=dp)
    base=np.load(charge.OUT/'wave-128-g8.npz');old=base['free_scalar'];norm=float(max(abs(old)))
    correction=waves['fine']['free_scalar'];new=old+correction;dyn=waves['dynamic']['free_scalar']
    cadence=float(max(abs(correction-waves['nine']['free_scalar']))/norm)
    quadrature=float(max(abs(correction-waves['g4']['free_scalar']))/norm)
    time_error=float(max(abs(correction[::2]-waves['coarse']['free_scalar']))/norm)
    # Apply through the original source owner as an independent linearity
    # check, rather than relying only on adding the two reported waveforms.
    d,_=correction_data(128,rows,m.model.bulk.volume);total=dict(np.load(charge.run.OUT/'source-128.npz'))
    for key in charge.run.KEYS+['metric_stress_erg']:total[key]=total[key]+d[key]
    applied=charge.read(m,total,8);linear=float(max(abs(applied['free_scalar']-new))/max(max(abs(new)),1e-300));assert linear<1e-10
    np.savez_compressed(OUT/'applied-charge.npz',**applied,previous=old,correction=correction)
    materiality=float(max(abs(correction))/norm);force=forces(m.model,rows)
    result=dict(classification='Counterexample candidate',passed=bool(cadence<.002 and quadrature<.002 and linear<1e-10),
        old_endpoint=float(old[-1]),native_pressure_endpoint=float(new[-1]),absolute_pressure_correction=float(correction[-1]),
        dynamic_pressure_correction=float(dyn[-1]),initial_pressure_offset_contribution=float(correction[-1]-dyn[-1]),
        maximum_source_change_relative_to_signal=materiality,correction_time_relative_to_signal=time_error,
        correction_cadence_relative_to_signal=cadence,correction_quadrature_relative_to_signal=quadrature,
        independent_application_relative=linear,interior_force=force,
        constitutive_repair_required=materiality>.02 or force['maximum_L1_relative']>.02,
        pressure_source_actually_applied=True,initial_offset_retained=True,
        energy_and_mass_sources_unchanged=True,atmosphere_EOS_audited=False,full_source_error_enclosed=False,
        coupled_EOS_evolution=False,coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,
        seconds=time.monotonic()-started)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1]
    if action in ['pilot','production']:native_samples(action=='pilot')
    else:globals()[action]()
