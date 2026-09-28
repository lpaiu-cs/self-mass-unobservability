"""Read the accepted pressure-work response through the same GR sources.

Counterexample candidate. Every increment comes from the175 coupled paths.
No diagnostic work is added to mass; physical closure is not inferred from
the conditional linear combination with the earlier selected charge.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,resource,shutil,sys,time
import numpy as np
import mpmath as mp
import sympy as sp
import complete_native_pressure_reciprocity as coupled
import read_native_incident_response as compact_owner
import read_native_incident_infinity as infinity_owner
import read_native_matched_mass as matched

OUT=Path('native-pressure-charge176-work');BEFORE=coupled.OUT
PHOTON=coupled.paths(1)[0];MATERIAL=OUT/'material';GR=OUT/'gr'
read,write,sha=coupled.read,coupled.write,coupled.sha
LD=np.longdouble;C=infinity_owner.C;G=infinity_owner.G
CAPS=dict(prepare=20,compact=180,packets=240,mass=60,audit=60)
TOTAL=570
KERNELS=Path('native-incident-infinity158-work/packet-kernels.npz')
GEOMETRY="""fields=np.array([m.fields(t)[2] for t in m.t])
        phi=fields[:,0]*AMP;lam=fields[:,2]*AMP;s=3*phi+lam"""


def active_geometry_source(source):
    replacements=[("phi=m.metric['delta_u'];lam=m.metric['delta_lambda'];s=3*phi+lam",
        GEOMETRY+"\n        assert not np.any(fields), 'Current175 response requires zero additional geometry'\n        np.savez_compressed(OUT/f'active-geometry-{steps}.npz',t=m.t,phi=phi,lam=lam)"),
        ('point=m.point(k);field=m.fields(t)[2];p0,delta,err=',
         'point=m.point(k);field=fields[k];p0,delta,err=')]
    for old,new in replacements:
        assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
    return source


def prepare():
    assert not OUT.exists();r=read(BEFORE/'result.json')
    assert r['passed'] and r['reciprocal_block_accepted'] and r['new_free_material_return_completed']
    block_receipt=BEFORE/('closure-block_resume-receipt.json' if (BEFORE/'block-budget-reassessment.json').exists() else 'closure-block-receipt.json')
    assert read(block_receipt)['error'] is None
    assert read(BEFORE/'normalization.json')['factor']==1.
    OUT.mkdir();MATERIAL.mkdir();GR.mkdir()
    for old,name in [('.phase176-metric-preflight.json','metric-preflight.json'),
                     ('.phase176-pre-metric-producer.py','pre-metric-producer.py')]:
        shutil.copyfile(old,OUT/name)
    preflight=read(OUT/'metric-preflight.json')
    assert preflight['no_evolution_steps'] and preflight['seconds']+CAPS['prepare']<30
    write(OUT/'metric-preflight-receipt.json',dict(action='metric_preflight',seconds=preflight['seconds'],error=None))
    # Exercise the collector with a time-varying active field and a stale cache.
    fake=SimpleNamespace(t=np.array([0.,1.,2.]),metric={'delta_u':np.full((3,2),99.)},
        fields=lambda t:(None,None,np.array([[t,t+1],[0.,0.],[2*t,3*t]])))
    ns=dict(np=np,m=fake,AMP=.25)
    exec(GEOMETRY.replace('\n        ','\n'),ns)
    assert np.array_equal(ns['phi'],.25*np.array([[0,1],[1,2],[2,3]]))
    assert np.array_equal(ns['lam'],.25*np.array([[0,0],[2,3],[4,6]]))
    files=[Path(__file__),Path(coupled.__file__),Path(coupled.prior.__file__),Path(compact_owner.__file__),
        Path(infinity_owner.__file__),Path(matched.__file__),BEFORE/'result.json',BEFORE/'block-result.json',
        BEFORE/'material-block-plan.json',block_receipt,BEFORE/'normalization.json',
        KERNELS,infinity_owner.BACKGROUND,matched.OUT/'final-result.json',matched.OUT/'charge-parts.npz',matched.OUT/'mass-fine.npz',
        OUT/'metric-preflight.json',OUT/'pre-metric-producer.py']
    for p in [coupled.paths(1)[1]/'production.json']+[coupled.paths(1)[1]/f'steps-{n}-reference-128.npz' for n in [64,128]]:
        q=MATERIAL/p.name;shutil.copyfile(p,q);assert sha(p)==sha(q);files.extend([p,q])
    files.extend(PHOTON/f'steps-{n}-reference-128.npz' for n in [64,128])
    write(OUT/'readout-plan.json',dict(classification='Counterexample candidate',checkpoint='8e6ce3b83',
        question='How does the actual pressure-work-corrected coupled response change the selected charge and its same-interface mass residual?',
        input='Both accepted175 photons and returned free-material histories, with their actual Radau stage times, weights, collision transfers and ports. No evolution replay or diagnostic0.38538erg addition.',
        method='Reuse the conservative pressure/trace source export, characteristic compact GR, fixed-ray signed packet kernel and exact center-J constraint. Combine source, arrived energy and homogeneous mass separately in the same high-precision normalization.',
        geometry_repair='The source collector must call the active fields at EACH canonical time, then use those same arrays for geometry subtraction and pressure. Constructor metric arrays are stale after the driver changes. A stage-local driver view is not a full history. Preserve the old producer and actual-object preflight; assert the175 additional geometry is zero.',
        mass='J is independent of the scalar argument in the existing constraint. Evaluate its source-dependent boundary and independently reconstruct it from cell lapse-weighted energies plus the same radial ports; zero scalar placeholders are not used as metric solutions.',
        stages='Read all saved local Radau substeps. At each canonical output endpoint only the matching final closing-stage age is exactly zero; verify the stored roundoff mismatch below1e-18s. Earlier positive ages remain unchanged.',
        comparison='Report the new computed increment and conditional combined selected sign separately. The latter requires the declared first-variation superposition; generated-GR feedback, complete boundary matching and full physical error control remain open. Do not inherit the165 sign enclosure for the changed transport.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,independent_GR=1e-9,source_identity=1e-12,port=1e-12,mass_identity=1e-12,normalization=1e-12),
        budgets=CAPS,total_action_seconds=TOTAL,preflight_seconds=preflight['seconds'],CPU_threads=1,virtual_GiB=3,
        forecast='Previous same-owner compact readout21.21s and packet readout15.55s. New Radau ages may miss the cache; measure three actual uncached ages and require2x worst age-group cost times remaining groups plus5s to fit240s. No extra rays/orders/clocks on failure.',
        stop='Any failed source, pressure, time, spatial, mass-identity or cost gate stops this readout. Preserve the result and do not expand the run automatically.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,bindings={str(p):sha(p) for p in files}))
    ds,de,dk,alpha,q,D=sp.symbols('ds de dk alpha q D',nonzero=True)
    x=sp.symbols('x');changed=(q*D+x*ds+alpha*x*(de-dk))/(D+x*(dk-de))
    assert sp.simplify(sp.diff(changed,x).subs(x,0)-(ds+(alpha+q)*(de-dk))/D)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='First variation of the same mass-normalized scalar charge; does not prove response superposition at branch corners or physical GR closure.'))


def compact():
    def setup():
        coupled.prior.initialize();coupled.run.c.base.Material=coupled.run.c.Material
    fn=FunctionType(compact_owner.compact.__code__,dict(compact_owner.compact.__globals__,
        OUT=OUT,PHOTON=PHOTON,MATERIAL=MATERIAL,GR=GR,fixed=SimpleNamespace(initialize=setup,write=write)))
    owner=compact_owner.run.base.old.aligned.base;original=owner.source
    owner.source=active_geometry_source(original)
    try:fn()
    finally:owner.source=original


def emission(n):
    path=PHOTON/f'steps-{n}-reference-128.npz';tt,energy,error=coupled.run.packets(path)
    with np.load(path) as p:
        clock=p['t'].copy();edges=p['actual_step_edges'];ids=np.array([np.argmin(abs(edges-t)) for t in clock])
        assert np.max(abs(edges[ids]-clock))<1e-18
        ages=[]
        for now,index in zip(clock,ids):
            age=now-tt[:2*index]
            if index:
                assert abs(age[-1])<1e-18 and np.all(age[:-1]>0)
                age[-1]=0.
            ages.append(age)
    return clock,energy,ages,error


def packets():
    assert read(GR/'result.json')['passed'];started=time.monotonic()
    infinity_owner.prior.initialize();infinity_owner.OUT=OUT
    packet,model,kernels,inverse_error=infinity_owner.packet_factory()
    with np.load(KERNELS) as cache:kernels.update({tuple(k):v for k,v in zip(cache['keys'],cache['values'])})
    cached=set(kernels);emitted={n:emission(n) for n in [64,128]};paths=[(128,8,8),(128,4,8),(128,8,4),(64,8,8)]
    required=set((float(age),a,r) for n,a,r in paths for group in emitted[n][2] for age in group if age>0)
    groups={age for age,a,r in required if (age,a,r) not in kernels};costs=[]
    if groups:
        ordered=sorted(groups)
        for age in dict.fromkeys([ordered[0],ordered[len(ordered)//2],ordered[-1]]):
            mark=time.monotonic()
            for a,r in [(8,8),(4,8),(8,4)]:packet(age,a,r)
            costs.append(time.monotonic()-mark)
    remaining={age for age,a,r in required if (age,a,r) not in kernels}
    upper=2*len(remaining)*max(costs,default=0.)+5
    write(OUT/'packet-admission.json',dict(classification='Counterexample candidate',point_seconds=costs,
        actual_distinct_required_keys=len(required),cached_required_keys=len(required&cached),remaining_age_groups=len(remaining),
        upper_remaining_seconds=upper,eligible=time.monotonic()-started+upper<CAPS['packets']))
    assert time.monotonic()-started+upper<CAPS['packets'],('Actual packet cost admission',upper)
    results={};rows=[]
    for n,a,r in paths:
        clock,energy,ages,port=emitted[n];assert np.max(abs(clock-model.t))<1e-18
        pieces=[]
        for group in ages:
            value=np.zeros((2,4),LD)
            for age,e in zip(group,energy):value+=packet(float(age),a,r)*e
            pieces.append(np.sum(value,axis=1,dtype=LD))
        parts=np.array(pieces);results[n,a,r]=parts
        np.savez_compressed(OUT/f'exterior-{n}-a{a}-r{r}.npz',t=clock,scalar=parts[:,0],arrived_energy_erg=parts[:,1])
        rows.append(dict(steps=n,angular=a,radial=r,actual_stage_count=len(energy),port_relative=port))
    fine=results[128,8,8];norm=np.maximum(np.max(abs(fine),axis=0),LD('1e-290'))
    errors={name:(np.max(abs(results[key]-fine),axis=0)/norm).astype(float).tolist()
        for name,key in [('time',(64,8,8)),('angular',(128,4,8)),('radial',(128,8,4))]}
    result=dict(classification='Counterexample candidate',passed=max(errors['time'])<.02 and max(errors['angular']+errors['radial'])<.002,
        controls=errors,rows=rows,new_required_kernel_keys=len(required-cached),delay_inverse_seconds=inverse_error(),
        seconds=time.monotonic()-started,final_charge_conclusion='unadjudicated')
    # Save the new reusable geometry only; the original cache stays frozen.
    new=required-cached;np.savez_compressed(OUT/'new-packet-kernels.npz',keys=np.array(sorted(new)),values=np.array([kernels[k] for k in sorted(new)]))
    write(OUT/'packets-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def mass():
    assert read(GR/'result.json')['passed'];coupled.prior.initialize()
    model=compact_owner.charge.gr.Response();rows=[];values={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));model.setup(d,q)
        zero=np.zeros((len(d['t']),len(model.tr)))
        z=compact_owner.charge.independent.centers(model,d,dict(delta_phi=zero,delta_Phi=zero),q)
        port=z['J'][:,-1].astype(LD)*LD(model.tz['lapse'][-1])/np.sqrt(LD(model.tz['b'][-1]))+LD(G)/LD(C)**4*d['outer_cumulative_energy_erg']
        mean=coupled.prior.original.previous.weights(model,d,q).astype(LD)
        energy=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(C)**2+d['gas_nonrest_energy_erg']+d['photon_energy_erg']
        direct=np.sum(energy*mean,axis=1,dtype=LD)-d['inner_cumulative_energy_erg']+d['outer_cumulative_energy_erg']
        identity=float(np.max(abs(port-direct*LD(G)/LD(C)**4))/max(np.max(abs(port)),LD('1e-290')))
        assert identity<1e-12,identity;values[n,q]=port
        np.savez_compressed(OUT/f'mass-{n}-g{q}.npz',t=d['t'],source_constraint_port_cm=port,energy_and_ports_erg=direct,mean_lapse=mean)
        rows.append(dict(steps=n,order=q,identity_relative=identity,endpoint_cm=float(port[-1]),endpoint_erg=float(direct[-1])))
    norm=max(np.max(abs(values[128,8])),LD('1e-290'))
    controls={name:float(np.max(abs(values[key]-values[128,8]))/norm) for name,key in [('time',(64,8)),('quadrature',(128,4))]}
    result=dict(classification='Counterexample candidate',passed=controls['time']<.02 and controls['quadrature']<.002,
        rows=rows,controls=controls,no_work_diagnostic_added=True,actual_mass_flux_balance_verified=False)
    write(OUT/'mass-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def audit():
    assert all(read(p)['passed'] for p in [GR/'result.json',OUT/'packets-result.json',OUT/'mass-result.json',OUT/'symbolic.json'])
    source_error=0.
    for n in [64,128]:
        d=np.load(GR/f'source-{n}-reference-128.npz');s=np.load(MATERIAL/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(C)**2;total=rest+d['gas_nonrest_energy_erg']
        checks=[total-s[:,0],d['nonrest_trace_erg']+rest-(s[:,0]-s[:,1]-2*s[:,3]),
            d['nonrest_stress_erg']+rest-(s[:,0]-s[:,1]),d['pressure_volume_erg']-s[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        source_error=max(source_error,float(max(np.max(abs(x)) for x in checks)/max(np.max(abs(s)),LD('1e-290'))))
    assert source_error<1e-12,source_error
    bg=np.load(matched.OUT/'charge-parts.npz');old=read(matched.OUT/'final-result.json')
    constructor=compact_owner.run.base.drive.METRIC/'corrected/metric-128-g8.npz'
    cache=np.load(constructor);geometry_rows=[];geometry_files=[constructor]
    for label,folder in [('self_GR157',infinity_owner.SELF),('mixed_GR162',matched.R)]:
        path=folder/'metric/metric-128-g8.npz';normal=folder/'normalization.json'
        g=np.load(path);factor=LD(read(normal)['factor']);ids=[np.argmin(abs(g['t']-t)) for t in cache['t']]
        assert np.max(abs(g['t'][ids]-cache['t']))<1e-18
        row=dict(component=label)
        for key in ['delta_u','delta_lambda']:
            active=g[key][ids].astype(LD)*factor;difference=cache[key].astype(LD)-active
            row[key]=dict(active_max=float(np.max(abs(active))),cached_max=float(np.max(abs(cache[key]))),
                difference_max=float(np.max(abs(difference))),relative_to_active=float(np.max(abs(difference))/max(np.max(abs(active)),LD('1e-290'))))
        geometry_rows.append(row);geometry_files.extend([path,normal])
    write(OUT/'previous-source-geometry.json',dict(classification='Counterexample candidate',rows=geometry_rows,
        scope='Stored active157/162 driver histories at canonical output times versus the155 constructor cache read by their unchanged shared source exporter. This identifies an upstream mismatch; it neither recomputes their charge nor validates a combined solution.',
        bindings={str(p):sha(p) for p in geometry_files}))
    base_mass=np.load(matched.OUT/'mass-fine.npz');clock=bg['t'];model=np.load(GR/'source-128-reference-128.npz')
    compact=np.load(GR/'wave-128-g8.npz');exterior=np.load(OUT/'exterior-128-a8-r8.npz');mass=np.load(OUT/'mass-128-g8.npz')
    for p in [compact,exterior,mass,base_mass]:assert np.max(abs(p['t']-clock))<1e-18
    mp.mp.dps=110;B=lambda v:mp.mpf(str(v));rows=[];maximum=0.;M=B(model['M_cm'])
    for j,t in enumerate(clock):
        alpha=B(bg['alpha']);D0=B(bg['old_denominator'][j]);kb=B(bg['background_kappa'][j]);q0=B(bg['background_charge'][j])
        D=D0+kb;qb=q0-(alpha+q0)*kb/D
        ds=B(compact['free_scalar'][j])+B(exterior['scalar'][j]);de=B(G)/B(C)**4/M*B(exterior['arrived_energy_erg'][j]);dk=B(mass['source_constraint_port_cm'][j])/M
        value=(ds+(alpha+qb)*(de-dk))/D
        independent=((qb*D+ds+alpha*(de-dk))/(D+dk-de)-qb)*(D+dk-de)/D
        error=float(abs(value-independent)/max(abs(value),mp.mpf('1e-290')));maximum=max(maximum,error)
        selected=B(old['selected_decimal'][j]);updated=selected+value
        oldport=B(base_mass['matched_homogeneous_port_cm'][j]);newport=oldport+B(mass['source_constraint_port_cm'][j])
        rows.append(dict(t=float(t),compact=str(B(compact['free_scalar'][j])/D),exterior=str(B(exterior['scalar'][j])/D),
            arrived_mass=str((alpha+qb)*de/D),homogeneous_mass=str(-(alpha+qb)*dk/D),
            computed_charge_increment=str(value),previous_selected=str(selected),conditional_selected=str(updated),
            previous_mass_port_erg=str(oldport*B(C)**4/B(G)),conditional_mass_port_erg=str(newport*B(C)**4/B(G))))
    assert maximum<1e-12,maximum
    end=rows[-1];selected=B(end['conditional_selected']);before=B(end['previous_selected'])
    result=dict(classification='Counterexample candidate',passed=True,source_identity_relative=source_error,
        normalization_independent_relative=maximum,rows=rows,endpoint=end,
        conditional_selected_sign_preserved=bool(selected*before>0),computed_increment_over_previous=float(B(end['computed_charge_increment'])/abs(before)),
        same_accepted_photon_material_sources_used=True,actual_Radau_packets_used=True,no_diagnostic_work_added=True,
        conditional_linear_readout_completed=True,generated_GR_field_returned_to_evolution=False,
        full_branch_superposition_certified=False,changed_transport_full_error_enclosed=False,
        source_geometry_matches_active_driver=True,previous_selected_source_geometry_revalidated=False,
        actual_mass_flux_balance_verified=False,exact_ADM_conservation_verified=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False,
        scope='Actual corrected finite response and its source/energy/mass readout. The combined selected number is conditional on first-variation superposition, revalidation of the earlier source geometry, and incomplete physical closure; the previous165 error interval is not transferred to it.')
    write(OUT/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));infinity_owner.incident.native.deadline(CAPS[action])
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'readout-plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
