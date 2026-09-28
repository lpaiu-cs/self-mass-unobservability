"""Extend the identical causal GR history without recomputing accepted times."""
from pathlib import Path
import gc,json,os,resource,sys,time
import numpy as np
from scipy.interpolate import PPoly
import read_full_captured_history as source
import return_complete_history_gr as previous

OUT=Path('native-retarded-extension248-work');INPUT=source.OUT;OLD=previous.OUT
base=source.prior;read,write,sha,bind,LD=source.read,source.write,source.sha,source.bind,source.LD
FIELDS=previous.FIELDS
CAPS=dict(check=1800,prepare=600,field1288=7200,field648=7200,field1284=7200,collect=300,audit=900)


def source_prefix(new,old):
    count=len(old['t']);geometry=len(old['geometry_times'])
    for key,value in old.items():
        v=new[key]
        if key=='t' or key in base.KEYS:v=v[:count]
        elif key.startswith('state_coeff_'):v=v[:,:count-1]
        elif key=='geometry_times':v=v[:geometry]
        elif key.startswith('geometry_coeff_'):v=v[:,:geometry-1]
        assert np.array_equal(v,value),('Source prefix changed',key)


class Response(base.Response):
    setup=bind(previous.prior.Response.setup,INPUT=INPUT)


def at(m,poly,times):
    saved=m.t;m.t=times
    try:return m.propagate(poly)
    finally:m.t=saved


def extend(m,d,old,keep,order):
    """Reuse full prefix values; evaluate the tail and two overlap samples."""
    start=time.monotonic();m.setup(d,order);times=m.t
    assert 2<=keep<len(times) and np.array_equal(times[:keep],old['t'][:keep])
    assert len(old['potential_iterations'])==1,'Prefix has a different potential iteration'
    overlap=2;query=times[keep-overlap:]
    free=at(m,m.source,query)
    free_U=np.concatenate([old['direct_and_mass_stress_U'][:keep],free[0][overlap:]])
    potential=-m.dx*m.z['V']*free_U[:,:-1][:,m.ids]
    potential_poly=base.base.gr.base.flow.green.polynomial(times,potential)
    add=at(m,potential_poly,query)
    potential_U=np.concatenate([old['potential_U'][:keep],add[0][overlap:]])
    U=free_U+potential_U
    Ut=np.concatenate([old['U_t'][:keep],(free[1]+add[1])[overlap:]])
    Ux=np.concatenate([old['U_x'][:keep],(free[2]+add[2])[overlap:]])
    relative=lambda v,reference:float(np.max(abs(v))/max(np.max(abs(reference)),1e-290))
    controls={k:relative(v[:overlap]-old[k][keep-overlap:keep],old[k][:keep]) for k,v in
        [('direct_and_mass_stress_U',free[0]),('potential_U',add[0]),('U',free[0]+add[0]),('U_t',free[1]+add[1]),('U_x',free[2]+add[2])]}
    assert max(controls.values())<1e-9,controls
    eta=float(base.C*times[-1]/2*np.sum(m.dx*abs(m.z['V']),dtype=LD));assert eta<.01
    err=relative(potential_U,U);assert err<1e-8,('Original potential iteration gate',err)
    z=m.tz;r=m.tr;f=U/r;fr=Ux/(z['lapse']*np.sqrt(z['b'])*r)-U/r**2
    J=np.array([np.interp(m.tx,m.x,row) for row in m.J]);dm=r*r*z['b']*z['Phi']*f+J
    dl=dm/(r*z['b']);volume=3*z['alpha']*f+dl
    result=dict(t=times,radius_E=r,U=U,U_t=Ut,U_x=Ux,delta_phi=f,delta_Phi=fr,delta_mass_cm=dm,
        delta_lambda=dl,delta_log_proper_volume=volume,
        delta_gas_energy_geom=-(z['Eg']+z['Pg'])*volume,delta_gas_radial_pressure_geom=-z['Kg']*volume,
        delta_photon_energy_geom=-4*z['alpha']*z['Er']*f-(z['Er']+z['Pr'])*dl,
        delta_photon_radial_pressure_geom=-4*z['alpha']*z['Pr']*f-(3*z['Pr']-z['R4'])*dl,
        direct_and_mass_stress_U=free_U,potential_U=potential_U,J_center_interpolated=J,potential_iterations=np.array([err]))
    # Every accepted prefix field must be reproduced, including mass/volume.
    for k,v in old.items():
        if k not in ['t','radius_E','potential_iterations']:
            assert relative(result[k][:keep]-v[:keep],v[:keep])<1e-12,('Prefix field changed',k)
    direct=at(m,m.direct,query)[0][:,-1];mass=float(d['M_cm'])
    row=dict(classification='Counterexample candidate',steps=int(d['original_clock']),order=order,
        seconds=time.monotonic()-start,reused_output_times=keep,new_output_times=len(times)-keep,
        independently_recomputed_overlap=overlap,overlap_relative=controls,
        compact_potential_contraction=eta,potential_iteration_relative=[err],
        endpoint_direct=-float(direct[-1])/mass,endpoint_compact_with_metric=-float(U[-1,-1])/mass,
        compact_potential_change=-float(potential_U[-1,-1])/mass,
        same_causal_prefix_reused=True,physical_steps=0,full_GR_feedback=False,final_charge_solved=False)
    return result,row


def check():
    assert read(OLD/'metric-result.json')['passed'];folder=OUT/'check';folder.mkdir(parents=True,exist_ok=True)
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(folder/part).mkdir(parents=True,exist_ok=True)
    for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=folder/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    bind(base.endpoint.initialize,OUT=folder)();rows=[]
    # Replay only the last three existing times against three independently
    # completed full fields. The other 532 times are reused, not rerun.
    for n,q in FIELDS.values():
        d=dict(np.load(OLD/f'gr/source-{n}.npz'));old=dict(np.load(OLD/f'gr/fields-{n}-g{q}.npz'))
        m=previous.Response();actual,row=extend(m,d,old,len(old['t'])-3,q)
        errors={k:float(np.max(abs(v-old[k]))/max(np.max(abs(old[k])),1e-290)) for k,v in actual.items() if k!='potential_iterations'}
        assert max(errors.values())<1e-9,errors;row['independent_full_field_reproduction']=errors;rows.append(row)
        del m,d,old,actual;gc.collect()
    raw=dict(np.load(base.OUT/'gr/source-64.npz'));source_prefix(raw,raw)
    bad=dict(raw);bad['state_coeff_baryon_g']=raw['state_coeff_baryon_g'].copy()
    bad['state_coeff_baryon_g'].flat[0]=np.nextafter(bad['state_coeff_baryon_g'].flat[0],LD(np.inf))
    try:source_prefix(bad,raw)
    except AssertionError:pass
    else:raise AssertionError('Changed causal source accepted')
    write(OUT/'check-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        changed_causal_source_rejected=True,scope='Exact saved finite operator extension, not a continuum error certificate.',
        bindings={str(p):sha(p) for p in [Path(__file__),OLD/'metric-result.json']+
                  [OLD/f'gr/fields-{n}-g{q}{ext}' for n,q in FIELDS.values() for ext in ['.npz','.json']]}))


def prepare():
    assert read(OUT/'check-result.json')['passed'];assert read(INPUT/'controller-status.json')['state']=='completed'
    r=read(INPUT/'result.json');assert r['representation_controls_passed'] and r['source_time_passed']
    assert read(INPUT/'driver-polynomial-audit.json')['passed'];assert not (OUT/'plan.json').exists()
    files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/part).mkdir(parents=True,exist_ok=True)
    for p in list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);files += [p,dst]
    candidates=[0.]
    for n in [64,128]:
        z=np.load(source.saved(n));candidates+=list(z['joint_stage_times'])+list(z['actual_step_edges']);files.append(source.saved(n))
    clock=[]
    for t in sorted(candidates):
        if not clock or t-clock[-1]>1e-18:clock.append(t)
    clock=np.array(clock);counts=[]
    for n in [64,128]:
        raw=dict(np.load(INPUT/f'gr/source-{n}.npz'));old=dict(np.load(base.OUT/f'gr/source-{n}.npz'))
        source_prefix(raw,old);knots,co=base.coefficients(raw)
        d=dict(raw,t=clock,original_clock=np.array(n));d.update({k:PPoly(np.asarray(v[::-1],float),knots)(clock) for k,v in co.items()})
        ids=np.array([np.argmin(abs(clock-t)) for t in raw['t']]);assert np.max(abs(clock[ids]-raw['t']))<1e-18
        for k in base.KEYS:d[k][ids]=raw[k]
        previous_source=np.load(OLD/f'gr/source-{n}.npz');count=len(previous_source['t'])
        assert np.array_equal(clock[:count],previous_source['t']);counts.append(count)
        for k in base.KEYS:assert np.array_equal(d[k][:count],previous_source[k]),('Stage source changed',n,k)
        np.savez_compressed(OUT/f'gr/source-{n}.npz',**d)
        files += [INPUT/f'gr/source-{n}.npz',base.OUT/f'gr/source-{n}.npz',OLD/f'gr/source-{n}.npz']
    assert counts[0]==counts[1]
    for action,(n,q) in FIELDS.items():
        worker=OUT/action
        for p in list((OUT/'sweep-0').rglob('*.npz'))+list((OUT/'gr').glob('source-*.npz'))+[OUT/v for v in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
            dst=worker/p.relative_to(OUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
        for part in ['sweep-1/photons','sweep-1/material']:(worker/part).mkdir(parents=True)
        files += [OLD/f'gr/fields-{n}-g{q}{ext}' for ext in ['.npz','.json']]
    files += [INPUT/v for v in ['result.json','sources.json','driver-polynomial-audit.json','source-receipt.json']]
    files += [base.OUT/'expanded-fields.py',OUT/'check-result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'sources.json',read(INPUT/'sources.json'))
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Complete the original full-period primary GR field at actual existing coupled stages, reusing exactly the accepted causal15/16prefix.',
        method='Require bit-identical source endpoints, state coefficients, geometry maps, drive and background in the common past. Keep all535existing output values, evaluate only added times plus two overlaps, and compare overlap U/U_t/U_x and potential/free fields below1e-9. Enforce original potential iteration, time, radial quadrature and independent integral gates.',
        decision='Admitted full-period primary GR becomes input for actual full-period return. It is not itself the returned-solution charge, selfGR closure, a spatial/physical error certificate or final charge.',
        output_times=len(clock),reused_output_times=counts[0],new_output_times=len(clock)-counts[0],
        physical_horizon_seconds=float(clock[-1]),full_declared_period=True,
        budgets=CAPS,CPU_affinities=[0,6,12],virtual_GiB_per_worker=16,
        forecast='Use the three measured five-query prefix-regression receipts, actual added-time count and full/prefix source-cut ratio. Each worker has2hours; final scaling must fit this before dispatch. This avoids repeating the common535time field solve.',
        stop='Any bit identity, overlap, original numerical gate or resource cap; preserve failure, no automatic reintegration, resolution change or weakened tolerance.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(INPUT/'symbolic.json'))


def field(action):
    n,q=FIELDS[action];worker=OUT/action;bind(base.endpoint.initialize,OUT=worker)()
    d=dict(np.load(OUT/f'gr/source-{n}.npz'));old=dict(np.load(OLD/f'gr/fields-{n}-g{q}.npz'))
    result,row=extend(Response(),d,old,len(old['t']),q)
    np.savez_compressed(worker/f'gr/fields-{n}-g{q}.npz',**result);write(worker/f'gr/fields-{n}-g{q}.json',row)


def collect():
    for action,(n,q) in FIELDS.items():
        assert read(OUT/f'{action}-receipt.json')['error'] is None
        for ext in ['.npz','.json']:os.link(OUT/action/f'gr/fields-{n}-g{q}{ext}',OUT/f'gr/fields-{n}-g{q}{ext}')


def audit():
    s=(base.OUT/'expanded-fields.py').read_text()
    s=base.replace(s,'rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]',"rows=[read(OUT/f'gr/fields-{n}-g{q}.json') for n,q in [(128,8),(64,8),(128,4)]]")
    coefficients=lambda d:base.coefficients(dict(np.load(INPUT/f"gr/source-{int(d['original_clock'])}.npz")))
    ns=dict(base.prior.prior.fields.__globals__,OUT=OUT,Response=Response,coefficients=coefficients,PPoly=PPoly)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-audit.py').write_text(s);ns['fields']()
    r=read(OUT/'result.json');r.update(full_declared_period=True,causal_prefix_reused=True,
        primary_field_only=True,actual_full_period_GR_return_closed=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False);write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;OUT.mkdir(exist_ok=True)
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        field(action) if action in FIELDS else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
