"""Apply the complete one-return outer photon lapse to the actual coupled solution.

Counterexample candidate. The phase-259 applied returned metric is kept except
for one missing boundary term: the accepted incident-metric background-photon
particular lapse (phase 261), evaluated at all 575 applied times (phase 265
clock). Lapse.run integrates inward from the outer value, so the same time
function is added to delta_nu, its faces, delta_log_lapse and delta_log_speed.
lambda, u, phi and every derivative used by the stages are unchanged; the mass
residual is unchanged because the particular solution already contains the
packets' own mass. The original 119/231 pair is re-evolved from t=0 with the
phase-259 solver and gates and read by the unchanged same-solution readers.
"""
from pathlib import Path
from types import SimpleNamespace
import fcntl,os,resource,sys,time
import numpy as np
import continue_returned_short_krylov as short

live=short.live
OUT=Path('native-photon-boundary265-work');OLD=short.OUT
GEOM=Path('native-geometric-clock265-work')
CHARGE=Path('native-photon-boundary-charge265-work');EXTERIOR=Path('native-photon-boundary-exterior265-work')
read,write,sha,bind=short.read,short.write,short.sha,short.bind
base=live.base;joint=short.joint;LD=short.LD
CAPS=dict(prepare=900,check=1800,coarse=21600,fine=28800,audit=900)
SHIFTED=['delta_nu','delta_nu_faces','delta_log_lapse','delta_log_speed']
METRICS=[(64,8),(128,4),(128,8)]
TAG=[None]


def boundary():
    g=np.load(GEOM/'boundary-575.npz');t=np.load(OLD/'metric/metric-128-g8.npz')['t']
    assert np.array_equal(g['t'],t)
    return t,np.asarray(g['photon_geometric_lapse'],float),np.asarray(g['photon_geometric_mass_cm'],float)


def shifted(old,db,mass):
    new=dict(old)
    for k in SHIFTED:new[k]=old[k]+db[:,None]
    new['delta_nu_interval_rate']=old['delta_nu_interval_rate']+(np.diff(db)/np.diff(old['t']))[:,None]
    new['outer_photon_geometric_lapse']=db;new['photon_geometric_mass_cm']=mass
    return new


def prepare():
    assert read(OLD/'result.json')['passed'] and read(OLD/'controller-status.json')['state']=='completed'
    geometric=read(GEOM/'result.json');assert geometric['passed'] and all(geometric['knot_values_exact'])
    assert not OUT.exists();OUT.mkdir()
    files=list((OLD/'sweep-0').rglob('*.npz'))+[p for p in (OLD/'gr').iterdir() if p.is_file()]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','endpoint-check.json']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    for part in ['metric','sweep-1/photons','sweep-1/material']:(OUT/part).mkdir(parents=True,exist_ok=True)
    t,db,mass=boundary()
    for n,q in METRICS:
        src=OLD/f'metric/metric-{n}-g{q}.npz';old=dict(np.load(src))
        np.savez_compressed(OUT/f'metric/metric-{n}-g{q}.npz',**shifted(old,db,mass));files.append(src)
    for name,pairs in [('charge-reader.py',[('native-short-return-charge259-work',str(CHARGE)),('native-short-return259-work',str(OUT)),
                                           ('actual259full-period corrected-knot return','actual265full-period complete-photon-boundary return')]),
                       ('exterior-reader.py',[('native-short-return-exterior259-work',str(EXTERIOR)),('native-short-return-charge259-work',str(CHARGE)),
                                             ('native-short-return259-work',str(OUT)),('Only259actual','Only265actual')])]:
        text=(OLD/name).read_text()
        for a,b in pairs:
            assert a in text,(name,a);text=text.replace(a,b)
        assert '259' not in text.replace('phase259','').replace('259actual',''),name
        (OUT/name).write_text(text);files += [OLD/name,OUT/name]
    files += [OLD/n for n in ['result.json','controller-status.json','plan.json','metric-result.json']]
    files += [GEOM/n for n in ['result.json','boundary-575.npz','plan.json','execution-plan.json']]+[Path(__file__)]
    files += [OUT/f'metric/metric-{n}-g{q}.npz' for n,q in METRICS]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Evolve the same original corrected 119/231 pair with the complete one-return outer photon lapse and read that solution with the unchanged compact, frozen-exterior and mass readers.',
        cause='The applied outer lapse contained the emitted-packet, scalar and matter-mass terms, but not the accepted incident-metric response of the background photons (phase261: up to2.77e-3 of the outer lapse).',
        repair='Add the phase-265 575-time particular lapse to the outer boundary. Because the lapse is integrated inward from that value, add the same time function to delta_nu, its faces, delta_log_lapse and delta_log_speed. No other metric array, source, field, primary history, EOS, grid, clock, horizon, solver or tolerance changes. The mass residual is not changed: the particular lapse already contains the packets\' mass.',
        reuse='Reuse the259GR sources and fields, primary history, and every unchanged metric array. The new lapse is nonzero from the first interval, so no259returned state is causally reusable; both clocks restart at t=0. Readers are the259adapters with only directories changed.',
        solver='Unchanged259solver: short right-preconditioned GMRES80/1 proposals with the257original solver as fallback; original vector1e-14, physical1e-13, material1e-13, three-Newton, constitutive, conservation, dense and2percent paired-time gates.',
        parallel='Coarse and fine run concurrently on CPUs4/6. Their shared generated-code initialization is serialized with a separate lock, and their solver diagnostics are written per clock.',
        budgets=CAPS,CPU_affinities=[4,6],virtual_GiB=16,
        forecast='259measured coarse111steps3615.5s and fine231steps6055.8s. Full coarse119steps from t=0 about65min and fine about101min if the per-step cost holds; unmeasured for the new metric. Keep the259caps6h/8h. Readers about31min (259:1880s).',
        decision='Continue the actual pair through its own final compact/frozen-exterior/mass audit. Any original gate failure stops and is preserved; no automatic finer clock, grid, path or tolerance change, and no superposition with the259solution.',
        scientific_gates_changed=False,discrete_equations_changed=False,source_knot_side_corrected=True,photon_geometric_boundary_applied=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def tagged(p,value):
    p=Path(p)
    if p.name in ('right-calls.json','last-short-exhaustion.json') and TAG[0]:p=p.with_name(f'{p.stem}-{TAG[0]}{p.suffix}')
    return write(p,value)


def initialize():
    with Path('.phase265-evolve-initialize.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return bind(short.initialize,OUT=OUT,write=tagged)()


def check():
    t,db,mass=boundary();rows=[]
    for n,q in METRICS:
        old=np.load(OLD/f'metric/metric-{n}-g{q}.npz');new=np.load(OUT/f'metric/metric-{n}-g{q}.npz')
        assert set(new.files)==set(old.files)|{'outer_photon_geometric_lapse','photon_geometric_mass_cm'}
        for k in old.files:
            if k in SHIFTED:assert np.array_equal(new[k],old[k]+db[:,None]),(n,q,k)
            elif k=='delta_nu_interval_rate':assert np.array_equal(new[k],old[k]+(np.diff(db)/np.diff(old['t']))[:,None])
            else:assert np.array_equal(new[k],old[k]),(n,q,k)
        rows.append(dict(clock=n,order=q,exact_shift=True,maximum_shift_over_outer_lapse=float(np.max(abs(db))/np.max(abs(old['delta_nu_faces'][:,-1])))))
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'))
    keys=[k for k in base.geometry.KEYS if k!='delta_lambda_rate']+['delta_u_t','actual_delta_lambda_rate']
    controls={label:{k:base.base.endpoint.aligned(dict(np.load(OUT/f'metric/metric-{n}-g{q}.npz')),fine,k) for k in keys}
              for label,n,q in [('time',64,8),('quadrature',128,4)]}
    passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002
    result=read(OLD/'metric-result.json');result.update(passed=bool(passed),controls_before_photon_boundary=result.get('controls'),
        controls=controls,photon_geometric_boundary_rows=rows,photon_geometric_boundary_applied=True)
    write(OUT/'metric-result.json',result);assert passed,controls
    Model=initialize();geometry=base.geometry
    assert geometry.OUT==OUT and geometry.ReturnOnly is base.actual.StageDriver
    driver=geometry.ReturnOnly(SimpleNamespace(n=len(fine['radius_E'])));assert 'outer_photon_geometric_lapse' in driver.g
    probe=[]
    for i in [1,len(t)//2,len(t)-1]:
        row=driver.returned(t[i])
        for k in ['delta_log_lapse','delta_nu','delta_lambda','delta_u']:assert np.array_equal(row[k],fine[k][i].astype(LD)),(i,k)
        probe.append(dict(index=i,time=float(t[i]),log_lapse_shift=float(db[i])))
    del Model
    write(OUT/'boundary-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,metric_controls=controls,
        stage_driver_reads_new_metric=True,stage_driver_probe=probe,
        scope='Exact construction of the shifted applied metric, its original time/quadrature gates and the live stage-driver lookup; not physical convergence.'))


def evolve(n):bind(live.evolve,OUT=OUT,initialize=initialize)(n)


def audit():
    bind(live.audit,OUT=OUT)();r=read(OUT/'result.json')
    r.update(photon_geometric_boundary_applied=True,complete_one_return_outer_photon_lapse=True,scientific_gates_changed=False,
        reciprocal_scalar_operator_variation_complete=False,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None;TAG[0]=action
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            assert not (OUT/f'{action}-receipt.json').exists()
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        evolve(64 if action=='coarse' else 128) if action in ['coarse','fine'] else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
