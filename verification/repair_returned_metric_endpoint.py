"""Use the source interval's left derivative at closing physical stages.

Keep the existing two clocks, original coupled equations and acceptance gates.
The changed input must be evolved; a diagnosed work correction is never added
to a saved state or to its final charge.
"""
from pathlib import Path
from types import SimpleNamespace
import inspect,os,resource,shutil,sys,time
import numpy as np
from scipy.interpolate import PPoly
import finish_returned_material_accuracy as prior

OUT=Path('native-left-metric258-work');OLD=prior.OUT
DIAG=Path('native-mass-time258-work')
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
base=prior.base;joint=prior.joint
CAPS=dict(prepare=600,check=600,coarse=21600,fine=28800,audit=900)


def derivative_at(poly,t,side):
    """Horner evaluation at the exact knot, without a nextafter time shift."""
    p=poly.derivative();ids=np.clip(np.searchsorted(p.x,t,side=side)-1,0,len(p.x)-2)
    dt=(t-p.x[ids]).reshape((-1,)+(1,)*(p.c.ndim-2))
    value=np.zeros((len(t),)+p.c.shape[2:])
    for co in p.c:value=value*dt+co[ids]
    return value


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'result.json')['passed']
    assert read(DIAG/'canonical-affine-result.json')['passed']
    assert read(DIAG/'pressure-receipt.json')['error'] is None
    assert not read(Path('native-material-exterior257-work/full/audit.json'))['passed']
    files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for p in (OLD/'gr').iterdir() if p.is_file()]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    for part in ['metric','sweep-1/photons','sweep-1/material']:(OUT/part).mkdir(parents=True,exist_ok=True)
    for n,q in [(64,8),(128,4),(128,8)]:
        src=DIAG/f'metric-left-{n}-g{q}.npz';dst=OUT/f'metric/metric-{n}-g{q}.npz'
        with np.load(src) as new,np.load(OLD/f'metric/metric-{n}-g{q}.npz') as old:
            assert set(new.files)==set(old.files)
            for k in old.files:
                if k!='actual_delta_lambda_rate':assert np.array_equal(new[k],old[k]),(n,q,k)
        shutil.copy2(src,dst);files += [src,dst,OLD/f'metric/metric-{n}-g{q}.npz']
    files += [DIAG/p for p in ['pressure-work.json','pressure-plan.json','pressure-receipt.json','canonical-affine-result.json']]
    for name,replacements in [
        ('charge-reader.py', [('native-material-charge257-work','native-left-metric-charge258-work'),('native-material-accuracy257-work','native-left-metric258-work'),('actual257full-period repaired return','actual258full-period corrected-knot return')]),
        ('exterior-reader.py', [('native-material-exterior257-work','native-left-metric-exterior258-work'),('native-material-accuracy257-work','native-left-metric258-work'),('native-material-charge257-work','native-left-metric-charge258-work'),('Only257actual','Only258actual')])]:
        original=Path('verification')/('read_material_return_charge.py' if name.startswith('charge') else 'read_material_return_exterior.py')
        source=original.read_text()
        for a,b in replacements:
            assert a in source;source=source.replace(a,b)
        dst=OUT/name;dst.write_text(source);files += [original,dst]
    files += [OLD/'result.json',Path('native-material-exterior257-work/full/audit.json'),Path(__file__)]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the corrected source-knot derivative to actual photon/material evolution and read the same completed solution through the original compact and mass gates.',
        cause='The applied constraint derivative used the right polynomial at closing Radau endpoints. All700actual native rates were reproduced exactly. The clock energy difference is2.3066239671e-19erg; explicit pressure work contributes2.3065576554e-19erg. Correct left knot values remove7.6736098142e-20erg in the read-only forcing comparison. Remaining coarse quadrature is not declared solved.',
        repair='Change only actual_delta_lambda_rate at source polynomial knots, using the exact left derivative. Preserve original causal metric values, fields, source polynomials, initial state, grid, clocks, horizon, EOS and all tolerances. No primitive-work correction or diagnostic energy is added. Actual joint states, photons, ports and final charge are recomputed together.',
        reuse='Reuse the completed primary239, full248GR fields and all unchanged249metric arrays. Changed forcing starts in the first canonical interval, so old257returned states cannot be causally reused. Re-evolve only the original119/231return paths; no new resolution/path/period. Keep every new canonical checkpoint.',
        solver='Existing257right preconditioner80/20,12linear refinements,3Newton proposals and separate material1e-13gate. Original whole/physical/native/constitutive/balance/dense and2percent time gates remain.',
        budgets=CAPS,CPU_affinity=3,virtual_GiB=16,
        forecast='236coarse111steps measured2258s;257late16fine steps3986s. Full changed RHS is unmeasured: assume40..100minutes coarse and90..240minutes fine, plus about50..90minutes unchanged charge route. Allow6h coarse/8h fine; use first canonical interval receipts to revise the estimate without interrupting accepted progress.',
        decision='Continue the actual original pair through its final mass/charge audit. Any original physical gate failure stops; preserve it. Do not automatically add finer paths, weaken the2percentmass gate or call this full physical exterior closure.',
        scientific_gates_changed=False,discrete_equations_changed=False,source_knot_side_corrected=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def initialize():return bind(prior.initialize,OUT=OUT)()


def check():
    # Continuous q has a derivative jump at1. The closing interval must use1,
    # not3; choosing3 makes Radau's exact integral of its linear piece1.5.
    p=PPoly(np.array([[1.,3.],[0.,1.]]),[0.,1.,2.]);t=np.array([1/3,1.])
    left=derivative_at(p,t,'left');right=derivative_at(p,t,'right')
    assert np.array_equal(left,[1.,1.]) and np.array_equal(right,[1.,3.])
    assert np.array([.75,.25])@left==p(1.)-p(0.)==1.
    assert np.array([.75,.25])@right==1.5
    bind(base.actual.bridge.endpoint.initialize,OUT=OUT)();m=base.Lapse();rows=[]
    for n,q in [(64,8),(128,4),(128,8)]:
        d=dict(np.load(OLD/f'gr/source-{n}.npz'));f=np.load(OLD/f'gr/fields-{n}-g{q}.npz')
        _,poly=base.representation(n);values=[]
        derivative=dict(delta_phi=f['U_t']/f['radius_E'],delta_Phi=np.zeros_like(f['U_t']))
        for side in ['right','left']:
            dot=dict(d);dot.update({k:derivative_at(v,d['t'],side) for k,v in poly.items()})
            values.append(base.geometry.constraints.centers(m.response,dot,derivative,q)['delta_lambda'][:,:-1])
        old=np.load(OLD/f'metric/metric-{n}-g{q}.npz')['actual_delta_lambda_rate']
        new=np.load(OUT/f'metric/metric-{n}-g{q}.npz')['actual_delta_lambda_rate']
        assert np.array_equal(new,old+(values[1]-values[0]))
        rows.append(dict(clock=n,order=q,independent_left_metric_exact=True))
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'))
    keys=[k for k in base.geometry.KEYS if k!='delta_lambda_rate']+['delta_u_t','actual_delta_lambda_rate']
    controls={}
    for label,n,q in [('time',64,8),('quadrature',128,4)]:
        other=dict(np.load(OUT/f'metric/metric-{n}-g{q}.npz'))
        controls[label]={k:base.base.endpoint.aligned(other,fine,k) for k in keys}
    assert max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002,controls
    result=read(OLD/'metric-result.json');result.update(passed=True,left_knot_derivative=True,controls=controls,rows_left=rows)
    write(OUT/'metric-result.json',result)
    write(OUT/'endpoint-check.json',dict(classification='Proven',passed=True,
        scope='Exact piecewise-linear knot-side counterexample and saved-polynomial metric reproduction only; not physical convergence.',
        rows=rows,metric_controls=controls))


def evolve(n):
    # Reuse the full capture/ledger producer, changing its loop to the existing
    # complete period and using the already verified257solver implementation.
    source=inspect.getsource(base.actual.evolve)
    a="Model=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()"
    assert source.count(a)==1 and source.count('for j in [1,2]:')==1
    source=source.replace(a,'Model=initialize()').replace('for j in [1,2]:','for j in range(1,17):')
    ns=dict(base.actual.evolve.__globals__,OUT=OUT,prior=SimpleNamespace(saved=base.saved),initialize=initialize)
    exec(compile(source,__file__,'exec'),ns);ns['evolve'](n)
    row=read(OUT/f'run-{n}.json');row.update(full_declared_period=True,source_knot_side_corrected=True)
    write(OUT/f'run-{n}.json',row)


def audit():
    bind(base.actual.audit,OUT=OUT)();r=read(OUT/'result.json')
    r.update(full_declared_period=True,scientific_gates_changed=False,source_knot_side_corrected=True,
        original257_mass_time_failure_preserved=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
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
