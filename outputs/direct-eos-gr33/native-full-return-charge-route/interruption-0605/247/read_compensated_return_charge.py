"""Read compact scalar charge from both components of the same GR-return solve."""
from pathlib import Path
from types import SimpleNamespace
from decimal import Decimal,localcontext
import json,os,resource,sys,time
import numpy as np
from scipy.interpolate import PPoly
import read_dense_returned_source as prior

CHECK=prior.CHECK;ROOT=Path('native-compensated-charge247-work');OUT=ROOT/('check' if CHECK else 'full')
INPUT=prior.OUT;HIGH=prior.prior.INPUT;base=prior.base
read,write,sha,bind,LD=prior.read,prior.write,prior.sha,prior.bind,prior.LD
KEYS=prior.KEYS;FIELDS=dict(field1288=(128,8),field648=(64,8),field1284=(128,4))
CAPS=dict(prepare=600,field1288=7200,field648=7200,field1284=7200,collect=300,audit=900,compare=300)


def coefficients(d):
    return prior.coefficients(dict(np.load(INPUT/f"gr/source-{int(d['original_clock'])}.npz")))


bridge=SimpleNamespace(**dict(vars(base),coefficients=prior.coefficients))
class Response(base.Response):
    setup=bind(prior.prior.returned.Response.setup,INPUT=INPUT,bridge=bridge)


def prepare():
    r=read(INPUT/'result.json');assert r['representation_controls_passed'] and r['source_time_passed']
    assert read(INPUT/'source-receipt.json')['error'] is None
    assert not OUT.exists();OUT.mkdir(parents=True);files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/part).mkdir(parents=True)
    for p in list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);files += [p,dst]
    clocks=[]
    for n in [64,128]:
        raw=dict(np.load(INPUT/f'gr/source-{n}.npz'));high=dict(np.load(HIGH/f'gr/source-{n}.npz'))
        times=np.load(HIGH/f'gr/fields-{n}-g8.npz')['t'];clocks.append(times)
        for k in ['edges','radius','volume','M_cm','K_cm','cx']:assert np.array_equal(raw[k],high[k]),k
        knots,co=prior.coefficients(raw);d=dict(raw,t=times,original_clock=np.array(n))
        d.update({k:PPoly(np.asarray(v[::-1],float),knots)(times) for k,v in co.items()})
        ids=np.array([np.argmin(abs(times-t)) for t in raw['t']]);assert np.max(abs(times[ids]-raw['t']))<1e-18
        for k in KEYS:d[k][ids]=raw[k]
        np.savez_compressed(OUT/f'gr/source-{n}.npz',**d)
        files += [INPUT/f'gr/source-{n}.npz',HIGH/f'gr/source-{n}.npz',HIGH/f'run-{n}.json']
        assert read(HIGH/f'run-{n}.json')['same_saved_stage_equation']
    assert np.array_equal(*clocks)
    for action,(n,q) in FIELDS.items():
        worker=OUT/action
        for p in list((OUT/'sweep-0').rglob('*.npz'))+list((OUT/'gr').glob('source-*.npz'))+[OUT/v for v in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
            dst=worker/p.relative_to(OUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
        for part in ['sweep-1/photons','sweep-1/material']:(worker/part).mkdir(parents=True)
        files += [HIGH/f'gr/fields-{n}-g{q}{ext}' for ext in ['.npz','.json']]
    files += [INPUT/v for v in ['result.json','plan.json','sources.json','geometry-result.json','source-receipt.json']]
    files += [HIGH/'result.json',base.OUT/'expanded-fields.py']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'sources.json',read(INPUT/'sources.json'))
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the actual dense same-return source to retarded GR and compare the componentwise compact scalar charge with its own high anchor.',
        method='Reuse the characteristic polynomial integrator and independent Jordan-radius integral. Use exactly the high GR output clock and same background/radial orders so the potential feedback representation is shared. Keep actual post-floor source endpoints. Store high and low fields/charges separately and sum final binary inputs at100decimal digits; no rounded-state addition or unrelated response.',
        gates=dict(source_time=.02,charge_field_time=.02,quadrature=.002,independent_GR=1e-9,potential_iteration=1e-8),
        budgets=CAPS,CPU_threads_per_process=1,virtual_GiB_per_process=16,CPU_affinities=[0,8,12],
        forecast='224short three GR fields cost120/62/76seconds on the original30/16output times and533source cuts. This low source has fewer polynomial cuts but73output times, with marginal cost unmeasured. Plan3..12minutes for the longest short worker; allow2hours EACH. Full535time scaling is unmeasured and must be reported from these receipts. No material/photon/EOS/background rerun.',
        decision='Only source-time, compact-U time, radial quadrature and independent integral gates admit the one-return compact charge comparison. Preserve all U_t/U_x differences. A surviving compact sign is conditional on this retained operator/time representation; it is not infinity-normalized charge, selfGR closure or a physical error certificate.',
        scope='Same actual227short or236common15/16high/low solve. Sampled Hermite controls only; EOS/derivative/spatial/boundary/nonlinear/static/observation/infinity conditions remain.',
        stop='Source provenance/representation failure, original numerical gate or generous resource cap. No new physical clock, automatic refinement, additional return iterate or tolerance relaxation.',
        input_directory=str(INPUT),high_directory=str(HIGH),output_times=len(clocks[0]),
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(INPUT/'symbolic.json'))


def field(action):
    n,q=FIELDS[action];worker=OUT/action
    bind(base.endpoint.initialize,OUT=worker)()
    fn=bind(base.base.gr.base.Response.run,OUT=worker/'gr');fn(Response(),n,q)


def collect():
    for action,(n,q) in FIELDS.items():
        assert read(OUT/f'{action}-receipt.json')['error'] is None
        for ext in ['.npz','.json']:os.link(OUT/action/f'gr/fields-{n}-g{q}{ext}',OUT/f'gr/fields-{n}-g{q}{ext}')


def audit():
    s=(base.OUT/'expanded-fields.py').read_text()
    s=base.replace(s,'rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]',"rows=[read(OUT/f'gr/fields-{n}-g{q}.json') for n,q in [(128,8),(64,8),(128,4)]]")
    ns=dict(base.prior.prior.fields.__globals__,OUT=OUT,Response=Response,coefficients=coefficients,PPoly=PPoly)
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-audit.py').write_text(s);ns['fields']()


def compare():
    r=read(OUT/'result.json');admitted=r['GR_return_admitted'];components=[]
    with localcontext() as ctx:
        ctx.prec=100
        for n in [64,128]:
            high=read(HIGH/f'gr/fields-{n}-g8.json');low=read(OUT/f'gr/fields-{n}-g8.json');values={}
            for name in ['endpoint_direct','endpoint_compact_with_metric']:
                a=Decimal.from_float(high[name]);b=Decimal.from_float(low[name]);total=a+b
                values[name]=dict(high=str(a),low=str(b),componentwise_total=str(total),
                    low_over_high=float(b/a) if a else None,nominal_sign_unchanged=bool(a*total>0))
            components.append(dict(clock=n,values=values))
    r.update(classification='Counterexample candidate',same_actual_returned_solution_read=True,
        shared_high_low_GR_output_clock=True,components=components,charge_comparison_admitted=admitted,
        conditional_compact_sign_survives_one_return=all(v['values']['endpoint_compact_with_metric']['nominal_sign_unchanged'] for v in components) if admitted else None,
        uniform_temporal_error_certificate=False,self_GR_return_closed=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1].removeprefix('check_');assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        field(action) if action in FIELDS else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
