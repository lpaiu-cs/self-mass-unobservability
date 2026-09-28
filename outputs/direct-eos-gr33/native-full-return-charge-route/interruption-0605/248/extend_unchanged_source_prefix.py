"""Reuse GR only before the horizon-dependent final source segment."""
from pathlib import Path
from types import SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import extend_retarded_history as prior

OUT=prior.OUT;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prior.CAPS,check=600)


def cutoff():
    raw=np.load(prior.base.OUT/'gr/source-64.npz')
    boundary=raw['geometry_times'][14]
    clock=np.load(prior.OLD/'gr/source-64.npz')['t']
    return float(clock[clock<boundary-1e-18][-1])


def prefix(new,old,cut):
    # Compare every coefficient whose interval intersects the reused past,
    # including its closing knot, as well as all endpoint and background data.
    state=int(np.searchsorted(old['t'],cut,side='left'))
    geometry=int(np.searchsorted(old['geometry_times'],cut,side='left'))
    for key,value in old.items():
        v=new[key]
        if key=='t' or key in prior.base.KEYS:v=v[:state+1];value=value[:state+1]
        elif key.startswith('state_coeff_'):v=v[:,:state];value=value[:,:state]
        elif key=='geometry_times':v=v[:geometry+1];value=value[:geometry+1]
        elif key.startswith('geometry_coeff_'):v=v[:,:geometry];value=value[:,:geometry]
        assert np.array_equal(v,value),('Reused source changed',key)


def check():
    assert read(OUT/'regression.json')['passed']
    assert read(OUT/'check-receipt.json')['source_sha256']==sha(prior.__file__)
    cut=cutoff();rows=[];files=[Path(__file__),Path(prior.__file__)]
    for n in [64,128]:
        paths=[prior.INPUT/f'gr/source-{n}.npz',prior.base.OUT/f'gr/source-{n}.npz']
        new,old=[dict(np.load(p)) for p in paths];prefix(new,old,cut)
        try:prior.source_prefix(new,old)
        except AssertionError as exc:original=repr(exc)
        else:raise AssertionError('The recorded horizon-dependent prefix mismatch disappeared')
        bad=dict(new);k='geometry_coeff_gas_nonrest_energy_erg';bad[k]=new[k].copy()
        bad[k][0,0,0]=np.nextafter(bad[k][0,0,0],prior.LD(np.inf))
        try:prefix(bad,old,cut)
        except AssertionError:pass
        else:raise AssertionError('Changed causal input accepted')
        rows.append(dict(clock=n,original_rejection=original,unchanged_past_exact=True,changed_past_rejected=True));files+=paths
    write(OUT/'prefix-check.json',dict(classification='Counterexample candidate',passed=True,
        reused_through_seconds=cut,rows=rows,
        method='Reuse only the bit-identical causal input through the last output before canonical14/16. Recompute all later GR outputs, including formerly cached points; do not alter244or224source arrays. Original source-prefix failure remains.',
        field_extension_algorithm_unchanged=True,original_field_regression_sha256=sha(OUT/'regression.json'),
        bindings={str(p):sha(p) for p in files},final_charge_conclusion='unadjudicated'))


def link(src,dst):
    if Path(dst).exists():assert sha(src)==sha(dst),(src,dst)
    else:os.link(src,dst)


def prepare():
    proof=read(OUT/'prefix-check.json');assert proof['passed']
    for p,h in proof['bindings'].items():assert sha(p)==h,p
    s=inspect.getsource(prior.prepare)
    replace=prior.base.replace
    s=replace(s,'source_prefix(raw,old);knots','source_prefix(raw,old,cut);knots')
    s=replace(s,"count=len(previous_source['t'])","count=int(np.searchsorted(previous_source['t'],cut,side='right'))")
    s=replace(s,"clock[:count],previous_source['t']","clock[:count],previous_source['t'][:count]")
    s=replace(s,"d[k][:count],previous_source[k]","d[k][:count],previous_source[k][:count]")
    ns=dict(prior.prepare.__globals__,cut=proof['reused_through_seconds'],source_prefix=prefix,
        os=SimpleNamespace(**dict(vars(os),link=link)))
    exec(compile(s,__file__,'exec'),ns);ns['prepare']()
    r=read(OUT/'plan.json');r.update(
        claim='Complete the identical full-period244GR source, reusing only its exactly unchanged causal prefix before14/16.',
        method=proof['method']+' Original extend(), all1e-9overlap/1e-12prefix and numerical gates unchanged. Two overlap queries independently retest cache validity.',
        forecast='Existing five-query regression scaled by actual recomputed output count and source cuts;2x allowance plus5minutes,2hours and16GiB perfield.',
        original_prefix_failure_preserved=True,reused_through_seconds=proof['reused_through_seconds'])
    r['bindings'].update({str(Path(__file__)):sha(__file__),str(OUT/'prefix-check.json'):sha(OUT/'prefix-check.json')})
    write(OUT/'plan.json',r)


def field(action):
    s=inspect.getsource(prior.field)
    s=prior.base.replace(s,"len(old['t']),q","read(OUT/'plan.json')['reused_output_times'],q")
    ns=dict(prior.field.__globals__);exec(compile(s,__file__,'exec'),ns);ns['field'](action)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/('prefix-check-receipt.json' if action=='check' else f'{action}-receipt.json');assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    prior.base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in prior.FIELDS:field(action)
        elif action in ['check','prepare']:globals()[action]()
        else:getattr(prior,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
