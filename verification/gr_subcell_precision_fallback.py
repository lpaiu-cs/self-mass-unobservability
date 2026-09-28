"""Recompute the failed block with a recorded precision fallback, unchanged gates."""
from types import FunctionType, SimpleNamespace
import inspect, json, sys
import mpmath as mp
import numpy as np
import gr_full_subcell_reference as original
import gr_eos_precision_runner as precision

g=original.g;OUT=g.OUT/'gr-subcell-precision-fallback';quad=precision.original


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();precision.verify()
    assert json.loads((precision.OUT/'result.json').read_text())['finite_passed']
    before=inspect.getsource(original.block)
    old='''    inverse=FunctionType(g.audit.strict_invert.__code__,dict(g.audit.strict_invert.__globals__,
        ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))'''
    new='    inverse=make_inverse(stats)'
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    (OUT/'candidate-block.py').write_text(after)
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='226f237',recomputed_block=[2688,2816],
        intervention='The original block algorithm, geometry, baryon coordinate, all nodes and all gates. Retain the original binary64 evaluator for roots that pass. Only a recorded energy-unit root failure invokes the separately validated extended-arithmetic version of the same original source model. Return a single binary64 lnT and rounded EOS vector only after reevaluating that rounded input with the extended evaluator and passing the unchanged point gate.',
        fallback='15 extended Newton evaluations, then test the rounded root; only if needed test up to three binary64 neighbors on each side. No averaging, pressure adjustment, entropy change, tolerance increase or old failure overwrite.',
        interpretation='A mixed-precision numerical approximation to the same declared original EOS formulas, with explicit evaluator-difference records. No error enclosure across precision switches or certification of the underlying continuous/physical EOS. The new molecular spectrum candidate is not used here.',
        failed_original_sha256=g.c.sha(original.OUT/'block-2688-failure.json'))
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_subcell_precision_fallback.py',OUT/'candidate-block.py',precision.OUT/'manifest.json']})
    save('plan.json',plan)


def bindings():
    return FunctionType(original.bindings.__code__,dict(original.bindings.__globals__,OUT=OUT))()


def make_inverse(stats):
    fn=g.audit.strict_invert
    basic=FunctionType(fn.__code__,dict(fn.__globals__,ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))
    extended=None
    def inverse(eos,lp,entropy,eps,guess):
        nonlocal extended
        oldmax=stats['maximum_score']
        try:return basic(eos,lp,entropy,eps,guess)
        except ValueError as error:
            assert error.args[0][0]=='energy-unit entropy inverse',error
            failed_score=stats['maximum_score'];failed_error=str(error)
        if extended is None:extended=quad.EOS()
        mp.mp.dps=80
        best,trace=quad.solve(extended,mp.mpf(lp),mp.mpf(entropy),eps,guess,15)
        centre=float(best[1]);samples=[];lower=upper=centre
        candidates=[centre]
        for _ in range(3):
            lower=np.nextafter(lower,-np.inf);upper=np.nextafter(upper,np.inf);candidates.extend([float(lower),float(upper)])
        accepted=None
        for t in candidates:
            values=extended(mp.mpf(lp),mp.mpf(t),eps);a=np.array([float(v) for v in values])
            budget=max(2.,32*np.spacing(abs(a[2]+a[1]/a[0])))
            score=np.exp(t)*abs(a[3]-entropy)/budget
            samples.append(dict(lnT=t,score=float(score),budget_erg_g=float(budget)))
            if score<=1:accepted=(a,t,float(score));break
        index=stats.get('extended_fallbacks',0);record=dict(classification='Counterexample candidate',
            logP=float(lp),entropy=float(entropy),X=eps.tolist(),guess=float(guess),
            original_error=failed_error,original_failed_score=float(failed_score),extended_trace=trace,
            rounded_input_samples=samples,passed=accepted is not None)
        if accepted is not None:
            a,t,score=accepted;old=eos(1,lp,t,eps)
            record.update(old_evaluator_at_accepted_input=old.tolist(),accepted_EOS=a.tolist(),
                new_minus_old_EOS=(a-old).tolist(),returned_lnT=t,accepted_score=score)
            stats['maximum_score']=max(oldmax,score)
        save(f'fallback-{index:04}.json',record)
        assert accepted is not None,'precision fallback failed at unchanged gate; raw failure retained'
        stats['extended_fallbacks']=index+1
        stats['extended_evaluations']=stats.get('extended_evaluations',0)+len(trace)+len(samples)
        print('SUBCELL PRECISION FALLBACK',index,failed_score,score,flush=True)
        return a,t,len(trace)+len(samples)
    return inverse


def run():
    plan=bindings();assert not (OUT/'block-2688.json').exists()
    namespace=dict(original.block.__globals__,OUT=OUT,save=save,bindings=bindings,make_inverse=make_inverse)
    exec(compile((OUT/'candidate-block.py').read_text(),str(OUT/'candidate-block.py'),'exec'),namespace)
    record=namespace['block'](('block-2688',list(range(*plan['recomputed_block']))))
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(record['rows']),
        finite_passed=record['all_passed'],extended_fallbacks=record['entropy_roots'].get('extended_fallbacks',0),
        original_failure_retained=True,continuous_native_or_switch_error_certified=False,
        physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    plan=bindings();precision.verify()
    assert g.c.sha(original.OUT/'block-2688-failure.json')==plan['failed_original_sha256']
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['cells']==128
    for path in OUT.glob('fallback-*.json'):
        row=json.loads(path.read_text());assert row['passed'] and row['accepted_score']<=1
    print('PASS precision fallback block bindings; consult finite_passed, original failure retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
