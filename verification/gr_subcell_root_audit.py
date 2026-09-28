"""Audit the preserved entropy-root failure without changing its acceptance gate."""
from types import FunctionType, SimpleNamespace
import json, shutil, sys
import mpmath as mp
import numpy as np
import gr_full_subcell_reference as original
from molecular_partition_data import interval_record

g=original.g;OUT=g.OUT/'gr-subcell-root-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    shutil.copy2(original.OUT/'block-2688-failure.json',OUT/'original-failure.json')
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3d70a3b',
        discovery='The original full-subcell block 2688 stopped at score 1.0014343486391484. A preliminary scan of 129 neighboring binary64 lnT values found no pass. This is a confirming audit, not a blind discovery.',
        neighbors_each_side=64,orders=['ascending','descending'],
        gate='Unchanged abs(exp(lnT)*(s_native-s_target)) <= max(2 erg/g,32 ulp(abs(u+P/rho))). No entropy/pressure/temperature target, EOS, native solver or tolerance is changed.',
        implication='If adjacent representable lnT inputs straddle the target entropy but neither meets the gate, increasing the number of Newton iterations cannot insert another binary64 input between them. This is a local representability obstruction for this native point representation, not a proof of physical nonexistence or global root nonexistence.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/gr_subcell_root_audit.py',g.ROOT/'verification/audit_structured_enthalpy.py',
            g.ROOT/'verification/direct_eos_gr.py',g.ROOT/'verification/direct_ion_eos.py',
            OUT/'original-failure.json',original.OUT/'plan.json',g.d.OUT/'full-integral-build.json']},
        runtime=json.loads((g.d.OUT/'full-integral-build.json').read_text())['sha256']))


def run():
    plan=json.loads((OUT/'plan.json').read_text());failure=json.loads((OUT/'original-failure.json').read_text())
    eos=g.EOS();X=np.array(failure['X']);lp=failure['logP'];target=failure['entropy']
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label='replayed-root')
    fn=g.audit.strict_invert;inverse=FunctionType(fn.__code__,dict(fn.__globals__,
        ROOT_STATS=stats,e=SimpleNamespace(OUT=OUT,save=save)))
    try:inverse(eos,lp,target,X,failure['guess'])
    except ValueError as error:save('replay.json',dict(classification='Counterexample candidate',error=str(error),statistics=stats))
    else:raise AssertionError('original failure did not replay')
    assert json.loads((OUT/'replayed-root-failure.json').read_text())==failure
    t=failure['best_lnT'];low=high=t;ts=[t]
    for _ in range(plan['neighbors_each_side']):
        low=np.nextafter(low,-np.inf);high=np.nextafter(high,np.inf);ts.extend([low,high])
    ts=np.array(sorted(ts));arrays=[];rows=[]
    for order in plan['orders']:
        indices=range(len(ts)) if order=='ascending' else reversed(range(len(ts)))
        values=np.zeros((len(ts),21))
        for i in indices:values[i]=eos(1,lp,float(ts[i]),X)
        arrays.append(values)
    assert np.array_equal(*arrays),'native history/order dependence'
    values=arrays[0];errors=values[:,3]-target
    budgets=np.maximum(2.,32*np.spacing(abs(values[:,2]+values[:,1]/values[:,0])))
    scores=np.exp(ts)*abs(errors)/budgets
    crossings=np.flatnonzero(errors[:-1]*errors[1:]<0);assert len(crossings)==1
    k=int(crossings[0]);assert np.nextafter(ts[k],np.inf)==ts[k+1]
    assert scores.min()>1 and set(np.flatnonzero(scores<1.01))=={k,k+1}
    mp.iv.dps=80
    def exact(value):
        a,b=float(value).as_integer_ratio();return mp.iv.mpf(a)/b
    bounds=[]
    for i in [k,k+1]:
        # Entropy subtraction is done as exact input rationals, not rounded first.
        error=exact(values[i,3])-exact(target)
        score=mp.iv.exp(exact(ts[i]))*abs(error)/exact(budgets[i]);assert score.a>1
        bounds.append(interval_record(score))
    np.savez_compressed(OUT/'neighbor-values.npz',lnT=ts,eos=values,entropy_error=errors,budget_erg_g=budgets,score=scores)
    save('result.json',dict(classification='Counterexample candidate',completed=True,original_failure_reproduced=True,
        sampled_temperatures=len(ts),native_evaluations=2*len(ts),all_outputs_order_bitwise_equal=True,
        minimum_sampled_score=float(scores.min()),any_sampled_pass=False,
        adjacent_input_hex=[float(ts[i]).hex() for i in [k,k+1]],
        adjacent_entropy_residual_erg_g_K=[float(errors[i]) for i in [k,k+1]],
        budget_erg_g=[float(budgets[i]) for i in [k,k+1]],
        no_physical_root_claim=False,original_full_grid_failure_retained=True))
    save('local-representability.json',dict(classification='Proven',adjacent_binary64_inputs=True,
        two_point_energy_unit_scores=bounds,
        statement='For the two consecutive binary64 lnT inputs and the recorded native entropy outputs interpreted as exact rational oracle values, both energy-unit residuals exceed the unchanged budget, even with interval evaluation of exp(lnT). There is no binary64 lnT value strictly between these endpoints. Hence this input interval has no admissible point under that oracle/gate representation.',
        limit='No enclosure of native EOS evaluation error or physical entropy is assumed or proven. Opposite signs alone are not a validated continuum root bracket without an EOS error enclosure and continuity. Wider searches do not change this local fact. Higher-precision input/evaluation or a justified interval primitive representation is required before calling a repaired full-grid run a pass.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('SUBCELL ROOT LOCAL REPRESENTABILITY',len(ts),float(scores.min()),flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in json.loads((OUT/'plan.json').read_text())['runtime'].items():assert g.c.sha(path)==digest,path
    result=json.loads((OUT/'result.json').read_text());assert result['completed'] and result['original_failure_reproduced']
    assert not result['any_sampled_pass'] and result['minimum_sampled_score']>1
    assert json.loads((OUT/'local-representability.json').read_text())['adjacent_binary64_inputs']
    print('PASS preserved root failure and conditional local point-representation obstruction',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
