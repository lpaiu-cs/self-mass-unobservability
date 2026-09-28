"""EOS electron inputs at the actual captured opacity states.

Counterexample candidate. The declared FreeEOS evaluator is electron-only.
A conditional ideal-pair bound is kept separate from physical plasma error.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, sys
import numpy as np
from mpmath import iv
import native_opacity as o

g=o.g;OUT=o.OUT/'electrons'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert json.loads((o.OUT/'identity-control.json').read_text())['passed']
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='7772717',
        inputs_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            o.OUT/'baseline-type1-captured.npz',o.OUT/'identity-control.json',
            g.d.OUT/'full-integral-build.json',g.ROOT/'verification/direct_eos_gr.py']},
        formula='lnfree_e=ln(rmue)-ln(rho_B), with derivatives at fixed composition. This is the sum of free electrons and the zero positron population in the declared electron-only EOS model.',
        states='Use the actual native Type1 log10 rho, log10 T and 26-species composition captured at entry, including the native input rounding.',
        finite_log_steps=[5e-5,2.5e-5],derivative_relative_tolerance=1e-3,
        block_cells=128,processes=2,
        scope='Input matching to the declared EOS, not a certified physical pair/plasma EOS or opacity derivative enclosure.',
        physical_EOS_certified=False,full_nonlinear_transport=False))
    pair_bound()


def pair_bound():
    iv.dps=60;minimum_eta=iv.mpf(-50);maximum_beta=iv.mpf('.005')
    ratio=iv.exp(-minimum_eta-2/maximum_beta)+iv.exp(-2*minimum_eta-2/maximum_beta)
    relative=2*ratio/(1-ratio);assert relative.b<iv.mpf('1.03e-130').a
    save('ideal-pair-bound.json',dict(classification='Proven',passed=True,
        model='An ideal relativistic electron/positron gas with the same nonnegative density-of-states weight for both charges, beta=kBT/(me c^2), eta_plus=-eta_minus-2/beta, and finite positive Boltzmann integral.',
        domain=dict(eta_minimum=-50,beta_maximum='.005',beta_positive=True),
        proof='For x>=0, q_plus<=exp(-x-eta-2/beta), whereas q_minus>=exp(eta-x)/(1+exp(eta)). Integrating against the same positive weight gives n_plus/n_minus<=exp(-eta-2/beta)+exp(-2eta-2/beta). This bound decreases with eta and increases with beta. The relative difference of total and net free charge-particle counts is 2r/(1-r).',
        maximum_positron_to_electron_ratio=str(ratio),maximum_total_vs_net_count_relative_difference=str(relative),
        scope='A conditional ideal-gas value bound. It does not bound interacting-plasma pair corrections, opacity derivatives, EOS implicit-root errors or a stellar trajectory. Actual sampled argument-domain checks are separate.'))


def block(start):
    plan=json.loads((OUT/'plan.json').read_text());source=dict(np.load(o.OUT/'baseline-type1-captured.npz'))
    stop=min(start+plan['block_cells'],len(source['X']));stem=f'block-{start}'
    path=OUT/(stem+'.npz');record=OUT/(stem+'.json');binding=g.c.sha(OUT/'plan.json')
    if path.exists() and record.exists():
        saved=json.loads(record.read_text());assert saved['plan_sha256']==binding and saved['output_sha256']==g.c.sha(path)
        return saved
    eos=g.EOS();rows=[];errors=[];arguments=[]
    electron_mass=json.loads((g.d.OUT/'isotope-no-go.json').read_text())['finite_source_constant_electron_mass_amu']
    beta_factor=g.c.KB*g.c.NA/(electron_mass*(g.c.gr.C*100)**2)
    for i in range(start,stop):
        r,t=source['parameters'][i,3:5]*np.log(10);x=source['X'][i];seen=[]
        def sample(lr,lt):
            a=eos(2,lr,lt,x);assert a[13]>0
            seen.append([a[12],beta_factor*np.exp(lt)])
            return np.log(a[13])-lr
        value=sample(r,t);slopes=[]
        for h in plan['finite_log_steps']:
            slopes.append([(sample(r+h,t)-sample(r-h,t))/(2*h),(sample(r,t+h)-sample(r,t-h))/(2*h)])
        slopes=np.array(slopes);errors.append(float(np.max(abs(slopes[1]-slopes[0])/np.maximum(1,abs(slopes[1])))))
        rows.append([value,*slopes[1]]);seen=np.array(seen)
        arguments.append([seen[:,0].min(),seen[:,0].max(),seen[:,1].min(),seen[:,1].max()])
    np.savez_compressed(path,values=np.array(rows),derivative_scores=np.array(errors),sampled_eta_beta_ranges=np.array(arguments))
    passed=max(errors)<plan['derivative_relative_tolerance']
    saved=dict(classification='Counterexample candidate',start=start,stop=stop,passed=passed,
        maximum_derivative_score=max(errors),plan_sha256=binding,output_sha256=g.c.sha(path))
    save(record.name,saved);assert passed,saved;return saved


def compute():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    n=len(np.load(o.OUT/'baseline-type1-captured.npz')['X']);starts=list(range(0,n,plan['block_cells']));records=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,k) for k in starts]):
            records.append(done.result());save('progress.json',dict(completed_cells=sum(r['stop']-r['start'] for r in records),total_cells=n))
            print('OPACITY EOS INPUTS',len(records),'/',len(starts),'blocks',flush=True)
    pieces=[dict(np.load(OUT/f'block-{k}.npz')) for k in starts]
    values=np.concatenate([p['values'] for p in pieces]);arguments=np.concatenate([p['sampled_eta_beta_ranges'] for p in pieces])
    np.save(OUT/'replacement.npy',values)
    in_box=bool(np.all(arguments[:,0]>=-50) and np.all(arguments[:,2]>0) and np.all(arguments[:,3]<=.005))
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=n,
        blocks=sorted(records,key=lambda r:r['start']),maximum_derivative_score=max(r['maximum_derivative_score'] for r in records),
        sampled_eta_range=[float(arguments[:,0].min()),float(arguments[:,1].max())],
        sampled_beta_range=[float(arguments[:,2].min()),float(arguments[:,3].max())],
        all_sampled_points_in_ideal_pair_bound_box=in_box,
        replacement_sha256=g.c.sha(OUT/'replacement.npy'),physical_EOS_or_opacity_certified=False))
    print('OPACITY EOS INPUTS COMPLETE',n,'sampled pair-domain check',in_box,flush=True)


def apply():
    assert json.loads((OUT/'result.json').read_text())['passed']
    replacement=np.load(OUT/'replacement.npy');o.setup('new-electrons');o.trace('new-electrons',replacement)
    actual=dict(np.load(o.OUT/'new-electrons-captured.npz'));previous=dict(np.load(o.OUT/'baseline-type1-captured.npz'))
    assert np.array_equal(actual['used'],replacement)
    change=actual['outputs']-previous['outputs'];relative=change[:,0]/previous['outputs'][:,0]
    save('applied.json',dict(classification='Counterexample candidate',completed=True,
        actual_new_EOS_inputs_confirmed_bitwise=True,maximum_opacity_relative_change=float(abs(relative).max()),
        changed_opacity_cells=int(np.count_nonzero(change[:,0])),maximum_log_derivative_changes=abs(change[:,1:]).max(0).tolist(),
        scope='A fixed-state input intervention in the unchanged Type1 opacity function. The opacity tables and their plasma/composition assumptions remain the original ones; physical transport is not certified.'))
    print('APPLIED EOS OPACITY INPUTS',int(np.count_nonzero(change[:,0])),float(abs(relative).max()),flush=True)


def expanded_pair_domain():
    previous=json.loads((OUT/'result.json').read_text());assert previous['passed']
    assert not previous['all_sampled_points_in_ideal_pair_bound_box']
    assert not (OUT/'expanded-pair-bound.json').exists()
    iv.dps=60;eta=iv.mpf(-50);beta=iv.mpf('.006')
    ratio=iv.exp(-eta-2/beta)+iv.exp(-2*eta-2/beta);relative=2*ratio/(1-ratio)
    assert relative.b<iv.mpf('1e-100').a
    save('expanded-pair-bound.json',dict(classification='Proven',passed=True,
        original_bound_sha256=g.c.sha(OUT/'ideal-pair-bound.json'),
        domain=dict(eta_minimum=-50,beta_positive=True,beta_maximum='.006'),
        maximum_total_vs_net_count_relative_difference=str(relative),
        derivation='Apply the already proved ideal-gas inequality 2r/(1-r), r<=exp(-eta-2/beta)+exp(-2eta-2/beta), to this separately declared larger box.',
        scope='Conditional ideal-pair number bound only. Preserve the original smaller-box failure and all original scientific gates. No interacting-plasma or derivative certificate is inferred.'))
    pieces=[dict(np.load(OUT/f'block-{i}.npz'))['sampled_eta_beta_ranges'] for i in range(0,previous['cells'],128)]
    arguments=np.concatenate(pieces)
    contained=bool(np.all(arguments[:,0]>=-50) and np.all(arguments[:,2]>0) and np.all(arguments[:,3]<=.006))
    save('expanded-pair-domain-census.json',dict(classification='Counterexample candidate',passed=contained,
        source_result_sha256=g.c.sha(OUT/'result.json'),declared_samples_per_cell=9,cells=len(arguments),
        original_smaller_box_passed=False,physical_or_continuous_domain_certificate=False))
    assert contained
    print('EXPANDED CONDITIONAL PAIR BOX',len(arguments),'cells; relative bound <1e-100',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
