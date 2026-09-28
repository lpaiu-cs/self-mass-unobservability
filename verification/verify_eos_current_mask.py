"""Audit final native masks and bound truncation of their frozen log weights."""
from fractions import Fraction as F
import json, sys
import numpy as np
import sympy as sp
from mpmath import iv
from interval_records import interval_text
import eos_current_mask as m

g=m.g;OUT=m.OUT


def read(path): return json.loads(path.read_text())


def symbolic():
    d,u,v,U,V=sp.symbols('d u v U V',real=True)
    mixture=(1-d)*U+d*V-((1-d)*u+d*v)**2
    decomposed=(1-d)*(U-u*u)+d*(V-v*v)+d*(1-d)*(u-v)**2
    assert sp.expand(mixture-decomposed)==0
    theta=sp.symbols('theta',real=True)
    w=[theta,2*theta+theta**2,-theta]
    Z=sum(sp.exp(x) for x in w);p=[sp.exp(x)/Z for x in w]
    mean=sum(pi*sp.diff(wi,theta) for pi,wi in zip(p,w))
    second=sum(pi*sp.diff(wi,theta,2) for pi,wi in zip(p,w))
    second+=sum(pi*sp.diff(wi,theta)**2 for pi,wi in zip(p,w))-mean**2
    assert sp.simplify(sp.diff(sp.log(Z),theta)-mean)==0
    assert sp.simplify(sp.diff(sp.log(Z),theta,2)-second)==0
    # Hard-cutoff left and right limits differ even though the jump is tiny.
    iv.dps=70;cut=iv.exp(-288);jump=cut/(1+cut)
    assert jump.a>0
    m.save('derivative-boundary.json',dict(classification='Proven',passed=True,
        fixed_mask_conditions='The same nonempty retained set on a neighborhood, twice differentiable current log weights w_i(theta), positive exact exponentials and finite derivative oscillations over all retained and omitted stages.',
        value='0 <= log Z - log Z_retained = log(1+rho) <= rho.',
        gradient='|d(log Z-log Z_retained)/dtheta| <= delta * osc(w_prime).',
        Hessian='|d2(log Z-log Z_retained)/dtheta2| <= delta*osc(w_second) + delta*(5/4-delta)*osc(w_prime)^2. Here delta is the actual omitted probability. The nondecreasing envelope 5*delta_upper/4 may replace the second coefficient for a bound valid for any delta<=delta_upper.',
        proof='First and second derivatives are the mean score and the mean score derivative plus score variance. The distribution splits into retained and omitted conditional distributions. Their variances are in [0,osc(score)^2/4]. The verified mixture identity adds delta*(1-delta) times the squared difference of conditional means. Total variation controls the first- and second-score mean differences.',
        hard_cutoff=dict(log_weights=['0','theta'],remove_second_when='theta<=-288',
            left_value='0',right_limit_jump=interval_text(jump),
            conclusion='The returned second-stage probability is discontinuous at the threshold. Its classical derivative and Hessian do not exist there. Small truncation amplitude does not supply a global C2 guarantee for the raw hard-masked implementation.'),
        scope='Conditional softmax calculus and a hard-cutoff no-go boundary. Actual EOS log-weight derivatives, their evaluation errors, coupled molecular/nonideal equilibrium and global root conditioning are not certified.'))


def run():
    assert not (OUT/'audit.json').exists();plan=read(OUT/'plan.json');symbolic()
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in read(OUT/'build.json')['sha256'].items(): assert g.c.sha(m.s.Path(path))==digest,path
    assert g.c.sha(m.s.LIB)==plan['original_library_sha256']
    table=dict(np.load(g.OUT/'initial-adiabats-17.npz'));count=0;zeros=0;groups=0;reused=0
    max_gap=None;max_count=0;witness=None;max_probability_error=0.;replay_rows=[]
    for start in range(0,plan['cells'],plan['block_cells']):
        record=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert record['start']==count and record['stop']==min(count+128,plan['cells'])
        assert record['plan_sha256']==g.c.sha(OUT/'plan.json') and record['output_sha256']==g.c.sha(path)
        a=dict(np.load(path));oldpath=m.s.OUT/f'block-{start}.npz';old=dict(np.load(oldpath))
        assert g.c.sha(oldpath)==read(m.s.OUT/f'block-{start}.json')['output_sha256']
        for k in old: assert np.array_equal(a[k],old[k]),(start,k)
        assert np.array_equal(a['eos'],table['reference'][start:record['stop']])
        for i in plan['control_cells']:
            if start<=i<record['stop']:
                control=dict(np.load(OUT/f'control-{i}.npz'))
                assert all(np.array_equal(a[k][i-start],v) for k,v in control.items())
        for offset,active_elements in enumerate(a['active']):
            reused+=int(a['reused'][offset]);assert a['ifnr'][offset] in [0,1,3]
            for element,active in enumerate(active_elements):
                if not active.any(): continue
                logs=a['logweights'][offset,element,active];mask=a['mask'][offset,element,active]
                populations=a['number_fractions'][offset,element,active]
                assert np.all(np.isfinite(logs)) and (~mask).any() and np.array_equal(mask,populations==0)
                retained=logs[~mask].astype(np.longdouble);weights=np.exp(retained-retained.max());weights/=weights.sum()
                saved=populations[~mask].astype(np.longdouble);saved/=saved.sum()
                max_probability_error=max(max_probability_error,float(abs(saved-weights).max()))
                if not mask.any(): continue
                gap=F(float(logs[mask].max()))-F(float(logs[~mask].max()));n=int(mask.sum())
                assert gap<0 and n<=28
                groups+=1;zeros+=n;max_count=max(max_count,n)
                if max_gap is None or gap>max_gap: max_gap=gap;witness=dict(cell=start+offset,element=element)
                replay_rows.append([start+offset,element,n,str(gap)])
        count=record['stop']
    assert count==5735 and max_probability_error<1e-12
    assert zeros==read(m.s.OUT/'audit.json')['exact_zero_defined_stage_entries']
    iv.dps=70;gap=iv.mpf(max_gap.numerator)/max_gap.denominator
    rho=max_count*iv.exp(gap);delta=rho/(1+rho)
    assert delta.b<iv.mpf('3e-124').a
    m.save('frozen-logweight-certificate.json',dict(classification='Proven',passed=True,
        cells=count,masked_element_groups=groups,omitted_stage_entries=zeros,
        maximum_current_log_gap=str(max_gap),gap_witness=witness,maximum_omitted_stages=max_count,
        omitted_retained_ratio_upper=interval_text(rho.b),omitted_probability_upper=interval_text(delta.b),
        L1_error_upper=interval_text((2*delta).b),
        proof='For every nonempty retained set, its weight sum is at least exp(max_retained_logweight). The omitted weight sum is at most omitted_count*exp(max_omitted_logweight). Every gap is computed as an exact rational difference of saved binary64 log weights. The maxima of the count and gap yield a conservative common bound, evaluated outward with interval arithmetic. No prior fresh-mask decision is assumed.',
        scope='Mathematical atomic softmax formed from the final saved binary64 log weights. This encloses the current reused-mask truncation only, not native log-weight evaluation, native exponentiation/normalization, molecular equilibrium, self-consistent EOS root or physical model error.'))
    m.save('gap-replay.json',dict(classification='Proven',columns=['cell','element','omitted_count','exact_log_gap'],rows=replay_rows))
    m.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,
        all_21_outputs_and_species_bitwise=True,serial_controls_bitwise=True,
        reused_mask_cells=reused,mask_zero_mismatches=0,maximum_retained_probability_error=max_probability_error,
        current_frozen_logweight_truncation_certified=True,physical_EOS_certified=False,full_GR_evolution=False))
    paths=[p for p in OUT.rglob('*') if p.is_file()]+[g.ROOT/'verification/eos_current_mask.py',g.ROOT/'verification/verify_eos_current_mask.py']
    m.save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        runtime=read(OUT/'build.json')['sha256'],original_library_sha256=plan['original_library_sha256']))
    print('PASS CURRENT MASK AUDIT',count,'cells;',zeros,'zeros; reused',reused,'probability upper',interval_text(delta.b),flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items(): assert g.c.sha(m.s.Path(path))==digest,path
    assert g.c.sha(m.s.LIB)==manifest['original_library_sha256']
    print('PASS CURRENT MASK',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
