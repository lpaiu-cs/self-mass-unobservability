"""Conditional cutoff bounds and a frozen-mask counterexample."""
from fractions import Fraction as F
import json, sys
from mpmath import iv
import sympy as sp
import eos_species_inventory as s
from interval_records import interval_text

OUT=s.g.OUT/'ionization-cutoff'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();iv.dps=70
    source=s.OUT/'ionize.f90';text=source.read_text()
    assert 'parameter :: arglim = 144._fp_kind' in text
    assert 'if(ifsame_under)' in text and 'if(arg.le.-arglim)' in text
    assert 'can (very rarely) zero what ordinarily would be the maximum' in text
    # At a fresh mask, every omitted weight is <= exp(-288)*the maximum.
    r=28;ratio=iv.exp(-288);raw=r*ratio;delta=raw/(1+raw)
    # For tiny delta, -(1-delta)log(1-delta)<=delta avoids cancellation.
    entropy=-delta*iv.ln(delta)+delta+delta*iv.ln(28)
    # Exact two-group identity: retained renormalization has L1 error 2 delta.
    R,O=sp.symbols('R O',positive=True);q=O/(R+O)
    assert sp.simplify((1-R/(R+O))+O/(R+O)-2*q)==0
    d,Hr,Ho=sp.symbols('d Hr Ho',positive=True)
    mixture=-(1-d)*sp.log(1-d)-d*sp.log(d)+(1-d)*Hr+d*Ho
    assert sp.expand(mixture-Hr-(-d*sp.log(d)-(1-d)*sp.log(1-d)+d*(Ho-Hr)))==0
    # A component removed at a previous call can dominate a later call.
    assert iv.exp(-300).b<ratio.a
    current=[F(1),F(1000)];true=[w/sum(current) for w in current];masked=[F(1),F(0)]
    tv=sum(abs(a-b) for a,b in zip(true,masked))/2;assert tv==F(1000,1001)
    save('result.json',dict(classification='Proven',passed=True,checkpoint='8f878bb',
        fresh_mask=dict(maximum_atomic_stages=29,maximum_omitted_stages=r,relative_weight_cutoff=interval_text(ratio),
            omitted_probability_upper=interval_text(delta.b),L1_error_upper=interval_text((2*delta).b),
            charge_expectation_error_upper=interval_text((28*delta).b),charge_square_expectation_error_upper=interval_text((28**2*delta).b),
            dimensionless_mixing_entropy_error_upper=interval_text(entropy.b),
            proof='At least one maximum weight remains. If the omitted/retained weight ratio is <=r*exp(-288), the omitted probability is <=r*exp(-288)/(1+r*exp(-288)). Retained renormalization changes the distribution in total variation by exactly that omitted probability. A bounded observable changes by at most its range times total variation. Fannes entropy bound is <=-delta*ln(delta)+delta+delta*ln(dimension-1).',
            conditions='The current mask must be freshly computed from current log weights, with the exact declared log-cutoff rule; the atomic softmax has at most 29 stages. This is not automatically a bound on coupled hydrogen molecular equilibrium or a self-consistent nonideal root.'),
        fixed_mask_counterexample=dict(previous_weights=['1','exp(-300)'],current_weights=['1','1000'],
            reused_mask=[False,True],current_total_variation=str(tv),
            conclusion='An old mask alone supplies no current small-error guarantee. The source explicitly permits reusing such a mask and notes that the current maximum may be removed.'),
        restoration_condition='For a reused mask, certify the current omitted-to-retained weight ratio rho directly; then delta<=rho/(1+rho). A previous threshold alone is insufficient. Population/entropy value bounds do not supply free-energy derivative or equilibrium-location bounds without the corresponding derivative/curvature assumptions.',
        source_sha256=s.g.c.sha(source),actual_EOS_mask_current_ratio_certified=False,physical_EOS_certified=False))
    paths=[OUT/'result.json',s.g.ROOT/'verification/ionization_cutoff.py',source]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(s.g.ROOT).as_posix():s.g.c.sha(p) for p in paths}))
    print('PASS FRESH CUTOFF BOUNDS; REUSED MASK COUNTEREXAMPLE',str(tv),'delta upper',interval_text(delta.b),flush=True)


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items(): assert s.g.c.sha(s.g.ROOT/rel)==digest,rel
    print('PASS IONIZATION CUTOFF BINDINGS',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
