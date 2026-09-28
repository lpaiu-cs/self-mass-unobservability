"""Proven: bound a frozen scalar readout change, without subtracting two IVPs.

Reuse exact polynomial/interval helpers. The target is the declared rounded
coefficient problem, not physical EOS or causal scalar evolution.
"""
import argparse
import json
import time
from fractions import Fraction as Q
from pathlib import Path

import mpmath as mp
import numpy as np

import gr_driven_readout as old
import gr_scalar_global_regular as algebra

ROOT = old.ROOT
OUT = ROOT/'outputs/direct-eos-gr33/gr-driven-readout-bound'


def save(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def prepare():
    assert not OUT.exists()
    original = json.loads((old.OUT/'plan.json').read_text())
    for rel, digest in original['bindings'].items():
        assert old.accepted.e.digest(ROOT/rel) == digest, rel
    assert not json.loads((old.OUT/'result.json').read_text())['numerical_controls_passed']
    star = old.accepted.context['initialize'](None)
    models = []
    for refinement, step in [(1, 0), (1, 70), (2, 140), (4, 280)]:
        data, _ = old.snapshot(star, refinement, step)
        models.append(old.eos.scalar_model('frozen-bound', data, 2e-12))
    base = models[0]
    for model in models[1:]:
        assert np.array_equal(model['p'].x, base['p'].x)
        assert np.array_equal(model['p'].c, base['p'].c)
        assert model['M'] == base['M'] and model['R'] == base['R']
    OUT.mkdir()
    np.savez_compressed(OUT/'coefficients.npz', x=base['p'].x,
        a=base['p'].c, b=np.array([model['v'].c for model in models]),
        R=base['R'], M=base['M'], mu=base['mu'])
    paths = [Path(__file__), Path(algebra.__file__), old.OUT/'plan.json',
        old.OUT/'manifest.json', old.OUT/'result.json', OUT/'coefficients.npz']
    save('plan.json', dict(classification='Proven',
        bindings={p.relative_to(ROOT).as_posix():old.accepted.e.digest(p) for p in paths},
        target='Upper bound on |delta(alpha/phi_infinity)| for exact binary64 PCHIP coefficients and fixed R/M and Schwarzschild mu.',
        method='Exact rational Bernstein enclosures and moments; outward 60-digit interval normalization. Identical a, R and M verified across snapshots. No subtraction of separately integrated charges.',
        decision_floor=original['absolute_alpha_readout_agreement'],
        decision='If every bound is below the preregistered 1e-9 readout floor, the saved 0.421-second thermal transient cannot meet that target in this coefficient model. This is not an observational exclusion or a bound on other drives.',
        resource_budget=dict(cpu_workers=1, gpu=False, native_EOS_calls=0,
            maximum_wall_seconds=300, repeats=1,
            estimate='Earlier readout used 8.78 seconds and 187 MiB. Exact-rational bound runtime is unmeasured; hard 300-second timeout, no grid or duration increase.'),
        physical_EOS_certified=False, coefficient_generation_error_certified=False,
        causal_scalar_evolution=False, observational_closure=False))
    print('PREPARED exact-coefficient change bound', flush=True)


def run():
    began = time.monotonic()
    mp.iv.dps = 60
    plan = json.loads((OUT/'plan.json').read_text())
    assert not (OUT/'result.json').exists()
    for rel, digest in plan['bindings'].items():
        assert old.accepted.e.digest(ROOT/rel) == digest, rel
    data = np.load(OUT/'coefficients.npz')
    rat, iv = algebra.rational, algebra.iv
    x = list(map(rat, data['x']))
    amin, B1, B2 = None, [Q(0)]*4, [Q(0)]*4
    variation = [Q(0)]*3
    for i, (left, right) in enumerate(zip(x, x[1:])):
        width = right-left
        a = list(map(rat, data['a'][::-1, i]))
        lower = min(algebra.bernstein(a, width))
        assert lower > 0
        amin = lower if amin is None else min(amin, lower)
        coefficients = []
        for k in range(4):
            b = list(map(rat, data['b'][k, ::-1, i]))
            assert max(algebra.bernstein(b, width)) <= 0
            B1[k] += algebra.moment([-v for v in b], left, width, 1)
            B2[k] += algebra.moment([-v for v in b], left, width, 2)
            coefficients.append(b)
        for k in range(3):
            difference = [a-b for a, b in zip(coefficients[k+1], coefficients[0])]
            magnitude = max(map(abs, algebra.bernstein(difference, width)))
            variation[k] += magnitude*(right**3-left**3)/3
        if i % 512 == 0:
            assert time.monotonic()-began < plan['resource_budget']['maximum_wall_seconds']
    mu = iv(rat(data['mu']))
    L = -mp.iv.ln(1-2*mu)/(2*mu)
    scale = iv(rat(data['R']))/iv(rat(data['M']))
    kappa = [(a-b)/amin for a, b in zip(B1, B2)]
    bounds = [1/(1-2*iv(k)-L*iv(b)) for k, b in zip(kappa, B2)]
    assert all(bound.a > 0 for bound in bounds)
    rows = []
    for k, refinement in enumerate([1, 2, 4]):
        upper = scale*bounds[0]*bounds[k+1]*iv(variation[k])
        rows.append(dict(refinement=refinement,
            normalized_field_supremum_bound=algebra.interval_text(bounds[k+1]),
            weighted_absolute_potential_change_exact=str(variation[k]),
            absolute_readout_change_upper_bound=algebra.interval_text(upper),
            approximate_upper_bound=float(upper.b),
            below_registered_readout_floor=bool(upper.b < mp.iv.mpf(str(plan['decision_floor'])).a)))
    # Analytic flat, uniform-potential positive control, with a nonzero change.
    a, b = mp.iv.mpf('.09'), mp.iv.mpf('.09001')
    expected = mp.iv.tan(mp.iv.sqrt(a))/mp.iv.sqrt(a)-mp.iv.tan(mp.iv.sqrt(b))/mp.iv.sqrt(b)
    control_bound = (b-a)/3/((1-2*a/3)*(1-2*b/3))
    assert abs(expected).b < control_bound.a and abs(expected).a > 0
    result = dict(classification='Proven', completed=True, rows=rows,
        coefficient_scope='Exact frozen binary64 cubic coefficients; no enclosure of generating GR/EOS data.',
        same_metric_coefficient_radius_mass=True,
        maximum_volterra_contraction=float(max(kappa)),
        all_changes_below_registered_floor=all(row['below_registered_readout_floor'] for row in rows),
        analytic_nonzero_control_passed=True, identical_input_bound_exactly_zero=True,
        proof='The Volterra operator norm is <=kappa=(B1-B2)/amin. Unit-central solution norm <=1/(1-kappa), and exterior normalization >=(1-2*kappa-L*B2)/(1-kappa)>0. Thus normalized field norm <=1/(1-2*kappa-L*B2). With identical a,mu,R/M, the exact two-solution Wronskian gives |delta alpha| <= (R/M)*F0*F1*integral x^2*|delta b| dx.',
        elapsed_seconds=time.monotonic()-began,
        physical_EOS_certified=False, causal_scalar_evolution=False,
        driven_thermal_response=False, observational_closure=False)
    save('result.json', result)
    save('manifest.json', dict(sha256={p.relative_to(ROOT).as_posix():old.accepted.e.digest(p)
        for p in OUT.iterdir() if p.is_file() and p.name != 'manifest.json'}))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run'])
    globals()[parser.parse_args().command]()
