"""Phase 292: compare the withdrawn (eta = e cos varpi) and corrected (ELL1, eta = e sin varpi) phase-dependent outputs.

Asserts that every phase-independent field is unchanged and prints the corrected values that the manuscript reports.
"""
import json, math
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
R = root/'outputs/research-completion'; W = R/'withdrawn-periastron-convention'
load = lambda p: json.loads(p.read_text(encoding='utf-8'))
out = {}
def same(a, b):
    """Equal up to floating-point rounding (BLAS order), relative 1e-9."""
    if isinstance(a, dict): return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list): return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float): return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-15)
    return a == b

# 1. Phase matching.
o, n = load(W/'physical-matching.json'), load(R/'physical-matching.json')
for k in ('checks', 'archived_auxiliary_phase_radians', 'archived_closure_radians', 'input_sha256', 'matching'):
    assert same(o[k], n[k]), k
out['matching'] = {k: n[k] for k in ('pericenter_radians', 'leading_physical_phase_radians', 'physical_minus_archived_radians', 'physical_closure_radians', 'archived_closure_radians')}
out['matching']['pericenter_degrees'] = {k: math.degrees(v) for k, v in n['pericenter_radians'].items()}

# 2. Corrected drive and fixed-array fits.
o, n = load(W/'corrected-physical-drive.json'), load(R/'corrected-physical-drive.json')
for k in ('masses_solar', 'semi_major_m', 'f', 'eccentricities', 'potential_factors', 'amplitudes', 'Ustar', 'normalized_amplitudes', 'mass_o_parameter_difference'):
    assert same(o['drive'][k], n['drive'][k]), k
out['drive'] = dict(varpi=n['drive']['varpi'], phases=n['drive']['phases'], torus=n['torus_checks'], torus_old=o['torus_checks'], fits=n['fits'])

# 3. Comparator audit: only rows with the physical phase stress may change.
o, n = load(W/'comparator-audit.json'), load(R/'comparator-audit.json')
assert [k for k in o] == [k for k in n]
for k in o:
    if k not in ('rows', 'positive_fast_spectrum_witnesses'): assert same(o[k], n[k]), k
changed = {r.get('phase') or r.get('phases') or r.get('phase_label') for r0, r in zip(o['rows'], n['rows']) if not same(r0, r)}
print('comparator rows', len(n['rows']), 'changed phase labels', changed)
rows = n['rows']
print('row keys', sorted(rows[0]))
out['comparator_rows'] = rows; out['witnesses'] = n['positive_fast_spectrum_witnesses']; out['witnesses_old'] = o['positive_fast_spectrum_witnesses']
out['comparator_rows_old'] = o['rows']

# 4. Simultaneous region: simulations identical, data sections new.
o, n = load(W/'simultaneous-validation.json'), load(R/'simultaneous-validation.json')
for k in ('status', 'seed', 'nsim', 'rows', 'scope', 'threshold', 'nominal_chi6', 'signal_translation_control', 'calibration_sha256'):
    assert same(o[k], n[k]), k
for k in ('a', 'scale', 'omnibus_null_statistic', 'exceeds_threshold'):
    assert same(o['data'][k], n['data'][k]), k
out['sections'] = n['data']['physical_lag_sections']; out['threshold'] = n['threshold']; out['omnibus'] = n['data']['omnibus_null_statistic']

# 5. Live analysis: derivative, basis, nonlinear and rejection records identical.
o, n = load(W/'runtime12-analysis.json'), load(R/'runtime12-analysis.json')
for k in ('status', 'scope', 'derivative_rows', 'halfstep_basis', 'nonlinear_gap', 'original_displacement_rejection'):
    assert same(o[k], n[k]), k
for a, b in zip(o['transient'], n['transient']):
    for k in ('tau', 'weighted_derivative_change', 'nuisance_projected_derivative_change', 'convergence_5percent_pass'): assert same(a[k], b[k]), k
ratios = [f['sigma_beta_ratio'] for t in n['transient'] for f in t['fits']]
old_ratios = [f['sigma_beta_ratio'] for t in o['transient'] for f in t['fits']]
proj = [r['new_sigma_over_archived'] for r in n['physical_diagonal_projection_comparison']]
old_proj = [r['new_sigma_over_archived'] for r in o['physical_diagonal_projection_comparison']]
out['live'] = dict(transient_ratio=[min(ratios), max(ratios), len(ratios)], transient_ratio_old=[min(old_ratios), max(old_ratios)],
                   halfstep_sigma_factor=[min(proj), max(proj), len(proj)], halfstep_sigma_factor_old=[min(old_proj), max(old_proj)],
                   halfstep_rows=n['physical_diagonal_projection_comparison'], transient=n['transient'])
(Path(__file__).resolve().parent/'compare292.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
print('phase-independent fields unchanged; corrected values in compare292.json')
