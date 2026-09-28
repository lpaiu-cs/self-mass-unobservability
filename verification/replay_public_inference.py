"""Replay the conditional six-coefficient inference with public JSON and NumPy.

--export rebuilds the small record from the private frozen arrays; it does not
run the timing engine, stellar evolution, or Monte Carlo calibration.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / 'outputs/research-completion/public-inference.json'
REFERENCE = ROOT / 'outputs/research-completion/simultaneous-validation.json'


def export():
    import simultaneous_inference as si
    from physical_drive_completion import drive, stencil

    g = si.setup()
    fit = si.scan(g, g['y'], g['y_rest'])
    precision = g['r6'].T @ fit['grams'][int(fit['index'][0])] @ g['r6'] / fit['scale'][0]**2
    d = drive()
    reference = json.loads(REFERENCE.read_text())
    lags = [row['tau'] for row in reference['data']['physical_lag_sections']]
    for tau in lags:
        assert np.allclose(drive_matrix(d['normalized_amplitudes'], d['phases'],
                                       [g['inp']['OMS'][k] for k in ['in', 'out', 'dif']], tau),
                           stencil(g, d, tau), rtol=1e-12, atol=1e-15)
    source_paths = ['outputs/research-completion/simultaneous-validation.json',
                    'outputs/research-completion/simultaneous-calibration.json',
                    'outputs/research-completion/physical-matching.json',
                    'verification/simultaneous_inference.py', 'verification/physical_drive_completion.py',
                    'verification/comparator_audit.py', 'verification/coverage_audit.py',
                    'verification/nuisance_audit.py', 'verification/estimated_covariance_audit.py',
                    'request10_external/scripts/sep_common.py']
    manifest = json.loads((ROOT / 'paper/revision-manifest.json').read_text())
    record = dict(
        classification='Imported from prior work',
        scope='Frozen conditional inference only; no derivative certification, coverage recalibration, or timing-engine replay.',
        basis_order=g['order'], basis_convention='response to cos(omega*t), then cos(omega*t+pi/2) = -sin(omega*t)',
        carrier_coefficients=np.linalg.solve(g['r6'], fit['coefficient'][:, 0]).tolist(),
        coefficient_precision=precision.tolist(), threshold=reference['threshold'],
        covariance_a=float(fit['a'][0]), residual_scale=float(fit['scale'][0]),
        drive=dict(normalized_amplitudes=d['normalized_amplitudes'], phases_radians=d['phases'],
                   angular_frequencies_per_day=[g['inp']['OMS'][k] for k in ['in', 'out', 'dif']],
                   Ustar=d['Ustar']), lags_days=lags,
        source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in source_paths},
        frozen_input_sha256={p: h for p, h in manifest['sha256'].items() if p.startswith('request10_external/')})
    RECORD.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def drive_matrix(amplitudes, phases, omegas, tau):
    instantaneous = np.asarray(amplitudes) * np.exp(1j * np.asarray(phases))
    delayed = instantaneous / (1 + 1j * np.asarray(omegas) * tau)
    return np.column_stack([np.column_stack([v.real, v.imag]).ravel()
                            for v in (instantaneous, delayed)])


def replay():
    record = json.loads(RECORD.read_text(encoding='utf-8'))
    for path, expected in record['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    reference = json.loads(REFERENCE.read_text())
    mean = np.asarray(record['carrier_coefficients'])
    precision = np.asarray(record['coefficient_precision'])
    assert mean.shape == (6,) and precision.shape == (6, 6)
    assert np.isfinite(mean).all() and np.isfinite(precision).all()
    assert np.allclose(precision, precision.T, rtol=1e-12, atol=0)
    root = np.linalg.cholesky(precision).T
    target = root @ mean
    omnibus = float(target @ target)
    assert abs(omnibus - reference['data']['omnibus_null_statistic']) < 1e-6
    assert record['threshold'] == reference['threshold']
    d = record['drive']
    rows = []
    assert record['lags_days'] == [r['tau'] for r in reference['data']['physical_lag_sections']]
    for tau, expected in zip(record['lags_days'], reference['data']['physical_lag_sections']):
        x = root @ drive_matrix(d['normalized_amplitudes'], d['phases_radians'],
                                d['angular_frequencies_per_day'], tau)
        coefficients = np.linalg.lstsq(x, target, rcond=None)[0]
        residual = target - x @ coefficients
        minimum = float(residual @ residual)
        static_residual = target - x[:, 0] * (x[:, 0] @ target) / (x[:, 0] @ x[:, 0])
        positive_minimum = minimum if coefficients[1] >= 0 else float(static_residual @ static_residual)
        assert abs(minimum - expected['minimum_statistic']) < 1e-6
        assert np.isclose(coefficients[1], expected['beta'], rtol=1e-5, atol=1e-16)
        assert (minimum > record['threshold']) == expected['empty']
        assert positive_minimum >= minimum - 1e-9
        rows.append(dict(tau_days=tau, minimum=minimum, beta=float(coefficients[1]),
                         beta_nonnegative_minimum=positive_minimum, empty=minimum > record['threshold']))
    print(json.dumps(dict(classification=record['classification'], scope=record['scope'],
                          threshold=record['threshold'], omnibus=omnibus, sections=rows), indent=2))
    print('PASS: public conditional-inference replay; six physical sections verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export', action='store_true', help='requires private frozen arrays; rebuild the public record')
    if parser.parse_args().export:
        export()
    replay()
