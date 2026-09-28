"""Counterexample candidate: exact-time comparison for the saved fluid control.

The matrix exponential solves the same finite affine control, not the full
nonlinear GR equations. The original failed GR verdict remains unchanged.
"""
from pathlib import Path
import json
import os

import numpy as np
from scipy.linalg import eig, expm, solve

import gr_time_convergence_cause as control


def main():
    root = control.ROOT
    source = root/'outputs/direct-eos-gr33/gr-time-cause-central-128'
    output = source/'exact-time'
    assert not output.exists()
    for path, digest in json.loads((source/'manifest.json').read_text())['sha256'].items():
        assert control.gr.e.digest(source/path) == digest, path
    plan = json.loads((control.SOURCE/'plan.json').read_text())
    previous = json.loads((source/'plan.json').read_text())
    control.NC = previous['cells']
    output.mkdir()
    control.gr.e.write(output/'plan.json', dict(classification='Counterexample candidate',
        assumptions='Exact finite affine fluid control; binary64 matrix exponential. No full-GR or continuum certification.',
        sha256={str(p.relative_to(root)): control.gr.e.digest(p) for p in
                [Path(__file__), source/'manifest.json', source/'plan.json']}))
    with np.load(source/'operator.npz') as data:
        M, K, b = [data[k] for k in ['M', 'K', 'b']]
    A, force = solve(M, -K), solve(M, -b)
    momentum = b.reshape(-1, 3).copy()
    momentum[:, :2] = 0
    acceleration = solve(M, -momentum.ravel())
    n = len(b)
    augmented = np.zeros((n+2, n+2))
    augmented[:n, :n] = A
    augmented[:n, n:] = np.column_stack([force, acceleration])
    duration = plan['duration_seconds']
    exponential = expm(duration*augmented)
    half = expm(duration/2*augmented)
    scale = np.asarray(control.SCALE, float)
    exact = exponential[:n, n].reshape(-1, 3)*scale
    only_momentum = exponential[:n, n+1].reshape(-1, 3)*scale
    doubling = (half@half)[:n, n].reshape(-1, 3)*scale
    numerical_check = np.max(abs(exact-doubling), axis=0)/np.maximum(np.max(abs(exact), axis=0), 1e-300)
    assert np.max(numerical_check) < 1e-8, numerical_check
    with np.load(source/'trajectories.npz') as data:
        endpoints = {r: data[f'path_{r}'][-1] for r in previous['refinements']}
    records = []
    for r in [64, 128]:
        final, record, _ = control.integrate(M, K, b, plan, r)
        endpoints[r] = final[:16]
        records.append(record)
    errors, pairs = [], []
    last = None
    for r, endpoint in endpoints.items():
        error = np.max(abs(endpoint-exact[:16]), axis=0)
        row = dict(refinement=r, exact_time_error=error.tolist())
        if last is not None:
            row['order'] = np.log2(last/error).tolist()
        errors.append(row)
        last = error
    refinements = list(endpoints)
    for left, right in zip(refinements[:-1], refinements[1:]):
        pairs.append(dict(refinements=[left, right], errors=np.max(
            abs(endpoints[right]-endpoints[left]), axis=0).tolist()))
    values, vectors = eig(A)
    coefficient = solve(vectors, force.astype(complex))
    response = np.zeros(len(values), dtype=complex)
    nonzero = abs(values) > 1e-12
    response[nonzero] = coefficient[nonzero]*np.expm1(values[nonzero]*duration)/values[nonzero]
    response[~nonzero] = coefficient[~nonzero]*duration
    eigen_endpoint = (vectors@response).real.reshape(-1, 3)*scale
    eigen_error = np.max(abs(eigen_endpoint[:16]-exact[:16]), axis=0)/np.maximum(np.max(abs(exact[:16]), axis=0), 1e-300)
    contribution = abs(vectors[2]*response)*scale[2]
    modes = [dict(real=float(values[j].real), imag=float(values[j].imag),
                  endpoint_center_velocity_contribution=float(contribution[j]))
             for j in np.argsort(contribution)[-12:][::-1]]
    result = dict(classification='Counterexample candidate', fields=control.FIELDS,
        exponential_doubling_relative=numerical_check.tolist(), eigen_endpoint_relative=eigen_error.tolist(),
        exact_time_errors=errors, pair_errors=pairs, additional_paths=records,
        initial_momentum_residual=float(b[2]*scale[2]),
        momentum_source_only_difference_relative=(np.max(abs(exact[:16]-only_momentum[:16]), axis=0)/np.max(abs(exact[:16]), axis=0)).tolist(),
        minimum_mode_real=float(values.real.min()), maximum_mode_real=float(values.real.max()),
        dominant_endpoint_center_velocity_modes=modes,
        scope=previous['assumptions'], full_GR_time_convergence_fixed=False)
    np.savez_compressed(output/'endpoints.npz', exact=exact, momentum_source_only=only_momentum,
                        **{f'path_{r}': endpoint for r, endpoint in endpoints.items()})
    control.gr.e.write(output/'result.json', result)
    control.gr.e.write(output/'manifest.json', dict(sha256={p.name: control.gr.e.digest(p)
        for p in output.iterdir() if p.is_file()}))
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    assert set(os.sched_getaffinity(0)) <= set(range(16))
    main()
