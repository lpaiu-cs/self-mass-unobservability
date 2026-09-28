"""Bounded reproduction of Phase86 algebra, verdicts and recorded sources."""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/direct-eos-gr33/def-orbital-radiative-transport'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local(path):
    prefix = '/mnt/e/lab/self-mass-unobservability/'
    return ROOT/path[len(prefix):] if path.startswith(prefix) else ROOT/path


def main():
    a, b, c, x, y, light, gradient = s.symbols('a b c x y light gradient', positive=True)
    C = s.Matrix([[a+b, -b], [-b, c+b]])
    w = s.Matrix([x, y])
    first = -light/s.sqrt(3)*C.inv()*w*gradient
    flux = (light/s.sqrt(3)*w.T*first)[0]
    K = (light**2/3*w.T*C.inv()*w)[0]
    assert s.simplify(flux+K*gradient) == 0
    positive = light**2/3*((c+b)*x*x+2*b*x*y+(a+b)*y*y)/(a*c+a*b+b*c)
    assert s.simplify(K-positive) == 0
    for name in ['plan.json', 'audit-plan.json']:
        for path, digest in json.loads((OUT/name).read_text())['bindings'].items():
            assert sha(local(path)) == digest, path
    for name, source in [('first-plan.json', 'first-pilot-source.py'),
                         ('second-plan.json', 'second-pilot-source.py')]:
        binding = json.loads((OUT/name).read_text())['bindings']
        digest = next(v for k, v in binding.items() if k.endswith('/def_orbital_radiative_transport.py'))
        assert sha(OUT/source) == digest
    result = json.loads((OUT/'result.json').read_text())
    audit = json.loads((OUT/'audit.json').read_text())
    assert result['passed'] and audit['passed'] and audit['total_compute_seconds_charged'] < 600
    assert not result['stationary_background_certified'] and not result['full_goal_complete']
    for rows in result['cases'].values():
        for row in rows:
            assert row['physical_GR_residual'] < 1e-9
            assert row['heat_law_residual'] < 1e-9 and row['heat_balance'] < 2e-13
            assert row['adjoint_readout_relative'] < 1e-5
    projected = json.loads((OUT/'combined-projection.json').read_text())
    rows = result['cases']['all-p4']
    total = sum(r['delta_alpha_radiative_over_phi0'] for r in rows)
    assert total == projected['cases']['all-p4']['summed_thermal_amplitudes']
    residual = projected['cases']['all-p4']['fits'][2]['residual_norm']
    assert residual < projected['previous_adiabatic_control_change']
    spectral = audit['spectral_transport']
    assert np.isclose(spectral['K_spectral_cgs']/spectral['K_native_Rosseland_cgs'],
                      spectral['spectral_over_native'], rtol=1e-14)
    record = dict(classification='Proven', symbolic_passed=True,
        identity='The l=1 steady collision equation gives K=c^2/3*w^T*C1^-1*w; positive two-frequency absorption/reciprocal exchange has K>0.',
        scope='Algebraic model identity and saved-result/source consistency, not a new continuum, opacity, background or observation certificate.')
    (OUT/'verification.json').write_text(json.dumps(record, indent=2)+'\n')
    print('PASS Phase86 spectral identity, current/historical source bindings, original gates and incomplete physical boundary')


if __name__ == '__main__':
    main()
