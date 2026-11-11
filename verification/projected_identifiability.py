"""Exact controls for the projected-separation method; no timing-engine inputs.

Run with SymPy. The dimensionless example is an analytic illustration, not a
new experiment on PSR J0337+1715. Existing empirical records are not modified.
"""
import json

import sympy as s


def main():
    z = s.symbols('z')
    epsilon = s.symbols('epsilon', positive=True)
    v = s.Matrix([s.Rational(1, 2), -s.Rational(1, 2),
                  s.Rational(1, 5), -s.Rational(2, 5)])
    # Coefficients multiply cos(w*t) and cos(w*t + pi/2), respectively.
    V = s.Matrix([[1, 0], [0, 1], [1, 0], [0, 2]])
    P = V * (V.T * V).inv() * V.T
    Q = s.eye(4) - P
    residual = Q * v
    checks = []

    def zero(name, expression):
        assert all(s.simplify(x) == 0 for x in s.Matrix(expression)), name
        checks.append(name)

    zero('orthogonal_comparator_projector', P * P - P)
    zero('symmetric_comparator_projector', P.T - P)
    zero('residual_orthogonal_to_comparator', V.T * residual)
    zero('exact_coefficient_distance', [residual.dot(residual) - s.Rational(117, 1000)])
    zero('identity_map_retains_separation', [v.dot(Q * v) - s.Rational(117, 1000)])
    zero('projected_map_erases_noninterpolable_difference', Q * (P * v))
    # This map is full rank for epsilon > 0 but attenuates only the difference.
    A = P + epsilon * Q
    zero('comparison_columns_unchanged', A * V - V)
    zero('ill_conditioned_information', [(Q * A * v).dot(Q * A * v)
                                        - s.Rational(117, 1000) * epsilon**2])
    zero('map_eigenvalues', [A.det() - epsilon**2])
    # A normalized weighted gap witness can still be unobservable after a map.
    zero('unit_witness_duality', [residual.dot(v) - residual.dot(residual)])
    # The exact degree-three interpolant for K=2, tau=1, beta=1, c_Y=0.
    polynomial = s.cancel((1 - (z*z + 1)*(z*z + 4)/10)/(1 + z))
    assert s.Poly(polynomial, z).degree() == 3
    checks.append('degree_three_interpolant')
    for frequency in (1, 2):
        for sign in (-1, 1):
            point = sign * s.I * frequency
            zero(f'interpolation_{sign * frequency}',
                 [polynomial.subs(z, point) - 1/(1 + point)])
    result = dict(
        classification='Proven',
        scope='Dimensionless exact algebra; no empirical inference or derivative certification.',
        frequencies=[1, 2], tau=1,
        comparator_powers=[0, 1],
        coefficient_distance_squared=str(residual.dot(residual)),
        information=dict(identity='117/1000', projection='0',
                         attenuated='117*epsilon**2/1000'),
        checks=checks,
    )
    print(json.dumps(result, indent=2))
    print(f'PASS: {len(checks)} exact projected-identifiability checks')


if __name__ == '__main__':
    main()
