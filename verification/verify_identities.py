"""Independent check of the three load-bearing algebraic reduction identities
used by the Delta<=4 collapse-theorem basis (docs/reduction-rules.md).

These are the ONLY genuinely algebraic (non-bookkeeping) reductions; the rest
are total-derivative/EOM rules. If these hold for generic STF 3x3 tensors, the
'E4', 'B4', 'EBEB' reductions are valid. I test them on random symmetric
trace-free 3x3 tensors, both numerically (many samples) and symbolically
(exact, one generic pair).
"""

import numpy as np
import sympy as sp

np.random.seed(11235)


def rand_stf():
    M = np.random.randn(3, 3)
    S = 0.5 * (M + M.T)
    S -= np.eye(3) * np.trace(S) / 3.0
    return S


def num_checks(n=20000):
    print("=" * 72)
    print("Numerical check over", n, "random STF 3x3 pairs (E,B):")
    print("=" * 72)
    d_e4 = d_b4 = d_ebeb = 0.0
    for _ in range(n):
        E, B = rand_stf(), rand_stf()
        # (1) Tr E^4 = 1/2 (Tr E^2)^2
        lhs = np.trace(E @ E @ E @ E)
        rhs = 0.5 * np.trace(E @ E) ** 2
        d_e4 = max(d_e4, abs(lhs - rhs))
        # (2) Tr B^4 = 1/2 (Tr B^2)^2
        lhs = np.trace(B @ B @ B @ B)
        rhs = 0.5 * np.trace(B @ B) ** 2
        d_b4 = max(d_b4, abs(lhs - rhs))
        # (3) Tr(EBEB) = (E:B)^2 + 1/2 Tr(E^2) Tr(B^2) - 2 Tr(E^2 B^2)
        lhs = np.trace(E @ B @ E @ B)
        EB = np.sum(E * B)  # E:B = Tr(EB)
        rhs = EB**2 + 0.5 * np.trace(E @ E) * np.trace(B @ B) - 2 * np.trace(E @ E @ B @ B)
        d_ebeb = max(d_ebeb, abs(lhs - rhs))
    print(f"  max|Tr E^4 - 1/2(Tr E^2)^2|          = {d_e4:.3e}")
    print(f"  max|Tr B^4 - 1/2(Tr B^2)^2|          = {d_b4:.3e}")
    print(f"  max|Tr(EBEB) - [ (E:B)^2 + ...     ]| = {d_ebeb:.3e}")
    ok = d_e4 < 1e-9 and d_b4 < 1e-9 and d_ebeb < 1e-9
    print(f"  all three identities hold numerically: {ok}")
    return ok


def sym_check_ebeb():
    print()
    print("=" * 72)
    print("Symbolic (exact) check of the mixed quartic identity:")
    print("=" * 72)
    # generic symmetric 3x3, then subtract trace to make STF
    def sym_stf(prefix):
        a, b, c, d, e, f = sp.symbols(f"{prefix}0:6", real=True)
        M = sp.Matrix([[a, b, c], [b, d, e], [c, e, f]])
        M = M - sp.eye(3) * M.trace() / 3
        return M
    E = sym_stf("e")
    B = sym_stf("b")
    lhs = (E * B * E * B).trace()
    EB = (E * B).trace()
    rhs = EB**2 + sp.Rational(1, 2) * (E * E).trace() * (B * B).trace() - 2 * (E * E * B * B).trace()
    diff = sp.simplify(lhs - rhs)
    print(f"  simplify( Tr(EBEB) - RHS ) = {diff}")
    print(f"  mixed quartic identity exact: {diff == 0}")
    return diff == 0


if __name__ == "__main__":
    a = num_checks()
    b = sym_check_ebeb()
    print()
    print("=" * 72)
    print(f"Tr E^4 / Tr B^4 identities (numeric)   : {a}")
    print(f"mixed quartic Tr(EBEB) identity (exact): {b}")
    print("=" * 72)
