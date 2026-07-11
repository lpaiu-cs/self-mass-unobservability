"""Independent verification of collapse-theorem checkable cores.

This is MY OWN reconstruction, not a rerun of their scripts. The point is to
reproduce the load-bearing quantitative claims from first principles and see
whether they agree.

Two things here:
  (A) The A5 non-analytic activation counterexample: e^{-1/Y^2} is C-infinity
      but non-analytic at 0, so its Taylor jet at 0 is trivial while the
      function is not. This is what "breaks the finite Taylor-jet collapse".
  (B) Method validation for the invariant counter: functionally-independent
      SO(3) scalar invariants of building-block tensors, counted by the
      generic rank of the Jacobian d(invariants)/d(components). For a single
      symmetric trace-free (STF) rank-2 tensor in 3D the answer MUST be 2
      (tr E^2, tr E^3), a textbook fact. If my counter returns 2 here, the
      method is trustworthy for the real Delta<=4 enumeration.
"""

import numpy as np
import sympy as sp

np.random.seed(20260708)


# ---------------------------------------------------------------------------
# (A) Non-analytic activation counterexample (A5 escape)
# ---------------------------------------------------------------------------
def check_nonanalytic():
    print("=" * 70)
    print("(A) A5 non-analytic activation:  f(Y) = exp(-1/Y^2),  f(0):=0")
    print("=" * 70)
    Y = sp.symbols("Y", positive=True)
    f = sp.exp(-1 / Y**2)
    # One-sided derivatives at 0^+ for n = 0..8
    print("  right-hand derivatives at Y -> 0+ :")
    all_zero = True
    fn = f
    for n in range(0, 9):
        lim = sp.limit(fn, Y, 0, dir="+")
        val = sp.nsimplify(lim)
        print(f"    f^({n})(0+) = {val}")
        if val != 0:
            all_zero = False
        fn = sp.diff(fn, Y)
    # The Taylor "jet" at 0 is therefore identically 0; but f(1)=e^{-1}!=0
    f_at_1 = float(sp.exp(-sp.Integer(1)))
    print(f"  Taylor jet at 0 is identically 0, yet f(1) = e^-1 = {f_at_1:.6f} != 0")
    print(f"  => all derivatives vanish at 0 : {all_zero}")
    print(f"  => finite Taylor-jet reconstruction FAILS (function != its jet): "
          f"{all_zero and f_at_1 != 0}")
    print("  VERDICT: A5 counterexample is genuine (C-inf, flat jet, nonzero fn).")
    return all_zero and f_at_1 != 0


# ---------------------------------------------------------------------------
# (B) Invariant-counter method validation on a single STF rank-2 tensor
# ---------------------------------------------------------------------------
def stf_from_params(p):
    """Symmetric trace-free 3x3 matrix from 5 real parameters."""
    a, b, c, d, e = p
    E = np.array([[a, b, c],
                  [b, d, e],
                  [c, e, -a - d]], dtype=float)  # traceless by construction
    return E


def random_rotation():
    # Haar-ish random rotation via QR of a Gaussian matrix
    A = np.random.randn(3, 3)
    Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def count_functionally_independent(invariant_fns, n_params, n_samples=40, tol=1e-8):
    """Generic rank of the Jacobian d(invariants)/d(params) via finite diff.

    invariant_fns: list of callables p (len n_params) -> scalar
    Returns the max rank over random sample points = number of functionally
    independent invariants.
    """
    h = 1e-6
    best = 0
    for _ in range(n_samples):
        p0 = np.random.randn(n_params)
        J = np.zeros((len(invariant_fns), n_params))
        for i, fn in enumerate(invariant_fns):
            for k in range(n_params):
                pp = p0.copy(); pp[k] += h
                pm = p0.copy(); pm[k] -= h
                J[i, k] = (fn(pp) - fn(pm)) / (2 * h)
        r = np.linalg.matrix_rank(J, tol=tol)
        best = max(best, r)
    return best


def check_single_E_rotinvariance():
    """tr(E^2), tr(E^3) must be rotation invariant; other contractions reduce."""
    print()
    print("=" * 70)
    print("(B) Method validation: SO(3) invariants of one STF rank-2 tensor E")
    print("=" * 70)
    # rotation invariance check
    p = np.random.randn(5)
    E = stf_from_params(p)
    max_dev = 0.0
    for _ in range(200):
        Q = random_rotation()
        Er = Q @ E @ Q.T
        for power in (2, 3):
            i0 = np.trace(np.linalg.matrix_power(E, power))
            i1 = np.trace(np.linalg.matrix_power(Er, power))
            max_dev = max(max_dev, abs(i0 - i1))
    print(f"  max |tr(E^k) - tr((QEQ^T)^k)| over 200 rotations, k=2,3: {max_dev:.2e}")

    # functional-independence count over the 5-param STF space
    inv = [
        lambda p: np.trace(np.linalg.matrix_power(stf_from_params(p), 2)),
        lambda p: np.trace(np.linalg.matrix_power(stf_from_params(p), 3)),
        # redundant candidates that MUST NOT add rank (tr E = 0, tr E^4 reducible):
        lambda p: np.trace(stf_from_params(p)),
        lambda p: np.trace(np.linalg.matrix_power(stf_from_params(p), 4)),
        lambda p: np.sum(stf_from_params(p) * stf_from_params(p)),  # = tr E^2
    ]
    n_indep = count_functionally_independent(inv, 5)
    print(f"  candidate scalars offered: tr E, tr E^2, tr E^3, tr E^4, E:E")
    print(f"  functionally independent count (Jacobian rank): {n_indep}")
    print(f"  EXPECTED for a single STF rank-2 tensor in 3D: 2  (tr E^2, tr E^3)")
    print(f"  method valid: {n_indep == 2}")
    # Cayley-Hamilton: for traceless 3x3, tr E^4 = (1/2)(tr E^2)^2  -- verify
    lhs = np.trace(np.linalg.matrix_power(E, 4))
    rhs = 0.5 * np.trace(np.linalg.matrix_power(E, 2)) ** 2
    print(f"  Cayley-Hamilton check tr E^4 = (1/2)(tr E^2)^2 : "
          f"{lhs:.6f} vs {rhs:.6f}  (match: {abs(lhs-rhs)<1e-9})")
    return n_indep == 2


if __name__ == "__main__":
    a = check_nonanalytic()
    b = check_single_E_rotinvariance()
    print()
    print("=" * 70)
    print(f"A5 counterexample genuine   : {a}")
    print(f"invariant-counter validated : {b}")
    print("=" * 70)
