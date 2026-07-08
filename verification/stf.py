"""General symmetric-trace-free (STF) rank-L tensor machinery in 3D.

Used to independently reproduce the collapse-theorem higher-rank family audits.
A rank-L STF tensor in 3D spans the spin-L irrep, dimension 2L+1. We build a
basis of the STF subspace as the nullspace of the trace map on symmetric
tensors, then generate random STF tensors as random combinations. The self-test
must return dims [1,3,5,7,9,11,13] for L=0..6 -- if it does, the generator is
correct (this is the load-bearing setup-correctness check).
"""

import itertools
import numpy as np


def _sym_projector_apply(A):
    """Fully symmetrize a rank-L array over all index permutations."""
    L = A.ndim
    out = np.zeros_like(A)
    perms = list(itertools.permutations(range(L)))
    for p in perms:
        out += np.transpose(A, p)
    return out / len(perms)


def stf_basis(L):
    """Return a list of full 3^L arrays forming a basis of the rank-L STF space."""
    if L == 0:
        return [np.array(1.0)]
    if L == 1:
        return [np.eye(3)[i] for i in range(3)]
    dim_full = 3 ** L
    # basis of symmetric tensors: symmetrize each standard basis array, dedup by rank
    sym_vecs = []
    seen = []
    for idx in itertools.product(range(3), repeat=L):
        A = np.zeros([3] * L)
        A[idx] = 1.0
        S = _sym_projector_apply(A)
        v = S.reshape(-1)
        if np.linalg.norm(v) < 1e-12:
            continue
        # keep only linearly independent ones
        if seen:
            M = np.array(seen)
            # residual after projecting onto span(seen)
            coef, *_ = np.linalg.lstsq(M.T, v, rcond=None)
            if np.linalg.norm(v - M.T @ coef) < 1e-9:
                continue
        seen.append(v)
        sym_vecs.append(S)
    # trace map on symmetric tensors: contract first two indices -> rank L-2
    # build matrix T: (dim sym) -> (dim of rank L-2 arrays), then nullspace
    if L >= 2:
        rows = []
        for S in sym_vecs:
            tr = np.einsum("ii...->...", S)  # contract axes 0,1
            rows.append(tr.reshape(-1))
        Tmat = np.array(rows).T  # (features, n_sym)
        # nullspace of Tmat (columns = coefficients in sym_vecs basis giving traceless)
        u, s, vh = np.linalg.svd(Tmat)
        tol = 1e-9 * max(Tmat.shape)
        null_mask = np.concatenate([s, np.zeros(vh.shape[0] - len(s))]) <= tol
        null_coefs = vh[null_mask]
    else:
        null_coefs = np.eye(len(sym_vecs))
    basis = []
    for c in null_coefs:
        arr = sum(c[k] * sym_vecs[k] for k in range(len(sym_vecs)))
        basis.append(arr)
    return basis


_BASIS_CACHE = {}


def rand_stf(L, rng):
    """Random STF rank-L tensor as a full 3^L numpy array."""
    if L not in _BASIS_CACHE:
        _BASIS_CACHE[L] = stf_basis(L)
    basis = _BASIS_CACHE[L]
    coefs = rng.standard_normal(len(basis))
    if L == 0:
        return float(coefs[0])
    return sum(coefs[k] * basis[k] for k in range(len(basis)))


def check_traceless(arr, L, tol=1e-9):
    if L < 2:
        return True
    tr = np.einsum("ii...->...", arr)
    return np.max(np.abs(tr)) < tol


if __name__ == "__main__":
    rng = np.random.default_rng(2024)
    print("STF rank-L generator self-test (expected dim = 2L+1):")
    ok = True
    for L in range(0, 7):
        basis = stf_basis(L)
        d = len(basis)
        expected = 2 * L + 1
        # verify a random STF tensor is symmetric and traceless
        A = rand_stf(L, rng)
        sym_ok = True
        if L >= 2:
            arr = A
            sym_ok = np.allclose(arr, _sym_projector_apply(arr), atol=1e-9)
        tl_ok = check_traceless(A, L)
        good = (d == expected) and sym_ok and tl_ok
        ok = ok and good
        print(f"  L={L}: dim={d:>2} (expect {expected:>2})  symmetric={sym_ok}  traceless={tl_ok}  {'OK' if good else 'FAIL'}")
    print(f"\nSTF generator valid: {ok}")
