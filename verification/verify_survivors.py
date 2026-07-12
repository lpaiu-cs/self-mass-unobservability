"""Independent check of Lemma 09: the five Delta<=4 survivors
   {E2, E3, E2^2, dotE2, gradE2}
are LINEARLY independent operators under the allowed reductions, with the
gradient block modeled correctly as an STF-3 octupole.

My reconstruction, not their script. Building blocks treated as independent
jet data at a point:
  E_ij         : STF (sym, traceless) rank-2                  -> 5 comps
  DtE_ij       : STF rank-2 (independent tensor)              -> 5 comps
  G_kij=grad E : totally symmetric trace-free rank-3 (STF-3)  -> 7 comps

Why STF-3: G_kij = partial_k partial_i partial_j Phi is symmetric in ALL
three indices by equality of mixed partials (Schwarz -- kinematics, not an
optional identity), and trace-free on every pair in the external vacuum
(nabla^2 Phi = 0: the SAME condition already used to make E itself
traceless).  Under this model the two extra gradient contractions of the
M4-era catalog are not independent:
  divE2       = (G_iij)(G_kkj) == 0        (vacuum trace),
  mixedGradE2 = G_kij G_ikj    == gradE2   (total symmetry).

Operators (as polynomials in the independent components):
  E2   = E_ij E_ij                      (E-sector, deg 2)
  E3   = E_ij E_jk E_ki                  (E-sector, deg 3)
  E2^2 = (E_ij E_ij)^2                   (E-sector, deg 4)
  dotE2= DtE_ij DtE_ij                   (DtE-sector, deg 2)
  gradE2 = G_kij G_kij                   (gradE-sector, deg 2)

I verify:
  (i)   full 5-operator linear independence (rank 5),
  (ii)  the STF-3 identities hold exactly on samples:
        divE2 == 0 and mixedGradE2 == gradE2 to machine precision,
  (iii) as a DIAGNOSTIC control: re-running the old generic model
        (G[k] an independent STF-2 slice for each k, 15 comps) reproduces
        the defective gradient-sector rank 3 -- demonstrating that the
        M4-era count 7 was an artifact of the generic ansatz, not a
        property of grad E.
"""

import itertools

import numpy as np

rng = np.random.default_rng(271828)


def kron(a, b):
    return 1.0 if a == b else 0.0


def rand_stf():
    M = rng.standard_normal((3, 3))
    S = 0.5 * (M + M.T)
    S -= np.eye(3) * np.trace(S) / 3.0
    return S


def rand_stf3():
    """Random totally symmetric trace-free rank-3 tensor (STF-3, 7 dof)."""
    T = rng.standard_normal((3, 3, 3))
    S = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        S += np.transpose(T, p)
    S /= 6.0
    t = np.einsum("aac->c", S)
    d = np.eye(3)
    S -= (
        np.einsum("ab,c->abc", d, t)
        + np.einsum("bc,a->abc", d, t)
        + np.einsum("ca,b->abc", d, t)
    ) / 5.0
    return S


def rand_gradE_generic():
    """The OLD (defective) model: G[k] an independent STF-2 slice, free k."""
    G = np.zeros((3, 3, 3))
    for k in range(3):
        G[k] = rand_stf()
    return G


def ops_from_sample():
    E = rand_stf()
    DtE = rand_stf()
    G = rand_stf3()
    E2 = np.einsum("ij,ij->", E, E)
    E3 = np.einsum("ij,jk,ki->", E, E, E)
    E2sq = E2**2
    dotE2 = np.einsum("ij,ij->", DtE, DtE)
    gradE2 = np.einsum("kij,kij->", G, G)
    return [E2, E3, E2sq, dotE2, gradE2]


def grad_triple(G):
    gradE2 = np.einsum("kij,kij->", G, G)
    divE = np.einsum("iij->j", G)
    divE2 = np.einsum("j,j->", divE, divE)
    mixedGradE2 = np.einsum("kij,ikj->", G, G)
    return [gradE2, divE2, mixedGradE2]


def linear_rank(sampler, n=400):
    """Rank of the value matrix [op_j at sample_s]; = # linearly independent
    homogeneous operators."""
    M = np.array([sampler() for _ in range(n)])
    return np.linalg.matrix_rank(M, tol=1e-9)


if __name__ == "__main__":
    print("=" * 72)
    print("Lemma 09 independent reconstruction (corrected STF-3 gradient block)")
    print("=" * 72)
    r_full = linear_rank(ops_from_sample)
    print(f"  full 5-operator value-matrix rank : {r_full}   (claim: 5)")

    # (ii) exact kinematic identities on STF-3 samples
    max_div, max_mixed_dev = 0.0, 0.0
    for _ in range(200):
        g2, d2, m2 = grad_triple(rand_stf3())
        max_div = max(max_div, abs(d2))
        max_mixed_dev = max(max_mixed_dev, abs(m2 - g2))
    print(f"  STF-3 identity check over 200 samples:")
    print(f"    max |divE2|                 : {max_div:.3e}   (claim: 0, vacuum trace)")
    print(f"    max |mixedGradE2 - gradE2|  : {max_mixed_dev:.3e}   (claim: 0, Schwarz)")

    r_grad = linear_rank(lambda: grad_triple(rand_stf3()))
    print(f"  gradient triple rank, STF-3 model   : {r_grad}   (claim: 1)")

    # (iii) diagnostic control: the old generic model reproduces the defect
    r_grad_generic = linear_rank(lambda: grad_triple(rand_gradE_generic()))
    print(f"  gradient triple rank, generic model : {r_grad_generic}   "
          f"(old defective count: 3)")

    # E-sector subset {E2,E3,E2^2} should be rank 3
    def esec():
        E = rand_stf()
        E2 = np.einsum("ij,ij->", E, E)
        E3 = np.einsum("ij,jk,ki->", E, E, E)
        return [E2, E3, E2**2]
    r_e = linear_rank(esec)
    print(f"  E-sector subset {{E2,E3,E2^2}} rank : {r_e}   (claim: 3)")

    print()
    ok = (
        r_full == 5
        and r_grad == 1
        and r_e == 3
        and max_div < 1e-12
        and max_mixed_dev < 1e-12
        and r_grad_generic == 3
    )
    print(f"  Lemma 09 five-scalar independence reproduced: {ok}")
    print(f"  (gradient sector is 1-dimensional; the old rank-3 count is an")
    print(f"   artifact of the generic (STF-2 x vector) ansatz, shown above)")
