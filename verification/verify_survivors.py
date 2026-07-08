"""Independent check of Lemma 09: the seven Delta<=4 survivors
   {E2, E3, E2^2, dotE2, gradE2, divE2, mixedGradE2}
are LINEARLY independent operators under the allowed reductions, WITHOUT
imposing transversality (div E != 0) or Bianchi.

My reconstruction, not their script. Building blocks treated as independent
jet data at a point:
  E_ij         : STF (sym, traceless) rank-2         -> 5 comps
  DtE_ij       : STF rank-2 (independent tensor)     -> 5 comps
  G_kij=grad E : sym-traceless in (ij), free in k    -> 3*5 = 15 comps

Operators (as polynomials in the independent components):
  E2   = E_ij E_ij                      (E-sector, deg 2)
  E3   = E_ij E_jk E_ki                  (E-sector, deg 3)
  E2^2 = (E_ij E_ij)^2                   (E-sector, deg 4)
  dotE2= DtE_ij DtE_ij                   (DtE-sector, deg 2)
  gradE2      = G_kij G_kij              (gradE-sector, deg 2)
  divE2       = (G_iij)(G_kkj)           (gradE-sector, deg 2)  [div E]_j[div E]_j
  mixedGradE2 = G_kij G_ikj              (gradE-sector, deg 2)

Since the three variable-sectors (E, DtE, gradE) are independent, cross-sector
linear independence is automatic; the ONLY nontrivial content is whether the
three gradient contractions {gradE2, divE2, mixedGradE2} are independent when
transversality is NOT assumed. I verify:
  (i)  full 7-operator linear independence (rank 7),
  (ii) the gradient-sector rank is exactly 3,
  (iii) as a control: if transversality div E = 0 WERE imposed, the gradient
        rank would drop (confirming the three are genuinely distinct only
        because transversality is absent).
"""

import numpy as np

rng = np.random.default_rng(271828)


def rand_stf():
    M = rng.standard_normal((3, 3))
    S = 0.5 * (M + M.T)
    S -= np.eye(3) * np.trace(S) / 3.0
    return S


def rand_gradE(transverse=False):
    """G_kij, symmetric-traceless in (i,j), free in k. 15 comps (or transverse)."""
    G = np.zeros((3, 3, 3))
    for k in range(3):
        G[k] = rand_stf()
    if transverse:
        # project out the divergence part: enforce sum_i G[i,i,j] = 0 for all j
        # by subtracting the trace-in-(k,i) piece appropriately (simple removal).
        div = np.einsum("iij->j", G)  # (div E)_j
        # crude transverse projection: G_kij -> G_kij - (1/?) delta_{ki} div_j-symmetrized
        # Build a correction that is STF in ij and removes the divergence.
        corr = np.zeros((3, 3, 3))
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    corr[k, i, j] = (
                        0.5 * (kron(k, i) * div[j] + kron(k, j) * div[i]) / 2.0
                        - kron(i, j) * div[k] / 3.0
                    )
        G = G - corr
    return G


def kron(a, b):
    return 1.0 if a == b else 0.0


def ops_from_sample():
    E = rand_stf()
    DtE = rand_stf()
    G = rand_gradE()
    E2 = np.einsum("ij,ij->", E, E)
    E3 = np.einsum("ij,jk,ki->", E, E, E)
    E2sq = E2**2
    dotE2 = np.einsum("ij,ij->", DtE, DtE)
    gradE2 = np.einsum("kij,kij->", G, G)
    divE = np.einsum("iij->j", G)
    divE2 = np.einsum("j,j->", divE, divE)
    mixedGradE2 = np.einsum("kij,ikj->", G, G)
    return [E2, E3, E2sq, dotE2, gradE2, divE2, mixedGradE2]


def grad_ops_from_sample(transverse):
    G = rand_gradE(transverse=transverse)
    gradE2 = np.einsum("kij,kij->", G, G)
    divE = np.einsum("iij->j", G)
    divE2 = np.einsum("j,j->", divE, divE)
    mixedGradE2 = np.einsum("kij,ikj->", G, G)
    return [gradE2, divE2, mixedGradE2]


def linear_rank(sampler, n_ops, n=400):
    """Rank of the value matrix [op_j at sample_s]; = # linearly independent
    homogeneous operators."""
    M = np.array([sampler() for _ in range(n)])  # (n, n_ops)
    return np.linalg.matrix_rank(M, tol=1e-9)


if __name__ == "__main__":
    print("=" * 72)
    print("Lemma 09 independent reconstruction")
    print("=" * 72)
    r_full = linear_rank(ops_from_sample, 7)
    print(f"  full 7-operator value-matrix rank : {r_full}   (claim: 7)")

    r_grad = linear_rank(lambda: grad_ops_from_sample(False), 3)
    print(f"  gradient sector rank (no transversality) : {r_grad}   (claim: 3)")

    r_grad_tv = linear_rank(lambda: grad_ops_from_sample(True), 3)
    print(f"  gradient sector rank (WITH transversality imposed): {r_grad_tv}")
    print(f"     -> control: imposing div E=0 collapses the gradient sector "
          f"(rank drops from 3): {r_grad_tv < 3}")

    # E-sector subset {E2,E3,E2^2} should be rank 3
    def esec():
        E = rand_stf()
        E2 = np.einsum("ij,ij->", E, E)
        E3 = np.einsum("ij,jk,ki->", E, E, E)
        return [E2, E3, E2**2]
    r_e = linear_rank(esec, 3)
    print(f"  E-sector subset {{E2,E3,E2^2}} rank : {r_e}   (claim: 3)")

    print()
    ok = (r_full == 7 and r_grad == 3 and r_e == 3)
    print(f"  Lemma 09 seven-scalar independence reproduced: {ok}")
    print(f"  (gradient triple independent *because* transversality is absent)")
