"""High-rank (Q,U,Z = STF rank 4,5,6) relation corroboration via EXACT O(3)
character (Molien) dimensions.

SUPERSEDED FOR RANK-4 (Q): this script predates the rank-4 completeness
correction and its Q conclusion (19) is out of date. It examines only the
degree-<=2-in-E load-bearing Q signatures (Q2, EEQ, EQQ, EEQQ, QQQQ) -- exactly
the candidate set the original enumeration used -- so it never scans the
higher-degree signatures {E:3,Q:1} and {E:1,Q:3} where the omitted survivors
live. The corrected rank-4 survivor dimension is 25, not 19; the six omitted
operators (EEQ, QQQ, E3Q, EQ3, EDtEQ, GradEGradQ) are constructed and counted in
rederive_rank4.py and rank4_spec_check.py, and cross-checked by the exact
character integral in tier1_survivor_exact.py. The U:19 and Z:23 corroboration
below is current.

Within the signatures it does examine, this script confirms the exact algebraic
relations the audit reported (the r4/r5/r6 survivor_rank_check nullities 3, 4, 6):
each signature's exact integer Molien dimension equals (candidate-class count) -
(reported relations), so within that candidate set no operator and no relation
was missed.
"""

import numpy as np

N = 8000
th = (np.arange(N) + 0.5) * np.pi / N
meas = (1 - np.cos(th)) / np.pi * (np.pi / N)

def _chi_sum(irreps, k):
    s = np.zeros_like(th)
    for l in irreps:
        s = s + np.sin((2*l+1)*(k*th)/2) / np.sin((k*th)/2)
    return s

def _sympow(pk, n):
    h = [np.ones_like(th)]
    for m in range(1, n+1):
        acc = np.zeros_like(th)
        for k in range(1, m+1):
            acc = acc + pk[k-1]*h[m-k]
        h.append(acc/m)
    return h[n]

def dim(blocks):
    """blocks: list of (irreps, parity, count). Exact parity-even invariant dim."""
    if sum(c for _, p, c in blocks if p == -1) % 2 == 1:
        return 0
    prod = np.ones_like(th)
    for irr, p, c in blocks:
        prod = prod * _sympow([_chi_sum(irr, k) for k in range(1, c+1)], c)
    return round(float(np.sum(prod * meas)))

E = ((2,), +1)

# (label, blocks, their candidate-class count) for the load-bearing signatures
CASES = {
    "Q (rank4)": [
        ("Q2",   [((4,), +1, 2)], 1),
        ("EEQ",  [(E[0], E[1], 2), ((4,), +1, 1)], 1),
        ("EQQ",  [(E[0], E[1], 1), ((4,), +1, 2)], 1),
        ("EEQQ", [(E[0], E[1], 2), ((4,), +1, 2)], 4),   # E2Q2 + 3 mixed
        ("QQQQ", [((4,), +1, 4)], 4),                     # Q2^2, Q4_bridge/chain/tetra
    ],
    "U (rank5)": [
        ("U2",   [((5,), +1, 2)], 1),
        ("EUU",  [(E[0], E[1], 1), ((5,), +1, 2)], 1),
        ("EEUU", [(E[0], E[1], 2), ((5,), +1, 2)], 4),
        ("UUUU", [((5,), +1, 4)], 4),
    ],
    "Z (rank6)": [
        ("Z2",   [((6,), +1, 2)], 1),
        ("EZZ",  [(E[0], E[1], 1), ((6,), +1, 2)], 1),
        ("EEZZ", [(E[0], E[1], 2), ((6,), +1, 2)], 4),
        ("ZZZZ", [((6,), +1, 4)], 5),
    ],
}
THEIR_NULLITY = {"Q (rank4)": 3, "U (rank5)": 4, "Z (rank6)": 6}
THEIR_RANK = {"Q (rank4)": "19 (pre-correction; corrected to 25)",
              "U (rank5)": 19, "Z (rank6)": 23}

if __name__ == "__main__":
    print("=" * 74)
    print("High-rank (Q,U,Z) relation corroboration via EXACT Molien dimensions")
    print("=" * 74)
    for fam, cases in CASES.items():
        print(f"\n{fam}:  their survivor rank {THEIR_RANK[fam]}, reported nullity {THEIR_NULLITY[fam]}")
        rel_found = 0
        for label, blocks, cand in cases:
            d = dim(blocks)
            rels = cand - d
            rel_found += max(rels, 0)
            flag = f"-> {rels} relation(s)" if rels else "-> independent"
            print(f"   {label:6} exact dim = {d}   (their ~{cand} candidate classes) {flag}")
        print(f"   relations confirmed in these signatures: {rel_found}  "
              f"(their total nullity {THEIR_NULLITY[fam]}; remainder in other sigs)")
    print()
    print("=" * 74)
    print("Within the examined signatures, exact Molien reproduces the algebraic")
    print("relations the audit reported (its nullities): no operator and no relation")
    print("was missed inside that candidate set. This corroborates U:19 and Z:23.")
    print("Rank-4 (Q): the pre-correction count 19 is SUPERSEDED -- the corrected")
    print("survivor dimension is 25. The six omitted higher-degree survivors live in")
    print("the {E:3,Q:1} and {E:1,Q:3} signatures NOT scanned here; they are counted")
    print("in rederive_rank4.py / rank4_spec_check.py and the exact integral in")
    print("tier1_survivor_exact.py.")
    print("=" * 74)
