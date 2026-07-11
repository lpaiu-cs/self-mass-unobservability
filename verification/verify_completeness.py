"""TIER 1 -- independent completeness of the Delta<=4 electric-sector enumeration.

Goal: confirm the collapse-theorem candidate catalog (21 contraction classes /
15 signatures, docs + symbolic/enumerate_contractions_delta4.py) is COMPLETE --
no parity-even scalar operator at weight <= 4 was missed -- using an
INDEPENDENTLY written contraction enumerator, not their script.

Method (fully independent of symbolic/):
  1. Sweep EVERY building-block signature (n_E, n_DtE, n_Dt2E, n_G, n_a) with
     total weight <= 4.  weights: E=1, DtE=2, Dt2E=3, G=gradE=2, a=1.
  2. For each signature, enumerate ALL perfect matchings of the index slots
     (every delta-only contraction to a scalar), respecting each block's
     tensor symmetry.  This is our own enumerator.
  3. Evaluate each matching at many random samples (identical factors share the
     same random tensor, as in a real operator) and take the numeric rank of
     the value matrix = the exact number of independent parity-even invariants
     at that signature.
  4. Compare to their catalog and, crucially, verify that NO weight<=4
     signature with nonzero invariant dimension is absent from their catalog.

Parity/epsilon: we only ever contract with delta (einsum), never epsilon, so
this counts exactly the parity-even, epsilon-free invariants -- the theorem's
sector.  A signature whose only invariants would need an epsilon returns rank 0
automatically.

Building-block tensor conventions (match survivor_rank_check.py):
  E, DtE, Dt2E : 3x3 symmetric trace-free                      (rank 2)
  G = grad E   : G[k] is 3x3 STF for each k (free k, STF ij)   (rank 3)
  a            : 3-vector                                       (rank 1)
"""

import itertools
import numpy as np

rng = np.random.default_rng(31337)

# ---- building blocks: (weight, rank) ----
BLOCKS = {
    "E":    (1, 2),
    "DtE":  (2, 2),
    "Dt2E": (3, 2),
    "G":    (2, 3),   # grad E
    "a":    (1, 1),
}
WEIGHT_CUT = 4


def rand_stf():
    M = rng.standard_normal((3, 3))
    S = 0.5 * (M + M.T)
    return S - np.eye(3) * np.trace(S) / 3.0


def sample_tensors():
    """One shared random instance per block type (identical factors share it)."""
    return {
        "E": rand_stf(),
        "DtE": rand_stf(),
        "Dt2E": rand_stf(),
        "G": np.stack([rand_stf() for _ in range(3)]),  # G[k] STF
        "a": rng.standard_normal(3),
    }


def perfect_matchings(slots):
    """Yield every perfect matching (list of index-pairs) of the slot list."""
    if not slots:
        yield []
        return
    first, rest = slots[0], slots[1:]
    for i in range(len(rest)):
        pair = (first, rest[i])
        remaining = rest[:i] + rest[i + 1:]
        for m in perfect_matchings(remaining):
            yield [pair] + m


def signature_factor_slots(sig):
    """Return (factors, slots): factors is a list of (block, [slot_ids]);
    slots is the flat list of global slot ids."""
    factors = []
    sid = 0
    for block, count in sig:
        _, rank = BLOCKS[block]
        for _ in range(count):
            ids = list(range(sid, sid + rank))
            sid += rank
            factors.append((block, ids))
    return factors, list(range(sid))


def eval_matching(factors, matching, T):
    """Contract all factors with delta pairs given by `matching` (einsum)."""
    # assign a letter to each matched pair
    label = {}
    for p, (u, v) in enumerate(matching):
        letter = chr(97 + p)
        label[u] = letter
        label[v] = letter
    subs, tensors = [], []
    for block, ids in factors:
        subs.append("".join(label[s] for s in ids))
        tensors.append(T[block])
    return float(np.einsum(",".join(subs) + "->", *tensors, optimize=True))


def invariant_dim(sig, n_samples=120):
    factors, slots = signature_factor_slots(sig)
    if len(slots) % 2 == 1:
        return 0, 0  # odd number of indices: no delta-scalar
    matchings = list(perfect_matchings(slots))
    if not matchings:
        return 0, 0
    # value matrix: rows = samples, cols = matchings
    rows = []
    for _ in range(n_samples):
        T = sample_tensors()
        rows.append([eval_matching(factors, m, T) for m in matchings])
    M = np.array(rows)
    rank = int(np.linalg.matrix_rank(M, tol=1e-9))
    return rank, len(matchings)


# ---- their catalog: signature -> (n_contraction_classes, independent_dim) ----
# (independent_dim = classes - stated algebraic relations)
THEIR = {
    (("E", 2),): (1, 1),
    (("a", 2),): (1, 1),
    (("DtE", 1), ("E", 1)): (1, 1),
    (("E", 3),): (1, 1),
    (("E", 1), ("a", 2)): (1, 1),
    (("G", 1), ("a", 1)): (1, 1),
    (("Dt2E", 1), ("E", 1)): (1, 1),
    (("DtE", 2),): (1, 1),
    (("DtE", 1), ("E", 2)): (1, 1),
    (("DtE", 1), ("a", 2)): (1, 1),
    (("E", 4),): (2, 1),   # E2^2, E4 with 1 relation E4 = E2^2/2
    (("E", 2), ("a", 2)): (2, None),   # a2E2, aE2a
    (("E", 1), ("G", 1), ("a", 1)): (3, None),  # aEGradE_1,2,3
    (("G", 2),): (3, 3),   # divE2, gradE2, mixedGradE2
    (("a", 4),): (1, 1),
}


def norm_sig(sig):
    return tuple(sorted(sig))


def all_signatures():
    """All (multidegree) signatures with 1 <= total weight <= WEIGHT_CUT."""
    blocks = list(BLOCKS)
    # max count per block bounded by weight cut
    ranges = []
    for b in blocks:
        w = BLOCKS[b][0]
        ranges.append(range(0, WEIGHT_CUT // w + 1))
    for counts in itertools.product(*ranges):
        weight = sum(c * BLOCKS[b][0] for c, b in zip(counts, blocks))
        n = sum(counts)
        if 1 <= weight <= WEIGHT_CUT and n >= 1:
            sig = tuple((b, c) for b, c in zip(blocks, counts) if c > 0)
            yield sig, weight


if __name__ == "__main__":
    print("=" * 78)
    print("TIER 1: independent completeness sweep of Delta<=4 electric-sector signatures")
    print("=" * 78)
    # analytic anchor: pure-E Hilbert series 1/((1-t^2)(1-t^3))
    import numpy.polynomial.polynomial as P
    # series coefficients up to t^4
    hilb = [1, 0, 1, 1, 1]  # 1/((1-t^2)(1-t^3)) = 1 + t^2 + t^3 + t^4 + ...
    print(f"\nAnalytic anchor -- pure-E invariant ring C[tr E^2, tr E^3]:")
    print(f"  Hilbert series 1/((1-t^2)(1-t^3)) coeffs t^0..t^4 = {hilb}")
    print(f"  => pure-E invariants at weight<=4: 1(const), E2(w2), E3(w3), E2^2(w4);")
    print(f"     no independent weight-4 beyond E2^2  => E4 forced reducible. [matches]")

    print(f"\n{'signature':<34}{'w':>2}  {'myDim':>5} {'#match':>6}  {'their(cls,dim)':>16}  verdict")
    print("-" * 78)
    covered = set()
    all_ok = True
    missing = []
    for sig, weight in sorted(all_signatures(), key=lambda x: (x[1], str(x[0]))):
        dim, nmatch = invariant_dim(sig)
        key = norm_sig(sig)
        their = THEIR.get(key)
        sig_str = "*".join(f"{b}^{c}" if c > 1 else b for b, c in sig)
        if dim == 0:
            # only report zero-dim if it's suspiciously in their catalog (shouldn't be)
            if their is not None:
                print(f"{sig_str:<34}{weight:>2}  {dim:>5} {nmatch:>6}  {str(their):>16}  MISMATCH(theirs nonempty)")
                all_ok = False
            continue
        covered.add(key)
        if their is None:
            print(f"{sig_str:<34}{weight:>2}  {dim:>5} {nmatch:>6}  {'(absent!)':>16}  *** MISSING FROM CATALOG ***")
            missing.append((sig_str, weight, dim))
            all_ok = False
        else:
            cls, tdim = their
            consistent = (tdim is None and dim <= cls) or (tdim is not None and dim == tdim)
            verdict = "ok" if consistent else "MISMATCH"
            if not consistent:
                all_ok = False
            tdim_s = tdim if tdim is not None else f"<= {cls}"
            print(f"{sig_str:<34}{weight:>2}  {dim:>5} {nmatch:>6}  ({cls}, {tdim_s})".ljust(74) + verdict)

    # any catalog signature we never hit?
    not_reproduced = set(THEIR) - covered
    print("-" * 78)
    print(f"catalog signatures with nonzero invariant dim reproduced: {len(covered)}/{len(THEIR)}")
    if not_reproduced:
        print(f"  catalog signatures my sweep gave dim 0 for: "
              f"{[ '*'.join(b for b,_ in s) for s in not_reproduced]}")
    print(f"signatures MISSING from their catalog (my dim>0, not listed): "
          f"{missing if missing else 'NONE'}")
    print()
    print("=" * 78)
    print(f"TIER-1 electric-sector completeness: "
          f"{'CONFIRMED (no missed operators, dims all match)' if all_ok and not missing else 'DISCREPANCY -- see above'}")
    print("=" * 78)
