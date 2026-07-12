"""Correct survivor-dimension check per family, matching the collapse-theorem
sector scripts EXACTLY: their block set is {E,DtE,Dt2E,GradE}+{X,DtX,Dt2X,GradX}
(NO Dt3/GradDt/GradDt2 as primitive blocks), and 'survivors' = block-set
invariants modulo total derivatives (D_tau^k of block-set operators) and the
lower-order EOM (a=0, so a excluded). Algebraic identities are automatic (we
compute true invariant dimensions).

survivor(w) = dim(inv(w))  -  dim( inv(w) ∩ TotalDerivImage(w) )
            = dim( inv(w) + Dimg(w) ) - dim( Dimg(w) )
with Dimg(w) = span{ D_tau(inv(w-1)), D_tau^2(inv(w-2)), D_tau^3(inv(w-3)) }.
D_tau images may use promoted blocks (Dt3E, GradDtE, ...), kept in the ambient
sample so the image is represented faithfully; inv(w) itself uses ONLY the
truncated block set.

Validated on the electric sector (must give 7) before use on B,S,V,T.
"""

import itertools
import numpy as np
from stf import rand_stf   # validated STF rank-L generator (dim 2L+1)

rng = np.random.default_rng(9091)
N_SAMPLES = 160
MAX_MATCH = 12000  # skip brute force above this many matchings (use Molien instead)

# --- SO(3) character integral for skipped big pure-primitive weight-4 sigs ---
_Nq = 4000
_th = (np.arange(_Nq) + 0.5) * np.pi / _Nq
_meas = (1.0 - np.cos(_th)) / np.pi * (np.pi / _Nq)
def _chi(l, th):
    return np.sin((2*l+1)*th/2.0) / np.sin(th/2.0)
def _sympow(pk, n):
    h = [np.ones_like(_th)]
    for m in range(1, n+1):
        acc = np.zeros_like(_th)
        for k in range(1, m+1):
            acc = acc + pk[k-1]*h[m-k]
        h.append(acc/m)
    return h[n]
def molien_dim(sig, blocks):
    """Parity-even invariant dim of a signature via O(3)=SO(3)xZ2 character integral.
    Only used for skipped big signatures (all pure-primitive, weight 4)."""
    par_count = sum(c*(1 if blocks[n][2]==-1 else 0) for n,c in sig)
    if par_count % 2 == 1:
        return 0
    prod = np.ones_like(_th)
    for name, c in sig:
        rank, w, par, kind, dt = blocks[name]
        irr = (rank,) if kind == 'stf' else tuple(sorted({rank-2, rank-1, rank}))
        pk = [sum(_chi(l, ((k*_th) % (2*np.pi))) for l in irr) for k in range(1, c+1)]
        prod = prod * _sympow(pk, c)
    return round(float(np.sum(prod * _meas)))

# block: name -> (rank, weight, parity, kind, dt_next)
# kind: 'stf' rank r, or 'grad' of stf rank r (rank r+1, free first idx)
def make_sector_blocks(families):
    """families: list of (sym, rank, parity). E is always included."""
    B = {}
    def fam(sym, rank, parity):
        B[sym]         = (rank,     1, parity,        'stf',  "Dt"+sym)
        B["Dt"+sym]    = (rank,     2, parity,        'stf',  "Dt2"+sym)
        B["Dt2"+sym]   = (rank,     3, parity,        'stf',  "Dt3"+sym)
        B["Dt3"+sym]   = (rank,     4, parity,        'stf',  None)      # ambient only
        B["Grad"+sym]  = (rank+1,   2, -parity,       'grad', "GradDt"+sym)
        B["GradDt"+sym]= (rank+1,   3, -parity,       'grad', None)      # ambient only
    fam("E", 2, +1)
    # CORRECTION (2026-07-12): GradE = partial^3 Phi (and GradDtE = partial^3
    # of the harmonic dtau Phi) are STF-3 octupoles: totally symmetric by
    # Schwarz and trace-free in the external vacuum.  Model them as true
    # STF rank-3 blocks (7 comps), not generic 'grad' objects (15 comps).
    # Admitted families' own Grad blocks stay generic ('grad'): independent
    # primitives with no assumed potential structure.
    B["GradE"]   = (3, 2, -1, 'stf', "GradDtE")
    B["GradDtE"] = (3, 3, -1, 'stf', None)      # ambient only
    for sym, rank, parity in families:
        if sym != "E":
            fam(sym, rank, parity)
    return B

# primitive (truncated) blocks that inv(w) is built from: exclude Dt3*, GradDt*
def primitive_names(blocks):
    return [n for n in blocks if not (n.startswith("Dt3") or n.startswith("GradDt"))]

def sample_tensors(blocks):
    T = {}
    for name, (rank, w, par, kind, dt) in blocks.items():
        if kind == 'stf':
            T[name] = rand_stf(rank, rng)
        else:  # grad: free first index k, STF rank (rank-1) in the rest
            base_rank = rank - 1
            T[name] = np.stack([rand_stf(base_rank, rng) for _ in range(3)])
    return T

def perfect_matchings(slots):
    if not slots:
        yield []; return
    first, rest = slots[0], slots[1:]
    for i in range(len(rest)):
        for m in perfect_matchings(rest[:i] + rest[i+1:]):
            yield [(first, rest[i])] + m

def factors_slots(sig, blocks):
    factors, sid = [], 0
    for name, cnt in sig:
        rank = blocks[name][0]
        for _ in range(cnt):
            factors.append((name, list(range(sid, sid+rank)))); sid += rank
    return factors, sid

def eval_matching(factors, matching, T):
    label = {}
    for p,(u,v) in enumerate(matching):
        label[u]=label[v]=chr(97+p)
    subs = [ "".join(label[s] for s in ids) for _,ids in factors ]
    tensors = [ T[name] for name,_ in factors ]
    return float(np.einsum(",".join(subs)+"->", *tensors, optimize=True))

def promote_factor_name(name, blocks):
    return blocks[name][4]  # dt_next

def matching_value_vectors(sig, blocks, samples, dtau=0):
    """Return list of value-vectors (one per matching) for this signature.
    dtau: apply D_tau this many times (promoting one factor per derivative, summed)."""
    factors, nslots = factors_slots(sig, blocks)
    if nslots % 2 == 1: return []
    matchings = list(perfect_matchings(list(range(nslots))))
    if len(matchings) > MAX_MATCH: return None  # too big
    vecs = []
    for m in matchings:
        col = []
        for T in samples:
            col.append(_eval_dtau(factors, m, T, blocks, dtau))
        vecs.append(col)
    return vecs

def _eval_dtau(factors, matching, T, blocks, dtau):
    if dtau == 0:
        return eval_matching(factors, matching, T)
    # D_tau = sum over factors of (that factor promoted); apply recursively
    total = 0.0
    for i,(name,ids) in enumerate(factors):
        nxt = promote_factor_name(name, blocks)
        if nxt is None:
            continue
        newf = list(factors); newf[i] = (nxt, ids)
        total += _eval_dtau(newf, matching, T, blocks, dtau-1)
    return total

def enum_sigs(names, blocks, wmax=4):
    caps = [range(0, wmax//blocks[n][1]+1) for n in names]
    for counts in itertools.product(*caps):
        tw = sum(c*blocks[n][1] for c,n in zip(counts,names))
        if 1 <= tw <= wmax:
            par = sum(c*(1 if blocks[n][2]==-1 else 0) for c,n in zip(counts,names))
            if par % 2 == 0:  # parity-even only
                yield tuple((n,c) for n,c in zip(names,counts) if c>0), tw

def rank_of(vec_list, tol=1e-7):
    if not vec_list: return 0
    M = np.array(vec_list)  # rows=matchings, cols=samples
    return int(np.linalg.matrix_rank(M, tol=tol))

def basis_vectors(vec_lists):
    """flatten list of matching-vector-lists into rows."""
    rows = []
    for vl in vec_lists:
        if vl: rows.extend(vl)
    return rows

def survivor_dim(families, wmax=4):
    blocks = make_sector_blocks(families)
    prim = primitive_names(blocks)
    samples = [sample_tensors(blocks) for _ in range(N_SAMPLES)]
    # inv(w): primitive-block invariants per weight (brute force); big sigs -> Molien
    inv_rows = {w: [] for w in range(0, wmax+1)}
    skip_survivors = 0
    skip_list = []
    def is_primitive(sig):
        return all(not (n.startswith("Dt") or n.startswith("Grad")) for n, _ in sig)
    for sig, tw in enum_sigs(prim, blocks, wmax):
        vv = matching_value_vectors(sig, blocks, samples, dtau=0)
        if vv is None:
            if is_primitive(sig):
                # pure-primitive => never a total derivative => full survivor; count via Molien
                d = molien_dim(sig, blocks)
                skip_survivors += d
                skip_list.append((dict(sig), d))
            else:
                # derivative-bearing big signature needed for inv/D_tau: compute-bound
                return None, {}, {}, ("compute-bound", dict(sig))
            continue
        inv_rows[tw].extend(vv)
    inv_dim = {w: rank_of(inv_rows[w]) for w in range(1, wmax+1)}
    total = 0
    per_w = {}
    for w in range(1, wmax+1):
        dimg = []
        bail = False
        for k in range(1, w):
            for sig, tw in enum_sigs(prim, blocks, wmax):
                if tw != w-k: continue
                vv = matching_value_vectors(sig, blocks, samples, dtau=k)
                if vv is None:
                    if not is_primitive(sig):
                        bail = True; break
                    continue  # pure-primitive can't be a D_tau image target here
                if vv: dimg.extend(vv)
            if bail: break
        if bail:
            return None, {}, {}, ("compute-bound-dtau",)
        d_dimg = rank_of(dimg)
        d_union = rank_of(inv_rows[w] + dimg)
        per_w[w] = d_union - d_dimg
        total += per_w[w]
    total += skip_survivors
    return total, inv_dim, per_w, skip_list


# sector label -> (families list, their independent survivor rank)
# Targets: corrected family table (2026-07-12, STF-3 electric gradient block);
# cross-checked against the exact character method (tier1_survivor_exact.py).
SECTORS = [
    ("E (electric)",   [("E", 2, +1)],                            5),
    ("E/B (r2 magn.)", [("B", 2, -1)],                            16),
    ("E/B/S (es)",     [("B", 2, -1), ("S", 0, +1)],              30),
    ("E/V (r1)",       [("V", 1, +1)],                            15),
    ("E/T (r3)",       [("T", 3, +1)],                            17),
    ("E/Q (r4)",       [("Q", 4, +1)],                            23),
    ("E/U (r5)",       [("U", 5, +1)],                            17),
    ("E/Z (r6)",       [("Z", 6, +1)],                            21),
]

if __name__ == "__main__":
    print("="*78)
    print("Family survivor-dimension: independent D_tau-quotient vs their survivor rank")
    print("(big pure-primitive weight-4 signatures counted via O(3) character integral)")
    print("="*78)
    print(f"{'sector':<16}{'inv(w) 1..4':<20}{'survivorDim':>12}{'theirRank':>11}  verdict")
    print("-"*78)
    for label, fams, target in SECTORS:
        total, invd, per_w, skipped = survivor_dim(fams)
        if total is None:
            print(f"{label:<16}{'--':<20}{'compute-bound':>12}{target:>11}  (rank-{fams[-1][1]} brute force infeasible; their rank {target})")
            continue
        wv = [invd[w] for w in range(1,5)]
        verdict = "MATCH" if total==target else "*** MISMATCH ***"
        molien_part = [s for s,_ in skipped] if isinstance(skipped, list) else []
        note = f"  [+{sum(d for _,d in skipped)} via Molien: {molien_part}]" if molien_part else ""
        print(f"{label:<16}{str(wv):<20}{total:>12}{target:>11}  {verdict}{note}")
