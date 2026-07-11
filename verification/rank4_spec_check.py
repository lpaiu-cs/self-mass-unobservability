"""Fast test (numpy einsum): 6 missing rank-4 survivor specs + their 15 ->
new-sector rank should be 18 (=> total 7+18 = 25). Each (signature, matching)
spec is turned into an einsum; identical block instances share one tensor."""
import numpy as np
from stf import rand_stf

rng = np.random.default_rng(7)
BRANK = {"E": 2, "DtE": 2, "Q": 4, "DtQ": 4, "GradQ": 5, "GradE": 3}

def sample():
    return {
        "E": rand_stf(2, rng), "DtE": rand_stf(2, rng),
        "Q": rand_stf(4, rng), "DtQ": rand_stf(4, rng),
        "GradE": np.stack([rand_stf(2, rng) for _ in range(3)]),   # (3,3,3)
        "GradQ": np.stack([rand_stf(4, rng) for _ in range(3)]),   # (3,3,3,3,3)
    }

def eval_spec(signature, rep, S):
    # assign a letter per edge; build each factor's subscript from its slot letters
    label = {}
    for e, (u, v) in enumerate(rep):
        label[u] = label[v] = chr(97 + e)
    subs, arrs = [], []
    for inst, name in enumerate(signature):
        r = BRANK[name]
        subs.append("".join(label[(inst, s)] for s in range(r)))
        arrs.append(S[name])
    return float(np.einsum(",".join(subs) + "->", *arrs, optimize=True))

# their 15 new-sector specs (copied from r4_survivor_rank_check._r4_new_specs)
THEIR = {
 "Q2": (("Q","Q"), (((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2)),((0,3),(1,3)))),
 "EQQ": (("E","Q","Q"), (((0,0),(1,0)),((0,1),(2,0)),((1,1),(2,1)),((1,2),(2,2)),((1,3),(2,3)))),
 "dotQ2": (("DtQ","DtQ"), (((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2)),((0,3),(1,3)))),
 "EQDtQ": (("DtQ","E","Q"), (((0,0),(1,0)),((0,1),(2,0)),((0,2),(2,1)),((0,3),(2,2)),((1,1),(2,3)))),
 "E2Q2": (("E","E","Q","Q"), (((0,0),(1,0)),((0,1),(1,1)),((2,0),(3,0)),((2,1),(3,1)),((2,2),(3,2)),((2,3),(3,3)))),
 "E2Q2_m1": (("E","E","Q","Q"), (((0,0),(1,0)),((0,1),(2,0)),((1,1),(3,0)),((2,1),(3,1)),((2,2),(3,2)),((2,3),(3,3)))),
 "E2Q2_m2": (("E","E","Q","Q"), (((0,0),(2,0)),((0,1),(2,1)),((1,0),(3,0)),((1,1),(3,1)),((2,2),(3,2)),((2,3),(3,3)))),
 "E2Q2_m3": (("E","E","Q","Q"), (((0,0),(2,0)),((0,1),(3,0)),((1,0),(2,1)),((1,1),(3,1)),((2,2),(3,2)),((2,3),(3,3)))),
 "divQ2": (("GradQ","GradQ"), (((0,0),(0,1)),((0,2),(1,1)),((0,3),(1,2)),((0,4),(1,3)),((1,0),(1,4)))),
 "gradQ2": (("GradQ","GradQ"), (((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2)),((0,3),(1,3)),((0,4),(1,4)))),
 "mixedGradQ2": (("GradQ","GradQ"), (((0,0),(1,1)),((0,1),(1,0)),((0,2),(1,2)),((0,3),(1,3)),((0,4),(1,4)))),
 "Q2^2": (("Q","Q","Q","Q"), (((0,0),(3,0)),((0,1),(3,1)),((0,2),(3,2)),((0,3),(3,3)),((1,0),(2,0)),((1,1),(2,1)),((1,2),(2,2)),((1,3),(2,3)))),
 "Q4_chain": (("Q","Q","Q","Q"), (((0,0),(2,0)),((0,1),(3,0)),((0,2),(3,1)),((0,3),(3,2)),((1,0),(2,1)),((1,1),(2,2)),((1,2),(2,3)),((1,3),(3,3)))),
 "Q4_bridge": (("Q","Q","Q","Q"), (((0,0),(2,0)),((0,1),(2,1)),((0,2),(3,0)),((0,3),(3,1)),((1,0),(2,2)),((1,1),(2,3)),((1,2),(3,2)),((1,3),(3,3)))),
 "Q4_tetra": (("Q","Q","Q","Q"), (((0,0),(1,0)),((0,1),(2,0)),((0,2),(3,0)),((0,3),(3,1)),((1,1),(2,1)),((1,2),(2,2)),((1,3),(3,2)),((2,3),(3,3)))),
}
NEW = {
 "EEQ": (("E","E","Q"), (((0,0),(2,0)),((0,1),(2,1)),((1,0),(2,2)),((1,1),(2,3)))),
 "QQQ": (("Q","Q","Q"), (((0,0),(1,0)),((0,1),(1,1)),((0,2),(2,0)),((0,3),(2,1)),((1,2),(2,2)),((1,3),(2,3)))),
 "E3Q": (("Q","E","E","E"), (((0,0),(1,0)),((0,1),(2,0)),((1,1),(2,1)),((0,2),(3,0)),((0,3),(3,1)))),
 "EQ3": (("E","Q","Q","Q"), (((0,0),(1,0)),((0,1),(2,0)),((1,1),(2,1)),((1,2),(3,0)),((1,3),(3,1)),((2,2),(3,2)),((2,3),(3,3)))),
 "EDtEQ": (("E","DtE","Q"), (((0,0),(2,0)),((0,1),(2,1)),((1,0),(2,2)),((1,1),(2,3)))),
 "GradEGradQ": (("GradE","GradQ"), (((1,0),(1,1)),((0,0),(1,2)),((0,1),(1,3)),((0,2),(1,4)))),
}

def rank_of(specs, n=80):
    rows = []
    for _ in range(n):
        S = sample()                      # ONE shared sample per row
        rows.append([eval_spec(*specs[l], S) for l in specs])
    sv = np.linalg.svd(np.array(rows), compute_uv=False)
    return int(np.sum(sv > 1e-9 * sv[0]))   # relative threshold on singular values

if __name__ == "__main__":
    print("nonzero per new spec:")
    for name, (sig, rep) in NEW.items():
        vals = [eval_spec(sig, rep, sample()) for _ in range(5)]
        print(f"  {name:11} nonzero={any(abs(v)>1e-9 for v in vals)}")
    r_their = rank_of(THEIR)
    r_all = rank_of({**THEIR, **NEW})
    print(f"\ntheir 15 new specs rank        = {r_their}   (expected 12)")
    print(f"their 15 + my 6 specs rank     = {r_all}   (target 18)")
    print(f"=> corrected total survivor dim = 7 + {r_all} = {7 + r_all}")
