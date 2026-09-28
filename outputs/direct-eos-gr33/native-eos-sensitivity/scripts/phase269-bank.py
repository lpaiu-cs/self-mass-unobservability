"""Phase267 step A: native interior banks on a declared edge list (the phase-119 bank, generalized).

Counterexample candidate. The phase-119 bank() is reused verbatim except for its grid line and output folder:
 identity : n=19, keep=15, phase-119 edges  -> must reproduce the phase-119 outputs bitwise;
 refined  : n=27, keep=8, cells 8-15 of the 19-cell grid split in two, thin cells 16-18 recomputed natively.
Kept cells copy the original 16-cell arrays; every other cell gets its own native snapshot, probes, thermal/rate
table and mechanics/spectral derivatives exactly as in phase 119.

 quad     : n=43, keep=8, cells 8-15 split in four by nested midpoints (contains the 2x edges bitwise),
            caps scaled to 35 new cells (phase 268).

 refined-all: the refined grid with every cell's EOS-derived arrays from the running native EOS (phase 269, PL-off
            libraries); kept-cell geometry unchanged and asserted.

Usage: python3 .phase269-bank.py <identity|refined|quad|refined-all> <output folder> [<audit json>]
"""
import inspect, json, os, sys, textwrap
from pathlib import Path
mode, folder = sys.argv[1], Path(sys.argv[2]); audit = sys.argv[3] if len(sys.argv) > 3 else None
opened = {}
if audit:
    def hook(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            p = os.path.abspath(os.fsdecode(args[0])); m = args[1] if len(args) > 1 else None
            e = opened.setdefault(p, dict(reads=0, writes=0)); e['writes' if isinstance(m, str) and any(c in m for c in 'wax+') else 'reads'] += 1
    sys.addaudithook(hook)
sys.path.insert(0, 'verification')
import numpy as np
import def_native_boundary_layer as bl
assert mode in ('identity', 'refined', 'quad', 'refined-all')
source = textwrap.dedent(inspect.getsource(bl.bank))
grid_old = "bg=chem.prior.Background();baseline=prior.Coupled();bd=baseline.bulk.d;n=19;keep=15"
edges_old = "edges=np.r_[original['edges'][:-1],bg.R+np.array([-400000.,-200000.,-100000.,0.])]"
assert source.count(grid_old) == 1 and source.count(edges_old) == 1
if mode in ('refined', 'refined-all'):
    source = source.replace(grid_old, grid_old.replace('n=19;keep=15', 'n=27;keep=8'))
    source = source.replace(edges_old, "e=original['edges'];sub=[v for k in range(8,15) for v in ((e[k]+e[k+1])/2,e[k+1])]+[(e[15]+bg.R-400000.)/2]\n"
                            "    edges=np.r_[e[:9],sub,bg.R+np.array([-400000.,-200000.,-100000.,0.])];assert len(edges)==n+1")
    # Resource caps scale with the number of new native cells (4 in phase 119, 19 here); no accuracy gate changes.
    # The first refined run stopped on the phase-119 cap of 1800 native calls (preserved in the phase-267 log).
    for a, b in [('native=chem.old.Native(cap=1800)', 'native=chem.old.Native(cap=9000)'), ('signal.alarm(75)', 'signal.alarm(400)')]:
        assert source.count(a) == 1, a; source = source.replace(a, b)
if mode == 'refined-all':
    # Phase 269 (PL-off EOS): the refined 2x grid, but every cell's EOS-derived arrays (target, y0, coefficients, thermo,
    # raw, thermal/rate table, mechanics and spectral derivatives) come from the running native EOS, so the thermal step's
    # native-y0 consistency check holds in every cell. The geometry of the kept cells 0-7 (r, a, B, rho, T, phi) stays the
    # original bank's and is still asserted; only the keep-equality asserts on EOS-derived keys are dropped.
    loops = [("    for j in range(keep,n):\n        rho,T,a=d['rho'][j],d['T'][j],d['a'][j];lt=np.log(T)",
              "    for j in range(n):\n        rho,T,a=d['rho'][j],d['T'][j],d['a'][j];lt=np.log(T)"),
             ("    for j in range(keep,n):\n        chem.setup(native,d,j);theta=.09", "    for j in range(n):\n        chem.setup(native,d,j);theta=.09"),
             ("    for key in ['r','a','B','rho','T','phi','raw','y0','thermo','target']:assert np.array_equal(d[key][:keep],original[key][:keep]),key\n"
              "    assert np.array_equal(raw[:keep],oldthermal['raw'][:keep]) and np.array_equal(rates[:keep],oldthermal['rates'][:keep])",
              "    for key in ['r','a','B','rho','T','phi']:assert np.array_equal(d[key][:keep],original[key][:keep]),key"),
             # resource caps scaled by 27/19 native cells (no accuracy gate changes)
             ('native=chem.old.Native(cap=9000)', 'native=chem.old.Native(cap=14000)'), ('signal.alarm(400)', 'signal.alarm(700)')]
    for a, b in loops:
        assert source.count(a) == 1, a[:60]; source = source.replace(a, b)
if mode == 'quad':
    source = source.replace(grid_old, grid_old.replace('n=19;keep=15', 'n=43;keep=8'))
    source = source.replace(edges_old, "e=original['edges'];quarter=lambda a,b:[(a+(a+b)/2)/2,(a+b)/2,((a+b)/2+b)/2,b]\n"
                            "    sub=[v for k in range(8,15) for v in quarter(e[k],e[k+1])]+quarter(e[15],bg.R-400000.)[:3]\n"
                            "    edges=np.r_[e[:9],sub,bg.R+np.array([-400000.,-200000.,-100000.,0.])];assert len(edges)==n+1")
    for a, b in [('native=chem.old.Native(cap=1800)', 'native=chem.old.Native(cap=18000)'), ('signal.alarm(75)', 'signal.alarm(800)')]:
        assert source.count(a) == 1, a; source = source.replace(a, b)
folder.mkdir(parents=True, exist_ok=False); (folder/'thermal-support').mkdir()
ns = dict(vars(bl), OUT=folder)
bl.tns['OUT'] = folder  # Thermal() at the end of bank() reads the geometry it has just written
exec(compile(source, bl.__file__ + '#phase269-' + mode, 'exec'), ns)
try:
    ns['bank']()
finally:
    if audit:
        root = os.path.abspath('.') + '/'
        rows = [dict(path=p.replace(root, ''), **e) for p, e in sorted(opened.items())
                if not ('/verification/' in p or p.endswith('.pyc') or '/usr/lib/' in p or 'site-packages' in p or p.startswith('/proc') or p.startswith('/dev'))]
        Path(audit).write_text(json.dumps(dict(mode=mode, folder=str(folder), files=rows), indent=1) + '\n')
