"""Phase269 stage 1: same-state EOS comparison of the PL/MHD and MHD-only FreeEOS builds in the charge-dominating layers.

Conjectural. Each build is evaluated in its own process through the phase-66 check harness (GasEOS bound to the build's
gas.so bridge, InventoryEOS.snapshot at fixed composition X). States are the original-grid cells 8-15 of the driven
chain's interior bank (initial background). Equilibrium ionization; the chain's kinetic hydrogen makes this a proxy.
Stage 2 adds the chain-lineage builds (orig, identity, ploff) with the H I excitation trace L and, through the matching
level library, the chain's H level fractions (def_native_hydrogen_exchange.Native.state).
Usage: python3 phase269-eos.py probe <build>
       python3 phase269-eos.py eval <build> <out npz>           build = shared|mhd|orig|identity|ploff
       python3 phase269-eos.py compare <a npz> <b npz>
       python3 phase269-eos.py identity <a npz> <b npz>
"""
import json, signal, sys
from types import FunctionType
import numpy as np
sys.path.insert(0, 'verification')
BRIDGES = {'shared': '/home/lpaiu/work/direct-eos-gr33/photon-shared-atomic/gas.so', 'mhd': '/home/lpaiu/work/direct-eos-gr33/photon-mhd-spectrum/gas.so',
           'orig': '/home/lpaiu/work/direct-eos-gr33/native-cold-population/stable/gas.so',
           'identity': '/home/lpaiu/work/direct-eos-gr33/native-cold-population-identity/stable/gas.so',
           'ploff': '/home/lpaiu/work/direct-eos-gr33/native-cold-population-ploff/stable/gas.so'}
LEVELS = {'orig': '/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-repaired/levels.so', 'identity': '/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-identity/levels.so',
          'ploff': '/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-ploff/levels.so'}
GEOMETRY = 'outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz'
CELLS = list(range(8, 16))


def gas_for(bridge):
    import def_photon_shared_atomic as a
    prior = a.p.a.previous; gas_module = prior.old.previous.matter.old
    gas = object.__new__(gas_module.GasEOS); init = gas_module.GasEOS.__init__
    FunctionType(init.__code__, dict(init.__globals__, BRIDGE=bridge), closure=init.__closure__)(gas)
    gas.inventory_lib = gas.gas_lib; native = gas.gas_lib.ionization_inventory
    def call(mode, value, t, eps, out, info):  # as def_photon_shared_atomic_check.run
        raw = np.full(24, np.nan); native(mode, value, t, eps, raw, info)
        out[:] = raw[:22]; out[20] = 0.; gas.molecules = raw[22:].copy()
    gas.call = call
    d, _ = prior.base.inputs()
    return gas, prior, d


def derived(v):
    """rho, P, u, s, chi_rho, chi_T, du/dlnrho, du/dlnT, Gamma1, cv*T."""
    chi_r, chi_t, cvT = v[5], v[6], v[10]
    return np.array([v[0], v[1], v[2], v[3], chi_r, chi_t, v[9], cvT, chi_r + chi_t**2*v[1]/(v[0]*cvT), cvT])


mode = sys.argv[1]
if mode in ('probe', 'eval'):
    signal.alarm(600)
    gas, prior, d = gas_for(BRIDGES[sys.argv[2]]); X = d['X'][0]
    g = np.load(GEOMETRY)
    if mode == 'probe':
        print('X', X.shape, np.round(X[:6], 6), 'sum', X.sum(), 'lnd0', float(d['lnd'][0]))
        snap = prior.inventory_reader.InventoryEOS.snapshot(gas, float(np.log(g['rho'][11])), float(np.log(g['T'][11])), X)
        print({k: np.shape(v) for k, v in snap.items()}); print('eos', np.asarray(snap['eos'])[:22])
    else:
        import ctypes
        lib = gas.gas_lib; levels = None
        if sys.argv[2] in LEVELS:  # the chain's H level split (def_native_hydrogen_exchange.Native.state)
            levels = ctypes.CDLL(LEVELS[sys.argv[2]]).photon_levels; arr_t = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
            levels.argtypes = [ctypes.c_int, ctypes.c_double, arr_t, arr_t]; levels.restype = None
        def block(name, n, dtype=ctypes.c_double, module='mod_excitation_block'):
            return np.ctypeslib.as_array((dtype*n).in_dll(lib, '__' + module + '_MOD_' + name)).copy()
        rows = {}
        for c in CELLS:
            lt = float(np.log(g['T'][c]))
            snap = prior.inventory_reader.InventoryEOS.snapshot(gas, float(np.log(g['rho'][c])), lt, X)
            rows[c] = {k: np.asarray(v) for k, v in snap.items()}
            try:  # optical ground-level log weights exist only in the shared-atomic lineage
                rows[c]['ground_logw'] = np.ctypeslib.as_array((ctypes.c_double*318).in_dll(lib, '__mod_excitation_MOD_shared_ground_logw')).copy()
            except ValueError:
                rows[c]['ground_logw'] = np.full(318, np.nan)
            count = int(block('extrace_count', 1, ctypes.c_int)[0]); ids = block('extrace_ids', 636, ctypes.c_int).reshape(318, 2)[:count]
            row = np.flatnonzero(np.all(ids == [1, 0], axis=1)); x = block('x', 5)
            L = float((block('extrace_value', 318*6).reshape(318, 6)[:count, 3]/block('extrace_scale', 318)[:count])[row[0]]) if len(row) == 1 else np.nan
            rows[c].update(L=np.array(L), x=x)
            if levels is not None:
                terms = np.zeros(30); levels(1, lt, np.ascontiguousarray(x), terms); tail = terms.reshape(3, 10)[0]; t = tail - np.r_[tail[1:], 0.]
                fraction = np.zeros(10); fraction[0] = np.exp(-L)
                if L > 0: fraction[1:] = -np.expm1(-L)*t[1:]/t[1:].sum()
                rows[c].update(terms=terms, fraction=fraction)
            print(c, 'eos', np.round(np.asarray(snap['eos'])[:12], 8).tolist(), flush=True)
        np.savez_compressed(sys.argv[3], cells=np.array(CELLS), X=X, rho=g['rho'][CELLS], T=g['T'][CELLS],
                            **{f'{k}_{c}': v for c, r in rows.items() for k, v in r.items()})
elif mode == 'identity':
    a, b = np.load(sys.argv[2]), np.load(sys.argv[3]); bad = []
    for k in a.files:
        same = np.array_equal(a[k], b[k], equal_nan=True) if a[k].dtype.kind == 'f' else np.array_equal(a[k], b[k])
        if not same: bad.append(k)
    print(json.dumps(dict(identical=not bad and set(a.files) == set(b.files), differing=bad[:20], compared=len(a.files))))
elif mode == 'compare':
    a, b = np.load(sys.argv[2]), np.load(sys.argv[3]); names = ['rho', 'P', 'u', 's', 'chi_rho', 'chi_T', 'du_dlnrho', 'du_dlnT', 'Gamma1', 'cvT']
    out = {}
    for c in a['cells']:
        va, vb = derived(a[f'eos_{c}']), derived(b[f'eos_{c}'])
        rel = (vb - va)/np.where(abs(va) > 0, abs(va), 1.)
        out[int(c)] = {n: float(r) for n, r in zip(names, rel)}
        print(int(c), ' '.join('%s %+.2e' % (n, r) for n, r in zip(names[1:], rel[1:])))
    print(json.dumps(dict(max_abs_relative_cells_9_13={n: max(abs(out[c][n]) for c in range(9, 14)) for n in names[1:]})))
    # (P/rho - du/dlnrho) = chi_T P/rho is the combination the mechanics uses (adiabatic modulus); du/dlnrho = (P/rho)(1-chi_T)
    for c in a['cells']:
        va, vb = a[f'eos_{c}'], b[f'eos_{c}']
        comb = lambda v: v[1]/v[0] - v[9]
        ident = va[9]/((va[1]/va[0])*(1 - va[6]))
        fa, fb = a[f'number_fractions_{c}'], b[f'number_fractions_{c}']
        yh = lambda f: f[0, 0]/f[0, :2].sum()  # neutral hydrogen fraction
        yhe = lambda f: f[1, :3]/f[1, :3].sum()
        ga, gb = a[f'ground_logw_{c}'], b[f'ground_logw_{c}']; act = a[f'active_{c}'].ravel()[:318].astype(bool) if a[f'active_{c}'].size >= 318 else np.ones(318, bool)
        print(int(c), 'identity du/dlnrho=(P/rho)(1-chi_T) ratio %.12f' % ident, ' (P/rho-du/dlnrho) rel %+.2e' % ((comb(vb) - comb(va))/comb(va)),
              ' neutral H %.4e rel %+.2e' % (yh(fa), (yh(fb) - yh(fa))/yh(fa)), ' He stages', np.round(yhe(fa), 5).tolist(), 'rel', np.round((yhe(fb) - yhe(fa))/np.maximum(yhe(fa), 1e-300), 6).tolist(),
              ' ground_logw max |diff| %.3e (first 3: %s vs %s)' % (np.max(abs(gb - ga)), np.round(ga[:3], 6).tolist(), np.round(gb[:3], 6).tolist()))
    if f'L_{a["cells"][0]}' in a.files and f'L_{a["cells"][0]}' in b.files:
        for c in a['cells']:
            line = '%d  L %.6e -> %.6e (%+.2e)' % (c, a[f'L_{c}'], b[f'L_{c}'], (b[f'L_{c}'] - a[f'L_{c}'])/a[f'L_{c}'])
            if f'fraction_{c}' in a.files and f'fraction_{c}' in b.files:
                fa, fb = a[f'fraction_{c}'], b[f'fraction_{c}']
                line += '  H level fractions n=1..4 %s -> %s' % (np.array2string(fa[:4], precision=4), np.array2string(fb[:4], precision=4))
            print(line)
    # proxy check: the chain's own native state rows at the same cells (initial kinetic H = equilibrium)
    g = np.load(GEOMETRY)
    for c in a['cells']:
        raw, v = g['raw'][c], a[f'eos_{c}']
        print(int(c), 'chain raw[:4]', np.array2string(raw[:4], precision=10), 'shared eos[:4]', np.array2string(v[:4], precision=10))
