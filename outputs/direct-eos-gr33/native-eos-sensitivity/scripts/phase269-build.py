"""Phase269 stage 2: PL-off (MDH-only) variants of the chain's native EOS library and hydrogen level library.

Counterexample candidate. The recorded phase-121-era stable link is replayed with only mod_free_eos and
mod_free_eos_detailed recompiled from source (CMake flags of the base build); the excitation object is reused.
 identity: unmodified sources; must reproduce the original library's outputs bitwise.
 ploff   : ifpi 3 -> -3 in the (ifoption 3, ifmodified 11) branch the chain calls (PL off, MDH pressure ionization,
           Coulomb, molecules and every other option kept) and the PL arrays zeroed before first use (the phase-68 fix
           for tainted allocations consumed by the excitation/molecular code when PL is off).
Levels: the recorded levels bridge rebuilt unmodified (identity) or with qstar_calc's PL flag off (ploff), as the
phase-68 MHD-only hydrogenic provider did.
Usage: python3 phase269-build.py native|levels identity|ploff
"""
import hashlib, json, shutil, subprocess, sys, time
from pathlib import Path
ROOT = Path('/home/lpaiu/work/direct-eos-gr33')
SOURCE = ROOT/'molecular-spectral/source/src'; BUILD = ROOT/'molecular-spectral/build/src'
COLD = ROOT/'native-cold-population'; LEV = ROOT/'photon-eos-levels-repaired'
RUNTIME = Path('/home/lpaiu/work/native-retained-tail-runtime/outputs/direct-eos-gr33')
kind, variant = sys.argv[1], sys.argv[2]; assert kind in ('native', 'levels') and variant in ('identity', 'ploff')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
log = []; patches = {}; began = time.monotonic()
def run(cmd, cwd):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=600)
    log.append(dict(command=cmd, returncode=p.returncode, output=(p.stdout + p.stderr)[-3000:])); assert p.returncode == 0, p.stderr[-3000:]
def patch(path, anchor, new, note):
    t = path.read_text(); assert t.count(anchor) == 1, (path, anchor[:60]); path.write_text(t.replace(anchor, new)); patches[path.name] = note

if kind == 'native':
    D = ROOT/f'native-cold-population-{variant}'; D.mkdir(); (D/'stable').mkdir()
    shutil.copyfile(COLD/'mod_excitation-stable.o', D/'mod_excitation-stable.o')
    for name in ['mod_free_eos.f90', 'mod_free_eos_detailed.f90', 'free_eos_detailed.f90']: shutil.copyfile(SOURCE/name, D/name)
    if variant == 'ploff':
        anchor = "       elseif(ifmodified.eq.11) then\n          ! EOS1 without radiation pressure\n          ifcoulomb = 5\n          ifpi = 3\n          ifrad = 0"
        patch(D/'mod_free_eos.f90', anchor, anchor.replace('ifpi = 3', 'ifpi = -3'), 'ifpi 3 -> -3 in (ifoption 3, ifmodified 11): Planck-Larkin off, MDH kept')
        anchor = "  tc2 = c2/t\n  ! calculate planck-larkin occupation probabilities and equilibrium\n"
        patch(D/'free_eos_detailed.f90', anchor, "  tc2 = c2/t\n  ! PL off still passes these arrays to excitation and molecular consumers (phase-68 fix).\n"
              "  plop=0._fp_kind;plopt=0._fp_kind;plopt2=0._fp_kind\n  dv_pl=0._fp_kind;dv_plt=0._fp_kind\n"
              "  ! calculate planck-larkin occupation probabilities and equilibrium\n", 'PL arrays zeroed before first use')
    flags = ['gfortran', '-cpp', '-DUSINGDLL', '-Dfree_eos_EXPORTS', '-I' + str(D), '-I' + str(BUILD), '-I' + str(SOURCE), '-O3', '-DNDEBUG', '-O3', '-fPIC']
    run(flags + ['-c', 'mod_free_eos_detailed.f90', '-o', 'mod_free_eos_detailed.o'], D)
    run(flags + ['-c', 'mod_free_eos.f90', '-o', 'mod_free_eos.o'], D)
    receipts = json.load(open(RUNTIME/'def-native-cold-population/stable-build.json'))['receipts']
    swap = {str(BUILD/'CMakeFiles/free_eos.dir/mod_free_eos.f90.o'): str(D/'mod_free_eos.o'),
            str(BUILD/'CMakeFiles/free_eos.dir/mod_free_eos_detailed.f90.o'): str(D/'mod_free_eos_detailed.o')}
    def move(w):  # the original cache path, matched at a path boundary (D itself starts with the same prefix)
        return w[:-len(str(COLD))] + str(D) if w.endswith(str(COLD)) else w.replace(str(COLD) + '/', str(D) + '/')
    for i in (1, 2, 3):  # 1: library link, 2: stable/gas.so bridge, 3: verbose diagnostic bridge
        cmd = [swap.get(w, w) for w in map(move, receipts[i]['command'])]
        if i == 1: assert all(v in cmd for v in swap.values()), 'both recompiled objects must replace the recorded ones'
        assert not any(w.endswith(str(COLD)) or str(COLD) + '/' in w for w in cmd), cmd
        run(cmd, D)
    outputs = [D/'libfree_eos_native_cold_stable.so', D/'stable/gas.so', D/'gas-stable-verbose.so']
    objects = {n: dict(new=sha(D/n), recorded=sha(BUILD/'CMakeFiles/free_eos.dir'/(n.replace('.o', '.f90.o')))) for n in ['mod_free_eos.o', 'mod_free_eos_detailed.o']}
else:
    D = ROOT/f'photon-eos-levels-{variant}'; D.mkdir()
    for p in LEV.glob('*.f90'): shutil.copyfile(p, D/p.name)
    shutil.copyfile(LEV/'bridge.F90', D/'bridge.F90')
    if variant == 'ploff':
        patch(D/'bridge.F90', 'call qstar_calc(1,.true.,.true.,z.eq.1,', 'call qstar_calc(1,.false.,.true.,z.eq.1,', 'qstar_calc PL flag off (MDH kept)')
    cmd = json.load(open(RUNTIME/'def-photon-eos-populations/levels-repaired/build.json'))['command']
    run(cmd, D); outputs = [D/'levels.so']; objects = {}
receipt = dict(classification='Counterexample candidate', kind=kind, variant=variant, directory=str(D), patches=patches, seconds=time.monotonic() - began,
               outputs={str(p): sha(p) for p in outputs}, objects=objects, log=log)
(D/'phase269-build.json').write_text(json.dumps(receipt, indent=1) + '\n')
print(json.dumps({k: v for k, v in receipt.items() if k != 'log'}))
