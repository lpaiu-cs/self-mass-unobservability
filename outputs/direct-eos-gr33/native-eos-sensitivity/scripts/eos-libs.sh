# Locate the FreeEOS builds and their gas bridges (chain native, phase-66 shared atomic PL/MHD, phase-68 MHD-only).
find /home/lpaiu -xdev \( -name 'libfree_eos_native_cold_stable.so' -o -name 'libfree_eos_mhd_spectrum.so' -o -name 'libfree_eos*shared*atomic*.so' -o -name 'gas-bridge.f90' \) 2>/dev/null | head -20
cd /home/lpaiu/work/native-retained-tail-runtime/verification || exit 1
python3 - <<'EOF'
import sys; sys.path.insert(0, '.')
import def_photon_shared_atomic as a, def_photon_mhd_spectrum as m, def_native_cold_population as c
print('shared atomic CACHE', a.CACHE, 'OUT', a.OUT)
print('mhd CACHE', m.CACHE, 'LIB', m.LIB)
print('cold CACHE', c.CACHE)
EOF
