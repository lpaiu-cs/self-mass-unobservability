# Which FreeEOS option (ifmodified / ifoption) do the chain's native bridge and the phase-66/68 bridges pass?
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
for f in $(find outputs/direct-eos-gr33 -maxdepth 3 -name '*bridge*.f90' 2>/dev/null | head -20); do
  echo "== $f"; grep -n -i 'call free_eos\|ifoption\|ifmodified' "$f" | head -6 | cut -c1-220
done
echo "== verification references to option numbers"
grep -n -o 'ifoption[^,;]*\|ifmodified[^,;]*' verification/def_native_cold_population.py verification/def_native_hydrogen_exchange.py verification/def_photon_shared_atomic.py 2>/dev/null | head -10
