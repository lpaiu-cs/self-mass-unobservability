# FreeEOS option semantics in the chain's source (read-only): free_eos signature and the option dispatch.
S=/home/lpaiu/work/direct-eos-gr33/molecular-spectral/source/src
grep -n "subroutine free_eos(" -A3 $S/mod_free_eos.f90 | head -8
grep -n "ifoption.eq\|ifmodified.eq\|ifmodified.lt\|ifmodified.gt" $S/mod_free_eos.f90 | head -40
echo "=== ifpi meaning"
grep -n "ifpi" $S/mod_free_eos.f90 | head -30
echo "=== bridge call"
grep -n "call free_eos" /mnt/e/lab/self-mass-unobservability/outputs/direct-eos-gr33/gr-radiation-eos-split/gas-bridge.f90
