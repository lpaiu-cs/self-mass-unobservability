# FreeEOS dispatch for the chain's call (read-only).
S=/home/lpaiu/work/direct-eos-gr33/molecular-spectral/source/src
grep -n "ifoption" $S/mod_free_eos.f90 | head -30
echo "=== free_eos header"
grep -n "subroutine free_eos" $S/*.f90 | head
awk 'NR>=560 && NR<=600' $S/mod_free_eos.f90
