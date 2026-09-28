# Phase269 stage 2a: build the identity and PL-off libraries, evaluate every build at the charge-layer states, compare.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
tr -d '\r' < $S/phase269-build.py > .phase269-build.py; tr -d '\r' < $S/phase269-eos.py > .phase269-eos.py
for k in native levels; do for v in identity ploff; do
  python3 .phase269-build.py $k $v > .phase269-build-$k-$v.log 2>&1 || { echo "build $k $v failed"; tail -30 .phase269-build-$k-$v.log; exit 1; }
  tail -1 .phase269-build-$k-$v.log | cut -c1-700
done; done
for b in orig identity ploff mhd shared; do
  taskset -c 7 python3 .phase269-eos.py eval $b .phase269-s2-$b.npz > .phase269-s2-$b.log 2>&1 || { echo "eval $b failed"; tail -30 .phase269-s2-$b.log; exit 1; }
done
echo "== identity vs original (must be bitwise)"; python3 .phase269-eos.py identity .phase269-s2-orig.npz .phase269-s2-identity.npz
echo "== PL-off patch vs phase-68 MHD-only build"; python3 .phase269-eos.py compare .phase269-s2-mhd.npz .phase269-s2-ploff.npz 2>&1 | grep -v 'chain raw'
echo "== lineage baseline: phase-66 PL/MHD vs chain PL/MHD"; python3 .phase269-eos.py compare .phase269-s2-shared.npz .phase269-s2-orig.npz 2>&1 | grep -v 'chain raw' | head -12
echo "== PL effect in the chain lineage: original vs PL-off"; python3 .phase269-eos.py compare .phase269-s2-orig.npz .phase269-s2-ploff.npz 2>&1 | grep -v 'chain raw'
