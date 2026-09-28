# Record the phase-244 source component comparison (PL-off vs PL/MHD 2x) at the given macros. Usage: bash run-source-diff269.sh 61 62 ...
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
for n in "$@"; do python3 $S/source-diff269.py $n > /dev/null || { echo "source diff $n failed"; exit 1; }; done
ls -la /home/lpaiu/work/native-eos269-runtime/.phase269-source-diff-*.json | cut -c30-
