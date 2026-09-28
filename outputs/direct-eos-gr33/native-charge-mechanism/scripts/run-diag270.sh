# Phase270: record the baryon-conservation and frozen-memory diagnostics as files in the 4x runtime.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
cd /home/lpaiu/work/native-refined268-runtime || exit 1
python3 $S/phase270-baryon.py > .phase270-baryon.txt 2>&1 && python3 $S/phase270-memory.py > .phase270-memory.txt 2>&1 && wc -l .phase270-baryon.txt .phase270-memory.txt
