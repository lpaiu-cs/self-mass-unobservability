#!/bin/bash
# Phase 287: run phase287-relax.py in the phase-268 runtime (same environment as run274.sh) and copy the outputs back.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/bee163a0-6378-4e27-b375-db8eac6a2770/scratchpad
cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp "$S/phase287-relax.py" .phase287-relax.py
timeout 900 python3 .phase287-relax.py .phase287-relax.json > .phase287-relax.log 2>&1; echo "exit $?" >> .phase287-relax.log
cp .phase287-relax.json .phase287-relax.log "$S/"
cat .phase287-relax.log
