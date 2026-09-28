# Phase268: rerun the smoke existence check with the documented self-created alias, then start the 4x primary run.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
NEW=/home/lpaiu/work/native-refined268-runtime; LOG=$NEW/.phase268-prepare.log
step() { echo "$(date '+%F %T') $*" >> $LOG; }
cd $NEW || exit 1
tr -d '\r' < $S/phase268-checks.py > .phase268-checks.py
step "exists check rerun (the first run failed only on the self-created phase-151 alias)"
python3 .phase268-checks.py exists >> $LOG 2>&1 || { step "FAILED: exists-check rerun"; echo failed > .phase268-prepare.done; exit 1; }
tail -1 .phase268-smoke.stdout.log | cut -c1-400 >> $LOG
step "prepared"; echo ok > .phase268-prepare.done
tr -d '\r' < $S/run-quad64.sh > .run-quad64.sh && bash -n .run-quad64.sh || { step "FAILED: run script syntax"; exit 1; }
setsid nohup bash .run-quad64.sh > .run-quad64.out 2>&1 < /dev/null &
step "4x primary run launched (pid $!)"
