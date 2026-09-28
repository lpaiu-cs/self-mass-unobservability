# Copy a scratchpad launcher into WSL (CR stripped) and start it detached. Usage: bash launch268.sh <script name>
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
cd /home/lpaiu/work || exit 1
tr -d '\r' < $S/$1 > .$1 && bash -n .$1 || { echo "syntax error in $1"; exit 1; }
setsid nohup bash .$1 > .$1.out 2>&1 < /dev/null &
sleep 3; echo "launched $1 pid $!"; ps -p $! -o pid,etime,cmd | cut -c1-120
