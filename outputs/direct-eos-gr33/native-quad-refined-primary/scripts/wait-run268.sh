# Wait until the 4x run stops, completes, or segment <label> passes; print the status log. Usage: bash wait-run268.sh <label>
cd /home/lpaiu/work/native-refined268-runtime || exit 1
L=.phase268-quad64-logs; W=primary268-quad64-work
until [ -f $L/stopped.txt ] || [ -f $L/complete.txt ] || { [ -f $W/$1-driver.json ] && grep -q '"passed": true' $W/$1-driver.json; }; do sleep 30; done
cat $L/status.txt; cat $L/stopped.txt 2>/dev/null; date '+%F %T'
