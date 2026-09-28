# Wait until a segmented run stops, completes, or segment <label> passes. Usage: bash wait-run.sh <runtime> <logs> <work> <label>
cd "$1" || exit 1
L=$2; W=$3
until [ -f $L/stopped.txt ] || [ -f $L/complete.txt ] || { [ -f $W/$4-driver.json ] && grep -q '"passed": true' $W/$4-driver.json; }; do sleep 30; done
cat $L/status.txt; cat $L/stopped.txt 2>/dev/null; date '+%F %T'
