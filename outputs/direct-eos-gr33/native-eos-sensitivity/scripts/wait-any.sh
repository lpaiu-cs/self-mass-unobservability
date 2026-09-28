# Wait until a marker file appears in a runtime, then print it and a log tail. Usage: bash wait-any.sh <runtime dir> <marker> <log>
cd "$1" || exit 1
until test -f "$2"; do sleep 20; done
cat "$2"; tail -12 "$3" | cut -c1-700
