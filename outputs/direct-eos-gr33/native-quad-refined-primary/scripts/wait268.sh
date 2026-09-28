# Wait until a marker file appears in the 4x runtime, then print the prepare log tail. Usage: bash wait268.sh <marker> [<log>]
cd /home/lpaiu/work/native-refined268-runtime || exit 1
until test -f "$1"; do sleep 20; done
cat "$1"; tail -6 "${2:-.phase268-prepare.log}" | cut -c1-600
