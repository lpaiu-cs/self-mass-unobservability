cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase278-transit.py .phase278-transit.py
rm -f .phase278-transit.done
/usr/bin/time -v python3 .phase278-transit.py .phase278-transit.json > .phase278-transit.log 2>&1; echo "exit $?" >> .phase278-transit.log
touch .phase278-transit.done
