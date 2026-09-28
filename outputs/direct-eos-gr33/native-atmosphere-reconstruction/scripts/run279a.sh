cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase279-gray.py .phase279-gray.py
rm -f .phase279-gray.done
timeout 2400 python3 .phase279-gray.py readout268-quad64-work/kernel-quad64.json .phase279-gray.json > .phase279-gray.log 2>&1; echo "exit $?" >> .phase279-gray.log; touch .phase279-gray.done
