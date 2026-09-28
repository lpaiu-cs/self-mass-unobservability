cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase279-gray-sph.py .phase279-gray-sph.py
rm -f .phase279-gray-sph.done
timeout 3000 python3 .phase279-gray-sph.py readout268-quad64-work/kernel-quad64.json .phase279-gray-sph.json > .phase279-gray-sph.log 2>&1; echo "exit $?" >> .phase279-gray-sph.log; touch .phase279-gray-sph.done
