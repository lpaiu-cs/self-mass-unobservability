cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase275-photosphere.py .phase275-photosphere.py
python3 .phase275-photosphere.py readout268-quad64-work/kernel-quad64.json .phase275-photosphere.json > .phase275-photosphere.log 2>&1; echo "exit $?" >> .phase275-photosphere.log
