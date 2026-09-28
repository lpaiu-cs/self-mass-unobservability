cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.atmos274.py .phase274-atmos.py
python3 .phase274-atmos.py .phase274-atmos.json > .phase274-atmos.log 2>&1; echo "exit $?" >> .phase274-atmos.log
