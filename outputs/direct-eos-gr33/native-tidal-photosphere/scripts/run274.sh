cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase274-tides.py .phase274-tides.py
python3 .phase274-tides.py .phase274-tides.json > .phase274-tides.log 2>&1; echo "exit $?" >> .phase274-tides.log
