cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.phase277-bounds.py .phase277-bounds.py
python3 .phase277-bounds.py .phase277-bounds.json > .phase277-bounds.log 2>&1; echo "exit $?" >> .phase277-bounds.log
