cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cp /home/lpaiu/work/.diag279b.py .diag279b.py
timeout 1500 python3 .diag279b.py > .diag279b.log 2>&1; echo "exit $?" >> .diag279b.log; touch .diag279b.done
