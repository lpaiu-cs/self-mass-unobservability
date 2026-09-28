# Phase271: clean per-mode records on the three grids (x arrival) and the rejected delay-arrival comparison on 4x.
R=/home/lpaiu/work/.run-freefall271.sh
for spec in "native-retained-tail-runtime readout267-identity2-work original" "native-refined267-runtime readout267-refined64-work refined64" "native-refined268-runtime readout268-quad64-work quad64"; do
  set -- $spec
  bash $R /home/lpaiu/work/$1 compare $2 $3 > /dev/null 2>&1 || echo "compare $3 failed"
  bash $R /home/lpaiu/work/$1 charge $2 $3 2>&1 | tail -1
done
bash $R /home/lpaiu/work/native-refined268-runtime compare readout268-quad64-work quad64-delay delay > /dev/null 2>&1 || echo "delay compare failed"
ls /home/lpaiu/work/native-*-runtime/readout26*-work/freefall-*.json
