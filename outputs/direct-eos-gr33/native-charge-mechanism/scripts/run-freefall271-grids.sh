# Phase271: free-fall compare + charge on the 1x (original) and 2x grids (same field solver per grid, T readouts).
R=/home/lpaiu/work/.run-freefall271.sh
echo "=== 2x"; bash $R /home/lpaiu/work/native-refined267-runtime compare readout267-refined64-work refined64 2>&1 | grep -E 'cell (1[2-9]|2[0-6]) ' | head -12
bash $R /home/lpaiu/work/native-refined267-runtime charge readout267-refined64-work refined64 2>&1 | tail -1
echo "=== 1x"; ls /home/lpaiu/work/native-retained-tail-runtime/outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz
bash $R /home/lpaiu/work/native-retained-tail-runtime compare readout267-identity2-work original 2>&1 | grep -E 'cell ([89]|1[0-8]) ' | head -12
bash $R /home/lpaiu/work/native-retained-tail-runtime charge readout267-identity2-work original 2>&1 | tail -1
