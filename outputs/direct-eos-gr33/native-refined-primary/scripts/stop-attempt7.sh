cd /home/lpaiu/work/native-refined267-runtime
W=primary267-refined64-work; L=.phase267-refined64-logs
mkdir -p $W/attempt7-seg-64/captures
for f in seg-64-fallbacks.json integer-64.json polish.json seg-64-driver.json seg-64-nonlinear-exceptions.json rejected-joint-stage.npz rejected-joint-stage.json; do [ -f $W/$f ] && mv $W/$f $W/attempt7-seg-64/; done
for f in $W/sweep-1/photons/seg-64*; do [ -f "$f" ] && mv "$f" $W/attempt7-seg-64/; done
for i in $(seq 244 299); do f=$(printf '%s/captures/captured-64-%03d.npz' $W $i); [ -f $f ] && mv $f $W/attempt7-seg-64/captures/; done
for f in seg-64.stdout.log seg-64.stderr.log exists-seg-64.json stopped.txt; do [ -f $L/$f ] && mv $L/$f $L/attempt7-$f; done
echo "$(date '+%F %T') attempt7 of seg-64 (v5, approved linear+nonlinear exceptions 1e-10) accepted 244 (original), 246 and 248 (exceptions: nonlinear defects 4.7e-12, 4.1e-11) and stopped at stage 250: long-double best vector residual 8.7e-10 > approved 1e-10 (moments 2e-19, material 2e-17), production polish 3.5e-7; preserved in $W/attempt7-seg-64" >> $L/status.txt
cp .phase267-driver-run.py .phase267-driver-run-v5.py
ls $W/attempt7-seg-64 $W/attempt7-seg-64/captures | head -40; ls $W/captures | tail -2; tail -3 $L/status.txt; ls $L
