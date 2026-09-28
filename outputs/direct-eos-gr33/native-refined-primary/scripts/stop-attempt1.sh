cd /home/lpaiu/work/native-refined267-runtime
pkill -f 'phase267-driver-run.py primary267-refined64-work seg-64' ; sleep 3
while pgrep -f 'run-refined64.sh' >/dev/null || pgrep -f 'phase267-driver-run.py' >/dev/null; do sleep 1; done
W=primary267-refined64-work; L=.phase267-refined64-logs
mkdir -p $W/attempt1-seg-64/captures
for f in seg-64-fallbacks.json integer-64.json polish.json seg-64-driver.json; do [ -f $W/$f ] && mv $W/$f $W/attempt1-seg-64/; done
for f in $W/sweep-1/photons/seg-64*; do [ -f "$f" ] && mv "$f" $W/attempt1-seg-64/; done
for i in $(seq 224 299); do f=$(printf '%s/captures/captured-64-%03d.npz' $W $i); [ -f $f ] && mv $f $W/attempt1-seg-64/captures/; done
for f in seg-64.stdout.log seg-64.stderr.log exists-seg-64.json stopped.txt; do [ -f $L/$f ] && mv $L/$f $L/attempt1-$f; done
echo "$(date '+%F %T') attempt1 of seg-64 stopped by operator at stage 246 (production solver stalled near 5e-11); preserved in $W/attempt1-seg-64" >> $L/status.txt
cp .phase267-driver-run.py .phase267-driver-run-v2.py
tr -d '\r' < /mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad/phase267-driver.py > .phase267-driver-run.py
ls $W/attempt1-seg-64 $W/attempt1-seg-64/captures | head -40; ls $W/captures | tail -2; tail -3 $L/status.txt; ls $L
