# SHA-256 of the scripts installed in the 4x runtime vs the scratchpad sources (CR stripped as when installed).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
cd /home/lpaiu/work/native-refined268-runtime || exit 1
for f in phase267-template phase267-undriven phase267-p150-setup phase267-p150 phase267-born phase267-placeholders phase267-driver phase267-readout phase267-depth phase267-audited phase268-bank phase268-thermal phase268-initial phase268-initial-audit phase268-checks; do
  a=$(sha256sum .$f.py | cut -c1-16); b=$(tr -d '\r' < $S/$f.py | sha256sum | cut -c1-16)
  echo "$f $a $b $([ $a = $b ] && echo same || echo DIFF)"
done
