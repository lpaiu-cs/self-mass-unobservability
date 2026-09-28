# State of the charge-dominating layers (original grid cells 8-15) and available independent EOS codes.
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
python3 - <<'EOF'
import numpy as np
g = np.load('outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz')
print('keys', g.files)
R = g['edges'][-1]
for c in range(8, 16):
    print(c, 'depth km %.0f-%.0f' % ((R - g['edges'][c+1])/1e5, (R - g['edges'][c])/1e5), 'rho %.3e' % g['rho'][c], 'T %.3e' % g['T'][c])
EOF
ls /home/lpaiu/work | grep -i -E 'freeeos|helm|opal|eos' | head
find /home/lpaiu -maxdepth 4 -iname '*helm*' -o -maxdepth 4 -iname '*freeeos*' 2>/dev/null | head -8
