# Inputs of the free-fall model (read-only): background arrays, pulse parameters, state polynomial format.
cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import inspect, sys
import numpy as np
sys.path.insert(0, 'verification')
z = np.load('readout268-quad64-work/gr/source-64.npz'); g = np.load('outputs/direct-eos-gr33/def-native-boundary-layer/geometry.npz')
fs = np.load('readout268-quad64-work/gr/field-source-64.npz')
print('source keys', [k for k in z.files if not k.startswith(('state_coeff', 'geometry_coeff'))])
print('field-source extra keys', sorted(set(fs.files) - set(z.files)))
for k in ['drive_amplitude', 'drive_radius', 'drive_duration', 'M_cm', 'K_cm', 'RJ', 'inner_material_face_cm', 'cx']: print(k, z[k])
for c in [8, 16, 20, 24, 28, 42, 43, 100]:
    print('cell', c, 'r %.6e' % z['radius'][c], 'a %.9f' % z['a'][c], 'B %.9f' % z['B'][c], 're %.9f' % z['re'][c], 'x %.6e' % z['drive_x'][c], 'delay %.6e' % z['delay'][c])
print('geometry phi range', g['phi'].min(), g['phi'].max(), 'rho cells 8..', g['rho'][8:12], 'a', g['a'][8:10], 'B', g['B'][8:10])
import apply_driver_aware_radau_gr as ad
print('=== prior.poly'); print(inspect.getsource(ad.prior.poly)[:1500])
print('state_coeff shape', z['state_coeff_baryon_g'].shape, 't', z['t'].shape)
EOF
