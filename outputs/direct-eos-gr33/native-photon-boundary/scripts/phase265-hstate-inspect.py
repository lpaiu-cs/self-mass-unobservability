"""Read-only: compare the low H (neutral) state entering coarse step 3 in phase 259 and 265b."""
import numpy as np
z = np.load('native-short-return259-work/sweep-1/photons/interval-1-64.npz')
print('259 keys sample', [k for k in z.files if 'stage' in k or 'unit' in k][:20])
old = z['joint_stage_conserved_scaled']  # (stages, 4, n) order B,S,E,N
print('stage states', old.shape)
unit = z['material_neutral_units']
r = np.load('native-photon-boundary-newton265b-work/rejected-joint-stage.npz')
g = r['initial'][-531*4:].reshape(531, 4)
new_N = g[:, 1] * unit
for idx in [3]:
    o = old[idx, 3]
    print(f'259 stage {idx} N: sum|.|={np.sum(abs(o)):.6e}; 265b initial N: sum|.|={np.sum(abs(new_N)):.6e}; rel diff L1 {np.sum(abs(new_N - o))/np.sum(abs(o)):.3e}')
sol = r['solution'].reshape(2, -1)[:, -531*4:].reshape(2, 531, 4)
d = r['defect'].reshape(2, -1)[:, -531*4:].reshape(2, 531, 4)
for s in range(2):
    num = np.sum(abs(d[s, :, 1]) * unit); den = np.sum(abs(sol[s, :, 1]) * unit)
    top = np.argsort(-(abs(d[s, :, 1]) * unit))[:6]
    share = np.sum(np.sort(abs(d[s, :, 1]) * unit)[::-1][:6]) / num
    print(f'stage {s}: H relative {num/den:.5e}; top cells {top.tolist()} carry {share:.2%} of the defect; |sol N| there {[f"{abs(sol[s,i,1]*unit[i]):.2e}" for i in top]} vs median {np.median(abs(sol[s,:,1]*unit)):.2e}')
