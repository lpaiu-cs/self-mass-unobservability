"""Direct mechanical contrast on the identical full discrete operator.

The clamped control replaces only momentum equations by xi=0. Subtracting
those equations gives a momentum-only RHS, avoiding subtraction of unit fields.
"""
from pathlib import Path
import inspect
import json
import time
import numpy as np
from scipy.sparse.linalg import spsolve
import def_free_surface_response_normalized as repaired

old=repaired.old
h=old.h
OUT=repaired.OUT/'direct-contrast'

# Reuse the frozen assembly exactly; do not rerun or edit its saved solutions.
source=inspect.getsource(old.solve)
anchor='    scale=np.asarray(abs(matrix).sum(1)).ravel();'
assert source.count(anchor)==1
source=source[:source.index(anchor)]+'    return matrix,rhs,grid,Phi,hw,R\n'
ns=dict(repaired.namespace)
exec(compile(source,__file__,'exec'),ns)


def solve(k,j):
    M,b,grid,Phi,hw,R=ns['solve'](k,j);n=len(grid)-1
    momentum=np.r_[4*np.arange(n)+1,4*n+2]
    fixed=M.tolil(copy=True)
    for row,node in zip(momentum,range(n+1)):
        fixed.rows[row]=[4*node];fixed.data[row]=[1.+0j]
    fixed=fixed.tocsc()
    def linear(A,rhs):
        scale=np.asarray(abs(A).sum(1)).ravel()
        y=spsolve(A.multiply((1/scale)[:,None]).tocsc(),rhs/scale)
        residual=float(np.max(abs(A@y-rhs)/(np.asarray(abs(A)@abs(y)).ravel()+abs(rhs)+1e-100)))
        assert residual<1e-9,residual
        return y,residual
    clamped,rc=linear(fixed,b)
    force=np.zeros_like(b)
    # Exactly identical scalar/baryon equations cancel symbolically. Never form
    # b-M@clamped on those rows: their roundoff would swamp the tiny contrast.
    force[momentum]=b[momentum]-M[momentum,:]@clamped
    delta,rd=linear(M,force)
    total=clamped+delta
    original=np.load(repaired.OUT/f'harmonic-{k}-grid-{j}.npz')['response'].ravel()
    reconstruction=float(np.max(abs(total-original))/np.max(abs(original)))
    assert reconstruction<1e-10,reconstruction
    assert np.max(abs(clamped[::4]))==0
    assert np.count_nonzero(force[np.setdiff1d(np.arange(len(b)),momentum)])==0
    y=delta.reshape(n+1,4);outgoing=(y[-1,2]-grid[-1]*y[-1,0]*Phi)/hw
    Mstar=json.loads((h.OUT/'absolute-shoot/result-0.001.json').read_text())['ADM_geom_m']
    charge=-R*grid[-1]/Mstar*outgoing
    np.savez_compressed(OUT/f'contrast-{k}-{j}.npz',grid=grid,clamped=clamped.reshape(n+1,4),contrast=y,force=force)
    return dict(harmonic=k,refinement=j,clamped_residual=rc,contrast_residual=rd,
        reconstruction_relative=reconstruction,charge=[float(charge.real),float(charge.imag)],
        outgoing=[float(outgoing.real),float(outgoing.imag)],
        surface_cancellation_ratio=float((abs(y[-1,2])+abs(grid[-1]*y[-1,0]*Phi))/max(abs(outgoing*hw),1e-100)))


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    files=[Path(__file__),Path(repaired.__file__),Path(old.__file__),repaired.OUT/'background.npz',repaired.OUT/'charge-result.json']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        claim='Resolve the material-motion outgoing charge without subtracting near-unit scalar fields, with exactly the same discrete baryon/scalar rows in free and supported comparisons.',
        identity='M delta=(M_clamped-M) y_clamped; RHS is supported only on replaced momentum rows. M(y_clamped+delta)=b in exact arithmetic.',
        gates=dict(residual=1e-9,reconstruction_relative=1e-10,charge_grid_relative=.02),
        budget=dict(hard_timeout_seconds=60,native_calls=0,time_steps=0,solves=16,grids=[1,2],harmonics=[0,1,2,3]),
        stop='No additional grid or frequency if gates fail; preserve failed direct subtraction.'))
    rows=[solve(k,j) for k in range(4) for j in [1,2]];comparisons=[]
    benchmark=json.loads((old.exterior.s.OUT/'companion-benchmark.json').read_text())
    scale=benchmark['leading_drive']['maximum_scalar_excursion_from_inverse_semimajor_reference']/.001
    static=complex(*rows[1]['charge'])
    for k in range(4):
        a,b=[complex(*r['charge']) for r in rows if r['harmonic']==k]
        error=abs(a-b)/max(abs(b),1e-100)
        comparisons.append(dict(harmonic=k,charge=[b.real,b.imag],relative_grid_difference=error,passed=error<.02,
            maximum_drive_scaled_delta_alpha_over_phi0=abs(b)*scale,
            static_subtracted_maximum_drive_scaled=abs(b-static)*scale))
    result=dict(classification='Counterexample candidate',rows=rows,comparisons=comparisons,
        passed=all(r['passed'] for r in comparisons),seconds=time.monotonic()-start,
        discrete_control='External holding force, not a second freely evolving star; scalar and baryon rows exactly shared.',
        grid_scope='Response-grid comparison only on a frozen interpolated background, not a joint continuum/EOS/thermal certificate.',
        charge_scope='Outgoing scalar monopole coefficient divided by fixed background ADM mass; does not include a dynamic ADM-mass normalization term.',
        actual_orbital_harmonics_applied=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
