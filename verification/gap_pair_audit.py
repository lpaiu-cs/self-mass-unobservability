"""Registered two-gap integer lattice and a selected nuisance-displacement probe."""
import itertools
import json
import math
import numpy as np

import comparator_audit as comp
import nuisance_audit as audit


def main():
    g=comp.prepare(); t=g['inp']['t']; sw=g['inp']['sw']; n=len(t)
    z=np.load(audit.ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    p=dict(zip(z['names'],z['params'])); cycle=86400e6/p['spinfreq']
    ordered=np.sort(t)
    ix=np.argsort(np.diff(ordered))[-10:][::-1]
    cuts=(ordered[ix]+ordered[ix+1])/2
    raw=cycle*sw[:,None]*(t[:,None]>cuts)
    projected=audit.sc.proj_out(g['qf'],raw)
    candidates={tuple([0]*10)}
    for i,j in itertools.combinations(range(10),2):
        for a,b in itertools.product([-1,0,1],repeat=2):
            v=np.zeros(10,dtype=int); v[i]=a; v[j]=b; candidates.add(tuple(v))
    vectors=np.array(sorted(candidates)); assert len(vectors)==201
    rows=[]; selected=None
    for amplitude in comp.AMPLITUDES:
        m=comp.metric(g,amplitude); h=g['h']
        y=g['yp']+h@((m['root']-1)*g['y'])
        steps=projected+h@((m['root']-1)[:,None]*(h.T@projected))
        harmonic=h@m['u']
        y=audit.sc.proj_out(harmonic,y); steps=audit.sc.proj_out(harmonic,steps)
        gram=steps.T@steps; cross=steps.T@y
        changes=np.einsum('ki,ij,kj->k',vectors,gram,vectors)+2*vectors@cross
        nonzero=np.any(vectors!=0,axis=1)
        j=int(np.argmin(np.where(nonzero,changes,np.inf)))
        norm=math.sqrt(float(vectors[j]@gram@vectors[j]))
        rows.append(dict(a=amplitude,coefficients=vectors[j].tolist(),minimum_delta_chi2=float(changes[j]),remaining_step_norm=norm,
                         candidates=[dict(coefficients=v.tolist(),delta_chi2=float(c)) for v,c in zip(vectors,changes)]))
        if amplitude==1.:
            selected=vectors[j]
            # Fit the change induced by adding this step, using exactly the registered full basis.
            cscale=np.linalg.norm(g['c6'],axis=0)
            design=np.column_stack([g['b'],g['c6']/cscale])
            whitened=design+h@((m['root']-1)[:,None]*(h.T@design))
            step=raw@selected
            ws=step+h@((m['root']-1)*(h.T@step))
            shift=np.linalg.lstsq(whitened,-ws,rcond=1e-12)[0]
            delta=shift[:28]/g['norms'][:28]
            assert np.linalg.norm(whitened@shift+ws)**2 <= float(ws@ws)
            prediction=g['inp']['J']@delta
            np.savez(audit.OUT/'runtime12/nonlinear-gap-input.npz',parameter_delta=delta,prediction_us=prediction,
                     selected_step_us=step/sw,coefficients=selected,gap_indices=ix)
    result=dict(status='Imported from prior work',scope='201 finite assignments on the ten longest gaps; arbitrary multiple errors and nonlinear reconnection not certified',
                gap_indices=ix.tolist(),gap_days=np.diff(ordered)[ix].tolist(),rows=rows,
                nonlinear_probe='Incremental 28-parameter displacement for the weakest a=1 assignment; not a complete nonlinear joint timing-and-noise refit')
    (audit.OUT/'gap-pair-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in rows: print(r['a'],r['minimum_delta_chi2'],r['coefficients'])


if __name__=='__main__': main()
