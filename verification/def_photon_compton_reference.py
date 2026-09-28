"""Independent finite-time spectral reference; do not extend a failed run."""
from pathlib import Path
import argparse
import json
import time
import signal
import numpy as np
from scipy.linalg import eigh,qr
import def_photon_spatial_coupling as m

OUT=m.OUT


def prepare():
    assert not (OUT/'reference-plan.json').exists()
    paths=[Path(__file__),Path(m.__file__),OUT/'bank.npz',OUT/'controls.json',OUT/'compton-control.npz',OUT/'result.json']
    m.ex.write(OUT/'reference-plan.json',dict(classification='Counterexample candidate',
        reason='The 8/rate_C control is time-converged but differs from the infinite-time invariant projection by 0.006355. The old plan had no decay-gap evidence for equilibrium by that time. Also the dense expanded null residual was normalized by the physical coupling rate rather than the stiff matrix scale. Preserve both failed tests.',
        claim='Independently evaluate the same finite-time Kompaneets problem by symmetric eigendecomposition on logarithmic frequency cells. Do not extend the physical duration or relax the old gates.',
        method='Remove the known energy and photon-number null modes by orthogonal projection, diagonalize the symmetric positive dissipative block, and evaluate its exponential at y=8. Frequency counts 256,512,1024. Compare temperature and energy-norm spectrum with the previous unsplit native-grid result.',
        gates=dict(reference_last_difference_initial=.001,native_to_reference_initial=.001,temperature_difference_K=.001,stiff_null_backward_error=1e-12),
        original_asymptotic_gate_retained=False,original_failed_result_preserved=True,
        budget=dict(hard_seconds=90,frequency_counts=[256,512,1024],new_stellar_steps=0,new_network_queries=0,new_EOS_calls=0,automatic_expansion=False),
        bindings={p.relative_to(m.h.ROOT).as_posix():m.h.digest(p) for p in paths}))


def reference(b,n):
    original=m.Operator(dict(b,rate_a=np.zeros_like(b['rate_a']),rate_s=np.zeros_like(b['rate_s'])),1)
    hi=float(b['edges_u'][len(original.u)]);lo=float(b['edges_u'][1])
    edges=np.r_[0,np.geomspace(lo,hi,n)];u=(edges[1:]+edges[:-1])/2
    Ci=4*float(b['arad'])*float(b['T'])**3*m.previous.loss.weights(edges,1)[:,1]
    copied=dict(b,u=u,edges_u=edges,Ci=Ci,rate_a=np.zeros(n),rate_s=np.zeros(n),rate_C=np.array(1.))
    op=m.Operator(copied,1);L=op.L.toarray()
    matrix=np.block([[np.array([[op.aa]]),op.q[None,:]],[op.q[:,None],L]])
    invariants=np.column_stack([np.r_[np.sqrt(op.Cm),op.energy],np.r_[0,op.number]])
    Q,_=qr(invariants,mode='full');null=Q[:,:2];rest=Q[:,2:]
    reduced=rest.T@matrix@rest;lam,vec=eigh((reduced+reduced.T)/2)
    assert lam.min()>0,lam.min()
    initial=np.zeros(op.size+1);initial[0]=np.sqrt(op.Cm)
    coefficients=vec.T@(rest.T@initial)
    equilibrium=null@(null.T@initial)
    final=equilibrium+rest@(vec@(np.exp(-8*lam)*coefficients))
    # Reconstruct the natural photon chemical-potential perturbation u*E/C.
    eta=op.u*final[1:]/np.sqrt(op.Ci)
    return dict(u=op.u,Ci=op.Ci,T=final[0]/np.sqrt(op.Cm),eta=eta,
        equilibrium_distance=float(np.linalg.norm(final-equilibrium)/np.sqrt(op.Cm)),
        smallest_positive_rate=float(lam[0]),slowest_mode_initial=float(abs(coefficients[0])/np.sqrt(op.Cm)))


def run():
    assert not (OUT/'reference-result.json').exists();start=time.monotonic();signal.alarm(90)
    plan=json.loads((OUT/'reference-plan.json').read_text())
    for p,sha in plan['bindings'].items():assert m.h.digest(m.h.ROOT/p)==sha,p
    b=dict(np.load(OUT/'bank.npz'));old=np.load(OUT/'compton-control.npz');Cm=float(old['Cm'])
    oldT=float(old['T'].real/np.sqrt(Cm));oldq=old['E'].real/np.sqrt(old['Ci'])
    rows=[];answers=[]
    for n in plan['budget']['frequency_counts']:
        r=reference(b,n);predicted=np.interp(old['u'],r['u'],r['eta'])/old['u']
        norm=float(np.sqrt(Cm*(r['T']-oldT)**2+np.sum(old['Ci']*(predicted-oldq)**2))/np.sqrt(Cm))
        row=dict(frequency_cells=n,temperature_K=float(r['T']),native_temperature_difference_K=float(abs(r['T']-oldT)),
            native_energy_norm_difference_initial=norm,equilibrium_distance=r['equilibrium_distance'],
            smallest_positive_rate_over_rate_C=r['smallest_positive_rate'],slowest_mode_initial=r['slowest_mode_initial'])
        rows.append(row);answers.append((r,predicted));print('REFERENCE',row,flush=True)
    a,x=answers[-2];z,y=answers[-1]
    difference=float(np.sqrt(Cm*(a['T']-z['T'])**2+np.sum(old['Ci']*(x-y)**2))/np.sqrt(Cm))
    op=m.Operator(dict(b,rate_a=np.zeros_like(b['rate_a']),rate_s=np.zeros_like(b['rate_s'])),1)
    # Standard backward error for a matrix-vector null test on a stiff grid.
    matrix_bound=max(op.aa+np.sum(abs(op.q)),float(np.max(np.asarray(abs(op.L).sum(1)).ravel()+abs(op.q))))
    backward=[]
    for vec in [np.r_[np.sqrt(op.Cm),op.energy],np.r_[0,op.number]]:
        a,z=op.rhs(vec[0],vec[1:]);backward.append(float(max(abs(a),np.max(abs(z)))/(matrix_bound*np.max(abs(vec)))))
    result=dict(classification='Counterexample candidate',passed=bool(difference<.001 and rows[-1]['native_energy_norm_difference_initial']<.001 and rows[-1]['native_temperature_difference_K']<.001 and max(backward)<1e-12),
        rows=rows,reference_last_energy_norm_difference_initial=difference,stiff_null_backward_errors=backward,
        original_asymptotic_control_passed=False,original_rate_normalized_null_control_passed=False,
        interpretation='Independent finite-time exponential determines whether the residual is real relaxation. Numerical backward error on the stiff matrix is a separate diagnostic and does not retroactively pass the original physical-rate-normalized null test.',
        original_results_preserved=True,physical_duration_extended=False,seconds=time.monotonic()-start)
    m.ex.write(OUT/'reference-result.json',result);print('REFERENCE_RESULT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
