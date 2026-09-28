"""Audit saved nonlinear trajectories against native states; no new evolution."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_causal_nonlinear as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.alarm(30)
    result=json.loads((out/'result.json').read_text());assert result['passed']
    model=task.Model();d=model.d;saved=np.load(out/'steps-128.npz');native=task.old.Native(cap=200)
    theta=saved['theta'][-1];eta=saved['eta'][-1];p,u,*_=model.eos.gas(theta,eta);ab,em,*_=model.eos.radiation(theta,eta)
    ids=np.unique(np.r_[0,3,6,9,12,15,np.argmax(abs(theta)),np.argmax(abs(eta))]);checks=[]
    for j in ids:
        task.setup(native,d,j);s=native.state(0.,np.log(d['T'][j])+theta[j],native.y0*(1+eta[j]))
        chi,e=task.prior.coefficients(native,s,d['Einf']/d['a'][j]);a=chi+e;mask=(a>0)&(e>0)
        checks.append(dict(cell=int(j),logT_change=float(theta[j]),neutral_ratio=float(1+eta[j]),
            constitutive=float(max(abs(p[j]/s['raw'][1]-1),abs(u[j]/s['raw'][2]-1))),
            rate=float(max(np.max(abs(ab[j,mask]/a[mask]-1)),np.max(abs(em[j,mask]/e[mask]-1))))))
    # A directional finite difference checks the assembled analytic Jacobian
    # against the full nonlinear RHS, including EOS inversion and scattering.
    state=saved['state'][64];f,t,p,u,jac=model.rhs(state);direction=np.zeros_like(state)
    direction[model.it]=.01*np.sin(np.arange(model.n)+.3);direction[model.iy]=.01*np.cos(np.arange(model.n))
    direction[:model.nj]=(model.J/model.scale*.001).ravel();direction[model.nj:model.nj+model.nh]=(model.H/model.scale*.001).ravel()
    eps=1e-4;fd=(model.rhs(state+eps*direction,False)[0]-model.rhs(state-eps*direction,False)[0])/(2*eps);analytic=jac@direction
    relative=float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic));assert relative<2e-5,relative
    # The face flux and midpoint mean are distinct staggered samples. This
    # reconstruction is a warning about unresolved radiation structure, not
    # proof that a continuous distribution violates the moment inequality.
    h=saved['state'][:,model.nj:model.nj+model.nh].reshape(-1,model.n+1,model.m)*model.scale+model.H
    hc=(h[:,:-1]+h[:,1:])/2;j=saved['J'];weight=d['num']*d['Einf']
    ratio=abs((hc*weight).sum(2))/(j*weight).sum(2)
    moment_note=dict(classification='Counterexample candidate',
        maximum_reconstructed_bolometric_flux_ratio=float(ratio.max()),P1_positive_moment_limit=float(1/np.sqrt(3)),
        initial_maximum_reconstructed_ratio=float(ratio[0].max()),
        independently_sampled_initial_face_ratio=float(max(abs(model.initial_flux_fraction))),
        physical_realizability_certified=False,
        interpretation='Averaging adjacent face H and pairing with midpoint J exceeds the P1 positive-moment bound in this coarse discretization. Its initial Wien-tail reconstruction already fails severely although every initial face moment is realizable. Spatial/angle reconstruction must be resolved before treating this trajectory as a physical photon distribution; positive J and exact energy alone are insufficient.')
    result=dict(classification='Counterexample candidate',
        passed=bool(max(c['constitutive'] for c in checks)<.002 and max(c['rate'] for c in checks)<.002 and relative<2e-5),
        numerical_native_audit_only=True,actual_trajectory_native_checks=checks,analytic_jacobian_directional_relative=relative,
        radiation_reconstruction=moment_note,native_calls=native.ion.calls,seconds=time.monotonic()-start,
        physical_photon_transport_certified=False,final_charge_solved=False,source_sha256=task.prior.sha(__file__))
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':main()
