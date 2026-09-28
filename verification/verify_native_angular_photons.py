"""Independent native/root regression and saved angular-field checks."""
from pathlib import Path
from types import MethodType
import json
import signal
import time
import numpy as np
import def_native_angular_photons as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.alarm(15)
    d=np.load(task.prior.OUT/'bank-16-8.npz');native=task.old.old.Native(cap=27)
    j=10;task.old.setup(native,d,j);target=native.base['molecular_H_fractions'].copy()
    low=native.state(0.,9.74137215461128,native.y0*1e-4);actual=native.ion.states[-1]['molecular_H_fractions']
    relative=float(np.max(abs(actual/target-1)));assert low['population_error']<1e-12 and relative<1e-6
    fixed=native.ion.constrain
    source=task.prior.old.old.constrain_source
    before='(now>1e-18)|(target_molecules>1e-18)';assert source.count(before)==1
    legacy=source.replace(before,'now>1e-18');ns=dict(vars(task.prior.old.old.cold.old));exec(compile(legacy,'legacy-constraint','exec'),ns)
    controls=[]
    for old_path in [False,True]:
        native.ion.constrain=MethodType(ns['constrain'],native.ion) if old_path else fixed
        task.old.setup(native,d,15);state=native.state(0.,np.log(d['T'][15])+.017,native.y0*.6)
        controls.append(state['raw'])
    same=float(np.max(abs(controls[0]-controls[1])/np.maximum(abs(controls[0]),1e-100)));assert same<1e-10
    native.ion.constrain=fixed
    m=task.Model(16);data=np.load(out/'second-order-angles-16-steps-128.npz')
    theta=data['theta'][-1];eta=data['eta'][-1];j=int(np.argmax(abs(theta)))
    task.old.setup(native,d,j);state=native.state(0.,np.log(d['T'][j])+theta[j],native.y0*(1+eta[j]))
    p,u,*_=m.eos.gas(theta,eta);aa,ee,*_=m.eos.radiation(theta,eta);chi,em=task.prior.coefficients(native,state,d['Einf']/d['a'][j]);ab=chi+em;mask=(ab>0)&(em>0)
    constitutive=float(max(abs(p[j]/state['raw'][1]-1),abs(u[j]/state['raw'][2]-1)))
    rate=float(max(np.max(abs(aa[j,mask]/ab[mask]-1)),np.max(abs(ee[j,mask]/em[mask]-1))));assert constitutive<.002 and rate<.002
    angular=[]
    for q,n in [(8,64),(8,128),(16,128)]:
        mm=task.Model(q);v=np.load(out/f'second-order-angles-{q}-steps-{n}.npz');x=v['snapshots'];weights=v['angular_weights'];mu=v['mu'];mu2=mu**2+(2/q)**2/12
        J=np.einsum('tiqf,q->tif',x,weights);H=np.einsum('tiqf,q->tif',x,weights*mu);K=np.einsum('tiqf,q->tif',x,weights*mu2)
        indexes=np.rint(v['snapshot_times']/task.prior.END*n).astype(int)
        error=max(np.max(abs(w-v[key][indexes])/(J+1e-200)) for w,key in [(J,'J'),(H,'H'),(K,'K')])
        variance=np.divide(J*K-H*H,J*J,out=np.zeros_like(J),where=J>1e-140)
        assert error<1e-12 and x.min()>=0 and variance.min()>=-1e-12 and np.all(K<=J*(1+1e-14))
        flat=mm.A@np.ones(mm.n*q);incoming=np.zeros((mm.n,q))
        incoming[0]=task.C*mm.area[0]*np.maximum(mu,0)/mm.W[0]
        incoming[-1]=-task.C*mm.area[-1]*np.minimum(mu,0)/mm.W[-1]
        transport=float(max(abs(flat+incoming.ravel()))/max(abs(mm.A.data)));assert transport<1e-11
        angular.append(dict(angles=q,steps=n,minimum_occupation=float(x.min()),minimum_variance=float(variance.min()),saved_moment_relative=float(error),constant_occupation_transport_relative=transport))
    port=data['outgoing_occupation'];Hport=np.einsum('tqf,q->tf',port,data['angular_weights']*data['mu'])
    power=4*np.pi*task.C*m.area[-1]*Hport*d['num']*d['Einf']
    np.savez_compressed(out/'surface-port-candidate.npz',times=data['t'],occupation=port,mu=data['mu'],angular_weights=data['angular_weights'],Einf=d['Einf'],number_measure=d['num'],
        surface_radius_cm=d['edges'][-1],surface_lapse=d['face_a'][-1],Killing_luminosity_per_frequency=power,source_full_trace_accepted=False)
    comparison=json.loads((out/'second-order-result.json').read_text())
    assert not comparison['passed'] and comparison['comparisons']['time_trace']>.02
    assert comparison['comparisons']['time_surface_spectrum']<.02 and comparison['comparisons']['angle_surface_spectrum']<.02
    result=dict(classification='Counterexample candidate',passed=True,scope='Native constraint regression and saved positive angular moments; not a full physical or time-convergence acceptance.',
        restored_H2_relative=relative,restored_total_population_error=low['population_error'],old_supported_native_relative=same,
        actual_hot_cell=j,actual_hot_native_constitutive=constitutive,actual_hot_native_rate=rate,angular_checks=angular,
        surface_port_time_angle_comparison_passed=True,full_trace_time_comparison_passed=False,moving_atmosphere_connected=False,final_charge_solved=False,
        native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=task.sha(__file__))
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
