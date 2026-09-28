"""Independent native endpoint and saved-angular audit for the unsplit solve."""
import json
import signal
import time
import numpy as np
import def_native_unsplit_photons as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,task.timeout);signal.alarm(10)
    result=json.loads((out/'result.json').read_text());assert result['passed']
    m=task.Model(16);d=m.d;v=np.load(out/'angles-16-steps-128.npz');theta=v['theta'][-1];eta=v['eta'][-1]
    native=task.prior.old.old.Native(cap=24);p,u,*_=m.eos.gas(theta,eta);aa,ee,*_=m.eos.radiation(theta,eta);checks=[]
    for j in sorted(set([int(np.argmax(abs(theta))),m.n-1])):
        task.prior.old.setup(native,d,j);state=native.state(0.,np.log(d['T'][j])+theta[j],native.y0*(1+eta[j]))
        chi,em=task.prior.prior.coefficients(native,state,d['Einf']/d['a'][j]);ab=chi+em;mask=(ab>0)&(em>0)
        constitutive=float(max(abs(p[j]/state['raw'][1]-1),abs(u[j]/state['raw'][2]-1)))
        rate=float(max(np.max(abs(aa[j,mask]/ab[mask]-1)),np.max(abs(ee[j,mask]/em[mask]-1))))
        assert constitutive<.002 and rate<.002
        checks.append(dict(cell=j,constitutive=constitutive,rate=rate,population_error=state['population_error']))
    angular=[]
    for q,n in [(8,64),(8,128),(16,128)]:
        data=np.load(out/f'angles-{q}-steps-{n}.npz');I=data['snapshots'];w=data['angular_weights'];mu=data['mu'];mu2=mu*mu+(2/q)**2/12
        J=np.einsum('tiqf,q->tif',I,w);H=np.einsum('tiqf,q->tif',I,w*mu);K=np.einsum('tiqf,q->tif',I,w*mu2)
        index=np.rint(data['snapshot_times']/task.prior.prior.END*n).astype(int)
        mismatch=float(max(np.max(abs(a-data[key][index])/(J+1e-200)) for a,key in [(J,'J'),(H,'H'),(K,'K')]))
        variance=np.divide(J*K-H*H,J*J,out=np.zeros_like(J),where=J>1e-140)
        assert I.min()>=0 and mismatch<1e-12 and variance.min()>=-1e-12 and np.all(K<=J*(1+1e-14))
        angular.append(dict(angles=q,steps=n,minimum_occupation=float(I.min()),minimum_variance=float(variance.min()),saved_moment_relative=mismatch))
    port=v['outgoing_occupation'];Hport=np.einsum('tqf,q->tf',port,v['angular_weights']*v['mu'])
    power=4*np.pi*task.prior.C*m.area[-1]*Hport*d['num']*d['Einf']
    np.savez_compressed(out/'surface-port-candidate.npz',times=v['t'],occupation=port,mu=v['mu'],angular_weights=v['angular_weights'],Einf=d['Einf'],number_measure=d['num'],
        surface_radius_cm=d['edges'][-1],surface_lapse=d['face_a'][-1],Killing_luminosity_per_frequency=power,source_time_angle_accepted=True,spatial_frequency_accepted=False)
    audit=dict(classification='Counterexample candidate',passed=True,native_checks=checks,angular_checks=angular,
        source_time_angle_comparison_passed=True,spatial_frequency_accepted=False,moving_atmosphere_connected=False,final_charge_solved=False,
        native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=task.sha(__file__))
    task.write(out/'audit.json',audit);signal.alarm(0);print(json.dumps(audit),flush=True)


if __name__=='__main__':main()
