"""Counterexample candidate: reaction feasibility on saved material histories.

The radiative-recombination fit is an imported dilute-atom model, not a
rigorous bound on all physical neutralization channels. Ionization sources
are nonnegative and may be omitted only for this model's lower ion bound.
"""
import ctypes
import json
import signal
import subprocess
import sys
import time
import urllib.request
import numpy as np
from scipy.integrate import cumulative_trapezoid
import def_native_ion_closure as ions
import def_native_metric_release as release

OUT=ions.OUT
write=ions.write
URL='https://www.pa.uky.edu/~verner/dima/rec/rrfit.f'


def prepare():
    assert not (OUT/'rates-plan.json').exists()
    write(OUT/'rates-plan.json',dict(classification='Counterexample candidate',
        claim='Test whether the saved LTE outflow can supply its required hydrogen neutralization with imported spontaneous radiative recombination. Use material inventory labels, not fixed Eulerian cells.',
        model='Verner-Ferland1996 dilute-atom total radiative recombination. Electron density bounded by complete ionization of the actual24-element mixture. Omitted nonnegative ionization only increases the hydrogen ion fraction. Stimulated/three-body recombination, charge exchange and molecular channels are NOT included in this conditional bound.',
        decision='If required LTE recombination greatly exceeds the retained event budget, the existing LTE charge cannot serve as the final physical charge. The next actual evolution must transport populations and heat together.',
        material_labels='Six fixed baryon labels corresponding to10,25,50,75,90,99percent of final gas outside the original surface. At each time use cumulative resolved baryon minus cumulative inner inflow. Discarded tail uncertainty is recorded separately.',
        budget=dict(download_seconds=30,build_seconds=30,analysis_seconds=60,EOS_calls=30,new_fluid_steps=0,CPU_threads=1,memory_GB=1),
        stop='No rate multiplier, resolution increase or omitted-channel completion claim. Stop on material-label escape, failed inventory replay or negative rate.',
        source_url=URL,source_sha256=ions.sha(__file__),
        inputs={str(release.OUT/f'cells-{n}.npz'):ions.sha(release.OUT/f'cells-{n}.npz') for n in [896,1792]}))
    with urllib.request.urlopen(URL,timeout=30) as response:source=response.read()
    (OUT/'rrfit.f').write_bytes(source)
    # Keep the published single-precision routine unchanged. The C bridge
    # exposes its result as double; its fit uncertainty is not machine error.
    bridge='''subroutine hydrogen_rr(t,r) bind(C)
 use iso_c_binding
 implicit none
 real(c_double),value::t
 real(c_double),intent(out)::r
 real::a,b
 a=real(t)
 call rrfit(1,1,a,b)
 r=real(b,c_double)
end subroutine
'''
    (OUT/'rr-bridge.f90').write_text(bridge)
    cmd=['gfortran','-O2','-fPIC','-shared',str(OUT/'rrfit.f'),str(OUT/'rr-bridge.f90'),'-o',str(ions.CACHE/'rr.so')]
    begin=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
    (OUT/'rr-build.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    write(OUT/'rates-provider.json',dict(classification='Imported from prior work',url=URL,
        sha256=ions.sha(OUT/'rrfit.f'),library_sha256=ions.sha(ions.CACHE/'rr.so'),command=cmd,build_seconds=time.monotonic()-begin,
        attribution='D.A.Verner and G.J.Ferland1996,ApJS103,467; authors rrfit version4,1999-06-29. Spontaneous radiative recombination to all levels; do not call this the finite native partition inverse.'))


def provider():
    lib=ctypes.CDLL(str(ions.CACHE/'rr.so'));fn=lib.hydrogen_rr
    fn.argtypes=[ctypes.c_double,ctypes.POINTER(ctypes.c_double)];fn.restype=None
    def rr(T):
        values=[]
        for t in np.asarray(T).flat:
            result=ctypes.c_double();fn(float(t),ctypes.byref(result));values.append(result.value)
        return np.array(values).reshape(np.shape(T))
    return rr


def run():
    assert not (OUT/'history.json').exists();begin=time.monotonic();signal.alarm(60)
    rr=provider();ion=ions.Ions(cap=30);fan=ion.fan
    eps=(fan.X/ion.gas.c.A)@ion.gas.mapping if hasattr(ion.gas,'c') else (fan.X/ions.inventory.g.c.A)@ion.gas.mapping
    cx=float(eps@ion.gas.weights);eps/=cx
    nmax_per_rho=float(eps@ion.Z)*cx*6.02214076e23
    allrows=[]
    for n in [1792,896]:
        d=dict(np.load(release.OUT/f'cells-{n}.npz'));flow=release.Flow(n);m=flow.base
        t=np.r_[0.,d['history'][:,0]];U=np.concatenate([d['initial'][None],d['snapshots']])
        ports=np.r_[0.,d['history'][:,3]];mass=U[:,0]*d['volume'];edges=np.c_[np.zeros(len(t)),np.cumsum(mass,axis=1)]-ports[:,None]
        # Labels increase outwards. Positive left inflow is subtracted,
        # because it was not part of the initial material inventory.
        centers=(edges[:,:-1]+edges[:,1:])/2
        inside=d['x_cm']<0;surface=edges[-1,int(inside.sum())];total=edges[-1,-1]
        fractions=np.array([.10,.25,.50,.75,.90,.99]);labels=surface+fractions*(total-surface)
        rho=[];temp=[];where=[]
        for i,state in enumerate(U):
            r,v,lt=flow.primitive(state);active=r>=flow.eos.floor
            assert all((labels>=centers[i,active][0])&(labels<=centers[i,active][-1])),'Material label outside resolved support'
            rho.append(np.exp(np.interp(labels,centers[i,active],np.log(r[active]*fan.rho))))
            temp.append(np.exp(np.interp(labels,centers[i,active],lt[active])))
            where.append(np.interp(labels,centers[i,active],d['x_cm'][active]))
        rho=np.array(rho);T=np.array(temp);where=np.array(where)
        rate=rr(T)*rho*nmax_per_rho
        # a/W <=1 for this saved weak-field lapse. Coordinate time therefore
        # overestimates each parcel's proper time in the declared static model.
        assert max(flow.a)<=1+1e-14
        integrated=cumulative_trapezoid(rate,t,axis=0,initial=0)
        rows=[]
        for j,fraction in enumerate(fractions):
            first=ion.snapshot(np.log(rho[0,j]),np.log(T[0,j]),np.zeros(318))
            last=ion.snapshot(np.log(rho[-1,j]),np.log(T[-1,j]),np.zeros(318))
            x0=float(first['eos'][14]);x1=float(last['eos'][14]);bound=x0*np.exp(-integrated[-1,j])
            required=max(0.,np.log(x0/x1));available=float(integrated[-1,j])
            rows.append(dict(outside_baryon_quantile=float(fraction),initial_x_cm=float(where[0,j]),final_x_cm=float(where[-1,j]),
                final_density_ratio=float(rho[-1,j]/fan.rho),final_T=float(T[-1,j]),initial_H_ion=x0,LTE_final_H_ion=x1,
                retained_RR_lower_H_ion=float(bound),required_log_neutralization=float(required),available_log_neutralization=available,
                required_over_available=float(required/max(available,1e-300)),retained_RR_incompatible=bool(x1<bound)))
        np.savez_compressed(OUT/f'material-history-{n}.npz',t=t,rho=rho,T=T,x_cm=where,labels=labels,
            electron_ceiling=rho*nmax_per_rho,radiative_recombination_rate=rate,integrated_rate=integrated)
        allrows.append(dict(cells=n,rows=rows,discarded_baryon_label_width=float(d['conserved_discard'][0]),
            smallest_tracer_separation=float(min(np.diff(labels))),
            discarded_over_tracer_spacing=float(d['conserved_discard'][0]/min(np.diff(labels)))))
    result=dict(classification='Counterexample candidate',EOS_calls=ion.calls,seconds=time.monotonic()-begin,
        paths=allrows,required_rate_model='Spontaneous radiative recombination only; other neutralization channels need separate control.',
        instantaneous_LTE_supported_by_retained_RR=False,full_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'history.json',result);ion.save('history-native-states.npz');signal.alarm(0)
    print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
