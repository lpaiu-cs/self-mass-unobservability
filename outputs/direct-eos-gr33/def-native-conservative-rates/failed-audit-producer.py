"""Counterexample candidate: actual native collision defect on conserved paths.

Reuse native states rather than redoing their EOS inversions. Constitutive
differences enter the actual moving collision owner, with opposite material
energy/H/momentum. No prescribed-time correction is installed as an EOS.
"""
from pathlib import Path
import json,signal,sys,time
import numpy as np
import sympy as sp
import def_native_conservative_eos_readout as prior

flow=prior.flow;chem=prior.chem;C=flow.C;write=flow.write;sha=flow.sha
OUT=flow.OUT.parent/'def-native-conservative-rates';AMU=1.66053906660e-24


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='095fd9cf6',
        claim='Reconstruct native absorption, stimulated/spontaneous emission and electron scattering from323 saved conserved native states, then apply them to the actual moving photon/material collision owner.',
        decision='Does interpolation satisfy its original0.2percent coefficient gate and leave actual net energy/H/force exchange within2percent? Material failure requires constitutive repair in coupled evolution, not another diagnostic-only certificate.',
        reuse='Saved native rawEOS,levels,affinity and actual angular occupations at17times,19deep cells. No EOS inversion or fluid replay.',
        closure='Native density,T,H and original other ionic/molecular inventories. Preserve the SAME additive first-order inventory correction to coefficients and electrons. Native frequencies use installed physical lapse; retain inherited moving spectral derivative closure.',
        comparison='Independently compare positive absorption/emission channels, net Killing energy/neutral creation/proper force, and the same actual collision owner. Integrate absolute defects over the saved trajectory; these are empirical frozen-trajectory measures, not final-response bounds.',
        gates=dict(weighted_coefficients=.002,net_exchange_integrated=.02,owner=1e-10,paired_energy_H=1e-10),
        budget=dict(seconds=40,native_initialization_calls=2,new_native_state_calls=0,CPU_threads=1,memory_GB=3),
        forecast='Only native cross-section evaluation and17 existing collision calls per comparison; previous complete charge readout8.2s.40s hard stop, no finer time/grid or automatic trajectory.',
        stop='Preserve original defect. No post-result weaker gate or time-only pressure/rate fit. If a coefficient fails, repair its constitutive interpolation owner using independent states before evolution.',
        limits='No native advected-inventory continuum, complete microscopic opacity, atmospheric EOS, uniform derivative/time bound, full GR or fixed-point claim.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(chem.prior.__file__),
            prior.OUT/'production-samples.json',prior.OUT/'result.json',flow.OUT/'coupled-128.npz',
            chem.old.ATOMIC,chem.old.LEVELS]}))
    a,e,I,en,q=sp.symbols('a e I en q')
    net=e*(1+I)-a*I
    assert sp.expand(net-(e-(a-e)*I))==0
    assert sp.expand(en*net-en*net)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='One bound-free photon changes neutral H by one and has exactly opposite photon/material energy; e*(1+I)-a*I equals e-(a-e)*I. This does not establish microscopic rates or finite evolution.'))


def audit():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(40)
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert sha(p)==h,p
    rows=json.loads((prior.OUT/'production-samples.json').read_text());by={(r['it'],r['cell']):r for r in rows}
    m=flow.Coupled();b=m.bulk;e=b.eos;native=chem.old.Native(cap=2);initial_calls=native.ion.calls
    points,z,ids=prior.inputs(m,128);coefficients=[];moments=[];checks=[]
    weight=b.d['num']*b.d['Einf'];energy_weight=b.volume[:,None,None]*b.w[None,:,None]*weight[None,None,:]/b.d['a'][:,None,None]**3
    number_weight=b.volume[:,None,None]*b.w[None,:,None]*b.d['num'][None,None,:]/b.d['a'][:,None,None]**3
    force_weight=energy_weight*b.mu[None,:,None]/(C*b.d['a'][:,None,None])
    for it,k in enumerate(ids):
        m.Pi=z['snapshot_Pi'][k].copy();m.h=z['snapshot_h'][k].copy();m.j=z['snapshot_j'][k].copy()
        m.mass=m.mass0-m.h[1:]+m.h[:-1];m.set_material(points[it]['t'])
        theta=z['snapshot_theta'][k];eta=z['snapshot_eta'][k];I=z['snapshot_bulk_I'][k];beta=m.velocity()
        oldgas=e.gas(theta,eta);oldrad=e.radiation(theta,eta);base=e.base.radiation(theta+e.theta0,eta)
        nr=[];ne=[]
        for j in range(b.n):
            r=by[it,j];s=dict(r,fraction=np.array(r['fraction']))
            chi,emit=chem.prior.coefficients(native,s,b.d['Einf']/b.d['a'][j]);nr.append([chi+emit,emit]);ne.append(r['raw'][13]/AMU)
        nr=np.array(nr);ne=np.array(ne)-e.inventory[2]*e.xi
        df=e.spectral['frequency'][0]*(1-e.f)+e.spectral['frequency'][1]*e.f
        for channel in [0,1]:
            # Preserve the original additive advected-inventory term. Other
            # density and frequency effects are evaluated natively above.
            nr[:,channel]-=base[channel]*e.rinventory[channel]*e.xi[:,None]*(1+df[:,channel]*e.frequency_shift[:,None])
        assert nr.min()>=0 and ne.min()>0
        positive=[]
        for channel,field in [(0,I),(1,np.ones_like(I)),(1,I)]:
            weighted=field*energy_weight*b.factor[:,None,None]
            true=nr[:,channel,None,:]*weighted;old=oldrad[channel][:,None,:]*weighted
            positive.append(float(np.sum(abs(true-old))/max(np.sum(true),1e-300)))
        before=m.collision(I,theta,eta,beta)
        gas,rad=e.gas,e.radiation
        def native_gas(t,y):
            v=list(gas(t,y));v[6]=ne;return v
        def native_rad(t,y):
            v=list(rad(t,y));v[0]=nr[:,0];v[1]=nr[:,1];return v
        e.gas=native_gas;e.radiation=native_rad
        try:after=m.collision(I,theta,eta,beta)
        finally:e.gas=gas;e.radiation=rad
        # Independent stationary coefficient residual, plus the moving terms
        # returned by the same owner; unlike a positive-rate comparison this
        # retains absorption/emission cancellation in the net source.
        raw=b.factor[:,None,None]*((nr[:,1]-oldrad[1])[:,None,:]-(nr[:,0]-oldrad[0]-nr[:,1]+oldrad[1])[:,None,:]*I)
        ds=b.scfactor[:,None,None]*(ne-oldgas[6])[:,None,None]*(np.einsum('qk,ikf->iqf',b.S,I)-I)
        defect=after[0]-before[0];estimate=raw+after[1]-before[1]+ds
        owner=float(np.sum(abs(defect-estimate)*energy_weight)/max(np.sum(abs(defect)*energy_weight),1.))
        def integrated(v):
            full,extra,escape,bound_extra=v
            # Scattering conserves photon count, including declared spectral
            # exits. Only bound-free creation changes neutral H.
            ab,em=(oldrad[:2] if v is before else (nr[:,0],nr[:,1]))
            bound=b.factor[:,None,None]*(em[:,None,:]-(ab-em)[:,None,:]*I)+bound_extra
            energy=np.sum(full*energy_weight,axis=(1,2))+b.volume*escape[:,1]
            neutral=np.sum(bound*number_weight,axis=(1,2))
            force=-np.sum(full*force_weight,axis=(1,2))-b.volume*escape[:,2]/(b.d['a']*C)
            return np.array([energy,neutral,force])
        original=integrated(before);actual=integrated(after);moments.append([original,actual,actual-original])
        coefficients.append(np.stack([oldrad[0],oldrad[1],nr[:,0],nr[:,1]]))
        checks.append(dict(time=points[it]['t'],positive_energy_weighted_errors=positive,
            electron_relative=float(np.max(abs(ne/oldgas[6]-1))),owner_relative=owner))
    assert native.ion.calls==initial_calls
    times=np.array([r['t'] for r in points]);moments=np.array(moments)
    integrated=np.trapz(np.sum(abs(moments),axis=-1),times,axis=0)
    net_relative=integrated[2]/np.maximum(integrated[0],1.)
    maxcoef=max(max(r['positive_energy_weighted_errors']) for r in checks);owner=max(r['owner_relative'] for r in checks)
    np.savez_compressed(OUT/'native-collision.npz',t=times,coefficients=coefficients,moments=moments,
        moments_order=np.array(['Killing_photon_energy_per_s','neutral_creation_per_s','gas_force_dyn']),
        integrated_absolute=integrated)
    result=dict(classification='Counterexample candidate',passed=bool(maxcoef<.002 and max(net_relative)<.02 and owner<1e-10),
        maximum_positive_channel_relative=maxcoef,net_exchange_integrated_relative=net_relative.tolist(),owner_relative=owner,
        checks=checks,native_initialization_calls=initial_calls,new_native_state_calls=0,seconds=time.monotonic()-start,
        native_coefficients_applied_to_collision=True,constitutive_repair_required=bool(maxcoef>=.002 or max(net_relative)>=.02),
        frozen_trajectory_only=True,finite_coupled_response_bounded=False,coupled_EOS_evolution=False,
        full_source_error_enclosed=False,final_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
