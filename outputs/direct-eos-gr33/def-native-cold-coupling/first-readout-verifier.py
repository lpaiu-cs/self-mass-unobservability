"""Independent physical readout and native checks of the retained continuation."""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_cold_coupling as task

OUT=task.OUT;prior=task.prior;write=task.write


def readout():
    assert not (OUT/'readout-correction.json').exists()
    write(OUT/'readout-correction-plan.json',dict(classification='Counterexample candidate',
        discovery='Source inspection during the continuation found two reporting-only references inherited from run(): p0 used restarted theta/eta instead of zeros, and the initial atmosphere reference evaluated initial T at restarted U[0]. Neither reference enters evolution or conservation numerators.',
        correction='Subtract the known constant reference offset from history indices104 onward. Recompute response normalizations and all unchanged overlap/space gates. Preserve original output and producer. No trajectory replay or gate change.',
        acceptance='Independently recompute the final physical trace directly from final U,theta,eta and original initial fields; compare the104-108 overlap and original full-interval2percent spatial gates.',
        source_sha256=task.sha(__file__)))
    assert (OUT/'result.json').exists()
    model=prior.Coupled(896,8);f=model.flow;m=model.m;b=model.bulk;f.eos=task.ColdEOS()
    d=dict(np.load(OUT/'resumed-896-128.npz'));r=np.load(OUT/'restart104.npz');f.eos.y=f.eos.y0
    wrongp,wrongu,*_=f.eos(r['U'][0],f.initial_temperature);p0,u0,*_=f.eos(f.initial[0],f.initial_temperature)
    delta_at=float(((f.initial[0]*(wrongu-u0)-3*(wrongp-p0))*m.vol).sum()*model.gas_scale)
    pb=b.eos.gas(np.zeros(b.n),np.zeros(b.n))[0];p104=b.eos.gas(r['theta'],r['eta'])[0]
    delta_b=float(-3*((p104-pb)*b.volume*b.d['a']).sum())
    oldresponse=max(abs(d['bulk_trace']).max(),abs(d['atmosphere_trace']).max(),abs(d['boundary_energy'][0]),1.)
    oldgas=max(abs(d['atmosphere_trace']).max(),abs(d['ledger'][3]*model.gas_scale),1.)
    d['bulk_trace'][104:]+=delta_b;d['atmosphere_trace'][104:]+=delta_at
    rho,v,lt,y=f.primitive(d['U']);p,u,*_=f.eos(rho,lt)
    trace=-f.eos.cx*d['U'][0]*v*v/(1+np.sqrt(1-v*v))+rho*u-3*p-(f.initial[0]*u0-3*p0)
    direct_at=float(trace@m.vol*model.gas_scale)
    pp,uu,*_=b.eos.gas(d['theta'],d['eta']);direct_b=float((b.d['rho']*(uu-b.u0)-3*(pp-pb))@(b.volume*b.d['a']))
    direct=max(abs(d['bulk_trace'][-1]/direct_b-1),abs(d['atmosphere_trace'][-1]/direct_at-1));assert direct<1e-10
    old=json.loads((OUT/'result.json').read_text());row=old['continuation'].copy()
    response=max(abs(d['bulk_trace']).max(),abs(d['atmosphere_trace']).max(),abs(d['boundary_energy'][0]),1.)
    gas=max(abs(d['atmosphere_trace']).max(),abs(d['ledger'][3]*model.gas_scale),1.)
    row['total_energy_relative']*=oldresponse/response
    for k in ['gas_energy_relative','atmosphere_photon_energy_relative']:row[k]*=oldgas/gas
    row.update(bulk_trace_erg=float(d['bulk_trace'][-1]),atmosphere_trace_erg=float(d['atmosphere_trace'][-1]))
    row['passed']=bool(row['failure'] is None and row['completed_steps']==128 and max(row[k] for k in ['total_energy_relative','gas_energy_relative','atmosphere_photon_energy_relative'])<1e-8 and row['baryon_relative']<1e-10 and row['species_relative']<1e-9)
    saved=np.load(prior.OUT/'cells-896-steps-128.npz');coarse=np.load(prior.OUT/'cells-448-steps-128.npz');overlap={};errors={}
    for k in ['bulk_trace','atmosphere_trace','outside_mass']:overlap[k]=float(np.max(abs(d[k][104:109]-saved[k][104:109]))/max(float(np.max(abs(saved[k]))),1.))
    for k in ['atmosphere_trace','outside_mass']:errors[k]=float(np.max(abs(d[k]-coarse[k]))/max(float(np.max(abs(d[k]))),1.))
    result=dict(old,passed=bool(row['passed'] and max(overlap.values())<1e-8 and max(errors.values())<.02),continuation=row,overlap=overlap,space_comparison=errors)
    for name in ['result.json','resumed-896-128.json','resumed-896-128.npz']:(OUT/name).rename(OUT/('first-readout-'+name))
    np.savez_compressed(OUT/'resumed-896-128.npz',**d);write(OUT/'resumed-896-128.json',row);write(OUT/'result.json',result)
    receipt=dict(classification='Counterexample candidate',bulk_reference_offset_erg=delta_b,atmosphere_reference_offset_erg=delta_at,independent_endpoint_relative=direct,
        state_arrays_unchanged=True,original_verdict_preserved=True,passed=result['passed'],overlap=overlap,space_comparison=errors)
    write(OUT/'readout-correction.json',receipt);print(json.dumps(receipt),flush=True)


def audit():
    assert not (OUT/'audit.json').exists();start=time.monotonic()
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',native_calls=80,seconds=30,fluid_replays=0,
        checks='Independent native final coldest/fastest cells and the original170K failed conservative state. Compare actual Doppler-shifted positive absorption and emission with logarithmic native coefficients. Symbolically verify scaled partition products and first/second derivatives.',
        constitutive_gate=.002,spectrum_gate=.002,source_sha256=task.sha(__file__)))
    signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(30)
    model=prior.Coupled(896,8);f=model.flow;f.eos=task.ColdEOS();tab=task.ColdSpectrum();native=task.logarithmic_native(80)
    d=np.load(OUT/'resumed-896-128.npz');rho,v,lt,y=f.primitive(d['U']);p,u,g,T,kap=f.eos(rho,lt);field=d['I'].sum(0);active=np.flatnonzero(rho>=f.eos.floor)
    ids=sorted(set([int(active[np.argmin(lt[active])]),int(active[np.argmax(abs(v[active]))])]))
    rows=[]
    for j in ids:
        z=native.state(float(np.log(rho[j])),float(lt[j]),float(y[j]));D=(1-v[j]*model.mu)/np.sqrt(1-v[j]**2);E=model.E[None,:]*D[:,None]/model.m.a[j]
        ab=np.zeros_like(E);em=np.zeros_like(E)
        for k in range(10):
            sigma=tab.cross(E,k+1);ab+=np.exp(np.log(y[j])+z['log_fraction'][k])*sigma;good=sigma>0
            em[good]+=np.exp(np.log(y[j])+z['log_fraction'][k]+z['affinity']-E[good]/(prior.optical.ex.K*T[j]))*sigma[good]
        aa,ee=tab.coefficients(np.array([rho[j]*f.eos.rho0]),lt[j:j+1],y[j:j+1],E[None]);errors=[]
        for true,estimate,I in [(ab,aa[0],field[j]),(em,ee[0],1+field[j])]:
            for factor in [model.number,model.number*model.E]:
                weight=model.w[:,None]*D[:,None]*factor[None,:]*I;errors.append(float(np.sum(abs(estimate-true)*weight)/max(float(np.sum(true*weight)),1e-250)))
        constitutive=float(max(abs(p[j]*f.eos.rho0*prior.C**2/z['raw'][1]-1),abs(u[j]*prior.C**2/z['raw'][2]-1)))
        rows.append(dict(cell=j,T_K=float(T[j]),rho=float(rho[j]*f.eos.rho0),y=float(y[j]),native_population=z['population_error'],constitutive=constitutive,spectrum=max(errors)))
    a,n,ell,s=sp.symbols('a n ell s',positive=True);term=2*n*n*(sp.exp(a/(n*n))-1-a/(n*n))*sp.exp(ell)
    first=2*sp.exp(a/(n*n)+ell-s)*(1-sp.exp(-a/(n*n)));second=2*sp.exp(a/(n*n)+ell-s)/(n*n)
    assert sp.simplify(sp.diff(term,a)*sp.exp(-s)-first)==0
    assert sp.simplify(sp.diff(term,a,2)*sp.exp(-s)-second)==0
    Q,L=sp.symbols('Q L',positive=True);assert sp.simplify(sp.exp(L+s)*(sp.exp(-s)*Q)-sp.exp(L)*Q)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='At each evaluation point a common exponential scaling of the original partition and its jets preserves physical products and derivative ratios in exact real arithmetic. The Planck-Larkin first/second a derivatives used in scaled_summand agree symbolically.',
        exclusion='This is not a uniform floating-point, tail-truncation, table-derivative, complete microphysics or final-charge error theorem.'))
    result=dict(classification='Counterexample candidate',passed=bool(field.min()>=0 and max(max(r['constitutive'],r['spectrum']) for r in rows)<.002),checks=rows,
        minimum_photon_occupation=float(field.min()),native_calls=native.ion.calls+1,constructor_native_calls=3,seconds=time.monotonic()-start,source_sha256=task.sha(__file__),final_charge_solved=False)
    write(OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
