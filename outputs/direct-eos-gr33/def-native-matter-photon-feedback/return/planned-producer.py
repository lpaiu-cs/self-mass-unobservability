"""Apply returned photon transfers to matter, then to represented GR.

Counterexample candidate: preserve separate waveform-sweep and GR residuals.
Do not infer closed exterior/deep physics from a finite feedback correction.
"""
from pathlib import Path
from types import FunctionType
import json,signal,sys,time
import numpy as np
import def_native_matter_photon_feedback as photons
import def_native_material_branch_response as old
import verify_native_material_response as stress
import def_native_characteristic_gr as wave
import verify_native_anisotropic_gr as independent

OUT=photons.OUT/'return';GR=OUT/'gr';write=old.write;sha=old.sha;AMP=old.AMP;C=old.base.C;LD=np.longdouble


class Material(old.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference);self.steps=steps
        p=np.load(photons.OUT/f'steps-{steps}-reference-{reference}.npz');ids=[int(np.argmin(abs(p['t']-t))) for t in self.t];assert np.max(abs(p['t'][ids]-self.t))<1e-18
        collision=p['collision_transfer'][ids];impulse=p['moments'][ids,3]
        self.transfer=np.stack([np.zeros_like(collision[:,:,0]),impulse/self.a,collision[:,:,0],collision[:,:,1]],axis=1)/AMP
    run=FunctionType(old.base.Material.run.__code__,dict(vars(old.base),OUT=OUT),argdefs=old.base.Material.run.__defaults__)


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d6ccddff9',
        claim='Apply newly returned photon collision energy/H/momentum to the actual shared material response, then apply its rest-mass,pressure,trace and updated photons to the characteristic GR/scalar solver.',
        decision='Measure the finite feedback correction to the represented scalar charge and the unresolved material/photon waveform residual. Preserve incomplete closure if those residuals are not controlled.',
        reuse='Same original531 cells,3.434ms,EOS banks,64/128 clocks and saved GR fields. Reuse exact accepted material SSP and primitive readout owners; only paired photon transfer histories change.',
        budgets=dict(pilot_seconds=40,material_production_seconds=500,source_seconds=45,GR_seconds=90,CPU_threads=1,memory_GB=3,new_native_bank_calls=0),
        forecast='Use the prior complete material paths to capture late CFL work, plus new prefix seconds/raw-call. Require2x forecast plus setup within500s; do not dispatch from an early-step extrapolation alone.',
        gates=dict(conservation=1e-8,directional=.002,time=.02,background=.02,pressure=.002,GR_quadrature=.002,GR_independent=1e-9),
        stop='Stop if the photon paths fail, a new material path fails, or a cap is exceeded. No automatic extra mesh, clocks, horizon or unbounded waveform iteration.',
        limits='This applies the return sweep to represented compact GR sources. New GR fields still need transport feedback; external generated/scattered scalar and deep response remain explicit. No final-charge flag.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(photons.__file__),Path(old.__file__),Path(old.base.__file__),Path(stress.__file__),Path(wave.__file__),old.OUT/'production.json']}))


def pilot():
    assert json.loads((photons.OUT/'result.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    start=time.monotonic();signal.signal(signal.SIGALRM,old.base.flow.old.optical.timeout);signal.alarm(40);rows=[]
    for steps,ref in [[64,128],[128,128],[128,64]]:
        m=Material(ref,steps);r=m.run(steps,f'pilot-{steps}-{ref}',2);r.update(physical_branch_ratio=m.physical_branch_ratio);rows.append(r)
        if not r['passed']:break
    reference=json.loads((old.OUT/'production.json').read_text())['paths'];calls=sum(r['raw_owner_calls'] for r in reference)
    rate=sum(r['seconds'] for r in rows)/sum(r['raw_owner_calls'] for r in rows);forecast=calls*rate+9;upper=2*forecast+10
    result=dict(classification='Counterexample candidate',rows=rows,forecast_seconds=forecast,upper_seconds=upper,eligible=len(rows)==3 and all(r['passed'] for r in rows) and upper<500,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);print(json.dumps(result));signal.alarm(0)
    if result['eligible']:
        write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,production_seconds=500,forecast_seconds=forecast,upper_seconds=upper,
            bindings={str(p):sha(p) for p in [Path(__file__),Path(photons.__file__),Path(old.__file__),Path(old.base.__file__),OUT/'pilot.json',photons.OUT/'result.json']}))


def production():
    p=json.loads((OUT/'execution-plan.json').read_text());assert p['eligible'];assert not (OUT/'production.json').exists()
    for path,h in p['bindings'].items():assert sha(path)==h,path
    start=time.monotonic();signal.signal(signal.SIGALRM,old.base.flow.old.optical.timeout);signal.alarm(p['production_seconds']);rows=[]
    try:
        for steps,ref in [[64,128],[128,128],[128,64]]:
            m=Material(ref,steps);r=m.run(steps,f'steps-{steps}-reference-{ref}',restart=f'pilot-{steps}-{ref}');r.update(physical_branch_ratio=m.physical_branch_ratio);rows.append(r)
            if not r['passed']:break
        result=dict(classification='Counterexample candidate',passed=len(rows)==3 and all(r['passed'] for r in rows),paths=rows,seconds=time.monotonic()-start)
        write(OUT/'production.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(OUT/'production-failure.json',dict(error=repr(exc),completed_paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


def sources():
    assert json.loads((OUT/'production.json').read_text())['passed'];assert not (OUT/'sources.json').exists();start=time.monotonic()
    signal.signal(signal.SIGALRM,old.base.flow.old.optical.timeout);signal.alarm(45);gr=wave.Response();histories=[];allstress=[];rows=[]
    b=gr.model.bulk;model=gr.model
    weights=4*np.pi*np.r_[b.W,model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
    for steps,ref in [[64,128],[128,128],[128,64]]:
        m=Material(ref,steps);d=np.load(OUT/f'steps-{steps}-reference-{ref}.npz');ids=[int(np.argmin(abs(d['t']-t))) for t in m.t];z=d['history_scaled'][ids];histories.append(z)
        c=gr.coeff(m.rE);V=m.V;factor=C**4/wave.base.G
        Eg0,Pg0,Kg0,Er0,Pr0=[c[key]*factor*V for key in ['Eg','Pg','Kg','Er','Pr']]
        phi=m.metric['delta_u'];lam=m.metric['delta_lambda'];s=3*phi+lam
        p=np.load(photons.OUT/f'steps-{steps}-reference-{ref}.npz');pid=[int(np.argmin(abs(p['t']-t))) for t in m.t];mom=p['moments'][pid]
        background=np.load(old.base.flow.OUT/f'coupled-{ref}.npz');bid=m.ids
        I=np.concatenate([background['snapshot_bulk_I'][bid],background['snapshot_I'][bid].sum(1)],axis=1)
        Ebg=np.einsum('tnqf,nqf->tn',I,weights)/m.a;Pbg=np.einsum('tnqf,nqf,q->tn',I,weights,b.mu2)/m.a
        mu4=(b.edges_mu[1:]**5-b.edges_mu[:-1]**5)/(5*np.diff(b.edges_mu));I0=np.concatenate([b.initial,model.initial_I]);E0=np.sum(I0*weights,axis=(1,2));ratio4=np.sum(I0*weights*mu4[None,:,None],axis=(1,2))/E0;R40=ratio4*Er0
        photonE=mom[:,0]/m.a-Ebg*s+4*Er0*phi+(Er0+Pr0)*lam
        photonP=mom[:,5]/m.a-Pbg*s+4*Pr0*phi+(3*Pr0-R40)*lam
        source=[];pressure_errors=[];balance=float(np.max(abs(np.sum(d['history_scaled'],axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        for k,t in enumerate(m.t):
            point=m.point(k);field=m.fields(t)[2];p0,delta,err=stress.pressure(m,k,z[k],field)
            q=point['Q'];total=(z[k,2]+m.rest*z[k,0])/m.a*AMP;backgroundE=(q[2]+m.rest*q[0])/m.a
            eF=total+(Eg0+Pg0-backgroundE)*s[k];pF=delta*AMP+(Kg0[None]-p0)*s[k]
            source.append(np.array([eF,pF[1],eF-pF[1]-2*pF[0],pF[0]]));pressure_errors.append(err*AMP)
        source=np.array(source);allstress.append(source);pressure_error=float(np.max(np.sum(abs(pressure_errors),axis=2))/max(np.max(np.sum(abs(source[:,[3,1]]),axis=2)),1.))
        prior=np.load(old.OUT/f'steps-{steps}-reference-{ref}.npz');ii=[np.argmin(abs(prior['t']-t)) for t in m.t];oldz=prior['history_scaled'][ii]
        change=np.max(np.sum(abs(z-oldz),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(z),axis=2),axis=0),1.)
        # This is a measured waveform mismatch, not a contraction certificate.
        joint=mom[:,[1,2]];actual=z[:,[2,3]]*AMP;residual=np.max(np.sum(abs(joint-actual),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),1.)
        base=dict(np.load(wave.base.OUT/f'source-{ref}.npz'));rest=z[:,0].astype(LD)*LD(AMP)*LD(m.model.cx)*LD(C)**2
        base.update(baryon_g=z[:,0]*AMP,gas_nonrest_energy_erg=source[:,0].astype(LD)-rest,nonrest_trace_erg=source[:,2].astype(LD)-rest,
            photon_energy_erg=photonE,photon_radial_pressure_erg=photonP,metric_stress_erg=source[:,0]+photonE-source[:,1]-photonP,
            inner_cumulative_energy_erg=p['radial_ports'][pid,0,1],outer_cumulative_energy_erg=p['radial_ports'][pid,1,1])
        label=f'{steps}-reference-{ref}';np.savez_compressed(GR/f'source-{label}.npz',**base)
        np.savez_compressed(OUT/f'stress-{label}.npz',t=m.t,radius_E=m.rE,material=source,photon_energy=photonE,photon_radial_pressure=photonP)
        rows.append(dict(steps=steps,reference=ref,conservation=balance,pressure_probe=pressure_error,material_sweep_change=change.tolist(),energy_H_waveform_residual=residual.tolist()))
    def compare(a,b):return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    comparisons=dict(time=compare(histories[0],histories[1]),background=compare(histories[2],histories[1]),stress_time=compare(allstress[0],allstress[1]),stress_background=compare(allstress[2],allstress[1]))
    result=dict(classification='Counterexample candidate',passed=max(v for row in comparisons.values() for v in row)<.02 and max(max(r['conservation']/1e-8,r['pressure_probe']/.002) for r in rows)<1,
        comparisons=comparisons,paths=rows,seconds=time.monotonic()-start,returned_photon_transfer_applied_to_material=True,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(OUT/'sources.json',result);print(json.dumps(result));signal.alarm(0)


class GRResponse(wave.Response):
    run=FunctionType(wave.base.Response.run.__code__,dict(vars(wave.base),OUT=GR))


def fields():
    assert json.loads((OUT/'sources.json').read_text())['passed'];assert not (GR/'result.json').exists();start=time.monotonic()
    write(GR/'plan.json',dict(classification='Counterexample candidate',claim='Apply returned full material/photon forcing to the existing characteristic linear GR and scalar operator; preserve its canonical closure and actual inner photon debit.',
        budget_seconds=90,paths=[['128-reference-128',8],['128-reference-128',4],['64-reference-128',8],['128-reference-64',8]],
        gates=dict(quadrature=.002,time=.02,background=.02,independent=1e-9),
        limits='Additional compact component on the declared momentarily balanced initial operator. New lapse/exterior scalar and transport reapplication still remain; no final physical charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(wave.__file__),Path(wave.base.__file__),OUT/'sources.json']}))
    signal.signal(signal.SIGALRM,old.base.flow.old.optical.timeout);signal.alarm(90);m=GRResponse();rows=[]
    for label,order in [('128-reference-128',8),('128-reference-128',4),('64-reference-128',8),('128-reference-64',8)]:rows.append(m.run(label,order))
    fine=np.load(GR/'fields-128-reference-128-g8.npz');norm=max(np.max(abs(fine['U'])),1e-300);comparison={}
    for key,name in [('quadrature','128-reference-128-g4'),('time','64-reference-128-g8'),('background','128-reference-64-g8')]:comparison[key]=float(np.max(abs(np.load(GR/f'fields-{name}.npz')['U']-fine['U']))/norm)
    data=dict(np.load(GR/'source-128-reference-128.npz'));direct,error=independent.direct(m,data,8);agreement=abs(direct/rows[0]['endpoint_direct']-1)
    result=dict(classification='Counterexample candidate',passed=comparison['quadrature']<.002 and max(comparison['time'],comparison['background'])<.02 and agreement<1e-9,
        comparisons=comparison,independent_direct_relative=agreement,paths=rows,seconds=time.monotonic()-start,returned_source_applied_to_represented_GR=True,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(GR/'result.json',result);print(json.dumps(result));signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
