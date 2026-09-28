"""Counterexample candidate: GR sources from one saved joint solution.

Read the completed prefix without rerunning or modifying its producer. Recover
native radial stress in the stable energy coordinate; use its actual geometry,
four material variables, photons and boundary ledger together.
"""
from pathlib import Path
from types import FunctionType
import gc,inspect,json,os,resource,sys,time
import numpy as np
import sympy as sp
import complete_full_incident_horizon as run
import def_native_matter_photon_feedback as feedback
import verify_native_material_response as stress
import def_native_characteristic_gr as gr
import verify_native_stage_energy_charge as charge
import def_retained_native_return as retained

OUT=Path('native-joint-gr186-work');INPUT=run.OUT
read,write,sha=run.read,run.write,run.sha
LD,AMP,C=run.LD,run.AMP,feedback.C
CAPS=dict(prepare=25,source=90,fields=60)
INTERVAL=2


def saved(n):return INPUT/f'sweep-1/photons/interval-{INTERVAL:02d}-{n}.npz'


def pressure_function():
    # Reuse the stable primitive inverse and the existing radial kinetic stress.
    s=feedback.primitive_source.rsplit('    return dict(',1)[0]
    assert s.count('xi=np.zeros(nb)')==1
    s=s.replace('xi=np.zeros(nb)','xi=model.mech.xi@dh')
    tail=inspect.getsource(stress.pressure).split('    # Verify pressure',1)[1]
    s+='    # Verify pressure'+tail
    scope=dict(feedback.namespace);exec(compile(s,__file__,'exec'),scope)
    (OUT/'expanded-pressure.py').write_text(s)
    return scope['primitive']


def prepare():
    assert not OUT.exists();OUT.mkdir()
    pair=read(INPUT/f'comparison-{INTERVAL:02d}.json');assert pair['passed']
    for folder in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    inputs=list((INPUT/'sweep-0').rglob('*.npz'))
    inputs += [INPUT/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for src in inputs:
        dst=OUT/src.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True)
        os.link(src,dst);reused[str(dst)]=sha(src)
    files=inputs+[saved(n) for n in [64,128]]+[INPUT/f'comparison-{INTERVAL:02d}.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='97a3f86c9',
        claim='Connect the SAME accepted full-input joint solution to its actual GR stress and retarded field; preserve all four material coordinates, current geometry and its own energy/boundary history.',
        decision='A source or field failure blocks this adapter. Passing this prefix admits reuse of this mapping on the eventual accepted full horizon, not a final-charge verdict or an automatic new evolution.',
        scope='Saved interval02 only,0..T/8. Three canonical times; the extraT/32time stays in the original input. Retained directional input; no self-generated GR is returned to material yet.',
        precision='Read normalized material history directly to recover Etilde without subtracting rest-sized energies. Restore inventory displacement in the existing stable primitive inverse and retain radial kinetic stress. Remove only canonical geometry already in the GR operator.',
        reuse='No new physical steps, EOS bank, incident field, cadence or spatial resolution. Frozen185producer inputs and running code are untouched; no previous charge or different response is added.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,identity=1e-12,independent_GR=1e-9,conservation=1e-8),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        forecast='Previously measured source/charge work was9.3s and three full17-time GR fields12.43s. Current joint constructor and native pressure reuse add unmeasured overhead. This three-time prefix is capped at90s source plus60s fields; no automatic expansion.',
        stop='Any original gate, missing canonical state, identity or budget failure. Keep failures; no cadence/order/gate change or new producer dispatch.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))
    a,a0,c,b,e=sp.symbols('a a0 c b e')
    assert sp.simplify((e+(a-a0)*c*b+a0*c*b)/a-c*b-e/a)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Stable nonrest energy identity for Eref=Etilde+(a-a0)*cx*c^2*B only. No uniform EOS or complete GR certificate.'))


def source():
    run.prior.OUT=OUT;run.prior.initialize();pressure=pressure_function()
    model=gr.Response();rows=[];outputs=[]
    for n in [64,128]:
        start=time.monotonic();m=run.owner.Model(n);p=dict(np.load(saved(n)))
        audit,_,_,_=run.verify(saved(n),INPUT/f'sweep-1/photons/pilot-{n}.npz')
        material=m.material;times=m.t[:INTERVAL+1]
        ids=np.array([int(np.argmin(abs(p['t']-t))) for t in times]);assert np.max(abs(p['t'][ids]-times))<1e-18
        d=retained.template()
        for key,v in list(d.items()):
            if v.ndim and len(v)==17:d[key]=v[:len(times)].copy()
        d['t']=times.copy();coeff=model.coeff(material.rE)
        Eg0,Pg0,Kg0,Er0,Pr0=[coeff[key]*C**4/gr.base.G*material.V for key in ['Eg','Pg','Kg','Er','Pr']]
        b=model.model.bulk;weights=4*np.pi*np.r_[b.W,model.model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
        R40=model.ratio4*Er0;gas=[];photons=[];probes=[];pressure_match=[];identity=[];baryons=[]
        for k,j in enumerate(ids):
            t=float(times[k]);g=p['material_history'][j]/AMP
            z=np.array([g[:,2]*m.bu,g[:,3]*m.su,g[:,0]*m.eu,g[:,1]*m.nu],LD)
            field=m.geometry(t)[0];vol=(3*field[0]+field[2])*AMP
            bank=dict(np.load(feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz'))
            p0,delta,error=pressure(material,k,z,field,bank)
            point=material.point(k);q=point['Q'];backgroundE=(q[2]+material.rest*q[0])/material.a
            nonrest=z[2]*AMP/material.a+(Eg0+Pg0-backgroundE)*vol
            pg,pr=delta*AMP+(Kg0[None]-p0)*vol
            B=z[0]*AMP;rest=B*LD(m.model.cx)*LD(C)**2
            gas.append([nonrest,pg,pr]);baryons.append(B);probes.append(error*AMP)
            pressure_match.append(delta[0]*AMP-p['moments'][j,6])
            actual=m.conserved(g)*AMP
            identity.append(actual-p['conserved_material_history'][j])
            assert np.allclose(actual[[2,3]],p['moments'][j,[1,2]],rtol=1e-12,atol=0)
            phi=field[0]*AMP;lam=field[2]*AMP
            Ebg=np.einsum('nqf,nqf->n',m.I[k],weights)/material.a
            Pbg=np.einsum('nqf,nqf,q->n',m.I[k],weights,b.mu2)/material.a
            mom=p['moments'][j]
            photons.append([mom[0]/material.a-Ebg*vol+4*Er0*phi+(Er0+Pr0)*lam,
                mom[5]/material.a-Pbg*vol+4*Pr0*phi+(3*Pr0-R40)*lam])
        gas=np.array(gas);photons=np.array(photons);B=np.array(baryons);rest=B*LD(m.model.cx)*LD(C)**2
        d.update(baryon_g=B,gas_nonrest_energy_erg=gas[:,0],nonrest_trace_erg=gas[:,0]-gas[:,2]-2*gas[:,1],
            nonrest_stress_erg=gas[:,0]-gas[:,2],pressure_volume_erg=gas[:,1],photon_energy_erg=photons[:,0],
            photon_radial_pressure_erg=photons[:,1],metric_stress_erg=rest+gas[:,0]-gas[:,2]+photons[:,0]-photons[:,1],
            inner_cumulative_energy_erg=p['radial_ports'][ids,0,1],outer_cumulative_energy_erg=p['radial_ports'][ids,1,1])
        norm=max(float(np.max(np.sum(abs(p['moments'][ids,6]),axis=-1))),1e-290)
        probe=float(np.max(np.sum(abs(probes),axis=-1))/norm)
        mapping=float(np.max(np.sum(abs(pressure_match),axis=-1))/norm)
        errors=np.sum(abs(identity),axis=-1)/np.maximum(np.sum(abs(p['conserved_material_history'][ids]),axis=-1),LD('1e-290'))
        row=dict(clock=n,pressure_probe=probe,pressure_mapping=mapping,conserved_mapping=float(np.max(errors)),
            same_solution_audit=audit,seconds=time.monotonic()-start)
        rows.append(row);write(OUT/f'source-{n}-check.json',dict(classification='Counterexample candidate',**row))
        np.savez_compressed(OUT/'gr'/f'source-{n}.npz',**d)
        assert probe<.002 and mapping<1e-12 and np.max(errors)<1e-12,row
        outputs.append(np.concatenate([gas,photons],axis=1));del m,p;gc.collect()
    errors=run.owner.joint.previous.run.c.relative(*outputs)
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,source_time=errors,rows=rows,
        same_solution_energy_and_ports=True,current_nonzero_geometry=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'sources.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def fields():
    assert read(OUT/'sources.json')['passed']
    # Initialize background ownership once; no response constructor or evolution.
    run.prior.OUT=OUT;run.prior.initialize()
    model=gr.Response();model.run=FunctionType(gr.Response.run.__code__,dict(gr.Response.run.__globals__,OUT=OUT/'gr')).__get__(model,type(model))
    rows=[model.run(n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=dict(np.load(OUT/'gr/fields-128-g8.npz'));norm=max(np.max(abs(fine['U'])),1e-290)
    errors={name:float(np.max(abs(fine['U']-np.load(OUT/'gr'/file)['U']))/norm)
        for name,file in [('time','fields-64-g8.npz'),('quadrature','fields-128-g4.npz')]}
    d=dict(np.load(OUT/'gr/source-128.npz'));direct,coordinate=charge.independent.direct(model,d,8)
    errors['independent_GR']=abs(direct-rows[0]['endpoint_direct'])/max(abs(direct),1e-290)
    result=dict(classification='Counterexample candidate',passed=errors['time']<.02 and errors['quadrature']<.002 and errors['independent_GR']<1e-9,
        controls=errors,rows=rows,same_joint_solution_GR_fields_computed=True,physical_horizon_seconds=float(d['t'][-1]),
        inverse_radius_residual=float(coordinate),GR_return_to_matter_applied=False,full_horizon_completed=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'fields.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));run.owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            plan=read(OUT/'plan.json')
            for p,h in dict(plan['bindings'],**plan['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
