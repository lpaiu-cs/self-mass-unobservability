"""Request33: connect the repaired EOS to material GR and native sources.

Counterexample candidate. EOS replacement is a change of model, not evolution.
The same nuclear inventories and baryon cells are retained without a mass refit.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, sys
import numpy as np
import direct_ion_eos as d
import audit_structured_enthalpy as audit

e=d.e;s=e.s;c=e.c;ROOT=e.ROOT;OLD=e.OUT
OUT=ROOT/'outputs/direct-eos-gr33'
CACHE=Path('/home/lpaiu/work/direct-eos-gr33')


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


class EOS(d.EOS):
    """The 21 defined outputs of EOS1; gamma_e is not an available result."""
    def __init__(self):
        super().__init__('full');native=self.call
        def defined_call(mode,value,t,eps,out,info):
            native(mode,value,t,eps,out,info)
            # The disabled diffraction branch never assigns gamma_e. This is
            # private padding for the legacy finite check, discarded below;
            # zero is never exposed as a physical diffraction parameter.
            out[20]=0.
        self.call=defined_call

    def __call__(self,*args): return np.delete(super().__call__(*args),20)


def prepare():
    assert not (OUT/'plan.json').exists()
    control=json.loads((d.OUT/'full-integral-mixture-control.json').read_text())
    assert control['passed'] and control['cells']==5735
    OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
    base=dict(np.load(OLD/'initial-state.npz'))
    values=dict(np.load(d.OUT/'full-integral-mixture-control.npz'))['values']
    assert values.shape==(len(base['X']),22) and np.all(values[:,0]>0)
    assert np.all(abs(np.log(values[:,0])-base['lnd'])<1e-10)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='98baf6e',
        objective='Use the same repaired 24-element EOS in GR constraints and fresh native reaction inputs, then advance physical transport and the remaining full original closure goal.',
        inputs_sha256={str(p.relative_to(ROOT)):c.sha(p) for p in [OLD/'initial-state.npz',
            d.OUT/'full-integral-mixture-control.json',d.OUT/'full-integral-mixture-control.npz',
            d.OUT/'full-integral-build.json',ROOT/'verification/direct_eos_gr.py']},
        sequence=['Pressure/entropy inverse and call-history controls with the new EOS',
            'Recompute GR at fixed baryon cells and nuclear abundances using new reference entropies',
            'Reevaluate the native reaction vector with the same new EOS auxiliaries',
            'Recompute thermal coordinates and advance nonlinear transport and GR time paths'],
        initial_state='Keep old physical rho_B,T,X only to define the new EOS reference entropy. The inherited geometry is an initial guess, not already a new GR solution. No equality of old and new physical entropy is assumed.',
        GR=dict(points=17,subdivision=4,interface_tolerance=1e-8,baryon_masses_fixed=True,
            table_processes=8,table_block_cells=128,
            boundary='The same finite photospheric pressure. Radius and gravitational mass are outputs of centre/surface matching; baryon mass is not retuned.',
            entropy_inverse='abs(T*(s-s_target)) <= max(2 erg/g,32 ulp(abs(H)))'),
        controls=dict(inverse_log_density=1e-10,known_entropy_log_temperature=1e-10,
            call_history_bitwise=True,auxiliary_finite_log_steps=[5e-5,2.5e-5],auxiliary_derivative_relative=1e-3),
        physical_EOS_certified=False,continuous_EOS_certified=False,full_GR_evolution=False,
        policy='Preserve prior paths and failures. A model-replacement mass shift is not temporal energy release. Isotope partition physics, physical plasma uncertainty, transport/atmosphere, scalar driving and observation closure remain required.'))
    reference={**base,'logP':np.log(values[:,1]),'u_W':values[:,2],'s_B':values[:,3],
        'CX':(base['X']/c.A)@c.W}
    np.savez_compressed(OUT/'reference-state.npz',**reference)
    save('reference-state.json',dict(classification='Counterexample candidate',
        original_state_sha256=c.sha(OLD/'initial-state.npz'),reference_sha256=c.sha(OUT/'reference-state.npz'),
        same_baryon_and_nuclear_inventories=True,already_GR_solution=False,
        maximum_pressure_relative_change=float(abs(np.expm1(reference['logP']-base['logP'])).max())))
    print('PREPARED DIRECT EOS GR',len(base['dm']),'fixed baryon cells',flush=True)


def inverse_control():
    plan=json.loads((OUT/'plan.json').read_text());base=dict(np.load(OUT/'reference-state.npz'))
    e.OUT=OUT;audit.ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label='inverse-control')
    eos=EOS();indices=np.unique(np.r_[np.linspace(0,len(base['X'])-1,65).astype(int),
        np.argmax(base['X'],axis=0),2972]);rows=[]
    for k,i in enumerate(indices):
        x=base['X'][i];t=base['lnT'][i];rho=base['lnd'][i];lp=base['logP'][i]
        a=eos(2,rho,t,x);inverse=eos(1,lp,t,x)
        j=int(indices[(k+1)%len(indices)]);eos(2,base['lnd'][j],base['lnT'][j],base['X'][j])
        repeated=eos(2,rho,t,x)
        _,recovered,_=audit.strict_invert(eos,lp,inverse[3],x,t+(.01 if k%2 else -.01))
        row=dict(cell=int(i),log_density_error=float(abs(np.log(inverse[0])-rho)),
            log_temperature_error=float(abs(recovered-t)),call_history_bitwise=bool(np.array_equal(a,repeated)))
        row['passed']=bool(row['log_density_error']<plan['controls']['inverse_log_density'] and
            row['log_temperature_error']<plan['controls']['known_entropy_log_temperature'] and row['call_history_bitwise'])
        rows.append(row);print('DIRECT EOS INVERSE',row,flush=True)
    passed=all(row['passed'] for row in rows)
    save('inverse-control.json',dict(classification='Counterexample candidate',passed=passed,rows=rows,
        root_statistics=audit.ROOT_STATS,physical_or_continuous_certificate=False))
    assert passed


def table_block(start):
    """Disjoint cells, isolated serial Fortran state, the same strict roots."""
    data=dict(np.load(OUT/'reference-state.npz'));plan=json.loads((OUT/'plan.json').read_text())
    stop=min(start+plan['GR']['table_block_cells'],len(data['X']))
    stem=f'initial-table-{start}';path=OUT/(stem+'.npz');record=OUT/(stem+'.json')
    binding=c.sha(OUT/'reference-state.npz');plan_hash=c.sha(OUT/'plan.json')
    if path.exists() and record.exists():
        saved=json.loads(record.read_text())
        assert saved['reference_sha256']==binding and saved['plan_sha256']==plan_hash
        assert saved['output_sha256']==c.sha(path)
        return saved
    e.OUT=OUT;eos=EOS();audit.ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label=stem)
    offsets=np.linspace(-.04,.04,plan['GR']['points']);reference=[];values=[]
    for i in range(start,stop):
        x=data['X'][i];lt=data['lnT'][i];a=eos(2,data['lnd'][i],lt,x);reference.append(a);row=[]
        for offset in offsets:
            b,t,_=audit.strict_invert(eos,np.log(a[1])+offset,a[3],x,lt)
            row.append([np.log(b[0]),t,b[2]*1e-4/c.gr.C**2])
        values.append(row)
    np.savez_compressed(path,reference=np.array(reference),values=np.array(values))
    saved=dict(classification='Counterexample candidate',start=start,stop=stop,
        reference_sha256=binding,plan_sha256=plan_hash,output_sha256=c.sha(path),root_statistics=audit.ROOT_STATS.copy())
    save(record.name,saved);return saved


def table():
    assert json.loads((OUT/'inverse-control.json').read_text())['passed']
    data=dict(np.load(OUT/'reference-state.npz'));plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert c.sha(ROOT/rel)==digest,rel
    starts=list(range(0,len(data['X']),plan['GR']['table_block_cells']));records=[]
    with ProcessPoolExecutor(max_workers=plan['GR']['table_processes']) as pool:
        for completed in as_completed([pool.submit(table_block,k) for k in starts]):
            records.append(completed.result())
            save('table-progress.json',dict(completed_cells=sum(r['stop']-r['start'] for r in records),
                total_cells=len(data['X']),completed_blocks=len(records),total_blocks=len(starts)))
            print('DIRECT EOS TABLE',len(records),'/',len(starts),'blocks',flush=True)
    arrays=[dict(np.load(OUT/f'initial-table-{k}.npz')) for k in starts]
    ref=np.concatenate([a['reference'] for a in arrays]);values=np.concatenate([a['values'] for a in arrays])
    offset=np.linspace(-.04,.04,plan['GR']['points']);eos=EOS();e.OUT=OUT
    audit.ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label='table-serial-control')
    controls=[]
    for i in [0,2972,len(ref)-1]:
        a=eos(2,data['lnd'][i],data['lnT'][i],data['X'][i]);assert np.array_equal(a,ref[i])
        for j in [0,8,16]:
            b,t,_=audit.strict_invert(eos,np.log(a[1])+offset[j],a[3],data['X'][i],data['lnT'][i])
            answer=np.array([np.log(b[0]),t,b[2]*1e-4/c.gr.C**2])
            assert np.array_equal(answer,values[i,j]);controls.append([i,j])
    np.savez_compressed(OUT/'initial-adiabats-17.npz',reference=ref,values=values,offset=offset,
        X=data['X'],lnT=data['lnT'],lnd=data['lnd'])
    save('table-control.json',dict(classification='Counterexample candidate',passed=True,
        independent_serial_bitwise_controls=controls,blocks=sorted(records,key=lambda r:r['start']),
        output_sha256=c.sha(OUT/'initial-adiabats-17.npz'),physical_or_continuous_certificate=False))
    print('DIRECT EOS TABLE COMPLETE',len(ref),'cells; serial controls bitwise equal',flush=True)


def project():
    assert json.loads((OUT/'inverse-control.json').read_text())['passed']
    control=json.loads((OUT/'table-control.json').read_text());assert control['passed']
    assert control['output_sha256']==c.sha(OUT/'initial-adiabats-17.npz')
    data=dict(np.load(OUT/'reference-state.npz'));plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert c.sha(ROOT/rel)==digest,rel
    parameters=json.loads((s.v.OLD/'restored-structure-17-4.json').read_text())['parameters']
    # Reuse the existing GR/root algorithms in this isolated process.
    c.OUT=OUT;c.EOS=EOS;c.be.invert=audit.strict_invert;e.OUT=OUT
    audit.ROOT_STATS.update(calls=0,evaluations=0,maximum_score=0.,label='initial-GR')
    state,record=c.structure('initial',data,17,4,parameters)
    assert np.array_equal(state['X'],data['X']) and np.array_equal(state['dm'],data['dm'])
    original=dict(np.load(OLD/'initial-state.npz'))
    save('initial-GR.json',dict(classification='Counterexample candidate',completed=True,
        state_sha256=c.sha(OUT/'initial-state-17-4.npz'),record=record,root_statistics=audit.ROOT_STATS,
        same_baryon_and_nuclear_inventories=True,
        maximum_log_temperature_change=float(abs(state['lnT']-original['lnT']).max()),
        radius_change_m=float(state['radius_faces_m'][0]-original['radius_faces_m'][0]),
        gravitational_mass_change_geom_m=float(state['mass_faces_geom'][0]-original['mass_faces_geom'][0]),
        model_replacement_not_time_evolution=True,physical_EOS_certified=False,full_GR_evolution=False))
    print('DIRECT EOS GR COMPLETE',record,audit.ROOT_STATS,flush=True)


def reference_source():
    """Separate fixed-state source control; the GR solution is not required."""
    assert json.loads((OUT/'inverse-control.json').read_text())['passed']
    data=dict(np.load(OUT/'reference-state.npz'));plan=json.loads((OUT/'plan.json').read_text())
    eos=EOS();rows=[];checks=[]
    for i,(r,t,x) in enumerate(zip(data['lnd'],data['lnT'],data['X'])):
        a=eos(2,r,t,x);slopes=[]
        for h in plan['controls']['auxiliary_finite_log_steps']:
            slopes.append([(eos(2,r,t+h,x)[12]-eos(2,r,t-h,x)[12])/(2*h),
                (eos(2,r+h,t,x)[12]-eos(2,r-h,t,x)[12])/(2*h)])
        slopes=np.array(slopes);score=float(np.max(abs(slopes[1]-slopes[0])/np.maximum(1,abs(slopes[1]))))
        rows.append([a[13]/np.exp(r),a[12],*slopes[1]])
        checks.append(dict(cell=i,score=score,passed=score<plan['controls']['auxiliary_derivative_relative']))
        if i%500==0: print('DIRECT EOS SOURCE AUXILIARY',i,flush=True)
    values=np.array(rows);np.save(OUT/'reference-auxiliary.npy',values)
    passed=all(r['passed'] for r in checks)
    save('reference-auxiliary.json',dict(classification='Counterexample candidate',passed=passed,
        rows=checks,maximum_score=max(r['score'] for r in checks),
        state_sha256=c.sha(OUT/'reference-state.npz'),physical_or_continuous_certificate=False))
    assert passed
    s.OUT=OUT;s.CACHE=CACHE;s.auxiliary=lambda state:values
    raw,corrected=s.evaluate('reference',data)
    assert np.array_equal(raw['aux_used'],values)
    previous=dict(np.load(OLD/'initial-corrected.npz'))
    changes={name:float(np.max(abs(corrected[name]-previous[name])/
        np.maximum(1e-30,np.maximum(abs(corrected[name]),abs(previous[name])))))
        for name in ['dxdt','heat','neutrino']}
    save('reference-source.json',dict(classification='Counterexample candidate',completed=True,
        new_EOS_inputs_confirmed_bitwise=True,fixed_rho_T_X=True,relative_model_changes=changes,
        native_source_sha256=c.sha(OUT/'reference-native.npz'),corrected_source_sha256=c.sha(OUT/'reference-corrected.npz'),
        common_physical_EOS_certified=False,GR_time_evolution=False))
    print('DIRECT EOS FIXED-STATE SOURCE',changes,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
