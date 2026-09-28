"""Independent stored-probe, flux, conservation and charge audit."""
from pathlib import Path
import json,sys,time
import numpy as np
import def_retained_native_acoustic as native
import apply_retained_native_acoustic as apply
read,write,sha=native.read,native.write,native.sha
OUT=native.OUT;R=apply.OUT;LD=np.longdouble


def audit():
    start=time.monotonic();native.deadline(60);counts=0;native_error=0.;law=0.;deep={}
    rows=read(OUT/'native.json')['rows']
    for row in rows:
        z=np.load(OUT/'hydro'/f"{row['label']}.npz");h=1e-4;s=z['samples']
        # Rebuild the derivative from raw constrained outputs independently.
        dx=((s[:,1,0]-s[:,1,1])*8-(s[:,0,0]-s[:,0,1]))/(6*h)
        dt=((s[:,1,2]-s[:,1,3])*8-(s[:,0,2]-s[:,0,3]))/(6*h)
        p=z['physical_pressure'];rho=z['raw'][:,0]
        K=dx[:,1]+dt[:,1]*(p/rho-dx[:,2])/dt[:,2]
        native_error=max(native_error,float(np.max(abs(K/z['K']-1))))
        law=max(law,float(np.max(abs(z['first_law']))));counts+=len(K)
        if row['kind']=='deep':deep[int(row['label'].split('-')[1])]={tuple(c):v for c,v in zip(z['coordinates'],K)}
    assert counts==6715 and native_error<1e-9 and law<.0001
    force_balance=0.;force_resolution=0.;gravity=0.
    for k in range(17):
        p=np.load(apply.FORCE/f'point-{k}.npz');f=p['flux'];g=p['gravity'];rate=-np.diff(f,axis=1);rate[1]+=g
        scale=np.maximum(np.sum(abs(rate),axis=1),1.)
        err=np.max(abs(np.sum(rate,axis=1)-p['ledger'])/scale);force_balance=max(force_balance,float(err))
        assert np.max(abs(rate-p['rate'])/np.maximum(abs(rate),1.))<1e-12
        gravity=max(gravity,float(np.max(abs(g))))
        force_resolution=max(force_resolution,read(apply.FORCE/f'point-{k}.json')['flux_resolution'])
    assert force_balance<1e-8 and force_resolution<.002 and gravity==0
    # Distinguish a stale initial modulus from the change during this history.
    native.prior.initialize();m=native.prior.Material(driven=False);roots=read(native.prior.EOS/'production-samples.json')
    K=[]
    for k in range(17):
        rr=sorted([r for r in roots if r['kind']=='deep' and r['it']==k],key=lambda r:r['cell'])
        K.append(np.interp(m.model.edge,m.model.bulk.d['r'],[deep[r['cell']][r['x'],r['lt'],r['y']] for r in rr]))
    K=np.array(K);old=m.model.face_K
    K_initial=float(np.max(abs(K[0]/old-1)));K_change=float(np.max(abs(K/K[0]-1)))
    K_discrepancy=float(np.max(abs(K/old-1)))
    apply.initialize_material(True);active=apply.Material(128,128);k=16;p=active.point(k)
    field=np.zeros((5,active.n));zero=np.zeros_like(p['Q']);probe=zero.copy();probe[1,10]=p['Q'][0,10]*native.prior.C**2*1e-8/native.prior.AMP
    f0=active.raw(k,zero,field,0.,p)[0];f1=active.raw(k,probe,field,native.prior.AMP,p)[0]
    saved=active.native_K[k].copy();active.native_K[k]=active.model.face_K.copy()
    try:
        g0=active.raw(k,zero,field,0.,p)[0];g1=active.raw(k,probe,field,native.prior.AMP,p)[0]
    finally:active.native_K[k]=saved
    positive=float(np.sum(abs((f1-f0)-(g1-g0)))/np.sum(abs(f1-f0)))
    assert positive==0 and np.any(f1!=f0) and np.array_equal(f1,g1)
    assert 'vf=(vl+vr)/2;ps=(pl+pr)/2' in active.extended_sources[-1]
    material_balance=0.;source_identity=0.;port=0.;prefix_resolution=0.;collision_round=0.;precision=0.;precision_points=0
    for n in [64,128]:
        for k in range(17):
            row=read(apply.SOURCE/str(n)/f'point-{k}.json');assert row['passed']
            for check in row.get('rows',[]):
                collision_round=max(collision_round,check['rounding_over_full_source'])
                if check.get('used_high_precision'):
                    precision_points+=1
                    p=np.load(apply.SOURCE/str(n)/f"precision-{k}-{check['factor']}.npz")
                    relative=float(np.max(abs(p['coarse']-p['fine']))/max(np.max(abs(p['fine'])),LD(1e-300)))
                    precision=max(precision,relative)
                    assert relative<1e-10 and check['precision_owner']<1e-9 and check['precision_root']<1e-40
    assert collision_round<.002 and precision_points>=4
    for n in [64,128]:
        p=np.load(apply.PHOTON/f'steps-{n}-reference-128.npz')
        q=np.load(apply.MATERIAL/f'steps-{n}-reference-128.npz')
        for folder in [apply.DIRECT,apply.MATERIAL]:
            for label in [f'pilot-{n}',f'steps-{n}-reference-128']:
                row=read(folder/f'{label}.json');assert row['passed']
                prefix_resolution=max(prefix_resolution,row['finite_resolution'])
                material_balance=max(material_balance,max(row['balance_relative']))
        gamma=1-1/np.sqrt(2);h=p['t'][-1]/n;w=np.arange(1,8,2,dtype=LD)/32
        packets=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*p['accepted_angular_luminosity'].astype(LD)
        cumulative=np.r_[LD(0),np.cumsum(packets@w)];ids=2*np.rint(p['t']/h).astype(int)
        port=max(port,float(np.max(abs(cumulative[ids]-p['radial_ports'][:,1,1]))/max(np.sum(abs(packets)@w),1.)))
        s=np.load(apply.GR/f'source-{n}-reference-128.npz');stress=np.load(apply.READOUT/f'stress-{n}-reference-128.npz')['material']
        rest=s['baryon_g'].astype(LD)*LD(s['cx'])*LD(native.prior.C)**2;energy=rest+s['gas_nonrest_energy_erg']
        errors=[energy-stress[:,0],s['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
                s['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),s['pressure_volume_erg']-stress[:,3],
                s['metric_stress_erg']-(energy+s['photon_energy_erg']-stress[:,1]-s['photon_radial_pressure_erg'])]
        source_identity=max(source_identity,float(max(np.max(abs(e)) for e in errors)/max(np.max(abs(stress)),1e-300)))
        assert np.any(q['delta_scaled']) and np.any(p['delta_packet_scaled_occupation'])
    assert material_balance<1e-8 and prefix_resolution<.002 and source_identity<1e-12 and port<1e-12
    pressure=read(apply.READOUT/'pressure-precision.json');assert pressure['passed'] and len(pressure['rows'])==34
    assert all(row['owner']<1e-9 and row['precision']<1e-10 and row['root']<1e-60 for row in pressure['rows'])
    assert not read(apply.MATERIAL/'sources.json')['passed'] and read(apply.READOUT/'sources.json')['passed']
    assert read(R/'readout-serialization-check.json')['passed']
    norm=0.;baseline=Path('retained-metric-return152-work/infinity');model=apply.prior.template();alpha=-float(model['K_cm']/model['M_cm'])
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        v=np.load(R/f'infinity/charge-{n}-a{a}-r{r}.npz');before=np.load(baseline/f'charge-{n}-a{a}-r{r}.npz')
        wave=np.load(apply.GR/f'wave-{n}-g8.npz');eps=v['epsilon'];de=eps-before['epsilon']
        # Use the separately saved packet increment to avoid subtracting two
        # large mass losses in the independent normalization identity.
        import def_native_global_scalar_closure as ext
        de=ext.G/ext.C**4*v['arrived_parts_erg'][:,1]/float(model['M_cm'])
        value=(wave['free_scalar']+v['normalized_exterior_parts'][:,1]+(alpha+before['normalized'])*de)/(1-before['epsilon']-de)
        norm=max(norm,float(np.max(abs(value-v['acoustic_return']))/max(np.max(abs(value)),1e-300)))
        assert np.array_equal(v['previous'],before['normalized'])
    assert norm<1e-10
    final=read(R/'infinity/result.json');assert final['passed'] and read(R/'result.json')['passed']
    bound=0
    for name,producer in [('plan.json','registered-pilot-producer.py'),('warm-plan.json','registered-warm-producer.py'),('native-execution-plan.json','registered-native-producer.py')]:
        for p,h in read(OUT/name)['bindings'].items():
            target=OUT/producer if Path(p).name==Path(native.__file__).name else Path(p)
            assert sha(target)==h,(name,p);bound+=1
    result=dict(classification='Counterexample candidate',passed=True,native_states=counts,bindings_checked=bound,
        native_derivative_rebuild=native_error,native_first_law=law,force_balance=force_balance,force_resolution=force_resolution,
        gravity_increment_exact_zero=gravity==0,initial_deep_K_difference=K_initial,history_deep_K_change=K_change,
        maximum_deep_K_difference=K_discrepancy,material_balance=material_balance,maximum_prefix_or_full_resolution=prefix_resolution,
        centered_operator_K_null_control=positive,
        source_identity=source_identity,angular_port=port,independent_normalization=norm,seconds=time.monotonic()-start,
        collision_arithmetic=collision_round,high_precision_saved_controls=precision_points,independent_precision_comparison=precision,
        finite_pressure_high_precision_controls=34,original_pressure_readout_rejected=True,
        physical_amplitude_unchanged=True,full_goal_complete=False,uniform_derivative_certificate=False,coupled_fixed_point_verified=False)
    write(OUT/'verification.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':audit()
