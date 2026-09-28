"""Counterexample candidate: source, ports and separated-response controls."""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import sympy as sp
import solve_native_incident_self_gr as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha;LD=np.longdouble


def prepare_block(sweep):
    r=read(OUT/f'sweep-{sweep}/residual.json')
    assert not r['finite_reciprocal_residual_passed'] and not r['next_sweep_allowed']
    files=[Path(__file__),Path(run.__file__),Path(run.coupled.__file__),Path(run.coupled.feedback.__file__),
           OUT/'common-input-plan.json',OUT/'common-input-producer.py',OUT/f'sweep-{sweep}/residual.json']
    files += [p/f'steps-{n}-reference-128.npz' for s in [sweep-1,sweep] for p in run.paths(s) for n in [64,128]]
    # E/H is solved within the photon block, not used as a lagged input except
    # through M=(Etilde,H)-C. A change in total E/H is not its equation defect.
    y,c,m,d=sp.symbols('y c m d')
    assert sp.expand((y+d)-(c+d)-(y-c))==0
    assert sp.expand((y-c)-m)==y-c-m
    write(OUT/'block-plan.json',dict(classification='Counterexample candidate',sweep=sweep,
        original_verdict='The registered whole-state change/monotonicity test failed. Stop iteration; no third response is authorized by this audit and no old verdict is changed.',
        question='Does the completed pair satisfy the actual lagged block inputs and all applied stage forcing, despite the large change of E/H already solved simultaneously with photons?',
        operator='L at fixed retained background/GR is independent of the lagged material history. Affine input depends only on B,S,xi(B),geometry and dM/dt, with M=(Etilde,H)-C. Test both64/128 actual SDIRK stages, holding the first-stage mechanical rate at the closing stage exactly as the owner does.',
        checks='All lagged B/S/xi/M histories and interval rates; physical energy/number weighted photon forcing and separate gasE/H rates; unchanged stage matrices; photon/free-material E/H match. Source defects are normalized by maximum stage L1 of the current corresponding source, never by mixed units.',
        gates=dict(block=.002,paired=.002,time=.02,operator=0.,source_identity=1e-12),
        scope='A separate finite block-equation consistency audit, not the failed whole-state-change convergence criterion, a continuum contraction, or a solution-error bound. Compact readout may be reported only under this explicitly narrower criterion.',
        symbolic=dict(classification='Proven',passed=True,identity='Adding the same collision transfer to total gas and C leaves M unchanged. Thus total gas iterate change is not by itself the residual of the lagged-input equations.'),
        budget_seconds=300,total_action_seconds=4600,new_evolution_steps=0,
        forecast='Use completed17-point assembly costs and a measured warmed stage-pair cost with2x margin. No Krylov solve, new state, clock, horizon or background.',
        bindings={str(p):sha(p) for p in files}))


def block(sweep):
    plan=read(OUT/'block-plan.json');assert plan['sweep']==sweep
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    run.initialize(sweep);rows=[];gamma=1-1/np.sqrt(2);started=time.monotonic()
    for n in [64,128]:
        setup=time.monotonic();m=run.coupled.Response(n);photon,material=run.paths(sweep)
        p=np.load(photon/f'steps-{n}-reference-128.npz');d=np.load(material/f'steps-{n}-reference-128.npz')
        ids=[int(np.argmin(abs(d['t']-t))) for t in m.t];motion=d['history_scaled'][ids].copy()
        offset=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*LD(run.C)**2*motion[:,0];motion[:,2]-=offset
        new=dict(motion=np.asarray(motion,float),energy_offset=np.asarray(offset,float),
            mechanical=np.asarray(motion[:,[2,3]]-p['collision_transfer'].transpose(0,2,1)/run.AMP,float),
            xi=np.array([m.material.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:m.nb])] for z in motion]))
        old={k:getattr(m,k) for k in new}
        def relative(a,b):return float(np.max(np.sum(abs(a-b),axis=-1))/max(np.max(np.sum(abs(b),axis=-1)),1e-290))
        inputs={k:relative(old[k],new[k]) for k in ['xi','energy_offset']}
        for j,name in enumerate(['B','S']):inputs[name]=relative(old['motion'][:,j],new['motion'][:,j])
        for j,name in enumerate(['Etilde','H']):
            inputs['M_'+name]=relative(old['mechanical'][:,j],new['mechanical'][:,j])
            inputs['dM_'+name]=relative(np.diff(old['mechanical'][:,j],axis=0),np.diff(new['mechanical'][:,j],axis=0))
        def select(v):
            for k,w in v.items():setattr(m,k,w)
        def pair(t,mechanical=None):
            select(old);a=m.local(t);select(new);b=m.local(t)
            for k in ['loss','sc','B','Bb','Be','pressure_map']:assert np.array_equal(a[k],b[k]),k
            for k in ['S','esc']:
                diff=a[k]-b[k]
                assert not (np.any(diff.data) if hasattr(diff,'nnz') else np.any(diff)),k
            if mechanical is not None:a['mechanical'],b['mechanical']=mechanical
            s=m.source(t)[0]/(m.scale*run.AMP)
            x=a['q']+s;y=b['q']+s
            ga=m.gas(a['q'],a['qb'],a['qe'])+a['mechanical']
            gb=m.gas(b['q'],b['qb'],b['qe'])+b['mechanical']
            difference=[np.sum(abs(x-y)*m.Eweight),np.sum(abs(x-y)*m.Nweight),
                np.sum(abs(ga[:,0]-gb[:,0])*m.eu),np.sum(abs(ga[:,1]-gb[:,1])*m.nu)]
            norm=[np.sum(abs(y)*m.Eweight),np.sum(abs(y)*m.Nweight),np.sum(abs(gb[:,0])*m.eu),np.sum(abs(gb[:,1])*m.nu)]
            return np.asarray(difference,LD),np.asarray(norm,LD),(a['mechanical'].copy(),b['mechanical'].copy())
        h=m.t[-1]/n;pair(gamma*h);mark=time.monotonic();pair(gamma*h);warm=time.monotonic()-mark
        point=m.point_seconds/m.point_count;forecast=2*(2*17*point+384*warm+2*(mark-setup-m.point_seconds)+10)
        assert forecast<300,('block audit forecast',forecast)
        maximum=np.zeros(4,LD);scale=np.zeros(4,LD)
        for k in range(n):
            a,b,mech=pair(k*h+gamma*h);maximum=np.maximum(maximum,a);scale=np.maximum(scale,b)
            a,b,_=pair(k*h+h,mech);maximum=np.maximum(maximum,a);scale=np.maximum(scale,b)
        source=(maximum/np.maximum(scale,LD('1e-290'))).astype(float).tolist()
        paired=read(OUT/f'sweep-{sweep}/residual.json')['rows'][[64,128].index(n)]['photon_material_energy_H_residual']
        rows.append(dict(steps=n,lagged_inputs=inputs,stage_forcing_defect=source,paired_E_H=paired,
            stage_count=2*n,stage_operator_identical=True,forecast_seconds=forecast))
        assert max(list(inputs.values())+source+paired)<.002,rows[-1]
    maximum=max(max(list(r['lagged_inputs'].values())+r['stage_forcing_defect']+r['paired_E_H']) for r in rows)
    write(OUT/'block-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        maximum_block_defect=maximum,original_waveform_change_test_passed=False,original_stop_preserved=True,
        no_new_evolution_steps=True,scope=plan['scope'],seconds=time.monotonic()-started))
    print(json.dumps(read(OUT/'block-result.json')),flush=True)


def prepare(sweep):
    out=OUT/f'sweep-{sweep}';files=[Path(__file__),Path(run.__file__),out/'result.json',OUT/'normalization.json']
    files += [p/f'steps-{n}-reference-128.npz' for p in run.paths(sweep) for n in [64,128]]
    files += [run.BEFORE/f'{folder}/steps-128-reference-128.npz' for folder in ['photons-precise','material-analytic']]
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',sweep=sweep,
        claim='Independently check conserved GR source identities, collision/mechanical partition, actual angular ports, preserved prefixes, and sampled additivity/homogeneity of the material directional operator before adding the separated charge.',
        reason='Positive amplitude normalization is homogeneous, but an HLL/donor directional derivative need not be globally additive at a branch corner. Test the actual saved response directions and report the finite sampled scope.',
        gates=dict(source=1e-12,partition=1e-8,port=1e-12,conservation=1e-8,material_linearity=.002),
        seconds=90,total_action_seconds=4600,new_evolution_steps=0,
        stop='No new response path or relaxed gate if a sampled branch-additivity check fails. Full uniform EOS/branch certification remains open.',
        bindings={str(p):sha(p) for p in files}))


def audit(sweep):
    plan=read(OUT/'audit-plan.json');assert plan['sweep']==sweep
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    photon,material=run.paths(sweep);pp,mm=run.paths(sweep-1);rows=[]
    for n in [64,128]:
        p=np.load(photon/f'steps-{n}-reference-128.npz');m=np.load(material/f'steps-{n}-reference-128.npz')
        oldp=np.load(pp/f'steps-{n}-reference-128.npz');oldm=np.load(mm/f'steps-{n}-reference-128.npz')
        d=np.load(OUT/f'sweep-{sweep}/gr/source-{n}-reference-128.npz')
        stress=np.load(material/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(run.C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-stress[:,0],d['nonrest_trace_erg']+rest-(stress[:,0]-stress[:,1]-2*stress[:,3]),
            d['nonrest_stress_erg']+rest-(stress[:,0]-stress[:,1]),d['pressure_volume_erg']-stress[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-stress[:,1]-d['photon_radial_pressure_erg'])]
        source=float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(stress)),1e-290))
        ids=[int(np.argmin(abs(oldm['t']-t))) for t in p['t']]
        actual=p['moments'][:,[1,2]];old=oldm['history_scaled'][ids][:,[2,3]]*run.AMP
        a=actual-p['collision_transfer'].transpose(0,2,1);b=old-oldp['collision_transfer'].transpose(0,2,1)
        norm=np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),np.max(np.sum(abs(old),axis=2),axis=0))
        partition=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(norm,1e-290)).astype(float).tolist()
        h=p['t'][-1]/n;gamma=1-1/np.sqrt(2)
        times=h*(np.arange(n)[:,None]+[gamma,1.]).ravel();assert np.max(abs(times-p['accepted_angular_times']))<1e-18
        packet=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*p['accepted_angular_luminosity']
        flux=packet@(np.arange(1,8,2)/32);cumulative=np.r_[0.,np.cumsum(flux,dtype=LD)]
        port=float(np.max(abs(cumulative[2*np.rint(p['t']/h).astype(int)]-p['radial_ports'][:,1,1]))/max(np.sum(abs(flux)),1e-290))
        balance=float(np.max(abs(np.sum(m['history_scaled'],axis=2,dtype=LD)+m['discards_scaled']-m['ledgers_scaled'])/np.maximum(m['norms_scaled'],1.)))
        prefix=np.load(photon/f'pilot-{n}.npz')
        for k in ['t','moments','material_history','collision_transfer','radial_ports','accepted_angular_times','accepted_angular_luminosity']:
            assert np.array_equal(p[k][:len(prefix[k])],prefix[k]),k
        assert source<1e-12 and max(partition)<1e-8 and port<1e-12 and balance<1e-8
        rows.append(dict(steps=n,source_identity=source,mechanical_collision_partition=partition,angular_port=port,material_balance=balance,accepted_prefix_preserved=True))
    run.initialize(sweep);m=run.coupled.Material(128,128)
    own=np.load(material/'steps-128-reference-128.npz');old=np.load(run.BEFORE/'material-analytic/steps-128-reference-128.npz')
    current_transfer=m.transfer.copy();run.coupled.transfers(m,run.BEFORE/'photons-precise',128);original_transfer=m.transfer.copy()
    primary=run.coupled.base.drive.Driver(8);saved=run.SavedMetric(128);linearity=[]
    class Blend:
        def __init__(self,a,b):self.a,self.b=a,b
        def view(self,t,clock,side='left'):
            g=primary.view(t,clock,side);h=saved.view(t,clock,side)
            return {k:self.a*g[k]+self.b*h.get(k,0.) for k in g}
    def rhs(t,z,a,b):
        m.driver=Blend(a,b);m.transfer=a*original_transfer+b*current_transfer
        return m.rhs(t,z)[0]
    for k in [1,8,16]:
        t=m.t[k];i=int(np.argmin(abs(own['t']-t)));j=int(np.argmin(abs(old['t']-t)))
        z1=own['history_scaled'][i];z0=old['history_scaled'][j]
        q=m.point(k)['Q'];units=np.maximum(abs(q),1.);units[1]=np.maximum(q[0]*run.C**2,1.)
        scale=float(np.max(abs(z1)/units)/max(np.max(abs(z0)/units),1e-290))
        r0=rhs(t,z0*scale,scale,0.);r1=rhs(t,z1,0.,1.)
        summed=rhs(t,z0*scale+z1,scale,1.);half=rhs(t,z1*.5,0.,.5)
        norm=np.maximum(np.sum(abs(r0),axis=1)+np.sum(abs(r1),axis=1),1.)
        add=(np.sum(abs(summed-r0-r1),axis=1)/norm).astype(float).tolist()
        hom=(np.sum(abs(2*half-r1),axis=1)/np.maximum(np.sum(abs(r1),axis=1),1.)).astype(float).tolist()
        assert max(add+hom)<.002,(k,add,hom)
        linearity.append(dict(k=k,time=float(t),balanced_original_direction_scale=scale,additivity=add,positive_homogeneity=hom))
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,material_directional_controls=linearity,
        scope='Sampled directions on the retained finite background only; not uniform branch/native derivative or continuum certification.',full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]);assert action in ['prepare','audit','prepare_block','block']
    cap=300 if action=='block' else 90
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.coupled.base.drive.native.deadline(cap)
    start=time.monotonic();cpu=time.process_time();error=None;receipt=OUT/f'audit-{action}-receipt.json';assert not receipt.exists()
    try:
        assert sum(read(p)['seconds'] for p in OUT.rglob('*-receipt.json'))+cap<=4600
        globals()[action](sweep)
    except Exception as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
