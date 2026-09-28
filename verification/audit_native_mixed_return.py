"""Counterexample candidate: audit the actual mixed-field transport return."""
from pathlib import Path
from types import SimpleNamespace
import inspect,json,resource,sys,time
import numpy as np
import sympy as sp
import return_native_mixed_gr as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha;LD=run.LD;CAP=90


def block_retry(sweep,source_only=False):
    folder=OUT/f'sweep-{sweep}';failed=read(OUT/f'block-{sweep}-receipt.json')
    spent=failed['seconds'];cap=285
    if source_only:
        spent+=read(OUT/f'audit-block_retry-{sweep}-receipt.json')['seconds'];cap=271
    assert 'block audit forecast' in failed['error'] and spent+cap<300
    files=[Path(__file__),Path(run.__file__),Path(run.audit_owner.__file__),
        Path(run.reuse.coupled.feedback.__file__),folder/'residual.json',OUT/f'block-{sweep}-receipt.json']
    files += [p/f'steps-{n}-reference-128.npz' for s in [sweep-1,sweep] for p in run.paths(s) for n in [64,128]]
    plan=folder/('block-source-plan.json' if source_only else 'block-retry-plan.json');assert not plan.exists()
    write(plan,dict(classification='Counterexample candidate',sweep=sweep,
        scope='Same actual lagged inputs, all SDIRK stage sources and paired E/H. Omit only collision(zero,zero,source=False) from the local audit, after exact comparisons to the unchanged owner. No new evolution or altered gate.',
        source_difference='When enabled, evaluate the unchanged new stage once and contract the exact D/Db/De coefficient maps with old-minus-new B/S/xi; mechanical differences are explicit. Compare against both original full stage evaluations at early/middle/late actual stages.',
        budget_seconds=cap,original_block_budget_seconds=300,failed_seconds=spent,
        bindings={str(p):sha(p) for p in files}))
    checks=[]
    def initialize(s):
        run.initialize(s);owner=run.reuse.coupled.Response
        class Response(owner):
            def local(self,t):
                self.set_stage(t);c=run.reuse.coupled.feedback.Response.local(self,t)
                k=int(np.clip(np.searchsorted(self.t,t)-1,0,16));seen=getattr(self,'zero_check_counts',{})
                if k in [0,8,15] and seen.get(k,0)<2:
                    reference=owner.local(self,t)
                    assert set(c)==set(reference)
                    for key,v in c.items():
                        delta=v-reference[key]
                        assert not (np.any(delta.data) if hasattr(delta,'nnz') else np.any(delta)),key
                    checks.append(dict(steps=self.audit_steps,k=k,exact=True));seen[k]=seen.get(k,0)+1
                    self.zero_check_counts=seen
                return c
            def __init__(self,n):super().__init__(n);self.audit_steps=n
        run.reuse.coupled.Response=Response
    source=inspect.getsource(run.audit_owner.block)
    source=source.replace("OUT/'block-plan.json'", repr(str(plan)))
    source=source.replace("OUT/'block-result.json'", "OUT/f'sweep-{sweep}/block-result.json'")
    source=source.replace("assert max(list(inputs.values())+source+paired)<.002,rows[-1]", '# Preserve a failed finite block rather than accepting it.')
    source=source.replace('passed=True,rows=rows,','passed=maximum<.002,rows=rows,')
    source=source.replace('assert forecast<300,',f'assert forecast<{cap},')
    source_checks=[]
    if source_only:
        anchor='        h=m.t[-1]/n;pair(gamma*h);'
        assert source.count(anchor)==1
        source=source.replace(anchor,'''        full_pair=pair;seen={}
        def pair(t,mechanical=None):
            select(new);b=m.local(t)
            k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));f=(t-m.t[k])/(m.t[k+1]-m.t[k])
            blend=lambda v:(1-f)*v[k]+f*v[k+1]
            motion=blend(old['motion']-new['motion']);xi=np.r_[blend(old['xi']-new['xi']),np.zeros(m.n-m.nb)]
            dq=np.zeros_like(b['q']);db=np.zeros_like(b['qb']);de=np.zeros_like(b['qe'])
            for w,r in [(1-f,m.point(k)),(f,m.point(k+1))]:
                drive=np.stack([motion[0]/r['units'][0],motion[1]/r['units'][1],np.zeros(m.n),xi],axis=-1)
                dq+=w*np.einsum('nqfj,nj->nqf',r['D'],drive)
                db+=w*np.einsum('nqfj,nj->nqf',r['Db'],drive)
                de+=w*np.einsum('knj,nj->kn',r['De'],drive)
            ma=np.diff(old['mechanical'][k:k+2],axis=0)[0].T/(m.t[k+1]-m.t[k])/np.stack([m.eu,m.nu],axis=-1)
            mb=b['mechanical']
            if mechanical is not None:ma,mb=mechanical
            s=m.source(t)[0]/(m.scale*run.AMP);y=b['q']+s
            gb=m.gas(b['q'],b['qb'],b['qe'])+mb;dg=m.gas(dq,db,de)+ma-mb
            difference=np.asarray([np.sum(abs(dq)*m.Eweight),np.sum(abs(dq)*m.Nweight),np.sum(abs(dg[:,0])*m.eu),np.sum(abs(dg[:,1])*m.nu)],LD)
            norm=np.asarray([np.sum(abs(y)*m.Eweight),np.sum(abs(y)*m.Nweight),np.sum(abs(gb[:,0])*m.eu),np.sum(abs(gb[:,1])*m.nu)],LD)
            if k in [0,8,15] and not seen.get(k):
                a,b,_=full_pair(t,mechanical)
                err=float(max(np.max(abs(a-difference)/np.maximum(b,1e-290)),np.max(abs(b-norm)/np.maximum(b,1e-290))))
                assert err<1e-12,(k,err);source_checks.append(dict(steps=n,k=k,relative=err));seen[k]=True
            return difference,norm,(ma.copy(),mb.copy())
        h=m.t[-1]/n;pair(gamma*h);''')
    namespace=dict(run.audit_owner.block.__globals__,OUT=OUT,
        source_checks=source_checks,run=SimpleNamespace(initialize=initialize,coupled=run.reuse.coupled,paths=run.paths,AMP=run.reuse.AMP,C=run.C))
    exec(compile(source,__file__,'exec'),namespace);namespace['block'](sweep)
    result=read(folder/'block-result.json');result.pop('original_stop_preserved');result.pop('original_waveform_change_test_passed')
    before=read(OUT/f'sweep-{sweep-1}/block-result.json')['maximum_block_defect']
    assert len(checks)==12,checks
    if source_only:assert len(source_checks)==6,source_checks
    result.update(criteria_preregistered=True,zero_collision_exact_checks=checks,
        source_difference_checks=source_checks,
        next_sweep_allowed=not result['passed'] and sweep<3 and result['maximum_block_defect']<before)
    write(folder/'block-result.json',result)


def block_sources(sweep):block_retry(sweep,True)


def block_finish(sweep):
    folder=OUT/f'sweep-{sweep}';plan=folder/'block-finish-plan.json';assert not plan.exists()
    failures=[OUT/f'block-{sweep}-receipt.json',OUT/f'audit-block_retry-{sweep}-receipt.json',OUT/f'audit-block_sources-{sweep}-receipt.json']
    assert all('block audit forecast' in read(p)['error'] for p in failures)
    files=[Path(__file__),Path(run.__file__),Path(run.audit_owner.__file__),folder/'residual.json']+failures
    files += [p/f'steps-{n}-reference-128.npz' for s in [sweep-1,sweep] for p in run.paths(s) for n in [64,128]]
    write(plan,dict(classification='Counterexample candidate',sweep=sweep,
        scope='Unchanged original full lagged-input, identical-operator and all384-stage source audit. No optimized source-difference path, no new evolution or threshold change.',
        budget_reassessment='Original300s block admission stopped at313.83s forecast. Two attempts to remove zero/duplicate local work also stopped on forecast. Preserve all failures and use the original exact evaluator with a400s read-only action cap; the original4200s total remains binding. This reallocates total budget, not a claim of meeting the original sub-budget.',
        budget_seconds=400,original_block_budget_seconds=300,total_action_seconds=run.TOTAL,
        recorded_prior_admission_seconds=sum(read(p)['seconds'] for p in failures),
        decision='Accept the finite coupling only if every original input/source/E-H defect is below0.002; otherwise preserve its actual failure and original maximum-three/nondecrease rule.',
        bindings={str(p):sha(p) for p in files}))
    source=inspect.getsource(run.audit_owner.block)
    source=source.replace("OUT/'block-plan.json'",repr(str(plan)))
    source=source.replace("OUT/'block-result.json'","OUT/f'sweep-{sweep}/block-result.json'")
    source=source.replace("assert max(list(inputs.values())+source+paired)<.002,rows[-1]",'# Preserve failed finite-block verdict.')
    source=source.replace('passed=True,rows=rows,','passed=maximum<.002,rows=rows,').replace('assert forecast<300,','assert forecast<400,')
    namespace=dict(run.audit_owner.block.__globals__,OUT=OUT,
        run=SimpleNamespace(initialize=run.initialize,coupled=run.reuse.coupled,paths=run.paths,AMP=run.reuse.AMP,C=run.C))
    exec(compile(source,__file__,'exec'),namespace);namespace['block'](sweep)
    result=read(folder/'block-result.json');result.pop('original_stop_preserved');result.pop('original_waveform_change_test_passed')
    before=read(OUT/f'sweep-{sweep-1}/block-result.json')['maximum_block_defect']
    result.update(criteria_preregistered=True,budget_reassessed=True,
        next_sweep_allowed=not result['passed'] and sweep<3 and result['maximum_block_defect']<before)
    write(folder/'block-result.json',result)


def prepare(sweep):
    folder=OUT/f'sweep-{sweep}';photon,material=run.paths(sweep)
    files=[Path(__file__),Path(run.__file__),Path(run.__file__).with_name('resolve_native_mixed_return.py'),
        OUT/'result.json',OUT/'metric/result.json',OUT/'normalization.json',folder/'block-result.json',
        OUT/'corrected-fields/field-128-g8.npz',run.mixed.OUT/'applied-charge.npz',run.inf.BACKGROUND]
    files += [p/f'steps-{n}-reference-128.npz' for s in [sweep-1,sweep] for p in run.paths(s) for n in [64,128]]
    files += [p/f'pilot-{n}.npz' for p in [photon,material] for n in [64,128]]
    files += [folder/f'gr/source-{n}-reference-128.npz' for n in [64,128]]
    files += [material/f'stress-{n}-reference-128.npz' for n in [64,128]]
    assert not (OUT/'audit-plan.json').exists()
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',sweep=sweep,
        claim='Verify conserved source export, signed angular ports, mechanical/collision partition and reused prefixes of the actual accepted mixed-field return. Preserve the corrected driving charge separately from its much smaller transport correction.',
        directional='Reuse the existing three saved-state material additivity and positive-homogeneity comparison against the original incident direction; no new evolution.',
        charge='For compact sources and zero additional initial data, the free retarded field at x=0,T equals its outgoing value at x=R,T+R/c. This applies also to the compact potential source already propagated. It is not full dynamical exterior closure.',
        gates=dict(source=1e-12,port=1e-12,partition=1e-8,balance=1e-8,directional=.002,old_observer=1e-10),
        budget_seconds=CAP,total_action_seconds=run.TOTAL,new_evolution_steps=0,
        scope='Retained finite model and sampled directions only. Additional exterior mixed source, physical mixed ADM normalization, full EOS/derivative/nonlinear error and static/observational comparison remain open.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def audit(sweep):
    plan=read(OUT/'audit-plan.json');assert plan['sweep']==sweep
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    assert read(OUT/'result.json')['passed'];assert read(OUT/f'sweep-{sweep}/block-result.json')['passed']
    photon,material=run.paths(sweep);previous_photon,previous_material=run.paths(sweep-1);rows=[]
    for n in [64,128]:
        p=np.load(photon/f'steps-{n}-reference-128.npz');m=np.load(material/f'steps-{n}-reference-128.npz')
        pp=np.load(previous_photon/f'steps-{n}-reference-128.npz');mm=np.load(previous_material/f'steps-{n}-reference-128.npz')
        d=np.load(OUT/f'sweep-{sweep}/gr/source-{n}-reference-128.npz');s=np.load(material/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(run.C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-s[:,0],d['nonrest_trace_erg']+rest-(s[:,0]-s[:,1]-2*s[:,3]),
            d['nonrest_stress_erg']+rest-(s[:,0]-s[:,1]),d['pressure_volume_erg']-s[:,3],
            d['metric_stress_erg']-(total+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        source=float(max(np.max(abs(v)) for v in errors)/max(np.max(abs(s)),1e-290))
        ids=[int(np.argmin(abs(mm['t']-t))) for t in p['t']]
        actual=p['moments'][:,[1,2]];old=mm['history_scaled'][ids][:,[2,3]]*run.reuse.AMP
        a=actual-p['collision_transfer'].transpose(0,2,1);b=old-pp['collision_transfer'].transpose(0,2,1)
        norm=np.maximum(np.max(np.sum(abs(actual),axis=2),axis=0),np.max(np.sum(abs(old),axis=2),axis=0))
        partition=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(norm,1e-290)).astype(float).tolist()
        _,port=run.inf.emitted(photon,n,SimpleNamespace(T=p['t'][-1]))
        balance=float(np.max(abs(np.sum(m['history_scaled'],axis=2,dtype=LD)+m['discards_scaled']-m['ledgers_scaled'])/np.maximum(m['norms_scaled'],1.)))
        pilot=np.load(photon/f'pilot-{n}.npz')
        for key in ['t','moments','material_history','collision_transfer','radial_ports','accepted_angular_times','accepted_angular_luminosity']:
            assert np.array_equal(p[key][:len(pilot[key])],pilot[key]),key
        pilot=np.load(material/f'pilot-{n}.npz')
        for key in ['t','history_scaled']:
            assert np.array_equal(m[key][:len(pilot[key])],pilot[key]),key
        assert source<1e-12 and port<1e-12 and balance<1e-8 and max(partition)<1e-8
        rows.append(dict(steps=n,source_identity=source,port=port,material_balance=balance,
            mechanical_collision_partition=partition,photon_and_material_prefixes_reused=True))
    run.initialize(sweep);m=run.reuse.coupled.Material(128,128)
    own=np.load(material/'steps-128-reference-128.npz');original_folder=run.inf.BEFORE
    old=np.load(original_folder/'material-analytic/steps-128-reference-128.npz')
    current=m.transfer.copy();run.reuse.coupled.transfers(m,original_folder/'photons-precise',128);original=m.transfer.copy()
    primary=run.mixed.inc.Driver(8);saved=run.StageMetric(128);directions=[]
    class Blend:
        def __init__(self,a,b):self.a,self.b=a,b
        def view(self,t,clock,side='left'):
            g=primary.view(t,clock,side);h=saved.view(t,clock,side)
            return {k:self.a*g.get(k,0.)+self.b*h.get(k,0.) for k in set(g)|set(h)}
    def rhs(t,z,a,b):
        m.driver=Blend(a,b);m.transfer=a*original+b*current
        return m.rhs(t,z)[0]
    for k in [1,8,16]:
        t=m.t[k];i=int(np.argmin(abs(own['t']-t)));j=int(np.argmin(abs(old['t']-t)))
        z1=own['history_scaled'][i];z0=old['history_scaled'][j]
        q=m.point(k)['Q'];units=np.maximum(abs(q),1.);units[1]=np.maximum(q[0]*run.C**2,1.)
        scale=float(np.max(abs(z1)/units)/max(np.max(abs(z0)/units),1e-290))
        r0=rhs(t,z0*scale,scale,0.);r1=rhs(t,z1,0.,1.)
        added=rhs(t,z0*scale+z1,scale,1.);half=rhs(t,z1*.5,0.,.5)
        add=(np.sum(abs(added-r0-r1),axis=1)/np.maximum(np.sum(abs(r0),axis=1)+np.sum(abs(r1),axis=1),1.)).astype(float).tolist()
        hom=(np.sum(abs(2*half-r1),axis=1)/np.maximum(np.sum(abs(r1),axis=1),1.)).astype(float).tolist()
        assert max(add+hom)<.002,(k,add,hom)
        directions.append(dict(k=k,balanced_original_scale=scale,additivity=add,positive_homogeneity=hom))
    # The previous independent exterior observers provide a check of the
    # boundary readout identity before using the corrected stored boundary.
    bg=np.load(run.inf.BACKGROUND);den=1-bg['epsilon'].astype(LD)
    original=np.load(run.mixed.OUT/'field-128-g8.npz');corrected=np.load(OUT/'corrected-fields/field-128-g8.npz')
    source=np.load(run.mixed.previous.FIELDS/'source-128.npz');M=LD(source['M_cm'])
    before=np.load(run.mixed.OUT/'applied-charge.npz');old_charge=-original['U'][::2,-1].astype(LD)/M/den
    reproduction=float(np.max(abs(old_charge-before['reciprocal_compact_return']))/np.max(abs(old_charge)));assert reproduction<1e-10
    actual=-corrected['U'][::2,-1].astype(LD)/M/den;transport=np.load(OUT/'return-128-a8-r8.npz')['charge']
    np.savez_compressed(OUT/'charge-parts.npz',t=bg['t'],previous_selected=before['previous_selected'],
        original_mixed=before['reciprocal_compact_return'],corrected_mixed=actual,transport_return=transport,
        selected_without_tiny_return=before['previous_selected']+actual,original_body=before['old_body'],direct=before['direct'])
    R,y,T,c=sp.symbols('R y T c',positive=True)
    assert sp.simplify((T+R/c)-(R+y)/c-(T-y/c))==0
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,directional= directions,
        old_independent_observer_reproduction=reproduction,corrected_mixed_endpoint=float(actual[-1]),
        physical_transport_return_endpoint=float(transport[-1]),transport_over_corrected_mixed=float(transport[-1]/actual[-1]),
        selected_endpoint_without_tiny_return=float(before['previous_selected'][-1]+actual[-1]),
        normalization=read(OUT/'normalization.json'),
        inherited_readout_metadata='The compact owner incident_amplitude and endpoint_over_incident_amplitude describe its original primary driver, not this mixed metric. Only the explicit current metric normalization and the separately divided physical return are used here.',
        symbolic=dict(classification='Proven',passed=True,identity='For source at x=-y<=0, (T+R/c)-(R+y)/c=T-y/c. The corrected compact retarded boundary field is its leading outgoing null readout.'),
        actual_mixed_field_returned=True,full_goal_complete=False,scope=plan['scope'])
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]);assert action in ['prepare','audit','block_retry','block_sources','block_finish']
    if action=='block_retry':CAP=285
    if action=='block_sources':CAP=271
    if action=='block_finish':CAP=400
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.mixed.inc.native.deadline(CAP)
    receipt=OUT/f'audit-{action}-{sweep}-receipt.json';assert not receipt.exists()
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        assert sum(read(p)['seconds'] for p in OUT.rglob('*-receipt.json'))+CAP<=run.TOTAL
        globals()[action](sweep)
    except Exception as exc:error=repr(exc);raise
    finally:write(receipt,dict(action='audit_'+action,seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
