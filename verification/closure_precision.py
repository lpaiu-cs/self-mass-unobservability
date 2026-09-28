"""Request30: close observed numerical prerequisites without changing verdicts."""
from pathlib import Path
import json, math, shutil, sys
import numpy as np
import common_eos as c
import native_eos_bridge as shim

ROOT=c.ROOT
OLD=c.OUT
OUT=ROOT/'outputs/closure-precision30'
CACHE=Path('/home/lpaiu/work/closure-precision30')


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='df49775',
        previous_manifest_sha256=c.sha(OLD/'manifest.json'),
        full_goal='Resolve the remaining physical EOS, GR evolution, scalar driving and full observational inference requirements; numerical prerequisite checks do not redefine completion.',
        steps=['Determine cause of independent EOS inverse reproducibility failure',
            'Attribute native O18 and heat temperature derivative discrepancy to individual reactions',
            'Correct a demonstrated root cause and rerun independent unchanged controls',
            'Advance conservative composition precision and transport closure using corrected prerequisites'],
        unchanged_controls=dict(native_vector_derivative_relative=1e-3,native_heat_derivative_relative=1e-3,
            composition_absolute=1e-16,composition_relative_to_total_change=1e-3,
            temperature_log_absolute=2e-6,global_energy_relative_to_release=1e-6),
        policy='Do not reinterpret earlier failures as passes, silently change an EOS or reaction model, or treat source signatures as proof of the preserved binary. '
            'Keep physical error, deterministic numerical evaluation, finite refinement and rigorous interval certification distinct.'))
    bindings={}
    for rel in ['docs/'+n+'.md' for n in c.cell.DOCS]+['paper/revision-manifest.json']:
        source=ROOT/rel;dest=OUT/'previous-notes'/source.name;dest.parent.mkdir(exist_ok=True)
        shutil.copy2(source,dest)
        bindings[rel]=dict(snapshot=dest.relative_to(ROOT).as_posix(),sha256=c.sha(source),historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)


def repeatability():
    initial=dict(np.load(OLD/'restored-state-17-4.npz'))
    final=dict(np.load(OLD/'P-evolve-4-0-source.npz'))
    scores=np.load(OLD/'P-evolve-4-0-independent-energy.npz')['scores']
    selected=np.unique(np.r_[np.argsort(scores)[-8:],0,2097,c.cell.ZONE]);eos=c.EOS();rows=[];raw=[]
    for i in selected:
        lp=initial['logP'][i];t0=initial['lnT'][i];t1=final['lnT'][i];x0=initial['X'][i];x1=final['X'][i]
        def a(): return eos(1,lp,t0,x0)
        def b(): return eos(1,lp,t1,x1)
        values=[a(),a(),b(),b()]
        for k in range(26):
            if k==c.NAMES.index('he4'): continue
            for h in [1e-5,5e-6]:
                xx=x0.copy();xx[k]+=h;xx[c.NAMES.index('he4')]-=h
                eos(1,lp,t0,xx)
        values.extend([b(),b(),a(),a()])
        eos(2,initial['lnd'][0],initial['lnT'][0],initial['X'][0]);values.extend([a(),b()])
        values=np.array(values);raw.append(values)
        u=values[:,2]+values[:,1]/values[:,0]
        rows.append(dict(cell=int(i),old_inverse_score=float(scores[i]),rho=float(np.exp(initial['lnd'][i])),
            T=float(np.exp(t0)),enthalpies=u.tolist(),
            initial_enthalpy_span_erg_g=float(np.ptp(u[[0,1,6,7,8]])),
            endpoint_enthalpy_span_erg_g=float(np.ptp(u[[2,3,4,5,9]])),
            max_endpoint_state_relative_span=float(np.max(np.ptp(values[[2,3,4,5,9]],axis=0)/np.maximum(1,abs(values[2]))))))
    np.savez_compressed(OUT/'EOS-repeatability.npz',cells=selected,outputs=np.array(raw))
    save('EOS-repeatability.json',dict(classification='Counterexample candidate',rows=rows,
        sequence=['initial','initial','endpoint','endpoint','after composition gradient: endpoint','endpoint','initial','initial','after remote state: initial','endpoint'],
        library_sha256=c.sha(c.gr.CACHE/'build/src/libfree_eos.so.1.0.0'),
        scope='Same physical inputs, varying prior call history. This tests the numerical evaluator, not physical EOS accuracy.'))
    print('EOS REPEATABILITY',rows,flush=True)


def context():
    c.native.OUT=OUT;c.native.CACHE=CACHE;c.native.context()
    shim.OUT=OUT;shim.CACHE=CACHE
    # shim.save is an imported function whose globals are in common_eos.
    shim.save=save


def flux_rounding_plan():
    save('flux-rounding-plan.json',dict(classification='Counterexample candidate',
        observation='EOS call-history tests show at most hundredths of erg/g in the worst outer inverse-score cells, far below the earlier residual. The independent audit recomputed redshifted flux using exp(2 nu), while the producer used exp(nu) squared.',
        hypothesis='Subtracting nearly equal luminosities magnifies the mathematically equivalent binary64 redshift evaluations into local energy differences.',
        test='Reconstruct both face-flux arrays and both N squared denominators. Compare the predicted difference in local integrated divergence with the saved inverse residual. Re-evaluate all ten source endpoints using the exact registered producer forcing.',
        policy='Retain both original audit and producer. Do not change energy thresholds. A diagnosed audit forcing mismatch corrects its EOS attribution, not the raw measured residual.'))


def flux_rounding():
    assert (OUT/'flux-rounding-plan.json').exists()
    original=dict(np.load(OLD/'restored-state-17-4.npz'));Nface=np.exp(original['nu_faces'])
    canonical=np.r_[original['L']*Nface[:-1]**2,0.]
    alternate=np.r_[original['L']*np.exp(2*original['nu_faces'][:-1]),0.]
    np.savez_compressed(OUT/'registered-face-flux.npz',canonical=canonical,alternate=alternate,
        initial_state_sha256=c.sha(OLD/'restored-state-17-4.npz'))
    eos=c.EOS();rest=(c.W/c.A-1)*(c.gr.C*100)**2;rows=[]
    for path in sorted(OLD.glob('*evolve-*-source.npz')):
        label=path.name.removesuffix('-source.npz');end=dict(np.load(path))
        report=json.loads((OLD/(label+'-source.json')).read_text());iso=report.get('isobaric',False)
        count,step=map(int,label.split('-')[-2:]);prefix='P-' if iso else ''
        start=original if step==0 else dict(np.load(OLD/f'{prefix}evolve-{count}-{step-1}-state-17-4.npz'))
        N=np.exp(start['nu']);h=report['dt']*N
        q0=np.diff(-canonical)/(start['dm']*N*N)
        q1=np.diff(-alternate)/(start['dm']*np.exp(2*start['nu']))
        old=np.load(OLD/(label+'-independent-energy.npz'));prediction=h*(q1-q0)
        values=[]
        for i in range(len(h)):
            coordinate=start['logP'][i] if iso else start['lnd'][i];mode=1 if iso else 2
            a=eos(mode,coordinate,start['lnT'][i],start['X'][i]);b=eos(mode,coordinate,end['lnT'][i],end['X'][i])
            e0=a[2]+(a[1]/a[0] if iso else 0);e1=b[2]+(b[1]/b[0] if iso else 0)
            target=e0-rest@(end['X'][i]-start['X'][i])-end['loss_per_baryon_gram'][i]-h[i]*q0[i]
            defect=e1-target;budget=max(2.,32*np.spacing(abs(target)),abs(target-e0)*1e-8)
            values.append([defect,defect/budget,prediction[i],old['defects'][i]-defect-prediction[i]])
        values=np.array(values);np.savez_compressed(OUT/(label+'-canonical-inverse.npz'),rows=values)
        worst=int(np.argmax(abs(values[:,1])))
        rows.append(dict(label=label,max_score=float(max(abs(values[:,1]))),worst_cell=worst,
            local_gate_passed=bool(max(abs(values[:,1]))<=1),
            prior_worst_cell=int(np.argmax(old['scores'])),
            predicted_flux_defect_max_erg_g=float(max(abs(prediction))),
            residual_prediction_error_max_erg_g=float(max(abs(values[:,3])))))
        save('flux-rounding.json',dict(classification='Counterexample candidate',rows=rows,
            prior_audit_preserved=True,physical_EOS_certified=False))
        print('CANONICAL INVERSE',rows[-1],flush=True)


class ColdEOS(c.EOS):
    """Same physics with a declared composition/degeneracy seed per query."""
    def __call__(self,mode,coordinate,T,x):
        seed=np.zeros(26)
        seed[c.NAMES.index('c12') if x[c.NAMES.index('c12')]<.5 else c.NAMES.index('h1')]=1
        # kif=0 fixes the seed degeneracy independently of the previous call.
        # A zero/nonzero elemental abundance change requests a cold start in
        # the unchanged FreeEOS source. Actual determinism is tested below.
        super().__call__(0,-20.,math.log(1e6),seed)
        return super().__call__(mode,coordinate,T,x)


def cold_plan():
    save('cold-EOS-plan.json',dict(classification='Counterexample candidate',
        method='Before each actual EOS query, evaluate kif=0 with ln f=-20 and T=1e6 K at pure C12 (or H1 if the requested C12 mass fraction is >=0.5). '
            'This sets a fixed degeneracy seed and changes zero/nonzero elemental abundances, triggering the source-defined cold-start path. '
            'The target EOS physics, flags and library remain unchanged.',
        gates='Test output equality after different prior calls. For energy inversion keep the original absolute/32-ulp/local-change tolerance; search actual representable lnT neighbors rather than widening that tolerance.',
        limits='Deterministic evaluation is neither physical EOS validation nor a rigorous continuum error enclosure.'))


def cold_test():
    assert (OUT/'cold-EOS-plan.json').exists()
    base=dict(np.load(OLD/'restored-state-17-4.npz'));end=dict(np.load(OLD/'P-evolve-4-0-source.npz'))
    cells=np.unique(np.r_[np.linspace(0,len(base['dm'])-1,40).astype(int),110,2097,2591,3023,3887,4760])
    eos=ColdEOS();warm=c.EOS();rows=[]
    for i in cells:
        values=[]
        for trial in range(3):
            k=int(cells[(int(i)+trial)%len(cells)])
            warm(2,base['lnd'][k],base['lnT'][k],base['X'][k])
            values.append(eos(1,base['logP'][i],end['lnT'][i],end['X'][i]))
        values=np.array(values)
        other=warm(1,base['logP'][i],end['lnT'][i],end['X'][i])
        rows.append(dict(cell=int(i),bitwise_repeatable=bool(np.array_equal(values[0],values[1]) and np.array_equal(values[0],values[2])),
            maximum_relative_spread=float(np.max(np.ptp(values,axis=0)/np.maximum(1,abs(values[0])))),
            warm_relative_difference=float(np.max(abs(other-values[0])/np.maximum(1,abs(values[0]))))))
    save('cold-EOS-test.json',dict(classification='Counterexample candidate',rows=rows,
        all_repeatable=all(r['bitwise_repeatable'] for r in rows),physical_EOS_certified=False))
    print('COLD EOS',len(rows),'states; repeatable',all(r['bitwise_repeatable'] for r in rows),
        'max spread',max(r['maximum_relative_spread'] for r in rows),flush=True)


def cold_inverse():
    assert (OUT/'cold-EOS-plan.json').exists()
    original=dict(np.load(OLD/'restored-state-17-4.npz'))
    flux=np.load(OUT/'registered-face-flux.npz')['canonical'];eos=ColdEOS()
    rest=(c.W/c.A-1)*(c.gr.C*100)**2;rows=[]
    for path in sorted(OLD.glob('*evolve-*-source.npz')):
        label=path.name.removesuffix('-source.npz');end=dict(np.load(path))
        report=json.loads((OLD/(label+'-source.json')).read_text());iso=report.get('isobaric',False)
        count,step=map(int,label.split('-')[-2:]);prefix='P-' if iso else ''
        start=original if step==0 else dict(np.load(OLD/f'{prefix}evolve-{count}-{step-1}-state-17-4.npz'))
        N=np.exp(start['nu']);h=report['dt']*N;q=np.diff(-flux)/(start['dm']*N*N);results=[]
        for i in range(len(h)):
            mode=1 if iso else 2;coordinate=start['logP'][i] if iso else start['lnd'][i]
            a=eos(mode,coordinate,start['lnT'][i],start['X'][i]);e0=a[2]+(a[1]/a[0] if iso else 0)
            target=e0-rest@(end['X'][i]-start['X'][i])-end['loss_per_baryon_gram'][i]-h[i]*q[i]
            budget=max(2.,32*np.spacing(abs(target)),abs(target-e0)*1e-8)
            def sample(t):
                b=eos(mode,coordinate,t,end['X'][i])
                return b[2]+(b[1]/b[0] if iso else 0)-target,b[10]-(b[1]/b[0]*b[8] if iso else 0)
            t=float(end['lnT'][i]);best=(float('inf'),t)
            for _ in range(10):
                defect,capacity=sample(t)
                if abs(defect)<abs(best[0]): best=(defect,t)
                if abs(defect)<budget*.25: break
                proposed=t-defect/capacity
                if proposed==t: break
                t=proposed
            if abs(best[0])>budget*.25:
                lo=hi=best[1]
                for _ in range(3):
                    lo=np.nextafter(lo,-np.inf);hi=np.nextafter(hi,np.inf)
                    for t in [lo,hi]:
                        defect,_=sample(t)
                        if abs(defect)<abs(best[0]): best=(defect,t)
            # Independent repeated evaluation follows a deliberately different
            # physical EOS query. The registered cold seed must erase history.
            c.EOS()(2,original['lnd'][0],original['lnT'][0],original['X'][0])
            checked,_=sample(best[1]);assert checked==best[0],(label,i,checked,best)
            results.append([best[1],checked,budget,checked/budget,best[1]-end['lnT'][i]])
        results=np.array(results);np.savez_compressed(OUT/(label+'-cold-inverse.npz'),rows=results,
            source_sha256=c.sha(path),forcing_sha256=c.sha(OUT/'registered-face-flux.npz'))
        worst=int(np.argmax(abs(results[:,3])))
        rows.append(dict(label=label,cells=len(results),max_score=float(max(abs(results[:,3]))),worst_cell=worst,
            passed=bool(max(abs(results[:,3]))<=1),max_temperature_correction=float(max(abs(results[:,4]))),
            repeated_endpoint_bitwise_equal=True))
        save('cold-inverse.json',dict(classification='Counterexample candidate',rows=rows,
            full_trajectory_recomputed=False,physical_EOS_certified=False,
            scope='Re-solved saved source endpoint energy equations with the original reaction increments/losses and canonical forcing; no new stellar trajectory.'))
        print('COLD INVERSE',rows[-1],flush=True)


def failed_cell():
    i=3587;initial=dict(np.load(OLD/'restored-state-17-4.npz'))
    data={k:v[i:i+1].copy() for k,v in initial.items() if v.ndim and len(v)==len(initial['dm'])}
    d={k:v[i:i+1].copy() for k,v in np.load(OLD/'evolution-initial-native.npz').items()}
    _,profile=c.mesa(OLD/'evolution-initial-profile.data.gz');p={k:v[i:i+1] for k,v in profile.items()}
    flux=np.load(OUT/'registered-face-flux.npz')['canonical'][i:i+2]
    oldout,oldeos=c.OUT,c.EOS
    try:
        c.OUT=OUT;c.EOS=ColdEOS
        # Original EOS class reads these declared map constants via c.OUT.
        for name in ['atomic-binding.json','freeeos-weights.json']: shutil.copy2(OLD/name,OUT/name)
        try:
            result,loss=c.source_step(data,d,p,140780.16/4,flux,'previously-failed-cell',False)
            outcome=dict(passed=True,max_delta_lnT=float(max(abs(result['lnT']-data['lnT']))))
        except (ValueError,AssertionError) as error:
            outcome=dict(passed=False,error=repr(error))
    finally: c.OUT,c.EOS=oldout,oldeos
    save('failed-cell-replay.json',dict(classification='Counterexample candidate',cell=i,steps=4,step_zero_based=0,
        original_error='energy inverse residual +19 erg/g',outcome=outcome,
        scope='Actual original local coupled source equation with a seeded EOS, same inputs and unchanged stopping criterion; no claim that the full four-step GR evolution ran.'))
    print('FAILED CELL REPLAY',outcome,flush=True)


def cold_full_audit():
    for name in ['atomic-binding.json','freeeos-weights.json','restored-state-17-4.npz']: shutil.copy2(OLD/name,OUT/name)
    oldout,oldeos=c.OUT,c.EOS
    try:
        c.OUT=OUT;c.EOS=ColdEOS;c.full_audit()
    finally: c.OUT,c.EOS=oldout,oldeos


def channel_plan():
    channels=['r_f18_wk_o18','r_o18_pa_n15','r_o18_pg_f19','r_n15_ap_o18','r_f19_gp_o18']
    save('channel-plan.json',dict(classification='Counterexample candidate',channels=channels,
        method='Retain the entire explicit 26-species network and set each named rate factor to 0.5 in turn. '
            'Use the same frozen four native EOS auxiliary inputs at T and T exp(plus/minus 5e-5). '
            'Twice the full-minus-halved difference isolates a rate direction while preserving the full network context.',
        checks='Reconstruct the O18 vector component and its analytic/finite-difference temperature derivative. '
            'Do not assume that channel heat or weak Q is additive until checked.',
        h=5e-5,reference_manifest_sha256=c.sha(OLD/'manifest.json')))


def channels():
    plan=json.loads((OUT/'channel-plan.json').read_text());context()
    initial=dict(np.load(OLD/'restored-state-17-4.npz'));aux=np.load(OLD/'EOS-common-replacement.npy');h=plan['h']
    for n,reaction in enumerate(plan['channels']):
        for sign in [0,-1,1]:
            label=f'channel-{n}-{sign}'
            if (OUT/(label+'-native.npz')).exists(): continue
            data={k:v.copy() for k,v in initial.items()};data['lnT']+=sign*h
            c.native.setup(label,data,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
            folder=CACHE/label;path=folder/'inlist1';text=path.read_text()
            text=text.replace('&star_job',"&star_job\n num_special_rate_factors=1\n reaction_for_special_factor(1)='"+reaction+"'\n special_rate_factor(1)=0.5\n",1)
            path.write_text(text);shutil.copy2(path,OUT/'inputs'/label/path.name)
            c.native.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:c.sha(p) for p in folder.iterdir() if p.is_file()}))
            np.save(OUT/(label+'-replacement.npy'),aux);shim.trace(label,aux)
    channel_analysis()


def channel_analysis():
    plan=json.loads((OUT/'channel-plan.json').read_text());h=plan['h'];j=c.NAMES.index('o18')
    refs=[dict(np.load(OLD/(name+'-native.npz'))) for name in ['EOS-common','EOS-attribution-frozen--1','EOS-attribution-frozen-1']]
    active=abs(refs[0]['heat'])>1;rows=[];totals=np.zeros((3,len(active)));arrays={}
    for n,name in enumerate(plan['channels']):
        if not all((OUT/f'channel-{n}-{sign}-native.npz').exists() for sign in [0,-1,1]): break
        d=[dict(np.load(OUT/f'channel-{n}-{sign}-native.npz')) for sign in [0,-1,1]]
        f=2*(refs[0]['dxdt'][:,j]-d[0]['dxdt'][:,j])
        dt=2*(refs[0]['dxdt_T'][:,j]-d[0]['dxdt_T'][:,j])*refs[0]['T']
        fd=((refs[2]['dxdt'][:,j]-d[2]['dxdt'][:,j])-(refs[1]['dxdt'][:,j]-d[1]['dxdt'][:,j]))/h
        totals+=np.array([f,dt,fd]);arrays[name]=np.array([f,dt,fd])
        scale=np.maximum(1e-30,np.maximum(abs(refs[0]['dxdt'][:,j]),abs(refs[0]['dxdt_T'][:,j]*refs[0]['T'])))
        score=abs(fd-dt)/scale;index=int(np.argmax(np.where(active,score,-1)))
        rows.append(dict(reaction=name,worst_active_cell=index,max_active_score=float(score[index]),
            at_2097=dict(f=float(f[2097]),reported_dlnT=float(dt[2097]),finite_dlnT=float(fd[2097]))))
    full=np.array([refs[0]['dxdt'][:,j],refs[0]['dxdt_T'][:,j]*refs[0]['T'],(refs[2]['dxdt'][:,j]-refs[1]['dxdt'][:,j])/(2*h)])
    arrays['full']=full;arrays['sum']=totals
    np.savez_compressed(OUT/'channel-analysis.npz',**arrays)
    save('channel-analysis.json',dict(classification='Counterexample candidate',rows=rows,
        complete=len(rows)==len(plan['channels']),
        at_2097_sum=totals[:,2097].tolist(),at_2097_full=full[:,2097].tolist(),
        scope='Actual native reaction-factor interventions; attribution precedes any correction or certification.'))
    print('CHANNEL ANALYSIS',rows,flush=True)


def weak_control_plan():
    save('weak-control-plan.json',dict(classification='Counterexample candidate',
        reason='The F18 named rate-factor intervention has no O18 effect; the archived network uses a separate weaklib path and global weak_rate_factor. Do not interpret the ineffective factor as zero physical F18 flux.',
        method='Set controls weak_rate_factor=0.5 for the full network at T and T exp(plus/minus 5e-5), retaining the same four frozen EOS inputs. '
            'Compare twice the full-minus-half difference of f, f_T, heat and heat_T with their temperature differences.',
        h=5e-5,original_failed_named_factor_preserved=True))


def weak_control():
    assert (OUT/'weak-control-plan.json').exists();context()
    initial=dict(np.load(OLD/'restored-state-17-4.npz'));aux=np.load(OLD/'EOS-common-replacement.npy')
    for sign in [0,-1,1]:
        label=f'weak-half-{sign}'
        if (OUT/(label+'-native.npz')).exists(): continue
        data={k:v.copy() for k,v in initial.items()};data['lnT']+=sign*5e-5
        c.native.setup(label,data,extra='weak_rate_factor=0.5',species=c.NAMES,
            network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        np.save(OUT/(label+'-replacement.npy'),aux);shim.trace(label,aux)
    d=[dict(np.load(OUT/f'weak-half-{sign}-native.npz')) for sign in [0,-1,1]]
    r=[dict(np.load(OLD/(name+'-native.npz'))) for name in ['EOS-common','EOS-attribution-frozen--1','EOS-attribution-frozen-1']]
    active=abs(r[0]['heat'])>1;rows=[];arrays={}
    for key in ['dxdt','heat']:
        factor=r[0]['T'][:,None] if key=='dxdt' else r[0]['T']
        fullfd=(r[2][key]-r[1][key])/1e-4;fullanalytic=r[0][key+'_T']*factor
        weakfd=2*((r[2][key]-d[2][key])-(r[1][key]-d[1][key]))/1e-4
        weakanalytic=2*(r[0][key+'_T']-d[0][key+'_T'])*factor
        residual=(fullfd-fullanalytic)-(weakfd-weakanalytic)
        scale=np.maximum(1e-30,np.maximum(abs(fullanalytic).max(1),abs(r[0][key]).max(1))) if key=='dxdt' else np.maximum(1,abs(fullanalytic))
        def score(x): return float(max(((abs(x).max(1) if key=='dxdt' else abs(x))/scale)[active]))
        rows.append(dict(output=key,full_discrepancy=score(fullfd-fullanalytic),
            remaining_after_weak_attribution=score(residual)))
        arrays.update({key+'_fullfd':fullfd,key+'_fullanalytic':fullanalytic,
            key+'_weakfd':weakfd,key+'_weakanalytic':weakanalytic,key+'_remaining':residual})
    np.savez_compressed(OUT/'weak-control.npz',**arrays)
    save('weak-control.json',dict(classification='Counterexample candidate',rows=rows,
        scope='Global weak-path attribution, with the unrelated named-factor null preserved. Attribution is not yet a corrected analytic derivative.'))
    print('WEAK CONTROL',rows,flush=True)


def blend_plan():
    save('blend-plan.json',dict(classification='Counterexample candidate',
        source_finding='The archived eval_weak.f90 blends lambda=a lambda_table+(1-a) lambda_low but returns a times the table temperature derivative, omitting da/dlnT times (lambda_table-lambda_low).',
        registered_default=dict(T9_off=.01,T9_on=.02,high_Z_threshold=26),
        controls=dict(pure_low=[2.,3.],pure_table=[1e-7,2e-7]),
        method='Run both endpoints at the identical full-star rho,T,X and frozen EOS inputs using existing star_job settings. '
            'Verify affine reconstruction of the actual native full source before adding a_prime times endpoint difference to the reported derivative.',
        gates='Use the original all-burning-cell 1e-3 vector and heat derivative criteria, against frozen-auxiliary and fresh-auxiliary FD controls. Retain every failure and identify any stencil crossing a blend boundary.',
        limits='A mathematical correction of the specified blended rate model does not certify physical weak rates, a full composition Jacobian, or an arbitrary continuous derivative domain.'))
    for rel in ['rates/private/eval_weak.f90','rates/public/rates_def_mic.f90','net/private/net_eval.f90',
        'star/defaults/star_job.defaults','star/job/run_star_support.f90']:
        dest=OUT/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(c.fresh.MESA/rel,dest)


def blend_endpoints():
    plan=json.loads((OUT/'blend-plan.json').read_text());context()
    initial=dict(np.load(OLD/'restored-state-17-4.npz'));aux=np.load(OLD/'EOS-common-replacement.npy')
    assert max(c.Z)<plan['registered_default']['high_Z_threshold']
    for kind,(off,on) in plan['controls'].items():
        label='blend-'+kind
        if (OUT/(label+'-native.npz')).exists(): continue
        c.native.setup(label,initial,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        folder=CACHE/label;path=folder/'inlist1'
        path.write_text(path.read_text().replace('&star_job',f'&star_job\n T9_weaklib_full_off={off:.17g}\n T9_weaklib_full_on={on:.17g}\n',1))
        shutil.copy2(path,OUT/'inputs'/label/path.name)
        c.native.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:c.sha(p) for p in folder.iterdir() if p.is_file()}))
        np.save(OUT/(label+'-replacement.npy'),aux);shim.trace(label,aux)
    blend_analysis()


def blend_analysis():
    plan=json.loads((OUT/'blend-plan.json').read_text())
    low=dict(np.load(OUT/'blend-pure_low-native.npz'));high=dict(np.load(OUT/'blend-pure_table-native.npz'))
    base=dict(np.load(OLD/'EOS-common-native.npz'));T9=base['T']*1e-9
    off,on=plan['registered_default']['T9_off'],plan['registered_default']['T9_on']
    alpha=np.clip((T9-off)/(on-off),0,1);prime=np.where((T9>off)&(T9<on),T9/(on-off),0.)
    active=abs(base['heat'])>1;rows=[];arrays={}
    frozen=[dict(np.load(OLD/(name+'-native.npz'))) for name in ['EOS-attribution-frozen--1','EOS-attribution-frozen-1']]
    chain=dict(np.load(OLD/'input-attribution.npz'))
    for key in ['dxdt','heat']:
        a=alpha[:,None] if key=='dxdt' else alpha;ap=prime[:,None] if key=='dxdt' else prime
        temp=base['T'][:,None] if key=='dxdt' else base['T']
        reconstructed=(1-a)*low[key]+a*high[key]
        correction=ap*(high[key]-low[key]);reported=temp*base[key+'_T'];fixed=reported+correction
        full=fixed+chain[key+'_Ye']+chain[key+'_eta']
        scale=np.maximum(1e-30,np.maximum(abs(reported).max(1),abs(base[key]).max(1))) if key=='dxdt' else np.maximum(1,abs(reported))
        def score(x): return float(max(((abs(x).max(1) if key=='dxdt' else abs(x))/scale)[active]))
        fd=(frozen[1][key]-frozen[0][key])/1e-4
        row=dict(output=key,affine_reconstruction_score=score(reconstructed-base[key]),
            original_frozen_derivative_score=score(fd-reported),corrected_frozen_derivative_score=score(fd-fixed),
            fresh_EOS=[])
        for h in [1e-4,5e-5]:
            minus=np.load(OLD/f'EOS-common-T-{h}--1-native.npz');plus=np.load(OLD/f'EOS-common-T-{h}-1-native.npz')
            error=score((plus[key]-minus[key])/(2*h)-full)
            crossing=((T9*np.exp(-h)<off)&(T9*np.exp(h)>off))|((T9*np.exp(-h)<on)&(T9*np.exp(h)>on))
            row['fresh_EOS'].append(dict(h=h,score=error,passed=error<1e-3,active_stencil_crossings=int(sum(crossing&active))))
        arrays.update({key+'_correction':correction,key+'_corrected_frozen':fixed,key+'_corrected_full':full})
        rows.append(row)
    np.savez_compressed(OUT/'blend-correction.npz',**arrays)
    save('blend-correction.json',dict(classification='Counterexample candidate',rows=rows,
        physical_weak_rates_certified=False,full_composition_Jacobian_certified=False,
        original_native_values_unchanged=True,original_failed_derivative_verdict_preserved=True))
    print('BLEND CORRECTION',rows,flush=True)


def direct_blend_plan():
    prior=OUT/'previous-notes/native_eos_bridge.py'
    assert c.sha(prior)==json.loads((OLD/'manifest.json').read_text())['sha256']['verification/native_eos_bridge.py']
    bindings=json.loads((OUT/'historical-note-bindings.json').read_text())
    bindings['verification/native_eos_bridge.py']=dict(snapshot=prior.relative_to(ROOT).as_posix(),sha256=c.sha(prior),code_revision=True)
    save('historical-note-bindings.json',bindings)
    save('direct-blend-plan.json',dict(classification='Counterexample candidate',
        reason='The star_job endpoint controls left the initial native source unchanged; their interpretation as endpoints is unverified. Keep that failed test and observe the actual module inputs at net_get entry.',
        controls=dict(pure_low=[2.,3.],pure_table=[1e-7,2e-7]),
        method='In our disposable child only, record and replace the two SHA-bound rates_def_mic blend inputs at each net_get entry, then restore them after return. Never overwrite returned rates or derivatives. '
            'Compare actual endpoint differences and affine reconstruction before correcting any derivative.',
        previous_configuration_only_results_preserved=True))


def direct_blend():
    plan=json.loads((OUT/'direct-blend-plan.json').read_text());context()
    initial=dict(np.load(OLD/'restored-state-17-4.npz'));aux=np.load(OLD/'EOS-common-replacement.npy')
    for kind,values in plan['controls'].items():
        label='direct-blend-'+kind
        if (OUT/(label+'-native.npz')).exists(): continue
        c.native.setup(label,initial,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        np.save(OUT/(label+'-replacement.npy'),aux);shim.trace(label,aux,weak_blend=values)
    direct_blend_analysis()


def direct_blend_analysis():
    low=dict(np.load(OUT/'direct-blend-pure_low-native.npz'));high=dict(np.load(OUT/'direct-blend-pure_table-native.npz'))
    base=dict(np.load(OLD/'EOS-common-native.npz'));T9=base['T']*1e-9
    assert np.array_equal(low['weak_blend_before'],high['weak_blend_before'])
    off,on=low['weak_blend_before'].T;assert np.all(on>off)
    alpha=np.clip((T9-off)/(on-off),0,1);prime=np.where((T9>off)&(T9<on),T9/(on-off),0.)
    active=abs(base['heat'])>1;rows=[];arrays={}
    frozen=[dict(np.load(OLD/(name+'-native.npz'))) for name in ['EOS-attribution-frozen--1','EOS-attribution-frozen-1']]
    chain=dict(np.load(OLD/'input-attribution.npz'))
    for key in ['dxdt','heat']:
        a=alpha[:,None] if key=='dxdt' else alpha;ap=prime[:,None] if key=='dxdt' else prime
        temp=base['T'][:,None] if key=='dxdt' else base['T']
        contrast=high[key]-low[key];reconstructed=(1-a)*low[key]+a*high[key]
        reported=temp*base[key+'_T'];fixed=reported+ap*contrast;full=fixed+chain[key+'_Ye']+chain[key+'_eta']
        scale=np.maximum(1e-30,np.maximum(abs(reported).max(1),abs(base[key]).max(1))) if key=='dxdt' else np.maximum(1,abs(reported))
        def score(x): return float(max(((abs(x).max(1) if key=='dxdt' else abs(x))/scale)[active]))
        fd=(frozen[1][key]-frozen[0][key])/1e-4
        row=dict(output=key,endpoint_contrast_score=score(contrast),affine_reconstruction_score=score(reconstructed-base[key]),
            original_score=score(fd-reported),corrected_frozen_score=score(fd-fixed),fresh_EOS=[])
        for h in [1e-4,5e-5]:
            minus=np.load(OLD/f'EOS-common-T-{h}--1-native.npz');plus=np.load(OLD/f'EOS-common-T-{h}-1-native.npz')
            error=score((plus[key]-minus[key])/(2*h)-full)
            row['fresh_EOS'].append(dict(h=h,score=error,passed=error<1e-3))
        arrays.update({key+'_correction':ap*contrast,key+'_corrected_frozen':fixed,key+'_corrected_full':full});rows.append(row)
    np.savez_compressed(OUT/'direct-blend-correction.npz',**arrays)
    save('direct-blend-correction.json',dict(classification='Counterexample candidate',rows=rows,
        actual_default_blend_limits=np.unique(low['weak_blend_before'],axis=0).tolist(),
        endpoint_outputs_overwritten=False,original_failures_preserved=True,physical_weak_rates_certified=False))
    print('DIRECT BLEND',rows,flush=True)


def symbolic():
    import sympy as s
    T,off,on,low=s.symbols('T off on low',positive=True);tab=s.Function('tab')(T)
    a=(T-off)/(on-off);rate=a*tab+(1-a)*low
    correction=T/(on-off)*(tab-low)
    assert s.simplify(T*s.diff(rate,T)-a*T*s.diff(tab,T)-correction)==0
    b=s.Function('b')(T);d=s.Function('d')(T);F=b+a*d
    assert s.simplify(T*s.diff(F,T)-T*(s.diff(b,T)+a*s.diff(d,T))-T*s.diff(a,T)*d)==0
    # Separate endpoint interpolation from the exact rate mixture derivative.
    rate0,rate1=s.symbols('rate0 rate1')
    assert s.simplify((1-a)*rate0+a*rate1-(rate0+a*(rate1-rate0)))==0
    D=lambda expr:T*s.diff(expr,T)
    assert s.simplify(D(D(D(rate)))-a*D(D(D(tab)))-T/(on-off)*(3*D(D(tab))+3*D(tab)+tab-low))==0
    A,B=s.symbols('A B',real=True);L=s.exp(A+B*T);z=B*T
    assert s.simplify(D(D(D(L)))-(z+3*z*z+z**3)*L)==0
    save('symbolic.json',dict(classification='Proven',blend_product_rule=True,
        affine_full_source_chain_rule=True,third_log_temperature_derivative_identity=True,
        formula='D_lnT F(T,a(T)) = partial_lnT F at fixed a + (da/dlnT) (F1-F0), when F is affine in a.',
        scope='Interior of a declared smooth blend branch. No derivative is assigned at a kink, and this algebra does not enclose numerical evaluator error.'))
    print('PASS symbolic weak blend chain rule',flush=True)


def table_plan():
    save('weak-table-plan.json',dict(classification='Counterexample candidate',
        prior_corrected_derivative_sha256=c.sha(OUT/'direct-blend-correction.json'),
        observation='The blend correction reduces the native vector error from 0.107 to 0.00214 but does not pass the unchanged 1e-3 gate. The smaller FD step increases the residual; the archived table interpolation uses single precision throughout.',
        method='Read the initialized spline coefficients, coordinates, nuclide IDs and low-temperature rates from our SHA-bound disposable native child. The table capture is read only. Compare this run with the original native values before using the captured coefficients for an independent polynomial derivative check.',
        limits='Keep both failed native derivative tests. A smooth polynomial made from stored coefficients is a declared numerical model, not the rounded executable function or physically certified weak rates.'))


def table_capture():
    assert (OUT/'weak-table-plan.json').exists();context()
    data=dict(np.load(OLD/'restored-state-17-4.npz'));aux=np.load(OLD/'EOS-common-replacement.npy');label='weak-table-capture'
    c.native.setup(label,data,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    np.save(OUT/(label+'-replacement.npy'),aux);shim.trace(label,aux,capture_weak_table=True)
    base=np.load(OLD/'EOS-common-native.npz');new=np.load(OUT/(label+'-native.npz'))
    errors={k:float(np.max(abs(new[k]-base[k])/np.maximum(abs(base[k]),1e-100))) for k in base.files}
    assert max(errors.values())<1e-10,errors
    table=np.load(OUT/(label+'-weak-table.npz'))
    save('weak-table-capture-audit.json',dict(classification='Counterexample candidate',
        same_inputs_relative_errors=errors,table_shapes={k:list(table[k].shape) for k in table.files},
        table_sha256=c.sha(OUT/(label+'-weak-table.npz')),passed=True,physical_rates_certified=False))
    print('WEAK TABLE CAPTURE',max(errors.values()),{k:table[k].shape for k in table.files},flush=True)


def linear_weak(table,T9,YeRho,indices,rounded,components=False):
    """Declared bilinear table, or the archived mixed-precision expressions."""
    dtype=np.float32 if rounded else np.float64
    xs=table['weak_reaction_t9s'].astype(dtype);ys=table['weak_reaction_lyerhos'].astype(dtype)
    x=np.clip(np.asarray(T9,dtype=dtype),xs[0],xs[-1]);y=np.clip(np.asarray(np.log10(YeRho),dtype=dtype),ys[0],ys[-1])
    ix=np.clip(np.searchsorted(xs,x,side='left')-1,0,len(xs)-2)
    jy=np.clip(((y-ys[0])/((ys[-1]-ys[0])/(len(ys)-1))).astype(int),0,len(ys)-2)
    xb=(x-xs[ix])/(xs[ix+1]-xs[ix]);xa=dtype(1)-xb
    yb=(y-ys[jy])/(ys[jy+1]-ys[jy]);ya=dtype(1)-yb
    fin=table['weak_reactions_data'][0].astype(dtype);values=[];derivatives=[]
    for ir in indices:
        count=3 if components else 2
        a=fin[ix,jy,:count,ir];b=fin[ix+1,jy,:count,ir]
        cc=fin[ix,jy+1,:count,ir];d=fin[ix+1,jy+1,:count,ir]
        # fx0/fx1 are double variables, but the RHS arithmetic is REAL32 in
        # the archived code. The final x interpolation is double arithmetic.
        fx0=(ya[:,None]*a+yb[:,None]*cc).astype(float)
        fx1=(ya[:,None]*b+yb[:,None]*d).astype(float)
        lograte=xa[:,None].astype(float)*fx0+xb[:,None].astype(float)*fx1
        dx=(fx1-fx0)/(xs[ix+1]-xs[ix])[:,None]
        rates=10.**lograte
        values.append(rates if components else rates.sum(1))
        derivatives.append(np.log(10)*x[:,None]*rates*dx if components else np.log(10)*x*np.sum(rates*dx,axis=1))
    return np.moveaxis(values,0,1),np.moveaxis(derivatives,0,1)


def linear_table_plan():
    save('linear-table-plan.json',dict(classification='Counterexample candidate',
        table_sha256=c.sha(OUT/'weak-table-capture-weak-table.npz'),
        observation='The actual native weak_bicubic flag is zero. Therefore audit the mixed-precision bilinear branch, not the available but inactive bicubic branch.',
        method='Reconstruct both actual blended endpoint contrasts from the captured table and compare free-electron versus fully-ionized electron density choices. First validate the mixed-precision expression against native source values; only then compare the same coefficient table evaluated in double precision. Retain the original blend-corrected failed FD verdict.',
        limits='A double-precision evaluation changes the numerical function while retaining table data and physical prescription. It does not retroactively certify native derivatives or provide physical table error bounds.'))


def linear_table_audit():
    assert (OUT/'linear-table-plan.json').exists()
    table=dict(np.load(OUT/'weak-table-capture-weak-table.npz'));assert table['weak_bicubic']==0
    names=list(c.ISO);indices=[];parents=[];daughters=[];rates=[]
    for i,(lhs,rhs,ir) in enumerate(zip(table['weak_lhs_nuclide_id'],table['weak_rhs_nuclide_id'],table['weak_reaclib_id'])):
        if ir>0 and names[lhs-1] in c.NAMES and names[rhs-1] in c.NAMES:
            indices.append(i);parents.append(c.NAMES.index(names[lhs-1]));daughters.append(c.NAMES.index(names[rhs-1]));rates.append(table['weak_lowt_rate'][ir-1])
    base=dict(np.load(OLD/'EOS-common-native.npz'));low=np.load(OUT/'direct-blend-pure_low-native.npz');high=np.load(OUT/'direct-blend-pure_table-native.npz')
    scale=np.maximum(1e-30,np.maximum(abs(base['T'][:,None]*base['dxdt_T']).max(1),abs(base['dxdt']).max(1)));active=abs(base['heat'])>1
    def vector(rate,data):
        out=np.zeros_like(data['X'])
        for k,(a,b) in enumerate(zip(parents,daughters)):
            flux=data['X'][:,a]*rate[:,k];out[:,a]-=flux;out[:,b]+=flux
        return out
    rows=[];arrays={}
    for electron in ['free','fully_ionized']:
        for rounded in [True,False]:
            ye=base['aux_used'][:,0] if electron=='free' else (base['X']/c.A)@c.Z
            rate,derivative=linear_weak(table,base['T']*1e-9,base['rho']*ye,indices,rounded)
            predicted=vector(rate-np.array(rates),base);actual=high['dxdt']-low['dxdt']
            score=abs(predicted-actual).max(1)/scale
            rows.append(dict(electron=electron,rounded=rounded,contrast_score=float(max(score[active])),worst_active_cell=int(np.argmax(np.where(active,score,-1)))))
            arrays[f'{electron}-{rounded}-contrast']=predicted
    np.savez_compressed(OUT/'linear-table-audit.npz',**arrays)
    save('linear-table-audit.json',dict(classification='Counterexample candidate',rows=rows,
        table_indices=indices,parent_species=[c.NAMES[i] for i in parents],daughter_species=[c.NAMES[i] for i in daughters],low_rates=rates,
        threshold_for_native_contrast_reconstruction=1e-10,physical_rates_certified=False))
    print('LINEAR TABLE',rows,flush=True)


def weak_vector(data,rounded):
    table=dict(np.load(OUT/'weak-table-capture-weak-table.npz'))
    spec=json.loads((OUT/'linear-table-audit.json').read_text())
    T9=data['T']*1e-9;rate,derivative=linear_weak(table,T9,data['rho']*data['aux_used'][:,0],spec['table_indices'],rounded)
    a=np.clip((T9-.01)/.01,0,1);ap=np.where((T9>.01)&(T9<.02),T9/.01,0.)
    value=a[:,None]*rate+(1-a[:,None])*np.array(spec['low_rates'])
    # The final net_derivs route zeros rate_screened_dT/dRho for weaklib.
    # Preserve that actual returned derivative when reconstructing native.
    dt=np.zeros_like(derivative) if rounded else a[:,None]*derivative+ap[:,None]*(rate-np.array(spec['low_rates']))
    out=np.zeros_like(data['X']);dout=np.zeros_like(out)
    for k,(parent,daughter) in enumerate(zip(spec['parent_species'],spec['daughter_species'])):
        i=c.NAMES.index(parent);j=c.NAMES.index(daughter)
        flux=data['X'][:,i]*value[:,k];dflux=data['X'][:,i]*dt[:,k]
        out[:,i]-=flux;out[:,j]+=flux;dout[:,i]-=dflux;dout[:,j]+=dflux
    return out,dout


def smooth_vector_plan():
    save('smooth-vector-plan.json',dict(classification='Counterexample candidate',
        prior_native_contrast_reconstruction_sha256=c.sha(OUT/'linear-table-audit.json'),
        reconstruction_gate=1e-10,derivative_gate=1e-3,
        model='Replace only the four identified tabulated weak vector contributions by bilinear interpolation in double precision, preserving the captured float32 table knots/data, low-T rates, linear blend, nuclide inventory and all other native terms. Add the blend product-rule term to the derivative.',
        validation='Compare analytic partial derivative plus the previously saved finite EOS auxiliary chain directions against the new function evaluated at all original frozen/fresh auxiliary temperature stencils. Preserve failed native and blend-only verdicts. The finite auxiliary chain is not an analytic composition Jacobian.',
        limits='This is a new declared numerical vector evaluator. Heat, neutrino, EOS physical accuracy, continuous-domain rigorous derivatives and GR evolution are separate gates.'))


def smooth_vector_audit():
    assert (OUT/'smooth-vector-plan.json').exists()
    prior=json.loads((OUT/'linear-table-audit.json').read_text())
    assert prior['rows'][0]['contrast_score']<1e-10
    data=dict(np.load(OLD/'EOS-common-native.npz'));w0,t0=weak_vector(data,True);w1,t1=weak_vector(data,False)
    f=data['dxdt']+w1-w0;partial=data['T'][:,None]*data['dxdt_T']+t1-t0
    chain=np.load(OLD/'input-attribution.npz');full=partial+chain['dxdt_Ye']+chain['dxdt_eta']
    scale=np.maximum(1e-30,np.maximum(abs(data['T'][:,None]*data['dxdt_T']).max(1),abs(data['dxdt']).max(1)))
    active=abs(data['heat'])>1;rows=[];arrays=dict(vector=f,partial_lnT=partial,full_lnT_finite_chain=full)
    for kind,h in [('frozen',5e-5),('fresh',1e-4),('fresh',5e-5)]:
        sides=[]
        for sign in [-1,1]:
            name=f'EOS-attribution-frozen-{sign}' if kind=='frozen' else f'EOS-common-T-{h}-{sign}'
            d=dict(np.load(OLD/(name+'-native.npz')));old,_=weak_vector(d,True);new,_=weak_vector(d,False)
            sides.append(d['dxdt']+new-old)
        fd=(sides[1]-sides[0])/(2*h);expected=partial if kind=='frozen' else full
        score=abs(fd-expected).max(1)/scale;worst=int(np.argmax(np.where(active,score,-1)))
        arrays[f'{kind}-{h}-FD']=fd
        rows.append(dict(auxiliary=kind,h=h,max_score=float(score[worst]),worst_active_cell=worst,
            failed_active_cells=int(sum((score>=1e-3)&active)),passed=bool(max(score[active])<1e-3)))
    np.savez_compressed(OUT/'smooth-vector-audit.npz',**arrays)
    save('smooth-vector-audit.json',dict(classification='Counterexample candidate',rows=rows,active_cells=int(sum(active)),
        max_value_change_score=float(max((abs(w1-w0).max(1)/scale)[active])),
        passed=all(row['passed'] for row in rows),native_failed_verdicts_preserved=True,
        native_derivative_certified=False,full_auxiliary_derivative_analytic=False,
        heat_and_neutrino_model_closed=False,physical_rates_certified=False))
    print('SMOOTH VECTOR',rows,flush=True)


def weak_thermal(data,rounded,alpha=None):
    table=dict(np.load(OUT/'weak-table-capture-weak-table.npz'));spec=json.loads((OUT/'linear-table-audit.json').read_text())
    T9=data['T']*1e-9;v,dv=linear_weak(table,T9,data['rho']*data['aux_used'][:,0],spec['table_indices'],rounded,True)
    total=v[:,:,:2].sum(2);dtotal=dv[:,:,:2].sum(2);qnu=v[:,:,2]/total
    dqnu=dv[:,:,2]/total-qnu*dtotal/total
    a=np.clip((T9-.01)/.01,0,1) if alpha is None else np.full(len(T9),alpha)
    ap=np.where((T9>.01)&(T9<.02),T9/.01,0.) if alpha is None else np.zeros(len(T9))
    rate=a[:,None]*total+(1-a[:,None])*np.array(spec['low_rates']);drate=a[:,None]*dtotal
    if not rounded: drate+=ap[:,None]*(total-np.array(spec['low_rates']))
    conv=np.float32(1.3806504e-16/1.602176487e-6) if rounded else 1.3806504e-16/1.602176487e-6
    temp=np.asarray(np.clip(T9,table['weak_reaction_t9s'][0],table['weak_reaction_t9s'][-1]),dtype=np.float32)*np.float32(1e9) if rounded else data['T']
    mu=data['aux_used'][:,1]*float(conv)*temp.astype(float)
    # net_derivs_support multiplies the rate-derivative vector by scalar Q
    # and Qnu, omitting their derivatives. Match that actual returned value
    # only in the rounded reconstruction; differentiate the new model fully.
    dmu=mu
    heat=np.zeros(len(T9));neu=heat.copy();ht=heat.copy();nt=heat.copy()
    for k,(parent,daughter) in enumerate(zip(spec['parent_species'],spec['daughter_species'])):
        i=c.NAMES.index(parent);qx=float(c.ISO[parent]['ex']-c.ISO[daughter]['ex'])
        weight=data['X'][:,i]/c.A[i]*c.reaction.QCONV
        heat+=weight*rate[:,k]*(qx+mu-qnu[:,k]);neu+=weight*rate[:,k]*qnu[:,k]
        ht+=weight*(drate[:,k]*(qx+mu-qnu[:,k])+(0 if rounded else rate[:,k]*(dmu-dqnu[:,k])))
        nt+=weight*(drate[:,k]*qnu[:,k]+(0 if rounded else rate[:,k]*dqnu[:,k]))
    if rounded: ht[:]=0;nt[:]=0
    return dict(heat=heat,neutrino=neu,heat_lnT=ht,neutrino_lnT=nt)


def thermal_table_plan():
    save('thermal-table-plan.json',dict(classification='Counterexample candidate',
        previous_vector_audit_sha256=c.sha(OUT/'smooth-vector-audit.json'),
        method='Independently reconstruct the four table weak heat/neutrino terms from the captured decay, capture, neutrino tables and archived Q/mu conventions. Compare the actual pure-table minus pure-low native controls before applying any energy correction.',
        reconstruction_gate=1e-10,derivative_gate=1e-3,
        limits='No replacement heat/neutrino function is accepted until native endpoint contrasts are reconstructed. Q/mu conventions are not a physical energy-reference certificate.'))


def thermal_table_audit():
    assert (OUT/'thermal-table-plan.json').exists()
    base=dict(np.load(OLD/'EOS-common-native.npz'));low=np.load(OUT/'direct-blend-pure_low-native.npz');high=np.load(OUT/'direct-blend-pure_table-native.npz')
    a=weak_thermal(base,True,0.);b=weak_thermal(base,True,1.);active=abs(base['heat'])>1;rows=[]
    for key in ['heat','neutrino']:
        scale=np.maximum(1,abs(base[key]));error=abs((b[key]-a[key])-(high[key]-low[key]))/scale
        rows.append(dict(output=key,max_score=float(max(error[active])),passed=bool(max(error[active])<1e-10),worst_active_cell=int(np.argmax(np.where(active,error,-1)))))
    save('thermal-table-audit.json',dict(classification='Counterexample candidate',rows=rows,
        passed=all(row['passed'] for row in rows),physical_energy_reference_certified=False))
    print('THERMAL TABLE',rows,flush=True)


def smooth_thermal_plan():
    save('smooth-thermal-plan.json',dict(classification='Counterexample candidate',
        prior_thermal_contrast_sha256=c.sha(OUT/'thermal-table-audit.json'),
        additional_source_finding='net_derivs_support passes Q and Qneu as scalar factors multiplying the rate derivative vector. Their temperature derivatives are omitted even though eval_weak computed them. Check this with actual endpoint derivative contrasts.',
        model='On the active weak-table branch T9>0.01, replace the same four heat/neutrino contributions with double-precision table rates, neutrino energies and electron chemical term. Add the derivatives of blend weights, Q and Qnu. Preserve other native contributions and the original low-T branch.',
        derivative_gate=1e-3,derivative_contrast_reconstruction_gate=1e-10,
        limits='Finite initial-state verification only. Physical atomic/nuclear data error, composition derivatives, EOS total-energy convention and a new GR trajectory remain independent obligations.'))


def smooth_thermal_audit():
    assert (OUT/'smooth-thermal-plan.json').exists()
    base=dict(np.load(OLD/'EOS-common-native.npz'));low=np.load(OUT/'direct-blend-pure_low-native.npz');high=np.load(OUT/'direct-blend-pure_table-native.npz')
    a=weak_thermal(base,True,0.);b=weak_thermal(base,True,1.);active=abs(base['heat'])>1
    contrast=(high['heat_T']-low['heat_T'])*base['T'];predicted=b['heat_lnT']-a['heat_lnT']
    scale=np.maximum(1,abs(base['heat_T']*base['T']));ec=float(max((abs(contrast-predicted)/scale)[active]))
    save('thermal-derivative-reconstruction.json',dict(classification='Counterexample candidate',
        score=ec,passed=ec<1e-10,criterion=1e-10,
        scope='Independent reconstruction of the returned endpoint derivative; a source-level omission argument alone does not prove the complete actual derivative formula.'))
    assert ec<1e-10,ec
    old=weak_thermal(base,True);new=weak_thermal(base,False);mask=base['T']*1e-9>.01
    def values(d):
        before=weak_thermal(d,True);after=weak_thermal(d,False);gate=d['T']*1e-9>.01
        return {k:d[k]+np.where(gate,after[k]-before[k],0.) for k in ['heat','neutrino']}
    corrected=values(base);partial=base['heat_T']*base['T']+np.where(mask,new['heat_lnT']-old['heat_lnT'],0.)
    chain=np.load(OLD/'input-attribution.npz');full=partial+chain['heat_Ye']+chain['heat_eta'];rows=[]
    arrays={**corrected,'heat_partial_lnT':partial,'heat_full_lnT_finite_chain':full}
    for kind,h in [('frozen',5e-5),('fresh',1e-4),('fresh',5e-5)]:
        sides=[]
        for sign in [-1,1]:
            name=f'EOS-attribution-frozen-{sign}' if kind=='frozen' else f'EOS-common-T-{h}-{sign}'
            sides.append(values(dict(np.load(OLD/(name+'-native.npz')))))
        fd=(sides[1]['heat']-sides[0]['heat'])/(2*h);expected=partial if kind=='frozen' else full
        error=abs(fd-expected)/scale;worst=int(np.argmax(np.where(active,error,-1)))
        rows.append(dict(auxiliary=kind,h=h,max_score=float(error[worst]),worst_active_cell=worst,
            failed_active_cells=int(sum((error>=1e-3)&active)),passed=bool(max(error[active])<1e-3)))
        arrays[f'{kind}-{h}-heat-FD']=fd;arrays[f'{kind}-{h}-neutrino-FD']=(sides[1]['neutrino']-sides[0]['neutrino'])/(2*h)
    np.savez_compressed(OUT/'smooth-thermal-audit.npz',**arrays)
    save('smooth-thermal-audit.json',dict(classification='Counterexample candidate',rows=rows,
        native_derivative_endpoint_contrast_error=ec,passed=all(row['passed'] for row in rows),
        native_derivative_certified=False,full_auxiliary_derivative_analytic=False,
        physical_energy_reference_certified=False,new_GR_trajectory_computed=False))
    print('SMOOTH THERMAL',ec,rows,flush=True)


def routing_plan():
    assert not (OUT/'weak-routing-plan.json').exists()
    for name in ['smooth-vector-audit.json','smooth-vector-audit.npz','thermal-derivative-reconstruction.json']:
        shutil.copy2(OUT/name,OUT/('pre-routing-'+name))
    for name in ['net_derivs.f90','net_derivs_support.f90']:
        source=c.reaction.OUT/'sources/net/private'/name;dest=OUT/'sources/net/private'/name
        shutil.copy2(source,dest)
    save('weak-routing-plan.json',dict(classification='Counterexample candidate',
        finding='The archived net_derivs weaklib branch explicitly assigns rate_screened_dT=0 and rate_screened_dRho=0. Therefore none of the tabulated weak-rate derivatives, including derivatives that eval_weak computed, reaches the final vector/heat derivative.',
        test='Verify that both native pure-table and pure-low controls have identical final dxdt_T and heat_T at every cell; reconstruct the old returned weak derivative as zero, and add the complete derivative of the declared double model.',
        previous_blend_only_and_pre_routing_results_preserved=True,
        derivative_gate=1e-3,physical_certification=False))


def routing_audit():
    assert (OUT/'weak-routing-plan.json').exists()
    low=np.load(OUT/'direct-blend-pure_low-native.npz');high=np.load(OUT/'direct-blend-pure_table-native.npz')
    rows={key:dict(differing_entries=int(np.count_nonzero(high[key]-low[key])),maximum_absolute_difference=float(np.max(abs(high[key]-low[key])))) for key in ['dxdt_T','heat_T','dxdt_rho','heat_rho']}
    assert all(row['differing_entries']==0 for row in rows.values()),rows
    save('weak-routing-audit.json',dict(classification='Counterexample candidate',rows=rows,
        actual_zero_temperature_and_density_route_confirmed=True,returned_values_not_overwritten=True))
    smooth_vector_audit();smooth_thermal_audit()


def corrected_values(data):
    before,_=weak_vector(data,True);after,_=weak_vector(data,False)
    a=weak_thermal(data,True);b=weak_thermal(data,False);mask=data['T']*1e-9>.01
    return dict(dxdt=data['dxdt']+(after-before),
        **{k:data[k]+np.where(mask,b[k]-a[k],0.) for k in ['heat','neutrino']})


def unified_source_audit():
    # Reevaluate the existing independently varied Ye/eta controls with the
    # same new function as the full-temperature control, avoiding mixed models.
    base=dict(np.load(OLD/'EOS-common-native.npz'));h=5e-5
    partials={'dxdt':np.load(OUT/'smooth-vector-audit.npz')['partial_lnT'],
        'heat':np.load(OUT/'smooth-thermal-audit.npz')['heat_partial_lnT']}
    chains={k:np.zeros_like(v) for k,v in partials.items()}
    for mode in ['Ye','eta']:
        sides=[corrected_values(dict(np.load(OLD/f'EOS-attribution-{mode}-{sign}-native.npz'))) for sign in [-1,1]]
        for key in chains: chains[key]+=(sides[1][key]-sides[0][key])/(2*h)
    active=abs(base['heat'])>1;rows=[];arrays=corrected_values(base)
    for key in partials:
        arrays[key+'_partial_lnT']=partials[key];arrays[key+'_total_lnT']=partials[key]+chains[key]
        scale=np.maximum(1e-30,np.maximum(abs(base['T'][:,None]*base['dxdt_T']).max(1),abs(base['dxdt']).max(1))) if key=='dxdt' else np.maximum(1,abs(base['T']*base['heat_T']))
        for step in [1e-4,5e-5]:
            sides=[corrected_values(dict(np.load(OLD/f'EOS-common-T-{step}-{sign}-native.npz'))) for sign in [-1,1]]
            fd=(sides[1][key]-sides[0][key])/(2*step);d=abs(fd-arrays[key+'_total_lnT'])
            score=(d.max(1) if key=='dxdt' else d)/scale;worst=int(np.argmax(np.where(active,score,-1)))
            rows.append(dict(output=key,h=step,score=float(score[worst]),worst_cell=worst,passed=bool(max(score[active])<1e-3)))
            arrays[f'{key}-{step}-FD']=fd
    np.savez_compressed(OUT/'unified-source-audit.npz',**arrays)
    save('unified-source-audit.json',dict(classification='Counterexample candidate',rows=rows,
        passed=all(row['passed'] for row in rows),cells=len(base['T']),burning_cells=int(sum(active)),
        new_function_used_in_all_auxiliary_controls=True,finite_chain_derivatives=True,
        full_composition_Jacobian_certified=False,continuous_domain_certified=False,
        physical_rates_and_EOS_certified=False,GR_trajectory_recomputed=False))
    print('UNIFIED SOURCE',rows,flush=True)


def weak_interval_plan():
    save('weak-interval-plan.json',dict(classification='Conjectural',
        target='Outward interval derivative and centered-difference truncation bound for the four declared double-table weak vector contributions at all initial burning cells, with rho, free Ye and X fixed.',
        step=5e-5,decimal_precision=40,
        exact_model='Treat the captured binary32 table coefficients/knots and binary64 initial T,rho,Ye,X, low-T rates as exact numbers. Temperature is T0 exp(s). Use exact bilinear log-rate interpolation, exp(log(10)*lograte), and the declared piecewise linear blend.',
        boundary='Require the whole temperature stencil to stay in one table and blend branch. Record any crossing as not certified rather than differentiating through a kink.',
        proof='On one branch L=exp(A+B*T9), D^3 L=(z+3*z^2+z^3)*L, z=B*T9. For a linear blend a, D^3[a*L+(1-a)*low]=a*D^3L+a_prime*(3*D^2L+3*DL+L-low). The centered derivative error is <=h^2 sup|D^3f|/6. Interval arithmetic also encloses the error of the saved floating derivative at the center.',
        excludes='Other native reactions, EOS auxiliary derivatives, composition derivatives, physical nuclear/EOS data, GR evolution and observations.'))


def weak_interval_audit():
    from mpmath import iv
    from fractions import Fraction as F
    plan=json.loads((OUT/'weak-interval-plan.json').read_text());iv.dps=plan['decimal_precision']
    def rational(endpoint):
        sign,man,exponent,_=endpoint._mpi_[0]
        return F((-1 if sign else 1)*man)*F(2)**exponent
    def upper(x): return max(abs(rational(x.a)),abs(rational(x.b)))
    def exact(x): return iv.mpf(float(x))
    table=dict(np.load(OUT/'weak-table-capture-weak-table.npz'));spec=json.loads((OUT/'linear-table-audit.json').read_text())
    data=dict(np.load(OLD/'EOS-common-native.npz'));_,analytic=weak_vector(data,False)
    active=np.where(abs(data['heat'])>1)[0];xs=table['weak_reaction_t9s'];ys=table['weak_reaction_lyerhos'];h=exact(plan['step'])
    rows=[]
    for i in active:
        t0=exact(data['T'][i])/10**9;lo=t0*iv.exp(-h);hi=t0*iv.exp(h)
        interval=iv.mpf([lo.a,hi.b]);off=iv.mpf('.01');on=iv.mpf('.02')
        if rational(hi.b)<F(1,100): branch='low'
        elif rational(lo.a)>F(1,50): branch='high'
        elif rational(lo.a)>F(1,100) and rational(hi.b)<F(1,50): branch='blend'
        else:
            rows.append(dict(cell=int(i),passed=False,reason='blend boundary crossing'));continue
        center=[iv.mpf(0) for _ in c.NAMES];third=[iv.mpf(0) for _ in c.NAMES]
        if branch!='low':
            j=int(np.clip(np.searchsorted(xs,float(t0.mid),side='left')-1,0,len(xs)-2))
            if rational(lo.a)<F(float(xs[j])) or rational(hi.b)>F(float(xs[j+1])):
                rows.append(dict(cell=int(i),passed=False,reason='table boundary crossing'));continue
            y=iv.ln(exact(data['rho'][i])*exact(data['aux_used'][i,0]))/iv.ln(10)
            if rational(y.b)<F(float(ys[0])): y=exact(ys[0])
            if rational(y.a)>F(float(ys[-1])): y=exact(ys[-1])
            k=int(np.clip(np.searchsorted(ys,float(y.mid),side='left')-1,0,len(ys)-2))
            assert F(float(ys[k]))<=rational(y.a)<=rational(y.b)<=F(float(ys[k+1]))
            yb=(y-exact(ys[k]))/(exact(ys[k+1])-exact(ys[k]));ya=1-yb
            for ir,parent,daughter,lowrate in zip(spec['table_indices'],spec['parent_species'],spec['daughter_species'],spec['low_rates']):
                sums=[]
                for temperature in [t0,interval]:
                    value=iv.mpf(0);d1=iv.mpf(0);d2=iv.mpf(0);d3=iv.mpf(0)
                    for channel in range(2):
                        grid=table['weak_reactions_data'][0,:,:,channel,ir]
                        fx0=ya*exact(grid[j,k])+yb*exact(grid[j,k+1]);fx1=ya*exact(grid[j+1,k])+yb*exact(grid[j+1,k+1])
                        B=iv.ln(10)*(fx1-fx0)/(exact(xs[j+1])-exact(xs[j]));z=B*temperature
                        val=iv.exp(iv.ln(10)*fx0+B*(temperature-exact(xs[j])))
                        value+=val;d1+=z*val;d2+=(z+z*z)*val;d3+=(z+3*z*z+z*z*z)*val
                    if branch=='blend':
                        a=(temperature-off)/(on-off);ap=temperature/(on-off)
                        d3=a*d3+ap*(3*d2+3*d1+value-exact(lowrate))
                        d1=a*d1+ap*(value-exact(lowrate))
                    sums.append((d1,d3))
                p=c.NAMES.index(parent);q=c.NAMES.index(daughter);weight=exact(data['X'][i,p])
                center[p]-=weight*sums[0][0];center[q]+=weight*sums[0][0]
                third[p]-=weight*sums[1][1];third[q]+=weight*sums[1][1]
        round_error=max(upper(center[j]-exact(analytic[i,j])) for j in range(len(c.NAMES)))
        truncation=max(upper(value)*F(plan['step'])**2/6 for value in third)
        scale=F(max(1e-30,float(max(abs(data['T'][i]*data['dxdt_T'][i]))),float(max(abs(data['dxdt'][i])))))
        rows.append(dict(cell=int(i),branch=branch,passed=bool((round_error+truncation)/scale<F(1,1000)),
            derivative_round_error_upper_exact=str(round_error),FD_truncation_upper_exact=str(truncation),
            total_score_upper_exact=str((round_error+truncation)/scale),total_score_upper=float((round_error+truncation)/scale)))
    result=dict(classification='Proven',rows=rows,passed=all(r['passed'] for r in rows),
        cells=len(rows),weak_contributions=4,table_sha256=c.sha(OUT/'weak-table-capture-weak-table.npz'),
        maximum_score_upper=max((r.get('total_score_upper',0) for r in rows)),
        interval_arithmetic='mpmath.iv, 40 decimal digits, exact rational output bounds',
        scope=plan['target'],excludes=plan['excludes'])
    save('weak-interval-audit.json',result)
    print('WEAK INTERVAL',result['passed'],len(rows),result['maximum_score_upper'],flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
