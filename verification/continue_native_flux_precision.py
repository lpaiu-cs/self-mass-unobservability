"""Counterexample candidate: precise native B flux in the actual coupled solve."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import continue_postfloor_native as prior
import precise_native_baryon as precision

OUT=Path('native-flux-precision202-work');OLD=prior.OUT
before,stable,owner,joint=prior.before,prior.stable,prior.owner,prior.joint
LD=joint.LD;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=120,check=180,controls=180,prefix=1200,coarse=7200,fine=10800,audit=120)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert 'True native joint Radau equation' in read(OLD/'coarse-receipt.json')['error']
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==114
    controls=read(Path('native-energy-coordinate201-work/fast-result.json'))
    assert controls['exact_match_to_full_bank_conversion'] and controls['precision_change']<1e-25
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for p in files:dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','coarse-receipt.json','failure-64.json']]
    files += [Path('native-energy-coordinate201-work')/n for n in ['result.json','roundoff-result.json','high-result.json','fast-result.json','fast-kernel.py']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='3da6f64b6',
        claim='Apply independent high-precision native mass flux to both the Newton RHS and the actual nonlinear residual of the remaining SAME coupled solution.',
        evidence='200failed the actual native stage after8proposals, despite passing all linear solves.199stored near-converged branches reproduce the actual native RHS exactly but its affine prediction differs by1.837e-11in the stage. Energy-coordinate round-trip change is zero. Independently reused primitive/reconstruction/HLL formulas at40and70digits agree to1.882e-41in that stage, while the saved proposal still fails. This admits a new solve, never acceptance of an old rejected state.',
        arithmetic='Use60digit true native B rates and shared face flux. Assemble only the B Newton RHS with the same high-precision native rate minus the stored sparse Jacobian action, and evaluate the actual B Radau defect before rounding that small defect. The independent native function chooses its own current branches; no fitted Jacobian is substituted for actual native acceptance.',
        covariance='Update raw B and its faces from the same high-precision flux. Retain raw reference-energy flux and adjust Etilde rate by minus kappa times the B-rate arithmetic change, preserving Eref=Etilde+kappa*B. Keep all other physical equations, photon maps, EOS banks, floor, period and clocks.',
        reuse='Restore114accepted states/histories exactly. Reuse200last rejected pair only as a Newton proposal with8fresh proposals for the repaired arithmetic. Re-evaluate all saved native B rates and the associated Etilde conversion, checking the complete prefix ledger before production. This does not recertify every old local vector equation or a uniform derivative bound.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='Independent high evaluation of two actual stages costs1.28s after setup.658saved stages imply about7minutes plus setup, with20minute prefix cap.200used1659s for8failed proposals of one step. Allow2/3hours for5coarse/16fine remaining substeps; later/fine costs remain uncertain, not an ETA. No accepted prefix reintegration.',
        stop='Any exact-restart, high-precision/constitutive, prefix ledger, original actual stage or paired-time failure;8newNewton/12linear limits or generous wall caps. No new grid,period or additional parameter path.',
        decision='Finish the actual two paths and original10channel time comparison; only then read their own complete GR fields. Self-GR and final infinity-normalized charge remain separate.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


def precise_native(m,t,g,native):
    with precision.mp.workdps(60):
        rate,flux=precision.native_B(m,t,g,m.precise_tangent,return_flux=True)
        new_raw=precision.cast(rate*precision.hp(m.bu));difference=new_raw-native[1][0]
        native[0][:,0]-=m.kappa*difference/m.eu
        native[0][:,2]=precision.cast(rate);native[1][0]=new_raw;native[3][0]=precision.cast(flux)
        m.precise_values[t]=(g.copy(),rate)
    return native


def precise_rhs(m,t,h,v,maps,guides,rhs):
    m.precise_values={}
    with precision.mp.workdps(60):
        affine=np.array([precision.native_B(m,t+c*h,g,m.precise_tangent)-precision.baryon_product(J,g)
            for c,g,(J,_) in zip(joint.C,guides,maps)])
        _,initial=m.unpack(v)
        value=precision.hp(initial[:,2])+precision.hp(h)*(precision.hp(joint.A)@affine)
        for row,new in zip(rhs.reshape(2,-1),value):m.unpack(row)[1][:,2]=precision.cast(new)
    return rhs


def precise_defect(m,t,h,v,sol,defect):
    pairs=[m.unpack(row)[1] for row in sol.reshape(2,-1)];rates=[]
    for c,g in zip(joint.C,pairs):
        saved,rate=m.precise_values[t+c*h];assert np.array_equal(saved,g);rates.append(rate)
    with precision.mp.workdps(60):
        _,initial=m.unpack(v)
        value=precision.hp(np.array(pairs)[:,:,2])-precision.hp(initial[:,2])-precision.hp(h)*(precision.hp(joint.A)@np.array(rates))
        for row,new in zip(defect.reshape(2,-1),value):m.unpack(row)[1][:,2]=precision.cast(new)
    return defect


def initialize(seed=True):
    init=FunctionType(prior.prior.initialize.__code__,dict(prior.prior.initialize.__globals__,OUT=OUT));init(False)
    run=owner.Model.run
    restore=FunctionType(prior.restore_substep.__code__,dict(prior.restore_substep.__globals__,OUT=OUT,restore_previous=run.__globals__['restore_substep']))
    def restored(m,previous,moment):
        values=restore(m,previous,moment)
        m.guide_g=np.load(OLD/'last-accepted-64.npz')['restart_guide'].copy()
        return values
    source=(OUT/'expanded-precise-native-stage.py').read_text()
    changes=[
        ('rhs=(np.tile(v,(2,1))+h*(A@src)).ravel();op=',
         'rhs=(np.tile(v,(2,1))+h*(A@src)).ravel();rhs=precise_rhs(m,t,h,v,maps,guides,rhs);op='),
        ('native=m.native(t+C[j]*h,gg,details=True)',
         'native=m.native(t+C[j]*h,gg,details=True);native=precise_native(m,t+C[j]*h,gg,native)'),
        ('defect=(sol.reshape(2,dim)-v-h*(A@np.array(rates))).ravel();relative=',
         'defect=(sol.reshape(2,dim)-v-h*(A@np.array(rates))).ravel();defect=precise_defect(m,t,h,v,sol,defect);relative=')]
    for old,new in changes:assert source.count(old)==1,old;source=source.replace(old,new)
    stage=run.__globals__['stages'];ns=dict(stage.__globals__,precise_rhs=precise_rhs,precise_native=precise_native,precise_defect=precise_defect)
    exec(compile(source,__file__,'exec'),ns)
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,restore_substep=restored,stages=ns['stages']),argdefs=run.__defaults__)
    (OUT/'expanded-high-native-stage.py').write_text(source)
    constructor=owner.Model.__init__
    def construct(m,n):
        constructor(m,n)
        with precision.mp.workdps(60):m.precise_tangent=precision.build(joint,owner)
        m.precise_values={}
        if seed and n==64:
            z=dict(np.load(OLD/'rejected-joint-stage.npz'));z['equations']=[];m.resume_seed=z
    owner.Model.__init__=construct


def check():
    fn=FunctionType(prior.check.__code__,dict(prior.check.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:initialize(False)));fn()


def controls():
    initialize(False);m=owner.Model(64);z=dict(np.load(OLD/'rejected-joint-stage.npz'))
    t,h=z['time'][()],z['step'][()];pairs=[m.unpack(row)[1] for row in z['solution'].reshape(2,-1)]
    _,initial=m.unpack(z['initial']);scale=np.linalg.norm(z['defect'])/read(OLD/'rejected-joint-stage.json')['equations'][-1]['relative']
    vectors=[];probes=[];mapping=[]
    for digits in [40,70]:
        with precision.mp.workdps(digits):
            tangent=precision.build(joint,owner)
            rates=np.array([precision.native_B(m,t+c*h,g,tangent) for c,g in zip(joint.C,pairs)])
            defect=precision.hp(np.array(pairs)[:,:,2])-precision.hp(initial[:,2])-precision.hp(h)*(precision.hp(joint.A)@rates)
            vectors.append(precision.cast(defect))
            if digits==70:
                for j,(c,g) in enumerate(zip(joint.C,pairs)):
                    ordinary=m.native(t+c*h,g)[:,2];nominal=precision.cast(rates[j])
                    mapping.append(float(np.sum(abs(nominal-ordinary))/max(np.sum(abs(nominal)),LD('1e-290'))))
                    for probe in [.5,2.]:
                        value=precision.cast(precision.native_B(m,t+c*h,g,tangent,probe))
                        probes.append(float(np.sum(abs(value-nominal))/max(np.sum(abs(nominal)),LD('1e-290'))))
    error=float(np.linalg.norm(vectors[0]-vectors[1])/scale)
    result=dict(classification='Counterexample candidate',passed=error<1e-25 and max(probes)<.002 and max(mapping)<1e-12,
        precision_change=error,constitutive_B=probes,ordinary_native_B_difference=mapping,
        saved_proposal_actual_B=[float(np.linalg.norm(v)/scale) for v in vectors],saved_proposal_accepted=False,new_physical_steps=0,
        final_charge_conclusion='unadjudicated')
    write(OUT/'controls.json',result);print(json.dumps(result),flush=True);assert result['passed'],result
    import sympy as sp
    e,k,b,d=sp.symbols('e k b d');assert sp.expand((e-k*d)+k*(b+d)-(e+k*b))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='The B-rate arithmetic correction preserves the reference-energy rate through Etilde_dot=Eref_dot-kappa*B_dot. No physical or uniform error certificate.'))


def prefix():
    initialize(False);rows=[]
    for n in [64,128]:
        m=owner.Model(n);p=OLD/'last-accepted-64.npz' if n==64 else OLD/'sweep-1/photons/interval-15-128.npz'
        z=dict(np.load(p));rates=[];changes=[]
        for j,t in enumerate(z['joint_stage_times']):
            q=z['joint_stage_conserved_scaled'][j]
            g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
            m.precise_values={};value=precise_native(m,t,g,m.native(t,g,details=True))[0]*m.units
            old=z['joint_native_rates_scaled'][j];rates.append(value)
            changes.append(np.sum(abs(value-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),LD('1e-290')))
            if (j+1)%32==0:write(OUT/f'prefix-progress-{n}.json',dict(completed=j+1,total=len(z['joint_stage_times'])))
        expected=np.sum(z['joint_stage_weights'][:,None,None].astype(LD)*(np.array(rates)+z['joint_collision_rates_scaled']),axis=0,dtype=LD)
        actual=z['g' if n==64 else 'restart_g']*m.units+z['material_floor_discard_scaled']
        error=np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))
        row=dict(clock=n,actual_stage_count=len(rates),native_rate_change=np.max(changes,axis=0).astype(float).tolist(),
            same_prefix_material_balance=error.astype(float).tolist(),passed=bool(max(error)<1e-8 and np.max(changes)<.002))
        rows.append(row);write(OUT/f'prefix-{n}.json',row);assert row['passed'],row
    write(OUT/'prefix-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        new_physical_steps=0,old_local_vector_audits_preserved=True,uniform_error_certificate=False,final_charge_conclusion='unadjudicated'))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'] and read(OUT/'controls.json')['passed'] and read(OUT/'prefix-result.json')['passed']
            fn=FunctionType(prior.evolve.__code__,dict(prior.evolve.__globals__,OUT=OUT,initialize=initialize));fn(64 if action=='coarse' else 128)
        elif action=='audit':
            fn=FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT));fn()
            r=read(OUT/'result.json');r.update(high_precision_native_B_applied=True,native_precision_digits=60);write(OUT/'result.json',r)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
