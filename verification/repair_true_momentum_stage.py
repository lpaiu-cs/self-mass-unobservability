"""Measure the true momentum arithmetic in the frozen rejected final stage."""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import continue_thermal_native as prior

OUT=Path('native-true-momentum237-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
joint=prior.joint;precision=prior.repair.precision;flux=prior.prior.flux
conserved=prior.prior.repair.conserved
CAPS=dict(prepare=180,check=1200)


def native_S(m,t,g,probe=1.,details=False):
    hp=precision.hp;material=m.material;field,rates=m.geometry(t);material.full_field=field
    z=conserved(m,g);k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k])
    face=np.zeros(m.n+1,dtype=object);gravity=np.zeros(m.n,dtype=object)
    for j,weight in [(k,1-w),(k+1,w)]:
        if weight:
            f,r,_=m.precise_tangent(material,j,z,probe);row=material.point(j)
            face+=hp(weight)*hp(f[1])
            # Preserve the original fixed metric-rate coefficient exactly.
            gravity+=hp(weight)*(hp(r)-hp(rates[0]+rates[1])*hp(row['Q'][1]))
    raw=-np.diff(face)+gravity;rate=raw/hp(m.su)
    return (rate,raw,face,gravity) if details else rate


def prepare():
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==118
    assert read(OLD/'controller-status.json')['state']=='failed';assert not OUT.exists();OUT.mkdir();files=[]
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for part in ['photons','material']:(OUT/'sweep-1'/part).mkdir(parents=True)
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','failure-64.json','coarse-receipt.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Correct the remaining actual119momentum arithmetic if a bit-exact reconstruction and independent precision checks demonstrate it, then apply to the actual final step.',
        evidence='235thermal repair accepted118at1.10532e-15, then119failed8actual Newton proposals at1.15246e-11. Full saved defect norm.01953125 is dominated by S(.01953125), with exact binary stage components-.01171875/.015625. B norm2.75e-7 is subdominant. Runtime1041.49s is below2hour cap; merely adding walltime does not remove this defect.',
        method='Reproduce the saved entire nonlinear defect. First hold original rounded native/collision rates fixed and promote only S Radau accumulation. Separately evaluate the identical S face/gravity source through234same high-precision native tangent and exact conserved map, preserving the original fixed metric coefficient. Compare40/80digits and0.5/2constitutive probes. No Jacobian surrogate for true native acceptance.',
        decision='A saved proposal remains unaccepted. Supported arithmetic must be applied consistently to actual Newton RHS, nonlinear defect, faces/gravity and momentum ledger before physical acceptance. Original118state/history and all gates remain.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=8,forecast='234similar native controls51.53s; allow20minutes. One rejected pair, no accepted-prefix replay, Krylov solve, new grid/clock or physical step.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    f0,f1,q0,q1,a,b,h,x,y=s.symbols('f0 f1 q0 q1 a b h x y')
    assert s.expand(x-y-h*(a*(f0+q0)+b*(f1+q1))-(x-y-h*(a*f0+b*f1)-h*(a*q0+b*q1)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Radau source and collision grouping identity only, not a physical error bound.'))


def check():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    m=prior.prior.prior.owner.Model(64);z=np.load(OLD/'rejected-joint-stage.npz');t,h=z['time'][()],z['step'][()]
    sol=z['solution'];v=z['initial'];cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    pairs=[m.unpack(row) for row in sol.reshape(2,-1)];norm=np.linalg.norm(z['defect'])/read(OLD/'rejected-joint-stage.json')['equations'][-1]['relative']
    rates=[];collisions=[];ordinary=[];m.precise_values={}
    for j,(x,g) in enumerate(pairs):
        ph,q,*_=m.collision(cs[j],x,g,True);native=flux.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
        collisions.append(q[:,3]);ordinary.append(native[0][:,3])
    original=(sol.reshape(2,-1)-v-h*(joint.A@np.array(rates))).ravel();original=flux.precise_defect(m,t,h,v,sol,original)
    assert np.array_equal(original,z['defect']),'Exact full actual defect reconstruction'
    def accumulate(native):
        hp=precision.hp;initial=m.unpack(v)[1][:,3];end=np.array([g[:,3] for _,g in pairs])
        value=hp(end)-hp(initial)-hp(h)*(hp(joint.A)@(native+hp(collisions)))
        d=original.copy()
        for row,new in zip(d.reshape(2,-1),value):m.unpack(row)[1][:,3]=precision.cast(new)
        return d
    with precision.mp.workdps(80):accumulator=accumulate(precision.hp(ordinary))
    vectors=[];native_values=[];probes=[]
    for digits in [40,80]:
        with precision.mp.workdps(digits):
            native=np.array([native_S(m,t+c*h,g) for c,(_,g) in zip(joint.C,pairs)])
            vectors.append(accumulate(native));native_values.append(precision.cast(native))
            if digits==80:
                for j,(c,(_,g)) in enumerate(zip(joint.C,pairs)):
                    for probe in [.5,2.]:
                        value=native_S(m,t+c*h,g,probe)
                        probes.append(float(np.sum(abs(value-native[j]))/max(np.sum(abs(native[j])),precision.mp.mpf('1e-290'))))
    new=vectors[-1];physical=joint.physical_norm(m,new)/z['physical_scales']
    result=dict(classification='Counterexample candidate',exact_saved_actual_defect=True,
        original_relative=float(np.linalg.norm(original)/norm),accumulation_only_relative=float(np.linalg.norm(accumulator)/norm),
        precise_native_S_relative=[float(np.linalg.norm(d)/norm) for d in vectors],
        precision_change=float(np.linalg.norm(vectors[0]-new)/norm),native_S_arithmetic_effect=float(np.linalg.norm(new-accumulator)/norm),
        physical=physical.astype(float).tolist(),constitutive_S=probes,
        controls_passed=bool(np.linalg.norm(vectors[0]-new)/norm<1e-25 and max(probes)<.002),
        saved_proposal_actual_stage_passed=bool(np.linalg.norm(new)/norm<1e-12 and max(physical)<1e-13),
        physical_state_accepted=False,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'comparison.npz',original=original,accumulation_only=accumulator,precise_native_S=new,native_S=native_values[-1])
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['controls_passed'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
