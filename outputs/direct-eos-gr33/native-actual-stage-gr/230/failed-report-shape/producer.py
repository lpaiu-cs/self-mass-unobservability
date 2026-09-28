"""Keep the gas-to-conserved map inside the true high-precision native flux."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import continue_integer_native as prior
import continue_native_flux_precision as flux

OUT=Path('native-conserved-precision230-work');OLD=prior.OUT
precision=flux.precision;read,write,sha=prior.read,prior.write,prior.sha
joint,LD=prior.joint,prior.LD
CAPS=dict(prepare=180,check=900)


def conserved(m,g):
    g=precision.hp(g)
    z=np.array([g[:,2]*precision.hp(m.bu),g[:,3]*precision.hp(m.su),g[:,0]*precision.hp(m.eu),g[:,1]*precision.hp(m.nu)],object)
    z[2]+=precision.hp(m.kappa)*z[0]
    return z


def native_function():
    return precision.compile_function(precision.native_B,[('z=m.conserved(g)','z=conserved(m,g)')],dict(conserved=conserved))


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[]
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==117
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for part in ['photons','material']:(OUT/'sweep-1'/part).mkdir(parents=True)
    files += [OLD/n for n in ['rejected-joint-stage.npz','rejected-joint-stage.json','last-accepted-64.npz','coarse-receipt.json','failure-64.json','defect-location.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Determine and repair the true118stage B residual from the gas-to-conserved conversion before the high-precision native flux.',
        evidence='229applied the passing228proposal exactly and passed8linear solves but failed actual nonlinear8.977e-10. Almost all residual is B at atmospheric cell261, with stage state1.34e19/7.70e19and long-double ULP1/8. The60digit native owner currently receives an already-rounded g-to-conserved array.',
        method='Reproduce the saved full actual defect exactly. Evaluate the IDENTICAL conserved map using exact binary inputs at40/80digits before the same primitive/reconstruction/HLL native function. Compare the actual full equation and half/double constitutive probes. Preserve all coefficients, branches, source, grid, clocks, floor and acceptance gates. No accepted step in this diagnosis; if arithmetic correction is supported, apply it to actual RHS and nonlinear acceptance together.',
        gates=dict(precision=1e-25,stage=1e-12,physical_stage=1e-13,constitutive=.002),budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='229actual8proposals took201.8s. Model setup plus a few native evaluations should be below2minutes; allow15minutes with no Krylov or prefix replay.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    e,k,b=s.symbols('e k b');assert s.expand((e+k*b)-k*b-e)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Exact conserved energy coordinate identity only; no numerical or physical error theorem.'))


def check():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    m=prior.owner.Model(64);z=np.load(OLD/'rejected-joint-stage.npz');t,h=z['time'][()],z['step'][()]
    sol=z['solution'];v=z['initial'];cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    norm=np.linalg.norm(z['defect'])/read(OLD/'rejected-joint-stage.json')['equations'][-1]['relative']
    def defect():
        rates=[];m.precise_values={}
        for j,row in enumerate(sol.reshape(2,-1)):
            x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True)
            native=flux.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
            rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
        d=(sol.reshape(2,-1)-v-h*(joint.A@np.array(rates))).ravel()
        return flux.precise_defect(m,t,h,v,sol,d)
    old=defect();assert np.array_equal(old,z['defect']),'Unchanged full nonlinear defect must reproduce exactly'
    native=native_function();vectors=[];probes=[]
    for digits in [40,80]:
        # Only this process changes the owner; running227/224 keep frozen sources.
        def promoted(*a,**kw):
            with precision.mp.workdps(digits):return native(*a,**kw)
        precision.native_B=promoted
        d=defect();vectors.append(d)
        if digits==80:
            with precision.mp.workdps(80):
                for c,row in zip(joint.C,sol.reshape(2,-1)):
                    _,g=m.unpack(row);now=t+c*h;base=native(m,now,g,m.precise_tangent)
                    for scale in [.5,2.]:
                        value=native(m,now,g,m.precise_tangent,scale)
                        probes.append(float(np.sum(abs(value-base))/max(np.sum(abs(base)),precision.mp.mpf('1e-290'))))
    last=vectors[-1];physical=joint.physical_norm(m,last)/joint.scales(m,v,sol)
    result=dict(classification='Counterexample candidate',exact_saved_actual_defect=True,old_relative=float(np.linalg.norm(old)/norm),
        promoted_relative=[float(np.linalg.norm(d)/norm) for d in vectors],precision_change=float(np.linalg.norm(vectors[0]-last)/norm),
        arithmetic_effect=float(np.linalg.norm(last-old)/norm),physical=physical.astype(float).tolist(),constitutive_B=probes,
        controls_passed=bool(np.linalg.norm(vectors[0]-last)/norm<1e-25 and max(probes)<.002),
        saved_proposal_actual_stage_passed=bool(np.linalg.norm(last)/norm<1e-12 and max(physical)<1e-13),
        physical_state_accepted=False,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'comparison.npz',old_defect=old,promoted_defect=last)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['controls_passed']


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
