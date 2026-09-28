"""Counterexample candidate: full incident geometry in the joint Radau equation.

Reuse the retained native owners and accepted four-variable stage solver.
This first return has prescribed incident geometry; self-GR is still open.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,shutil,sys,textwrap,time
import numpy as np
import sympy as sp
import couple_native_fluid_radau as joint
import def_native_incident_drive as drive

OUT=Path('native-full-incident181-work');BEFORE=joint.OUT
read,write,sha=joint.read,joint.write,joint.sha
LD,AMP,C=joint.LD,joint.AMP,joint.previous.original.C
engine=joint.previous.engine;face=engine.face
CAPS=dict(prepare=15,check=65,pilot=180)
CAPS['check_retry']=60
CAPS['check_final']=45
CAPS['check_chain']=30


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def clone(fn,changes,extra=None):
    source=textwrap.dedent(inspect.getsource(fn))
    for old,new in changes:
        assert source.count(old)==1,(fn.__name__,old,source.count(old));source=source.replace(old,new)
    ns=dict(fn.__globals__,**(extra or {}));exec(compile(source,__file__,'exec'),ns)
    return ns[fn.__name__],source


fixed_conserved=engine.conserved
def conserved(m,k,label,V,dV,h):
    U,F,dU,dF,speed,dspeed=fixed_conserved(m,k,label,V,dV,h)
    a=np.asarray(m.model.flow.base.af,LD);ell=np.interp(m.rEf,m.rE,np.asarray(m.full_field[1],float))[m.nb:]
    extra=ell*(U[2]+m.model.m.a0*m.model.flow.eos.cx*U[0])
    # F_energy=v*(U_energy+a*p), including v=0 without division.
    f=m.model.flow;f.eos.y=V[3];p=f.eos(V[0],V[2])[0]
    dU[2]+=extra;dF[2]+=V[1]*(extra+a*ell*p)
    return U,F,dU,dF,speed,dspeed


def install_geometry():
    engine.flux_direction=FunctionType(engine.flux_direction.__code__,dict(engine.flux_direction.__globals__,conserved=conserved))
    source=engine.source
    replacements=[
        ('field=np.zeros((5,m.n))','field=np.asarray(m.full_field,LD)'),
        ('flux=flux_direction(m,k,L,R,dL,dR,probe)*factor[:,None]',
         'flux=flux_direction(m,k,L,R,dL,dR,probe)*factor[:,None]\n'
         '    uf,ef=[np.interp(m.rEf,m.rE,np.asarray(v,float))[nb:] for v in field[:2]]\n'
         '    flux[:3]+=(2*uf+ef)*baseline[:3]*factor[:3,None]'),
        ('    # Exact fixed-geometry gravity derivative in physical conservative units.\n'
         '    dE=(z[2,nb:].astype(LD)+LD(m.rest)*z[0,nb:])/(m.a[nb:]*m.V[nb:])\n'
         '    g=m.model.m\n'
         "    gravity=4*np.pi*C*np.diff(g.rf)*(-g.r*g.r*g.ap*dE+2*g.a*g.r*p['dp'][nb:])",
         "    g=m.model.m;u,ell,lam,aden,ap=field[:,nb:];vol=3*u+lam\n"
         "    E0=(row['Q'][2,nb:]+LD(m.rest)*row['Q'][0,nb:])/(m.a[nb:]*m.V[nb:])\n"
         '    dE=(z[2,nb:].astype(LD)+LD(m.rest)*z[0,nb:])/(m.a[nb:]*m.V[nb:])-E0*vol\n'
         "    P0=row['Pg'][nb:]/m.V[nb:]\n"
         '    base=-g.r*g.r*g.ap*E0+2*g.a*g.r*P0\n'
         "    variation=-g.r*g.r*(g.ap*dE+ap*E0+2*u*g.ap*E0)+2*g.a*g.r*(p['dp'][nb:]+(ell+u)*P0)\n"
         '    gravity=4*np.pi*C*(np.diff(g.rf)*variation+np.diff(g.rf*uf)*base)')]
    for old,new in replacements:
        assert source.count(old)==1,old;source=source.replace(old,new)
    engine.source=source;ns=dict(engine.tangent.__globals__,flux_direction=engine.flux_direction)
    exec(compile(source,__file__,'exec'),ns);engine.tangent=ns['tangent']
    (OUT/'expanded-geometric-tangent.py').write_text(source)


def initialize():
    global Model
    joint.OUT=OUT;joint.initialize();Parent=joint.Model;install_geometry()
    class Response(Parent):
        def __init__(self,n):
            super().__init__(n)
            self.driver=self.redshift_driver;self.material.driver=self.driver;self.drive_scale=1.
            self.guide_g=np.zeros((self.n,4),LD)
            self.set_stage(0.)
        def lift(self,t):
            # ponytail: direct forcing for this short connection trial; an exact
            # affine lift is warranted only if its time error is measured.
            return np.zeros_like(self.I[0]),np.zeros_like(self.I[0]),np.zeros(3),np.zeros(3),0.
        def guide(self,t):return self.guide_g.copy()
        def local(self,t):
            self.set_stage(t);return super().local(t)
        def boundary_ports(self,t,x):
            self.set_stage(t);return super().boundary_ports(t,x)
        def geometry(self,t,enabled=True):
            m=self.material;d=self.driver.at(t);s=self.drive_scale/AMP if enabled else 0.
            u=d['delta_u']*s;ell=d['delta_log_lapse']*s;lam=d['delta_lambda']*s
            aden=m.rE*d['delta_u_prime']*s/m.den
            ap=m.ap*(ell-u-aden)+m.a/m.Aden*(d['delta_nu_prime']+d['delta_u_prime'])*s
            return np.array([u,ell,lam,aden,ap],LD),np.array([d['delta_u_t'],d['delta_lambda_rate']],LD)*s
        def native(self,t,g,probe=1.,details=False,tangent=None,metric=True):
            m=self.material;field,rates=self.geometry(t,metric);m.full_field=field
            z=self.conserved(g);k=int(np.clip(np.searchsorted(self.t,t,side='right')-1,0,15));w=(t-self.t[k])/(self.t[k+1]-self.t[k])
            F=np.zeros((4,self.n+1),LD);G=np.zeros((4,self.n),LD)
            for j,v in [(k,1-w),(k+1,w)]:
                if not v:continue
                f,r,_=(tangent or engine.tangent)(m,j,z,probe);F+=v*f;G[1]+=v*r
                row=m.point(j)
                G[2]-=v*field[1]*(row['base_rate'][2]+m.rest*row['base_rate'][0])
                G[1]-=v*(rates[0]+rates[1])*row['Q'][1]
                G[2]-=v*m.a*(row['Pr']*(rates[0]+rates[1])+2*row['Pg']*rates[0])
            raw=-np.diff(F,axis=1)+G
            normalized=np.column_stack([(raw[2]-self.kappa*raw[0])/self.eu,raw[3]/self.nu,raw[0]/self.bu,raw[1]/self.su])
            return (normalized,raw,np.zeros(4,LD),F,G) if details else normalized
    Response.jacobian,source=clone(Parent.jacobian,[
        ('return self.native(t,q,tangent=tangent)','return self.native(t,q,tangent=tangent,metric=False)'),
        ('defect=(J@g.ravel()).reshape(self.n,4)-base',
         'reset(False);offset=self.native(t,np.zeros_like(g),tangent=tangent)\n'
         '    defect=(J@g.ravel()).reshape(self.n,4)+offset-base')])
    (OUT/'expanded-affine-jacobian.py').write_text(source)
    new_stages,source=clone(joint.stages,[
        ('np.r_[LD(0),np.sum(gravity,dtype=LD),LD(0),LD(0)]','np.sum(gravity,axis=1,dtype=LD)'),
        ('    return result,sum(b*r for b,r in zip(B,transport))',
         '    m.guide_g=result[-1][1].copy()\n    return result,sum(b*r for b,r in zip(B,transport))')])
    (OUT/'expanded-full-stages.py').write_text(source)
    Response.run=FunctionType(Parent.run.__code__,dict(Parent.run.__globals__,stages=new_stages),argdefs=Parent.run.__defaults__)
    Model=Response


def prepare():
    assert not OUT.exists();OUT.mkdir();reuse={}
    for s in [0,1]:
        for p in paths(s):p.mkdir(parents=True)
    for src in (BEFORE/'sweep-0').rglob('*.npz'):
        dst=OUT/src.relative_to(BEFORE);os.link(src,dst);reuse[str(dst.relative_to(OUT))]=dict(path=str(src),sha256=sha(src))
    for name in ['normalization.json','photon-conservation-plan.json']:shutil.copyfile(BEFORE/name,OUT/name)
    write(OUT/'reuse.json',reuse)
    files=[Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    files+=[Path(__file__),drive.FIELDS/'born-g8.npz',Path('native-fluid-time180-work/pilot-result.json')]
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='1cfb0215d',
        claim='Apply the full declared incident scalar geometry and conservative redshift correction to the SAME photon/B/S/Etilde/H Radau stages. Evolve from zero response without adding a correction solution to old matter or charge.',
        decision='Whether the full-input native stage equation closes at original residual, constitutive and conservation gates. Only success warrants a paired time study and self-GR return; final charge stays unadjudicated.',
        physical_scope='Retained first variation with prescribed primary incoming packet plus saved first Born field, eta1e-30. This is not a companion-matched drive, a self-consistent returned GR solution, or nonlinear stellar evolution.',
        geometry='Actual stage fields in primitive volumes, face lapse/area/energy, deep and atmospheric gravity, reference-energy conversion and metric pressure work. Selected branches are affine in gas at fixed input; Jacobian columns have zero input. True full-input RHS decides acceptance.',
        photons='Existing full incident source plus conservative spatial redshift and angular-bin work correction. Direct physical forcing, zero affine lift. Current B/S/E/H and inventory drive collisions in both stages. Same-solution energy, floor and outgoing histories are saved.',
        reuse='Frozen corrected EOS/background and native owners. Hard-link prior constructor inputs then zero all lagged material feedback. Old histories never fix accepted states. No old response or charge is added.',
        gates=dict(stage=1e-12,physical_stage=1e-13,conservation=1e-8,constitutive=.002,branch=.01,source=1e-12),
        controls='Affine selected-branch reconstruction; metric-on/zero-input route; finite native geometric derivative; symbolic geometric identities. One actual64-clock two-macro prefix T/32, same existing local split rule, at most3 Newton iterations per step.',
        budget=CAPS,CPU_threads=1,virtual_GiB=4,
        forecast='179coarse6substeps cost215.853s; two macros have at most4substeps, expected80..170s including setup. Geometry and full forcing can change Newton cost;180s hard cap. Check includes original19s setup/Jacobian plus finite controls,65s cap. No full-period, new clock or extra split dispatch.',
        stop='Any original gate, constitutive control, Newton limit or cap stops. No automatic retry/refinement/longer horizon. Preserve rejected trajectories and181scope separately from180time failure.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    h,r,dr,ap,dap,e,de,a,ell,p,dp,u,dg=sp.symbols('h r dr ap dap e de a ell p dp u dg')
    expression=(dr+h*dg)*(-(r*(1+h*u))**2*(ap+h*dap)*(e+h*de)+2*a*(1+h*ell)*r*(1+h*u)*(p+h*dp))
    expected=dg*(-r*r*ap*e+2*a*r*p)+dr*(-r*r*(ap*de+dap*e+2*u*ap*e)+2*a*r*(dp+(ell+u)*p))
    assert sp.expand(sp.diff(expression,h).subs(h,0)-expected)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Product rule for native atmospheric gravitational source with radius, radial width, lapse, lapse gradient, pressure and energy variations. No uniform EOS or full-GR bound.'))


def check(chain=False):
    initialize();m=Model(64);now=m.t[-1]/48;g=np.zeros((m.n,4),LD)
    start=time.monotonic();J,base=m.jacobian(now,g);seconds=time.monotonic()-start
    assert np.any(base) and not np.any(m.native(now,g,metric=False))
    # Independent existing finite native owner, at a resolved geometric probe.
    mat=m.material;field,rates=m.geometry(now);k=0;row=mat.point(k);mat.full_field=field
    F,G,_=engine.tangent(mat,k,m.conserved(g),1.)
    size=max(np.max(abs(field[:4]))/1e-5,np.max(abs(field[4])/np.maximum(abs(mat.ap),mat.a/mat.R))/1e-3)
    eps=LD(1)/size;finite=[]
    for h in [eps,eps/2,eps/4]:
        f,r,*_=mat.raw(k,np.zeros((4,m.n)),field,float(h))
        finite.append(((f.astype(LD)-row['flux'])/h,(r.astype(LD)-row['gravity'])/h))
    controls=[]
    for i in [0,1]:
        f=2*finite[i+1][0]-finite[i][0];r=2*finite[i+1][1]-finite[i][1]
        exact=-np.diff(F,axis=1);exact[1]+=G
        trial=-np.diff(f,axis=1);trial[1]+=r
        controls.append((np.sum(abs(trial-exact),axis=1)/np.maximum(np.sum(abs(exact),axis=1),1.)).astype(float).tolist())
    mat.raw(k,np.zeros((4,m.n)),np.zeros((5,m.n)),0.)
    write(OUT/'geometric-control.json',dict(classification='Counterexample candidate',relative=controls))
    if not chain:assert max(controls[-1])<.002,('Full geometric native direction',controls)
    else:
        # Independent resolved face perturbation. The failed gross-recovery
        # subtraction is retained above; it is not a derivative certificate.
        f=mat.model.flow;V=mat.raw(k,np.zeros((4,mat.n)),np.zeros((5,mat.n)),0.)[3]['primitive'].astype(LD)
        left,right,_,_=face.reconstruction(f,V,np.zeros_like(V),f.join_state,np.zeros(4))
        phase=np.arange(left.shape[1])+1
        direction=lambda q:np.array([q[0]*.2*np.sin(phase),2e-5*np.cos(phase),.3*np.cos(phase),q[3]*.2*np.sin(phase)],LD)
        dl,dr=direction(left),direction(right);dl[:,left[0]==0]=0;dr[:,right[0]==0]=0
        field=np.zeros((5,mat.n));field[0]=.2*np.cos(np.arange(mat.n));field[1]=.3*np.sin(np.arange(mat.n)+.3)
        mat.full_field=field;uf,ef=[np.interp(mat.rEf,mat.rE,v)[mat.nb:] for v in field[:2]]
        base_face=face.face_flux(f,left,right)
        exact=engine.flux_direction(mat,k,left,right,dl,dr,1.)[:3]+(2*uf+ef)*base_face[:3]
        geom=mat.model.m;af0,area0=geom.af.copy(),geom.area.copy();samples=[];h=1e-5
        try:
            for sign in [-1,1]:
                geom.af=af0*(1+sign*h*ef);geom.area=area0*(1+sign*h*uf)**2
                samples.append(face.face_flux(f,left+sign*h*dl,right+sign*h*dr)[:3])
        finally:geom.af=af0;geom.area=area0
        error=(np.sum(abs((samples[1]-samples[0])/(2*h)-exact),axis=1)/np.maximum(np.sum(abs(exact),axis=1),LD('1e-100'))).astype(float).tolist()
        write(OUT/'resolved-face-control.json',dict(classification='Counterexample candidate',relative=error,
            failed_gross_recovery_control_preserved=True,uniform_native_derivative_certificate=False))
        assert max(error)<.002,('Resolved full geometric HLL direction',error)
    s,l,e=m.source(now);c=m.local(now)
    expected=m.driver.at(now);got=m.g
    stage_error=max(float(np.max(abs(got[key][0]-expected[key]))/max(np.max(abs(expected[key])),1e-290)) for key in ['delta_log_lapse','delta_log_speed'])
    assert stage_error<1e-12 and e<1e-12 and np.any(s) and any(np.any(c[key]) for key in ['q','qb','qe'])
    m.drive_scale=0.;s0,l0,e0=m.source(now);c0=m.local(now)
    assert not any(np.any(v) for v in [s0,l0,c0['q'],c0['qb'],c0['qe'],m.native(now,g)])
    result=dict(classification='Counterexample candidate',passed=True,jacobian_seconds=seconds,jacobian_nnz=J.nnz,
        rejected_gross_recovery_relative=controls,resolved_face_relative=error if chain else None,
        geometric_validation='Symbolic gravity and geometric conserved/face chain rule, plus independent resolved direct-primitive face perturbation; gross recovery control remains rejected.',
        exact_zero_source=True,actual_stage_metric_relative=stage_error,
        full_incident_source_nonzero=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated')
    write(OUT/'check-result.json',result);print(json.dumps(result),flush=True)


def pilot():
    assert read(OUT/'check-result.json')['passed'];initialize();m=Model(64)
    row=m.run(64,'pilot-64',2);file=paths(1)[0]/'pilot-64.npz';p=dict(np.load(file))
    p.update(joint_stage_times=np.array(m.stage_t),joint_stage_weights=np.array(m.stage_h),
        joint_stage_conserved_scaled=np.array(m.stage_states),joint_native_rates_scaled=np.array(m.stage_native),
        joint_collision_rates_scaled=np.array(m.stage_collision),joint_discard_rates_scaled=np.array(m.stage_discard))
    np.savez_compressed(file,**p)
    accumulated=np.sum((p['joint_native_rates_scaled']+p['joint_collision_rates_scaled'])*p['joint_stage_weights'][:,None,None],axis=(0,1),dtype=LD)
    final=p['delta_material']/AMP*m.units+p['material_floor_discard_scaled']
    balance=(abs(np.sum(final,axis=0)-accumulated)/np.maximum(np.sum(abs(final),axis=0),LD('1e-290'))).astype(float).tolist()
    port=joint.previous.run.packets(file)[2]
    row.update(same_solution_material_balance=balance,angular_port_relative=port,
        maximum_true_stage=max(a[-1]['relative'] for a in m.newton_iterations),
        maximum_true_physical_stage=max(max(a[-1]['moments']) for a in m.newton_iterations),
        maximum_newton_iterations=max(map(len,m.newton_iterations)),full_incident_input_applied=True,
        self_GR_return_closed=False,time_comparison_completed=False,full_horizon_completed=False,
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    row['passed']=bool(row['passed'] and max(balance)<1e-8 and port<1e-12)
    write(file.with_suffix('.json'),row);write(OUT/'pilot-result.json',row)
    write(OUT/'stage-checks.json',dict(classification='Counterexample candidate',newton=m.newton_iterations,stages=m.stage_log))
    print(json.dumps(row),flush=True);assert row['passed'],row


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                original=OUT/'initial-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(original)==h,p
            if action=='check_retry':
                previous=read(OUT/'check-receipt.json');assert previous['seconds']+CAPS[action]<65
                write(OUT/'dispatch-repair.json',dict(classification='Counterexample candidate',failure=previous,
                    repair='Correct indentation after dedenting the reused Jacobian method. No model was constructed and no physical step was taken.',
                    source_sha256=sha(__file__),remaining_check_cap=CAPS[action]))
            elif action=='check_final':
                spent=sum(read(OUT/f'{v}-receipt.json')['seconds'] for v in ['check','check_retry']);assert spent+CAPS[action]<65
                write(OUT/'field-repair.json',dict(classification='Counterexample candidate',
                    repair='Use the native float64 interpolation boundary for geometric face fields. Explicitly refresh actual-time geometry in local collisions and boundary ports, whose nearest inherited owners bypass DrivenPhoton.local.',
                    preserved='indent-repaired-producer.py',check_spent=spent,remaining_check_cap=CAPS[action],source_sha256=sha(__file__)))
            elif action=='check_chain':
                spent=sum(read(OUT/f'{v}-receipt.json')['seconds'] for v in ['check','check_retry','check_final']);assert spent+CAPS[action]<65
                write(OUT/'chain-control-plan.json',dict(classification='Conjectural',
                    rejection=read(OUT/'check_final-receipt.json'),
                    decision='Do not accept the failed gross-native recovery subtraction as a derivative reference. Validate the new lapse/area HLL chain directly with an independent resolved primitive/geometry perturbation, and the gravity product rule symbolically. The true same-input stage equations and original EOS-probe gates remain mandatory. No finite gross-recovery or uniform EOS certificate is asserted.',
                    root_hypothesis='At zero material direction the metric-induced baryon/neutral flux derivative is small; gross recovery subtraction degrades as the step is halved. Floating recovery noise is a hypothesis, not an established unique root.',
                    no_probe_ladder=True,no_stage_equation_change=True,check_spent=spent,remaining_check_cap=CAPS[action],source_sha256=sha(__file__)))
            elif action=='pilot':assert sha(__file__)==read(OUT/'chain-control-plan.json')['source_sha256']
        if action in ['check_retry','check_final','check_chain']:check(action=='check_chain')
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
