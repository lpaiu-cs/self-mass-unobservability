"""Differentiate HLL algebra directly, sampling only local constitutive jets.

Counterexample candidate. The failed finite-face path stays frozen. This uses
the same equation and minmod directions, with no subtraction of gross fluxes.
"""
from pathlib import Path
import inspect,json,resource,sys,time
import numpy as np
import resolve_native_material_face_response as face

run=face.run;OUT=face.OUT;read,write,sha=face.read,face.write,face.sha
LD,C,AMP=face.LD,face.C,face.AMP
CAPS=dict(check=40,pilot=60,production=450,residual=30)


def jets(m,k,label,V,h):
    key=(k,label,float(h))
    if not hasattr(m,'thermo_jets'):m.thermo_jets={}
    if key not in m.thermo_jets:
        f=m.model.flow;rho,v,lt,y=V;partials=[]
        for axis in [0,2,3]:
            samples=[]
            for sign in [-1,1]:
                q=V.copy()
                if axis==2:q[axis]+=sign*h
                else:q[axis]*=np.exp(sign*h)
                f.eos.y=q[3];samples.append(np.array(f.eos(q[0],q[2])[:3],dtype=LD))
            partials.append((samples[1]-samples[0])/(2*h))
        m.thermo_jets[key]=np.array(partials)
    return m.thermo_jets[key]


def conserved(m,k,label,V,dV,h):
    f=m.model.flow;a=np.asarray(f.base.af,LD);rho,v,lt,y=V;dr,dv,dt,dy=dV
    U,F,thermo=f.conserved(*V,a);p,u,gamma=map(lambda x:np.asarray(x,LD),thermo[:3])
    coords=np.array([np.divide(dr,rho,out=np.zeros_like(dr),where=rho>0),dt,
                     np.divide(dy,y,out=np.zeros_like(dy),where=y>0)])
    dp,du,dgamma=np.einsum('ajn,an->jn',jets(m,k,label,V,h),coords)
    root=np.sqrt(1-v*v);W=1/root;dW=W**3*v*dv;wm=v*v/(root*(1+root))
    D=rho*W;dD=dr*W+rho*dW
    den=np.maximum(rho,f.eos.floor);dden=np.where(rho>f.eos.floor,dr,0.)
    enthalpy=f.eos.cx+u+p/den
    H=rho*enthalpy;dH=dr*enthalpy+rho*(du+dp/den-p*dden/den**2)
    dS=dH*W*W*v+H*(2*W*dW*v+W*W*dv)
    dK=a*(f.eos.cx*(dD*wm+D*dW)+(dr*u+rho*du+dp)*W*W+(rho*u+p)*2*W*dW-dp)+(a-f.base.a0)*f.eos.cx*dD
    dU=np.array([dD,dS,dK])
    dF=np.array([dD*v+D*dv,dS*v+U[1]*dv+dp,(dK+a*dp)*v+(U[2]+a*p)*dv])
    cs=np.asarray(thermo[-1],LD);hd=np.maximum(H,LD('1e-100'))
    dcs=np.divide((dgamma*p+gamma*dp-gamma*p*np.where(H>1e-100,dH,0.)/hd)/hd,
                  2*cs,out=np.zeros_like(cs),where=cs>0)
    speed=[];dspeed=[]
    for sign in [-1,1]:
        num=v+sign*cs;den=1+sign*v*cs
        speed.append(num/den);dspeed.append(((dv+sign*dcs)*den-num*sign*(dv*cs+v*dcs))/den**2)
    return U[:3],F[:3],dU,dF,speed,dspeed


def extremum(values,derivatives,minimum):
    values=np.array(values);derivatives=np.array(derivatives)
    best=np.min(values,axis=0) if minimum else np.max(values,axis=0)
    masked=np.where(values==best,derivatives,np.inf if minimum else -np.inf)
    return best,(np.min(masked,axis=0) if minimum else np.max(masked,axis=0))


def flux_direction(m,k,L,R,dL,dR,probe):
    a=conserved(m,k,'left',L,dL,LD('1e-5')*probe)
    b=conserved(m,k,'right',R,dR,LD('1e-5')*probe)
    UL,FL,dUL,dFL,vl,dvl=a;UR,FR,dUR,dFR,vr,dvr=b
    zero=np.zeros(L.shape[1],dtype=LD)
    sl,dsl=extremum([zero,vl[0],vr[0]],[zero,dvl[0],dvr[0]],True)
    sr,dsr=extremum([zero,vl[1],vr[1]],[zero,dvl[1],dvr[1]],False)
    den=sr-sl
    base=np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),den,out=np.zeros_like(FL),where=den>0)
    num=dsr*FL+sr*dFL-dsl*FR-sl*dFR+(dsl*sr+sl*dsr)*(UR-UL)+sl*sr*(dUR-dUL)-base*(dsr-dsl)
    value=np.divide(num,den,out=np.zeros_like(num),where=den>0)*m.model.m.af*m.model.m.area
    return np.vstack([value,np.zeros(value.shape[1],dtype=LD)])


# Reuse the full primitive/join/gravity/donor implementation. Only the gross
# flux subtraction is replaced. The original failed producer remains intact.
source=inspect.getsource(face.tangent)
begin=source.index('    units=np.array(');end=source.index('    # Differentiate the true donor.',begin)
source=source[:begin]+'    flux=flux_direction(m,k,L,R,dL,dR,probe)*factor[:,None]\n'+source[end:]
namespace=dict(face.tangent.__globals__,flux_direction=flux_direction)
exec(compile(source,__file__,'exec'),namespace);tangent=namespace['tangent']


def paths(sweep):
    p,m=face.BASE_PATHS(sweep)
    return p,m.with_name('material-flux-tangent') if sweep else m


def initialize():
    face.paths=paths;face.tangent=tangent;face.initialize()


def check():
    failed=face.paths(1)[1]/'pilot-64.npz'
    assert not read(failed.with_suffix('.json'))['passed']
    write(OUT/'flux-tangent-plan.json',dict(classification='Counterexample candidate',
        failure='The local finite-face endpoint passed, but its actual4step path failed at the middle knot with0.927percent momentum derivative difference. Preserve that trajectory and reject production.',
        repair='Differentiate native HLL/conserved algebra analytically. Sample only local EOS p,u,gamma partials with independent normalized log-density,log-temperature,log-neutral steps. Keep the reused direct primitive map, exact directional minmod, shared face, true donor and fixed-geometry gravity.',
        admission='Check the failed middle and endpoint states at EOS derivative steps5e-6/1e-5/2e-5 and retain0.2percent total RHS gate. Check resolved finite HLL directions independently. Then actual4/8step prefixes and original2percent time/1e-8 conservation/1percent donor gates.',
        scope='Same retained response on the zero-additional-metric sweep only; no full EOS derivative bound, nonlinear star or final charge conclusion.',
        budgets=CAPS,original_total_seconds=640,CPU_threads=1,virtual_GiB=3,
        stop='No further probe,precision,grid or clock ladder. Any failed physical gate or forecast rejects production.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(face.__file__),Path(run.__file__),failed,failed.with_suffix('.json'),OUT/'local-face-plan.json',OUT/'local-face-check.json']}))
    initialize();m=run.c.Material(128,64);d=np.load(failed);rows=[]
    for i in [2,4]:
        t=d['t'][i];z=d['history_scaled'][i]
        rates=[m.rhs(t,z,h)[0] for h in [.5,1.,2.]]
        errors=[(np.sum(abs(r-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).tolist() for r in [rates[0],rates[2]]]
        rows.append(dict(index=i,time=float(t),half_nominal_double=errors))
        assert np.max(errors)<.002,rows[-1]
    # Independently compare the HLL chain rule on a resolved primitive direction.
    k=1;m.raw(k,np.zeros((4,m.n)),np.zeros((5,m.n)),0.)
    f=m.model.flow;V=m.raw(k,np.zeros((4,m.n)),np.zeros((5,m.n)),0.)[3]['primitive'].astype(LD)
    L,R,_,_=face.reconstruction(f,V,np.zeros_like(V),f.join_state,np.zeros(4))
    phase=np.arange(L.shape[1])+1.
    def direction(q):return np.array([q[0]*.2*np.sin(phase),2e-5*np.cos(phase),.3*np.cos(phase),q[3]*.2*np.sin(phase)],dtype=LD)
    dL,dR=direction(L),direction(R);dL[:,L[0]==0]=0;dR[:,R[0]==0]=0
    # Cache keys for this same knot carry the identical background L/R.
    exact=flux_direction(m,k,L,R,dL,dR,1.)[:3];h=1e-5
    fd=(face.face_flux(f,L+h*dL,R+h*dR)-face.face_flux(f,L-h*dL,R-h*dR))[:3]/(2*h)
    control=(np.sum(abs(fd-exact),axis=1)/np.maximum(np.sum(abs(exact),axis=1),LD('1e-100'))).astype(float).tolist()
    assert max(control)<.002,control
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,native_HLL_derivative_control=control,
                owner=m.face_owner_error,physical_branch_ratio=m.physical_branch_ratio,final_charge_conclusion='unadjudicated')
    write(OUT/'flux-tangent-check.json',result);print(json.dumps(result),flush=True)
    # Exact quotient chain rule for the two-wave HLL expression.
    import sympy as sp
    x=sp.symbols('x');s,r,Lv,Rv,Uv,Vv=sp.symbols('s r L R U V');ds,dr,dLx,dRx,dUx,dVx=sp.symbols('ds dr dL dR dU dV')
    N=r*Lv-s*Rv+s*r*(Vv-Uv);D=r-s
    changed=((r+x*dr)*(Lv+x*dLx)-(s+x*ds)*(Rv+x*dRx)+(s+x*ds)*(r+x*dr)*(Vv+x*dVx-Uv-x*dUx))/(D+x*(dr-ds))
    expected=(dr*Lv+r*dLx-ds*Rv-s*dRx+(ds*r+s*dr)*(Vv-Uv)+s*r*(dVx-dUx)-N/D*(dr-ds))/D
    assert sp.simplify(sp.diff(changed,x).subs(x,0)-expected)==0
    write(OUT/'flux-tangent-symbolic.json',dict(classification='Proven',passed=True,scope='HLL quotient chain rule on a selected wave-speed branch; exact ties use one-sided extrema. No uniform native-EOS or full coupled error theorem.'))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'flux-tangent-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+read(OUT/'probe-localization.json')['seconds']
        assert spent+CAPS[action]<=640,(spent,CAPS[action])
        if action=='check':check()
        else:
            assert read(OUT/'flux-tangent-check.json')['passed'] and read(OUT/'flux-tangent-symbolic.json')['passed']
            for p,h in read(OUT/'flux-tangent-plan.json')['bindings'].items():assert sha(p)==h,p
            paths(1)[1].mkdir(exist_ok=True);run.initialize=initialize;run.paths=paths
            if action in ['pilot','production']:run.material(action=='pilot')
            else:run.residual()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
