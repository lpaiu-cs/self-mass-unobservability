"""Local primitive directions through the same atmospheric HLL faces.

Counterexample candidate. Fixed additional geometry only. Preserve the native
EOS, exact directional minmod, shared face, conserved ledger and original gates.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import return_native_pressure_matter as run
import def_native_matter_photon_feedback as feedback

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha
LD=np.longdouble;C=LD(29979245800.);AMP=run.AMP
BASE_INIT=run.initialize;BASE_PATHS=run.paths
CAPS=dict(check=45,pilot=70,production=450,residual=40)
source=feedback.primitive_source.replace('raw=m.raw(k,np.zeros_like(z),np.zeros_like(field),0.);','')
assert source!=feedback.primitive_source
namespace=dict(feedback.namespace);exec(compile(source,__file__,'exec'),namespace)
primitive=namespace['primitive']


def minmod(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)


def minmod_direction(a,b,da,db):
    """One-sided derivative, including zero and equal-slope corners."""
    same=a*b>0
    tied=np.where(a>0,np.minimum(da,db),np.maximum(da,db))
    d=np.where(same,np.where(abs(a)<abs(b),da,np.where(abs(b)<abs(a),db,tied)),0.)
    d=np.where((a==0)&(b!=0)&(da*b>0),da,d)
    d=np.where((b==0)&(a!=0)&(db*a>0),db,d)
    return np.where((a==0)&(b==0),minmod(da,db),d)


def reconstruction(f,V,dV,join,djoin):
    # Same two ghost/background reconstructions as the existing owner.
    results=[]
    for background in [True,False]:
        v=V[:3]-f.background_cell if background else V[:3]
        ghost=join[:3]-f.background_left[:,0] if background else join[:3]
        x=np.column_stack([ghost,v,np.zeros(3)])
        dx=np.column_stack([djoin[:3],dV[:3],np.zeros(3)])
        a=x[:,1:-1]-x[:,:-2];b=x[:,2:]-x[:,1:-1]
        da=dx[:,1:-1]-dx[:,:-2];db=dx[:,2:]-dx[:,1:-1]
        s=np.zeros_like(x);ds=np.zeros_like(x)
        s[:,1:-1]=minmod(a,b);ds[:,1:-1]=minmod_direction(a,b,da,db)
        L=x[:,:-1]+s[:,:-1]/2;R=x[:,1:]-s[:,1:]/2
        dL=dx[:,:-1]+ds[:,:-1]/2;dR=dx[:,1:]-ds[:,1:]/2
        if background:L+=f.background_left;R+=f.background_right
        results.append((L,R,dL,dR))
    L,R,dL,dR=results[1]
    for a,b in zip(results[1],results[0]):a[:,:2]=b[:,:2]
    for v,d in [(L,dL),(R,dR)]:
        d[0]=np.where(v[0]>0,d[0],np.where(v[0]==0,np.maximum(d[0],0),0))
        v[0]=np.maximum(v[0],0)
    y=np.r_[join[3],V[3],f.eos.y0];dy=np.r_[djoin[3],dV[3],0.]
    s=np.zeros_like(y);ds=np.zeros_like(y)
    s[1:-1]=minmod(y[1:-1]-y[:-2],y[2:]-y[1:-1])
    ds[1:-1]=minmod_direction(y[1:-1]-y[:-2],y[2:]-y[1:-1],dy[1:-1]-dy[:-2],dy[2:]-dy[1:-1])
    L=np.vstack([L,y[:-1]+s[:-1]/2]);R=np.vstack([R,y[1:]-s[1:]/2])
    dL=np.vstack([dL,dy[:-1]+ds[:-1]/2]);dR=np.vstack([dR,dy[1:]-ds[1:]/2])
    # The common physical face uses the deep join state directly.
    L[:,0]=join;dL[:,0]=djoin
    return L,R,dL,dR


def face_flux(f,L,R):
    m=f.base
    UL,FL,tl=f.conserved(*L,m.af);UR,FR,tr=f.conserved(*R,m.af)
    sl=np.minimum(0,np.minimum((L[1]-tl[-1])/(1-L[1]*tl[-1]),(R[1]-tr[-1])/(1-R[1]*tr[-1])))
    sr=np.maximum(0,np.maximum((L[1]+tl[-1])/(1+L[1]*tl[-1]),(R[1]+tr[-1])/(1+R[1]*tr[-1])))
    return np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),sr-sl,out=np.zeros_like(FL),where=sr>sl)*m.af*m.area


def tangent(m,k,z,probe):
    row=m.point(k);zero=np.zeros_like(z);field=np.zeros((5,m.n));nb=m.nb
    raw=m.raw(k,zero,field,0.);f=m.model.flow;b=m.model.bulk
    V=raw[3]['primitive'].astype(LD);join=f.join_state.astype(LD).copy()
    if not hasattr(m,'face_banks'):m.face_banks={}
    if k not in m.face_banks:m.face_banks[k]=dict(np.load(feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz'))
    et=z.astype(LD).copy();et[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*et[0]
    p=primitive(m,k,et,field,m.face_banks[k])
    xi=m.model.mech.xi@np.r_[LD(0),-np.cumsum(z[0,:nb],dtype=LD)]
    p['dt'][:nb]+=b.eos.inventory[1]*xi/b.eos.gas(row['theta'],row['eta'])[2]
    dV=np.array([V[0]*p['dr'][nb:],p['dv'][nb:],p['dt'][nb:],V[3]*p['dy'][nb:]],dtype=LD)
    djoin=np.array([join[0]*p['dr'][nb-1],p['dv'][nb-1],p['dt'][nb-1],join[3]*p['dy'][nb-1]],dtype=LD)
    L,R,dL,dR=reconstruction(f,V,dV,join,djoin)
    baseline=face_flux(f,L,R)
    factor=4*np.pi*m.model.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH],dtype=LD)
    owner=float(np.max(np.sum(abs(baseline[:3]*factor[:3,None]-row['flux'][:3,nb:]),axis=1)/np.maximum(np.sum(abs(row['flux'][:3,nb:]),axis=1),1.)))
    assert owner<1e-8,('Reconstructed face owner',owner)
    m.face_owner_error=max(m.face_owner_error,owner)
    units=np.array([np.maximum.reduce([abs(L[0]),abs(R[0]),np.full(L.shape[1],f.eos.floor)]),
                    np.full(L.shape[1],1e-3),np.ones(L.shape[1]),np.maximum.reduce([abs(L[3]),abs(R[3]),np.full(L.shape[1],1e-30)])],dtype=LD)
    size=np.max(np.maximum(abs(dL),abs(dR))/units,axis=0)
    eps=np.divide(LD('1e-5')*probe,size,out=np.ones_like(size),where=size>0)
    values=[]
    for h in [eps,eps/2]:
        left=L+h*dL;right=R+h*dR
        assert np.all(left[0]>=0) and np.all(right[0]>=0)
        values.append((face_flux(f,left,right)-baseline)/h)
    flux=(2*values[1]-values[0])*factor[:,None]
    # Differentiate the true donor. The amplified probe cannot select it.
    base=row['flux'][0,nb:];donor=np.where(base==0,flux[0]>=0,base>=0)
    yy=np.where(donor,np.r_[join[3],V[3]],np.r_[V[3],f.eos.y0])
    dy=np.where(donor,np.r_[djoin[3],dV[3]],np.r_[dV[3],LD(0)])
    flux[3]=f.eos.nH*(flux[0]*yy+base*dy)
    df,dg=m.deep_tangent(k,z,field)
    F=np.c_[df,flux]
    # Exact fixed-geometry gravity derivative in physical conservative units.
    dE=(z[2,nb:].astype(LD)+LD(m.rest)*z[0,nb:])/(m.a[nb:]*m.V[nb:])
    g=m.model.m
    gravity=4*np.pi*C*np.diff(g.rf)*(-g.r*g.r*g.ap*dE+2*g.a*g.r*p['dp'][nb:])
    G=np.r_[dg,gravity]
    nonzero=row['flux'][0]!=0
    ratio=float(AMP*np.max(abs(F[0,nonzero]/row['flux'][0,nonzero]),initial=0))
    m.physical_branch_ratio=max(m.physical_branch_ratio,ratio);assert ratio<.01
    return F,G,row['dt']


def paths(sweep):
    p,m=BASE_PATHS(sweep)
    return p,m.with_name('material-local-face') if sweep else m


def initialize():
    run.original.paths=paths;run.paths=paths;BASE_INIT();Parent=run.c.Material
    class Material(Parent):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps);self.face_owner_error=0.
        def rhs(self,t,z,probe=1.):
            j,w,field,rates=self.fields(t)
            assert not np.any(field) and not np.any(rates),'This repair is scoped to the current zero-additional-metric sweep'
            F=np.zeros((4,self.n+1),dtype=LD);G=np.zeros((4,self.n),dtype=LD);dt=np.inf
            for k,weight in [(j,1-w),(j+1,w)]:
                if not weight:continue
                a,b,cfl=tangent(self,k,z,probe);F+=weight*a;G[1]+=weight*b;dt=min(dt,cfl)
            G+=(self.transfer[j+1]-self.transfer[j])/(self.t[j+1]-self.t[j])
            rate=-np.diff(F,axis=1)+G
            return np.asarray(rate,float),np.asarray(F[:,0]-F[:,-1]+np.sum(G,axis=1),float),dt
    run.c.Material=Material


def check():
    failed=BASE_PATHS(1)[1]/'pilot-64.npz'
    assert not read(failed.with_suffix('.json'))['passed']
    write(OUT/'local-face-plan.json',dict(classification='Counterexample candidate',checkpoint='0a4c0736c',
        claim='Resolve the failed same-input material derivative and then apply it to the actual4/8step prefix and full horizon.',
        change='Reuse the conservative-to-primitive map and deep tangent. Differentiate minmod analytically, including exact corners. Evaluate the same native HLL/EOS along locally normalized face directions with forward Richardson; direct fixed-geometry gravity variation and actual donor derivative. One shared face appears once.',
        scope='Current zero-additional-geometry material sweep only. This is a directional response of the retained operator, not full nonlinear evolution or a uniform native EOS certificate.',
        controls='Half/nominal/double local primitive probes5e-6/1e-5/2e-5 at the failed state, original0.2percent RHS threshold, reconstruction equality, primitive mapping control,1e-8 owner/conservation,1percent physical donor,2percent paired-time comparison. No global probe or precision ladder.',
        budgets=CAPS,total_original_seconds=640,CPU_threads=1,virtual_GiB=3,
        reuse='Byte-identical full173photon histories. Keep all failed174results. Resume new accepted material prefixes only.',
        stop='Any failed control stops this route; no automatic clock,grid,period,source or probe expansion.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),Path(feedback.__file__),failed,OUT/'primitive-precision-check.json',OUT/'plan.json']}))
    initialize();m=run.c.Material(128,64);d=np.load(failed);rows=[]
    # Exact corner controls for the minmod one-sided derivative.
    a=np.array([0,0,2,2,-2,-2,2,0.,1]);b=np.array([0,2,0,2,-2,0,-2,-2,3.])
    da=np.array([1,1,-1,3,-3,-1,1,1,2.]);db=np.array([2,3,-1,1,-1,-2,1,1,4.])
    exact=minmod_direction(a,b,da,db)
    assert np.max(abs((minmod(a+1e-6*da,b+1e-6*db)-minmod(a,b))/1e-6-exact))<1e-8
    # Owner reconstruction, including its special first two faces.
    for k in [0,1]:
        q=m.point(k);raw=m.raw(k,np.zeros((4,m.n)),np.zeros((5,m.n)),0.)
        f=m.model.flow;V=raw[3]['primitive'].astype(LD);join=f.join_state.copy()
        L,R,_,_=reconstruction(f,V,np.zeros_like(V),join,np.zeros(4))
        l,r=f.reconstruct(V,m.t[k]);l[:,0]=join
        err=float(max(np.max(abs(L-l)),np.max(abs(R-r))));assert err<1e-12
        # Independently push a resolved primitive direction through the native
        # conserved owner and invert its derivative using the reused map.
        active=q['active'][m.nb:];phase=np.arange(V.shape[1])+1.
        desired=np.array([.2*np.sin(phase),2e-5*np.cos(phase),.3*np.cos(phase),.2*np.sin(phase)])*active
        vals=[];h=1e-5
        for sign in [-1,1]:
            rr=V[0]*np.exp(sign*h*desired[0]);vv=V[1]+sign*h*desired[1]
            tt=V[2]+sign*h*desired[2];yy=V[3]*np.exp(sign*h*desired[3])
            vals.append(f.conserved(rr,vv,tt,yy,m.a[m.nb:])[0])
        dq=np.zeros((4,m.n),dtype=LD)
        dq[:,m.nb:]=(vals[1]-vals[0])/(2*h)*m.V[m.nb:]*f.eos.rho0*np.array([1,C*C,C*C,f.eos.nH])[:,None]
        dq[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*dq[0]
        bank=dict(np.load(feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz'))
        got=primitive(m,k,dq,np.zeros((5,m.n)),bank)
        actual=np.array([got[key][m.nb:] for key in ['dr','dv','dt','dy']])
        mapping=(np.sum(abs(actual-desired),axis=1)/np.maximum(np.sum(abs(desired),axis=1),1e-100)).tolist()
        assert max(mapping)<.002,('Conserved primitive map',k,mapping)
        rows.append(dict(k=k,reconstruction_absolute=err,owner=q['owner_error'],primitive_map=mapping))
    rates=[m.rhs(float(d['time']),d['delta_scaled'],h)[0] for h in [.5,1.,2.]]
    errors=[(np.sum(abs(r-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).tolist() for r in [rates[0],rates[2]]]
    np.savez_compressed(OUT/'local-face-rates.npz',rates=rates,state=d['delta_scaled'],t=d['time'])
    result=dict(classification='Counterexample candidate',passed=bool(np.max(errors)<.002),half_nominal_double=errors,
                reconstruction=rows,face_owner_error=m.face_owner_error,physical_branch_ratio=m.physical_branch_ratio,
                actual_material_prefix_completed=False,final_charge_conclusion='unadjudicated')
    write(OUT/'local-face-check.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'local-face-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+read(OUT/'probe-localization.json')['seconds']
        assert spent+CAPS[action]<=640,(spent,CAPS[action])
        if action=='check':check()
        else:
            assert read(OUT/'local-face-check.json')['passed']
            for p,h in read(OUT/'local-face-plan.json')['bindings'].items():assert sha(p)==h,p
            paths(1)[1].mkdir(exist_ok=True);run.initialize=initialize;run.paths=paths
            if action in ['pilot','production']:run.material(action=='pilot')
            else:run.residual()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
