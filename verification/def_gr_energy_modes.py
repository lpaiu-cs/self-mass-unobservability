"""Mass-orthogonal rational projection of the proved weak GR energy.

The reduced oscillators use the analytic response to every original heat pole.
No eigenvalue clipping, damping or fitting. Projection error is tested on the
same four actual response histories before any physical interpretation.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu,spsolve_triangular
from scipy.linalg import cholesky_banded,eigh
import def_gr_energy_evolution as evolution

fem=evolution.fem;OUT=fem.OUT/'modes';write=fem.write;digest=evolution.digest
TIMES=np.arange(65)/64


class Projection:
    def __init__(self,model,frequency=128.):
        self.model=model;M=model.M;K=model.K
        scale=1/np.sqrt(M.diagonal());scaled=(diags(scale)@M@diags(scale)).tocsc()
        n=model.size;band=np.zeros((4,n))
        for k in range(4):band[k,:n-k]=scaled.diagonal(-k)
        cb=cholesky_banded(band,lower=True);self.C=(diags(1/scale)@diags([cb[k,:n-k] for k in range(4)],-np.arange(4),shape=(n,n))).tocsc()
        assert np.max(abs((self.C@self.C.T-M).data),initial=0)<1e-12*np.max(abs(M.data))
        self.sigma=frequency**2;self.shift=K+self.sigma*M
        self.scale=np.sqrt(self.shift.diagonal());D=diags(1/self.scale)
        self.lu=splu((D@self.shift@D).tocsc());self.extended=self.shift.astype(np.longdouble);self.error=0.
        heat=model.heat;flux=np.zeros(len(heat.edges));flux[heat.face_ids]=np.sum(heat.amplitude,axis=1)*model.original.radiation.geometry.tc
        self.forcing=model.load@flux
        b=spsolve_triangular(self.C,self.forcing,lower=True)
        self.initial=b/np.linalg.norm(b)

    def solve(self,rhs):
        x=self.lu.solve(rhs/self.scale)/self.scale
        for _ in range(2):x+=self.lu.solve(np.asarray((rhs.astype(np.longdouble)-self.extended@x.astype(np.longdouble))/self.scale,float))/self.scale
        defect=rhs.astype(np.longdouble)-self.extended@x.astype(np.longdouble)
        self.error=max(self.error,float(np.max(abs(defect)/(abs(self.extended)@abs(x)+abs(rhs)+1e-100))))
        assert self.error<1e-9
        return x

    def basis(self,size):
        Q=np.empty((self.model.size,size),order='F');Q[:,0]=self.initial
        for j in range(1,size):
            q=self.sigma*(self.C.T@self.solve(self.C@Q[:,j-1]))
            for _ in range(2):q-=Q[:,:j]@(Q[:,:j].T@q)
            norm=np.linalg.norm(q);assert norm>1e-14,(j,norm)
            Q[:,j]=q/norm
        orth=float(np.max(abs(Q.T@Q-np.eye(size))));assert orth<1e-10
        physical=spsolve_triangular(self.C.T.tocsr(),Q,lower=False)
        return Q,physical,orth


def energy_matrix(model,Q):
    # Evaluate the weak quadratic form directly. Subtracting large assembled
    # K diagonal/off-diagonal terms near almost coincident nodes loses digits.
    dx=np.diff(model.grid);local=np.array([(1-1/np.sqrt(3))/2,(1+1/np.sqrt(3))/2])
    padded=np.vstack([Q,np.zeros((1,Q.shape[1]))]);values=padded[model.indices]
    sample=((1-local)[None,:,None,None]*values[:-1,None,:,:]+local[None,:,None,None]*values[1:,None,:,:]).reshape(-1,2,Q.shape[1])
    derivative=np.repeat((values[1:]-values[:-1])/dx[:,None,None],2,axis=0)
    data=model.canonical_data;A=data[:,:4].reshape(-1,2,2);B=data[:,4:8].reshape(-1,2,2);C=data[:,8:12].reshape(-1,2,2)
    cov=derivative-np.einsum('nij,njk->nik',A,sample);weights=np.repeat(dx/2,2)
    left=cov.reshape(-1,Q.shape[1]);right=(weights[:,None,None]*np.einsum('nij,njk->nik',B,cov)).reshape(-1,Q.shape[1])
    K=left.T@right
    K+=sample.reshape(-1,Q.shape[1]).T@(weights[:,None,None]*np.einsum('nij,njk->nik',C,sample)).reshape(-1,Q.shape[1])
    skew=float(np.max(abs(K-K.T))/max(np.max(abs(K)),1e-100));assert skew<1e-12,skew
    return K,skew


def functions(k,t):
    """Entire real functions S=sin(sqrt(k)t)/sqrt(k), C=(1-cos)/k, R=(t-S)/k."""
    x=k*t*t
    if abs(x)<.1:
        S=t;C=t*t/2;R=t**3/6
        termS=S;termC=C;termR=R
        for j in range(1,12):
            termS*=-x/((2*j)*(2*j+1));S+=termS
            termC*=-x/((2*j+1)*(2*j+2));C+=termC
            termR*=-x/((2*j+2)*(2*j+3));R+=termR
        cosine=1-k*C
    elif k>0:
        w=np.sqrt(k);S=np.sin(w*t)/w;C=2*np.sin(w*t/2)**2/k;R=(t-S)/k;cosine=np.cos(w*t)
    else:
        w=np.sqrt(-k);S=np.sinh(w*t)/w;C=2*np.sinh(w*t/2)**2/(-k);R=(t-S)/k;cosine=np.cosh(w*t)
    return S,C,R,cosine


def modal_response(model,K,Q):
    k,V=eigh(K,check_finite=True);heat=model.heat;tc=model.original.radiation.geometry.tc
    lam=heat.rates*tc;faces=V.T@(Q.T@model.load);face_source=faces[:,heat.face_ids]
    coefficients=[];decay_bound=0.
    assert k.min()>-1.,float(k.min())  # previously checked K+M positive, not eigenvalue clipping
    for i,ki in enumerate(k):
        amp=face_source[i,:,None]*heat.amplitude*tc
        denom=lam*lam+ki
        assert np.all(denom>0)
        coefficients.append([np.sum(amp),np.sum(amp/lam),np.sum(amp/denom),np.sum(amp/(lam*denom))])
        decay_bound=max(decay_bound,float(np.sum(abs(amp)/denom)))
    modal=np.zeros((65,len(k)));velocity=modal.copy()
    for j,t in enumerate(TIMES[1:],1):
        for i,(ki,(a,b,c,d)) in enumerate(zip(k,coefficients)):
            S,C,R,cosine=functions(ki,t)
            # Exact pole response except exp(-lambda*t) bounded below.
            modal[j,i]=a*R-b*C+c*S-d*cosine
            velocity[j,i]=a*C-b*S+c*cosine+d*ki*S
    # At all saved positive times lambda*t exceeds1e12; bound the tiny omitted
    # exponential using exp(-700), an explicitly looser representable bound.
    assert lam.min()*TIMES[1]>700
    bound=decay_bound*np.exp(-700)
    return (Q@V)@modal.T,(Q@V)@velocity.T,dict(minimum_stiffness=float(k.min()),maximum_stiffness=float(k.max()),
        omitted_exponential_modal_velocity_bound=float(bound),minimum_pole_times_first_sample=float(lam.min()*TIMES[1]))


def series(model,K,Q,label):
    start=time.monotonic();u,v,meta=modal_response(model,K,Q)
    history=[];native_v=[];native_psi=[]
    for j,t in enumerate(TIMES):
        row,q,qt,cv,cf=evolution.readout(model,u[:,j],v[:,j],t);history.append(row);native_v.append(cv);native_psi.append(cf)
    flux,energy=model.heat.faces(1.)
    balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    np.savez_compressed(OUT/(label+'.npz'),grid=model.grid,indices=model.indices,q=q,qt=qt,u=u[:,-1],ut=v[:,-1],
        native_radius=model.original.native,native_velocity=np.array(native_v),native_scalar=np.array(native_psi),
        weights=model.original.weights,masks=np.array(model.original.masks),heat_energy=energy,heat_flux=flux)
    result=dict(classification='Counterexample candidate',size=Q.shape[1],history=history,heat_telescoping=balance,seconds=time.monotonic()-start,**meta)
    write(OUT/(label+'.json'),result);print('MODAL',label,result['seconds'],history[-1],flush=True)
    return result


def control():
    # The closed expression must satisfy the forced ODE and zero initial data.
    import sympy as s
    t,w,lam=s.symbols('t w lam',positive=True)
    y=(t-s.sin(w*t)/w)/w**2-(1-s.cos(w*t))/(lam*w**2)+(s.sin(w*t)/w+(s.exp(-lam*t)-s.cos(w*t))/lam)/(lam**2+w**2)
    assert s.simplify(s.diff(y,t,2)+w*w*y-t-(s.exp(-lam*t)-1)/lam)==0
    assert s.simplify(y.subs(t,0))==0 and s.simplify(s.diff(y,t).subs(t,0))==0
    for k in [0.,1e-12,9.,-1e-4,1e12]:
        S,C,R,co=functions(k,1.)
        assert abs((1-k*C)-co)<1e-12 and abs((1-k*R)-S)<1e-12
    return dict(classification='Proven',passed=True,identity='Every projected oscillator satisfies y_tt+k*y=t+(exp(-lambda*t)-1)/lambda, with zero initial data. Entire functions retain the k=0 limit and negative k without clipping.',scope='Identity for a projected oscillator; not a bound on omitted spatial/projection modes.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(120)
    paths=[Path(__file__),Path(evolution.__file__),Path(fem.__file__),Path(fem.canonical.__file__),evolution.BANK/'fine-bank.npz',evolution.BANK/'coarse-bank.npz',evolution.OUT/'result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='29fe9acd',
        claim='Resolve remaining time-propagation error of the complete canonical GR energy by analytic heat-pole oscillator response in a mass-orthonormal rational space.',
        method='Shift128, nested dimensions128/256/512, exact symmetric weak-form projection and exact pole response up to an explicit negligible exponential bound. No eigenvalue clipping or damping.',
        decision='Only if all original four histories meet relative2% and order1.5, perform fixed coefficient, outer3R and coarser-space contrasts. No larger basis or refined spatial grid on failure.',
        gates=dict(relative=.02,order=1.5,coefficient=.02,outer=.002,spatial=.02,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(pilot_size=64,production_maximum_size=512,maximum_models=4,hard_seconds=450,CPU_threads=1,memory_GB=3,new_EOS_calls=0,automatic_expansion=False),
        bindings={str(p):digest(p) for p in paths}))
    write(OUT/'control.json',control());start=time.monotonic();model=fem.Model(evolution.BANK/'fine-bank.npz');P=Projection(model)
    setup=time.monotonic()-start;start=time.monotonic();Q,physical,orth=P.basis(64);basis=time.monotonic()-start
    K,skew=energy_matrix(model,physical);case=series(model,K,physical,'pilot-64')
    forecast=4*(setup+64*basis+12*case['seconds'])*1.5
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,basis64_seconds=basis,
        series_seconds=case['seconds'],orthogonality=orth,stiffness_skew=skew,linear_residual=P.error,
        forecast_seconds=forecast,assumption='Four models; basis orthogonalization scales quadratically from64 to512, modal work factor12,50 percent margin. High-dimensional timing and contrasts unmeasured.'))
    print('FORECAST',forecast,flush=True);signal.alarm(0)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<450
    signal.alarm(450);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    model=fem.Model(evolution.BANK/'fine-bank.npz');P=Projection(model);Q,physical,orth=P.basis(512);K,skew=energy_matrix(model,physical)
    cases={str(n):series(model,K[:n,:n],physical[:,:n],'fine-'+str(n)) for n in [128,256,512]};comparisons={}
    for field in evolution.FIELDS:
        a,b,c=[np.array([row[field] for row in cases[str(n)]['history']]) for n in [128,256,512]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        comparisons[field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
    time_passed=all(row['last']<.02 and row['order']>1.5 for row in comparisons.values())
    checks=[dict(orthogonality=orth,stiffness_skew=skew,linear_residual=P.error)]
    if time_passed:
        for label,bank,outer,coarse in [('coefficient','coarse-bank.npz',2,False),('outer','fine-bank.npz',3,False),('spatial','fine-bank.npz',2,True)]:
            other=fem.Model(evolution.BANK/bank,outer,coarse);proj=Projection(other);_,pq,orth=proj.basis(512);kr,skew=energy_matrix(other,pq)
            cases[label]=series(other,kr,pq,label+'-512');checks.append(dict(orthogonality=orth,stiffness_skew=skew,linear_residual=proj.error))
            for field in evolution.FIELDS:
                c=np.array([row[field] for row in cases['512']['history']]);d=np.array([row[field] for row in cases[label]['history']])
                comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
    passed=time_passed and all(row.get('coefficient',1)<.02 and row.get('outer',1)<.002 and row.get('spatial',1)<.02 for row in comparisons.values())
    result=dict(classification='Counterexample candidate',actual_GR_fluid_scalar_evolved=True,time_passed=time_passed,passed=passed,
        comparisons=comparisons,checks=checks,endpoint=cases['512']['history'][-1],paths=list(cases),seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,
        physical_pressure_metric_reconstruction_checked=False,full_dynamic_charge_solved=False,
        limitation='Acceptance of this different complete spatial scheme additionally requires physical pressure/metric and boundary consistency checks. Coarse-bank sensitivity is not an independent EOS error certificate.')
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
