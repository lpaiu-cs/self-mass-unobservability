"""Source-aligned polynomial elements for the same regular canonical GR problem.

No new heat cells, EOS states or time interval. Composite quadrature preserves
all existing input/background knots while polynomial displacement resolves
the short acoustic response inside each original cell.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from numpy.polynomial import Polynomial,legendre
from scipy.sparse import coo_matrix,diags,vstack
from scipy.sparse.linalg import splu,spsolve_triangular
from scipy.linalg import cholesky_banded
import def_gr_energy_modes as modes
import def_gr_canonical_regular as regular

fem=modes.fem;OUT=fem.OUT.parent/'def-gr-spatial-repair';write=fem.write;digest=modes.digest
FIELDS=modes.evolution.FIELDS;BANK=modes.evolution.BANK


def basis(degree):
    nodes=np.r_[0.,(legendre.Legendre.basis(degree).deriv().roots()+1)/2,1.]
    polynomials=[]
    for i,x in enumerate(nodes):
        p=Polynomial.fromroots(np.delete(nodes,i));polynomials.append(p/p(x))
    return nodes,polynomials


class Model:
    def __init__(self,degree,bank,outer=2,quadrature=6):
        fem.canonical=regular;self.original=fem.base.Problem(bank,outer)
        bg=self.original.bg;heat=self.original.heat;heat.d={k:heat.d[k] for k in heat.d.files}
        self.bg=bg;self.heat=heat;self.degree=degree;self.outer=outer
        self.cells=bg.grid;self.local,self.polynomials=basis(degree)
        dx=np.diff(self.cells);self.grid=np.r_[(self.cells[:-1,None]+dx[:,None]*self.local[:-1]).ravel(),outer]
        self.surface_index=bg.surface_index*degree;nf=self.surface_index+1;n=len(self.grid);self.size=n+nf-1
        self.indices=np.full((n,2),-1,int);self.indices[:nf,0]=2*np.arange(nf);self.indices[:nf,1]=2*np.arange(nf)+1
        self.indices[nf:,1]=2*nf+np.arange(n-nf);self.indices[-1,1]=-1
        # Keep heat faces, mechanical interpolation knots and native material
        # interpolation knots in every quadrature, regardless of trial degree.
        cuts=np.unique(np.r_[self.cells,bg.saved['r']/bg.Rs,heat.rnative,heat.edges])
        cuts=cuts[(cuts>=0)&(cuts<=outer)];gx,gw=legendre.leggauss(quadrature)
        points=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel()
        self.weights=(np.diff(cuts)[:,None]*gw/2).ravel()
        self.data,self.points=fem.coefficients(bg,points)
        self.V,self.D=self.evaluation(points)
        A=self.data[:,:4].reshape(-1,2,2);B=self.data[:,4:8].reshape(-1,2,2)
        C=self.data[:,8:12].reshape(-1,2,2);W=self.data[:,12:16].reshape(-1,2,2)
        self.cov=[self.D[i]-sum(diags(A[:,i,j])@self.V[j] for j in range(2)) for i in range(2)]
        def form(left,coeff,right):
            return sum(left[i].T@diags(self.weights*coeff[:,i,j])@right[j] for i in range(2) for j in range(2)).tocsc()
        self.K=form(self.cov,B,self.cov)+form(self.V,C,self.V);self.M=form(self.V,W,self.V)
        gs=self.data[:,16:22].reshape(-1,2,3);hs=self.data[:,22:28].reshape(-1,2,3)
        src=fem.source_points(heat,self.points);g=[sum(diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        h=[sum(diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
        F=sum(self.cov[i].T@diags(self.weights*B[:,i,j])@g[j] for i in range(2) for j in range(2))
        F-=sum(self.V[i].T@diags(self.weights)@h[i] for i in range(2))
        # Same exact nodal heat momentum lift as the original GR formulation.
        r=self.grid[:nf];node=bg.sample(r);N,a=self.original.radiation.geometry.metric(r)
        ids=np.clip(np.searchsorted(heat.edges,r,side='right')-1,0,len(heat.edges)-2)
        frac=(r-heat.edges[ids])/(heat.edges[ids+1]-heat.edges[ids])
        interp=coo_matrix((np.c_[1-frac,frac].ravel(),(np.repeat(np.arange(nf),2),np.c_[ids,ids+1].ravel())),shape=(nf,len(heat.edges))).tolil()
        inv=np.divide(1.,r**3,out=np.zeros_like(r),where=r>0)
        for i in np.flatnonzero(r<heat.edges[1]):interp.rows[i]=[1];interp.data[i]=[1.];inv[i]=heat.edges[1]**-3
        inside=(r<1)&(node['p']>0);factor=np.zeros(nf);gr=fem.base.task.h.gr
        factor[inside]=gr.G*1e-7*inv[inside]/(4*np.pi*gr.C**4*self.original.radiation.geometry.R*N[inside]*a[inside]*np.exp(-8*node['phi'][inside]**2)*(node['e'][inside]+node['p'][inside]))
        hz=(diags(factor)@interp.tocsr()).tocoo()
        self.H=coo_matrix((hz.data,(self.indices[hz.row,0],hz.col)),shape=(self.size,len(heat.edges))).tocsc()
        self.load=(F+self.K@self.H).tocsc();self.nativeV,_=self.evaluation(self.original.native)

    def evaluation(self,r):
        ids=np.clip(np.searchsorted(self.cells,r,side='right')-1,0,len(self.cells)-2)
        width=np.diff(self.cells)[ids];x=(r-self.cells[ids])/width
        shapes=np.column_stack([p(x) for p in self.polynomials]);derivatives=np.column_stack([p.deriv()(x) for p in self.polynomials])/width[:,None]
        nodes=ids[:,None]*self.degree+np.arange(self.degree+1)[None,:];rows=np.repeat(np.arange(len(r)),self.degree+1)
        values=[];slopes=[]
        for field in range(2):
            columns=self.indices[nodes,field].ravel();valid=columns>=0
            values.append(coo_matrix((shapes.ravel()[valid],(rows[valid],columns[valid])),shape=(len(r),self.size)).tocsr())
            slopes.append(coo_matrix((derivatives.ravel()[valid],(rows[valid],columns[valid])),shape=(len(r),self.size)).tocsr())
        return values,slopes

    def energy_matrix(self,Q):
        # Stream quadrature blocks: p4 must not create several full nq-by512 arrays.
        K=np.zeros((Q.shape[1],Q.shape[1]));n=len(self.weights)
        for lo in range(0,n,1024):
            ix=slice(lo,min(lo+1024,n));v=[a[ix]@Q for a in self.V];c=[a[ix]@Q for a in self.cov]
            B=self.data[ix,4:8].reshape(-1,2,2);C=self.data[ix,8:12].reshape(-1,2,2);w=self.weights[ix]
            for i in range(2):
                for j in range(2):
                    if np.any(B[:,i,j]):K+=c[i].T@((w*B[:,i,j])[:,None]*c[j])
                    if np.any(C[:,i,j]):K+=v[i].T@((w*C[:,i,j])[:,None]*v[j])
        skew=float(np.max(abs(K-K.T))/max(np.max(abs(K)),1e-100));assert skew<1e-12
        return K,skew


class Projection(modes.Projection):
    def __init__(self,model):
        self.model=model;M=model.M;K=model.K;n=model.size
        scale=1/np.sqrt(M.diagonal());scaled=(diags(scale)@M@diags(scale)).tocsc()
        bandwidth=int(max(abs(M.tocoo().row-M.tocoo().col)));band=np.zeros((bandwidth+1,n))
        for j in range(bandwidth+1):band[j,:n-j]=scaled.diagonal(-j)
        cb=cholesky_banded(band,lower=True)
        self.C=(diags(1/scale)@diags([cb[j,:n-j] for j in range(bandwidth+1)],-np.arange(bandwidth+1),shape=(n,n))).tocsc()
        assert np.max(abs((self.C@self.C.T-M).data),initial=0)<1e-12*np.max(abs(M.data))
        self.sigma=128.**2;self.shift=K+self.sigma*M;self.scale=np.sqrt(self.shift.diagonal());D=diags(1/self.scale)
        self.lu=splu((D@self.shift@D).tocsc());self.extended=self.shift.astype(np.longdouble);self.error=0.
        heat=model.heat;flux=np.zeros(len(heat.edges));flux[heat.face_ids]=np.sum(heat.amplitude,axis=1)*model.original.radiation.geometry.tc
        self.forcing=model.load@flux;b=spsolve_triangular(self.C,self.forcing,lower=True);self.initial=b/np.linalg.norm(b)


def series(model,K,Q,label):
    u,v,meta=modes.modal_response(model,K,Q);history=[];velocities=[];scalars=[]
    r=model.original.native;N,a=model.original.radiation.geometry.metric(r);speed=a/N*r*fem.base.task.h.gr.C
    w=model.original.weights;masks=model.original.masks
    for j,t in enumerate(modes.TIMES):
        flux,energy=model.heat.faces(t);q=u[:,j]-model.H@energy;qt=v[:,j]-model.H@(flux*model.original.radiation.geometry.tc)
        cv=speed*(model.nativeV[0]@qt);cf=model.nativeV[1]@q
        history.append(dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(w@(cv*cv))),scalar_mass_RMS=float(np.sqrt(w@(cf*cf))),
            **{name:float(np.sqrt(w[m]@cv[m]**2/w[m].sum())) for name,m in zip(FIELDS[2:],masks)}))
        velocities.append(cv);scalars.append(cf)
    balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    np.savez_compressed(OUT/(label+'.npz'),grid=model.grid,cells=model.cells,indices=model.indices,q=q,qt=qt,
        native_radius=r,native_velocity=velocities,native_scalar=scalars,weights=w,masks=masks,heat_energy=energy,heat_flux=flux)
    result=dict(classification='Counterexample candidate',degree=model.degree,projection=Q.shape[1],history=history,heat_telescoping=balance,**meta)
    write(OUT/(label+'.json'),result);print('SPATIAL',label,history[-1],flush=True);return result


def control():
    errors=[]
    for degree in [1,2,4]:
        nodes,polys=basis(degree);x=np.linspace(0,1,31)
        value=np.array([p(x) for p in polys]);slope=np.array([p.deriv()(x) for p in polys])
        error=max(np.max(abs(value.sum(0)-1)),np.max(abs(slope.sum(0))))
        for k in range(degree+1):error=max(error,np.max(abs(nodes**k@value-x**k)))
        assert error<1e-12;errors.append(float(error))
    return dict(classification='Counterexample candidate',passed=True,polynomial_reproduction_errors=errors,scope='Basis partition, derivative and degree reproduction; no claim of actual GR response accuracy.')


def prepare():
    assert not (OUT/'plan.json').exists();signal.alarm(90);start=time.monotonic()
    paths=[Path(__file__),Path(regular.__file__),Path(modes.__file__),Path(fem.__file__),OUT/'inspection.json',BANK/'fine-bank.npz',BANK/'coarse-bank.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='95d99987',
        claim='Resolve the measured short acoustic response by polynomial displacement within the original source-aligned cells, then apply the same input to actual GR evolution.',
        fixed='Same cell faces, source poles, physical GR equations, zero initial data,0.230804955685s horizon and65 readout times. Composite6-point integration keeps all existing input/background knots. No new EOS or input fitting.',
        method='Compare polynomial degrees1,2,4; each actual path has nested128/256/512 mass-orthogonal propagation spaces. Polynomial evaluation at the original native radii, no smoothing or changed readout definition.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial_relative=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='First measure degree2/64 pilot. Run only the fixed degrees1,2,4. If any final projection or spatial gate fails, stop; otherwise run exactly three degree4 contrasts: coefficient, outer3R, quadrature8. Never add degree, cells or projection vectors automatically.',
        budget=dict(main_models=3,conditional_models=3,maximum_degree=4,maximum_projection=512,hard_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):digest(p) for p in paths}))
    write(OUT/'control.json',control());model=Model(2,BANK/'fine-bank.npz');P=Projection(model);setup=time.monotonic()-start
    start=time.monotonic();_,Q,orth=P.basis(64);basis_seconds=time.monotonic()-start
    start=time.monotonic();K,skew=model.energy_matrix(Q);series(model,K,Q,'pilot-p2-64');response_seconds=time.monotonic()-start
    forecast=1.5*(6*setup+10*64*basis_seconds+6*64*response_seconds)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,basis64_seconds=basis_seconds,
        response64_seconds=response_seconds,dofs=model.size,quadrature_points=len(model.weights),forecast_seconds=forecast,
        orthogonality=orth,stiffness_skew=skew,linear_residual=P.error,
        assumption='Quadratic basis/projection cost from64 to512; weighted ten degree2-equivalent basis builds, six preparations and response builds,50 percent allowance. Degree4 and exterior cost unmeasured; cap300s.'))
    signal.alarm(0);print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<300
    signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();cases={};propagation={};checks=[]
    for degree in [1,2,4]:
        model=Model(degree,BANK/'fine-bank.npz');P=Projection(model);_,Q,orth=P.basis(512);K,skew=model.energy_matrix(Q)
        rows={n:series(model,K[:n,:n],Q[:,:n],f'p{degree}-{n}') for n in [128,256,512]};cases[str(degree)]=rows[512];propagation[str(degree)]={}
        for field in FIELDS:
            a,b,c=[np.array([h[field] for h in rows[n]['history']]) for n in [128,256,512]]
            norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
            propagation[str(degree)][field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
        checks.append(dict(degree=degree,dofs=model.size,quadrature_points=len(model.weights),orthogonality=orth,skew=skew,linear_residual=P.error))
        del model,P,Q,K
    comparisons={}
    for field in FIELDS:
        a,b,c=[np.array([h[field] for h in cases[str(n)]['history']]) for n in [1,2,4]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        comparisons[field]=dict(previous=float(d1),last=float(d2),decreased=bool(d2<d1))
    prop_pass=all(v['last']<.02 and v['order']>1.5 for p in propagation.values() for v in p.values())
    spatial_pass=all(v['last']<.02 and v['decreased'] for v in comparisons.values())
    if prop_pass and spatial_pass:
        for label,bank,outer,nquad in [('coefficient','coarse-bank.npz',2,6),('outer','fine-bank.npz',3,6),('quadrature','fine-bank.npz',2,8)]:
            model=Model(4,BANK/bank,outer,nquad);P=Projection(model);_,Q,orth=P.basis(512);K,skew=model.energy_matrix(Q)
            cases[label]=series(model,K,Q,label+'-p4-512')
            for field in FIELDS:
                c=np.array([h[field] for h in cases['4']['history']]);d=np.array([h[field] for h in cases[label]['history']])
                comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
            checks.append(dict(label=label,orthogonality=orth,skew=skew,linear_residual=P.error));del model,P,Q,K
    passed=prop_pass and spatial_pass and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('quadrature',1)<.002 for v in comparisons.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,passed=passed,
        propagation_passed=prop_pass,spatial_passed=spatial_pass,propagation=propagation,comparisons=comparisons,checks=checks,
        paths=list(cases),seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
